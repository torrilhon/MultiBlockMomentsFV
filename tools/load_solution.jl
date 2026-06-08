# load_solution.jl
# ----------------------------------------------------------------------------
# Loader and dense-output for the multi-block BINARY (.hdr/.dat) output of
# MultiBlockMomentsFV.
#
# Milestone 1: load the entire result (all blocks, all time instances) into a
#              single Julia data structure -- no plotting.
# Milestone 2: a "dense-output" function  value(x, y, t, variable)  that locates
#              the relevant block and grid cell for the point (x,y), and
#              interpolates both in time and in space.
#
# ---- Output file format (two-file sidecar) ---------------------------------
# Each run writes a pair of files sharing one basename:
#
#   <name>.hdr   pure ASCII run manifest (human-readable / greppable)
#   <name>.dat   pure binary, contiguous Float64 payloads, no embedded text
#
# .hdr structure:
#   # MultiBlockMomentsFV binary output header
#   # FORMAT_VERSION = 1
#   # CREATED_UTC = ... / GIT_COMMIT = ...
#   [RUN]      solver=.. Kn=.. R13=.. Type=.. FullGridSize=.. refine=..
#              t_output=0.0 t1 t2 ..   (the output times; 0.0 is the IC, then the post-IC times)
#   [BINARY]   BINARY_FILE=<name>.dat NVARS=10 ORDER=variable-fastest
#              PRECISION=Float32|Float64 ENDIAN=little
#   [VARIABLES]
#   rho vx vy T p Pxx Pxy Pyy qx qy        (10 PHYSICAL fields; x,y NOT stored)
#   [ZONES]
#   ZONE block=<n> I=<ni> J=<nj> NVARS=10 BX=.. BY=.. DX=.. DY=.. LFB_I=..
#        LFB_J=.. ISTART=.. JSTART=.. NGHOST=1 SOLUTIONTIME=<t> STRANDID=<id>
#        OFFSET=<o> NBYTES=<b>
#   ZONE ...
#
#   * Each ZONE is one (time, block) pair.
#   * STRANDID groups all blocks belonging to the *same* output time instance.
#   * Each block carries ONE layer of ghost cells on every side, so for an
#     interior grid of nstx x nsty cells the zone has  I = nsty+2,  J = nstx+2.
#   * I = points in y-direction (fastest), J = points in x-direction (slowest).
#   * OFFSET/NBYTES locate the zone payload inside the .dat file, so the reader
#     seeks directly to any block of any time with no scanning.
#   * Payload is VARIABLE-FASTEST: for each point the 10 PHYSICAL variables are
#     stored contiguously, points ordered x-slowest / y-fastest (matching I,J).
#     A zone thus reads as reshape(read, NVARS, I*J): columns = points.
#   * Coordinates x,y are NOT stored. They are reconstructed on load from the
#     per-zone uniform-grid geometry: for grid index (i along x, j along y),
#       x = BX + (i - LFB_I)*DX,   y = BY + (j - LFB_J)*DY
#     with i = ISTART + p÷I, j = JSTART + p%I for point p (0-based). These are
#     cell *centres*; the first centre sits at xmin - h/2 (the ghost) and the
#     last at xmax + h/2. Reconstructed x,y are prepended as columns 1,2 so the
#     in-memory layout still exposes the full 12 columns x,y,<10 physical>.
#   * Payload elements are Float32 or Float64 per PRECISION; they are read at
#     the stored precision and widened to Float64 in memory.
#
# Because the blocks have different grid sizes, the result is stored ragged
# (one array per block) rather than as one padded dense array; the conceptual
# 5-D access  data[time, block, x, y, variable]  is provided via getindex.
# ----------------------------------------------------------------------------

module LoadSolution

export SolutionData, load_solution, value, point_values, variable_index,
       block_bounds, gridsize, find_block, domain_extent

# ----------------------------------------------------------------------------
# Data structure  (ragged, with conceptual 5-D access)
# ----------------------------------------------------------------------------
"""
    SolutionData

Full multi-block, multi-time result of one `.dat` output file.

Metadata
- `title::String`                    dataset title
- `varnames::Vector{String}`         column names, incl. "x","y"
- `times::Vector{Float64}`           absolute solution times (length `nt`)
- `strandids::Vector{Int}`           STRANDID per time instance (length `nt`)
- `blockids::Vector{Int}`            block numbers as written in the file
- `nt, nblk, nvar::Int`              counts
- `nx, ny::Vector{Int}`              INTERIOR cells per block (no ghosts)
- `nxg, nyg::Vector{Int}`            grid points per block INCLUDING ghosts
- `xc, yc::Vector{Vector{Float64}}`  cell-centre coordinates per block (incl. ghosts)
- `bbox::Matrix{Float64}`            interior bounding box per block,
                                     row = block, cols = [xmin xmax ymin ymax]

Field data (ragged)
- `blocks::Vector{Vector{Array{Float64,3}}}`
        `blocks[it][ib]` is an array of size `(nxg[ib], nyg[ib], nvar)`,
        indexed `[ix, iy, ivar]` with ghost cells included.

Conceptual 5-D access is supported:
    sd[it, ib, ix, iy, iv]     # scalar
    sd[it, ib]                 # the (nxg,nyg,nvar) array for that time/block
both with 1-based indices (block position `ib`, not the file's block number).
"""
struct SolutionData
    title::String
    varnames::Vector{String}
    varsyms::Tuple{Vararg{Symbol}}   # varnames as Symbols, for NamedTuple keys
    times::Vector{Float64}
    strandids::Vector{Int}
    blockids::Vector{Int}
    nt::Int
    nblk::Int
    nvar::Int
    nx::Vector{Int}
    ny::Vector{Int}
    nxg::Vector{Int}
    nyg::Vector{Int}
    xc::Vector{Vector{Float64}}
    yc::Vector{Vector{Float64}}
    bbox::Matrix{Float64}
    blocks::Vector{Vector{Array{Float64,3}}}
end

# conceptual 5-D / 2-D indexing
Base.getindex(sd::SolutionData, it::Integer, ib::Integer) = sd.blocks[it][ib]
Base.getindex(sd::SolutionData, it::Integer, ib::Integer,
              ix::Integer, iy::Integer, iv::Integer) = sd.blocks[it][ib][ix, iy, iv]

function Base.show(io::IO, sd::SolutionData)
    print(io, "SolutionData(\"", sd.title, "\"): ",
          sd.nt, " times, ", sd.nblk, " blocks, ", sd.nvar, " variables\n")
    print(io, "  times      = ", sd.times, "\n")
    print(io, "  variables  = ", sd.varnames, "\n")
    for k in 1:sd.nblk
        print(io, "  block ", sd.blockids[k], ": interior ", sd.nx[k], "x", sd.ny[k],
              " (+ghost ", sd.nxg[k], "x", sd.nyg[k], "), x∈[",
              round(sd.bbox[k,1]; digits=4), ",", round(sd.bbox[k,2]; digits=4),
              "], y∈[", round(sd.bbox[k,3]; digits=4), ",", round(sd.bbox[k,4]; digits=4), "]\n")
    end
end

# ----------------------------------------------------------------------------
# Internal: parse one zone header
# ----------------------------------------------------------------------------
struct Zone
    block::Int
    time::Float64
    strandid::Int
    I::Int            # points in y-direction (fastest index)
    J::Int            # points in x-direction (slowest index)
    nvars::Int        # physical variables per point stored in .dat (no x,y)
    offset::Int       # byte offset of the payload inside the .dat file
    nbytes::Int       # payload byte count
    # block geometry for on-the-fly x,y reconstruction (uniform Cartesian):
    #   x = BX + (i - LFB_I)*DX,  y = BY + (j - LFB_J)*DY
    bx::Float64; by::Float64; dx::Float64; dy::Float64
    lfb_i::Int; lfb_j::Int; istart::Int; jstart::Int
end

# Parse the [VARIABLES] section, the binary PRECISION, and a TITLE-like label
# from the .hdr manifest. The manifest has no TITLE field; we synthesise one
# from the solver/Type lines. PRECISION defaults to Float64 if absent (older
# headers). Returned vnames always start with "x","y" (reconstructed) followed
# by the physical fields stored in .dat.
function _parse_header(lines)
    title  = "Unknown"
    solver = ""; typ = ""
    fvars  = String[]
    prec   = "Float64"
    invars = false
    for ln in lines
        s = strip(ln)
        isempty(s) && continue
        if startswith(s, "[VARIABLES]")
            invars = true; continue
        elseif startswith(s, "[")
            invars = false
        end
        if invars && !startswith(s, "#")
            append!(fvars, split(s))
            invars = false        # single variable line
        end
        m = match(r"^solver\s*=\s*(\S+)",    s); m !== nothing && (solver = m.captures[1])
        m = match(r"^Type\s*=\s*(\S+)",      s); m !== nothing && (typ    = m.captures[1])
        m = match(r"^PRECISION\s*=\s*(\S+)", s); m !== nothing && (prec   = m.captures[1])
    end
    if !isempty(solver)
        title = isempty(typ) ? solver : "$solver ($typ)"
    end
    # x,y are not stored; prepend them so downstream code sees the full 12 columns.
    vnames = ("x" in fvars || "y" in fvars) ? fvars : vcat(["x", "y"], fvars)
    eltyp = prec == "Float32" ? Float32 : Float64
    return title, vnames, eltyp
end

# Parse the per-zone records from the [ZONES] section.
function _parse_zones(lines)
    zones = Zone[]
    geti(line, key) = (m = match(Regex("\\b$key\\s*=\\s*(-?\\d+)"), line);
                       m === nothing ? 0 : parse(Int, m.captures[1]))
    getf(line, key) = (m = match(Regex("\\b$key\\s*=\\s*(-?[0-9.eE+]+)"), line);
                       m === nothing ? 0.0 : parse(Float64, m.captures[1]))
    for line in lines
        startswith(strip(line), "ZONE") || continue
        tm = match(r"SOLUTIONTIME\s*=\s*([0-9Ee+\-.]+)", line)
        push!(zones, Zone(
            geti(line, "block"),
            tm === nothing ? 0.0 : parse(Float64, tm.captures[1]),
            geti(line, "STRANDID"),
            geti(line, "I"),
            geti(line, "J"),
            geti(line, "NVARS"),
            geti(line, "OFFSET"),
            geti(line, "NBYTES"),
            getf(line, "BX"), getf(line, "BY"),
            getf(line, "DX"), getf(line, "DY"),
            geti(line, "LFB_I"), geti(line, "LFB_J"),
            geti(line, "ISTART"), geti(line, "JSTART")))
    end
    return zones
end

# Resolve the .dat companion path from a given path (accepts .hdr, .dat, or base).
function _companion_paths(path::AbstractString)
    base = endswith(path, ".hdr") || endswith(path, ".dat") ?
           path[1:end-4] : path
    return base * ".hdr", base * ".dat"
end

# ----------------------------------------------------------------------------
# load_solution
# ----------------------------------------------------------------------------
"""
    load_solution(path) -> SolutionData

Read a complete binary output result (all blocks, all times) into a
`SolutionData`. `path` may be the `.hdr`, the `.dat`, or the shared basename;
the companion file is located automatically. Metadata comes from the `.hdr`
manifest and the heavy field data is read directly from the binary `.dat`
using the per-zone byte offsets, so no scanning of the payload is needed.
"""
function load_solution(path::AbstractString)
    hdrpath, datpath = _companion_paths(path)
    isfile(hdrpath) || error("header file not found: $hdrpath")
    isfile(datpath) || error("binary file not found: $datpath")

    lines = readlines(hdrpath)
    title, vnames, eltyp = _parse_header(lines)
    nvar = length(vnames)

    zones = _parse_zones(lines)
    isempty(zones) && error("No ZONE records found in $hdrpath")

    # Time instances, grouped by STRANDID, in order of first appearance.
    strand_order = Int[]
    for z in zones
        z.strandid in strand_order || push!(strand_order, z.strandid)
    end
    nt = length(strand_order)
    sidx = Dict(s => k for (k, s) in enumerate(strand_order))

    # Blocks (assumed identical set across times), keyed by file block number.
    blockids = sort!(unique(z.block for z in zones))
    nblk = length(blockids)
    bidx = Dict(b => k for (k, b) in enumerate(blockids))

    zmap = Dict{Tuple{Int,Int},Zone}()
    times_vec = fill(NaN, nt)
    for z in zones
        zmap[(z.strandid, z.block)] = z
        times_vec[sidx[z.strandid]] = z.time
    end

    # Per-block grid sizes from the first time instance. I = y-pts, J = x-pts.
    nxg = zeros(Int, nblk); nyg = zeros(Int, nblk)
    for b in blockids
        z = zmap[(strand_order[1], b)]
        nxg[bidx[b]] = z.J
        nyg[bidx[b]] = z.I
    end
    nx = nxg .- 2
    ny = nyg .- 2

    xc = [fill(NaN, nxg[k]) for k in 1:nblk]
    yc = [fill(NaN, nyg[k]) for k in 1:nblk]
    bbox = fill(NaN, nblk, 4)

    blocks = [[fill(NaN, nxg[k], nyg[k], nvar) for k in 1:nblk] for _ in 1:nt]

    ix_col = findfirst(==("x"), vnames); ix_col === nothing && (ix_col = 1)
    iy_col = findfirst(==("y"), vnames); iy_col === nothing && (iy_col = 2)

    # Read the binary payloads. Each zone stores only the PHYSICAL fields (no
    # x,y), VARIABLE-FASTEST, with points ordered x-slowest (jx) / y-fastest
    # (iy). Reshaping the flat stream to (nvf, ni, nj) gives [v, iy, jx]. We
    # widen to Float64, place the physical fields into columns 3.. and
    # RECONSTRUCT x,y (columns 1,2) from the zone's uniform-grid geometry.
    # Grid indices: along x  i = ISTART + (jx-1),  along y  j = JSTART + (iy-1)
    #   x = BX + (i - LFB_I)*DX,   y = BY + (j - LFB_J)*DY.
    open(datpath, "r") do io
        for s in strand_order, b in blockids
            z  = zmap[(s, b)]
            it = sidx[s]; k = bidx[b]
            ni = z.I; nj = z.J            # ni = y-pts (fast), nj = x-pts (slow)
            nvf = z.nvars == 0 ? nvar - 2 : z.nvars   # physical fields in .dat
            arr = blocks[it][k]
            seek(io, z.offset)
            raw = Vector{eltyp}(undef, nvf * ni * nj)
            read!(io, raw)
            flat = reshape(raw, nvf, ni, nj)          # [v, iy, jx]
            @inbounds for jx in 1:nj, iy in 1:ni
                xij = z.bx + ((z.istart + (jx-1)) - z.lfb_i) * z.dx
                yij = z.by + ((z.jstart + (iy-1)) - z.lfb_j) * z.dy
                arr[jx, iy, ix_col] = xij
                arr[jx, iy, iy_col] = yij
                for v in 1:nvf
                    arr[jx, iy, 2 + v] = Float64(flat[v, iy, jx])
                end
                if it == 1
                    iy == 1 && (xc[k][jx] = xij)
                    jx == 1 && (yc[k][iy] = yij)
                end
            end
        end
    end

    # TRUE geometric bounding box = physical cell-FACE extent of the block.
    # Stored coords are cell centres with one ghost layer each side:
    #   xc[1]=ghost centre, xc[2]=first interior centre,
    #   xc[end-1]=last interior centre, xc[end]=ghost centre.
    # The block physically covers from the left face of the first interior cell
    # to the right face of the last interior cell, i.e. half a cell INSIDE each
    # ghost centre -- exactly the midpoint between the ghost centre and the
    # first interior centre. These true boxes TILE the domain: neighbouring
    # blocks share a face with no overlap and no gap, so any (x,y) lands in a
    # unique block. Interpolation inside that block still uses xc/yc (ghosts
    # included), so a point in the half-cell margin near a face naturally
    # interpolates between the ghost value and the first interior value --
    # correct both at a physical boundary and across an interface halo.
    for k in 1:nblk
        dx = nxg[k] > 1 ? xc[k][2] - xc[k][1] : 0.0
        dy = nyg[k] > 1 ? yc[k][2] - yc[k][1] : 0.0
        bbox[k, 1] = xc[k][2]     - 0.5dx   # left  face of first interior cell
        bbox[k, 2] = xc[k][end-1] + 0.5dx   # right face of last  interior cell
        bbox[k, 3] = yc[k][2]     - 0.5dy
        bbox[k, 4] = yc[k][end-1] + 0.5dy
    end

    varsyms = Tuple(Symbol(v) for v in vnames)
    return SolutionData(title, vnames, varsyms, times_vec, strand_order, blockids,
                        nt, nblk, nvar, nx, ny, nxg, nyg, xc, yc, bbox, blocks)
end

# ----------------------------------------------------------------------------
# Convenience accessors
# ----------------------------------------------------------------------------
variable_index(sd::SolutionData, name::AbstractString) =
    (i = findfirst(==(name), sd.varnames);
     i === nothing ? error("variable \"$name\" not in $(sd.varnames)") : i)
variable_index(::SolutionData, i::Integer) = Int(i)

block_bounds(sd::SolutionData, b::Integer) =
    (xmin=sd.bbox[b,1], xmax=sd.bbox[b,2], ymin=sd.bbox[b,3], ymax=sd.bbox[b,4])
gridsize(sd::SolutionData, b::Integer) =
    (nx=sd.nx[b], ny=sd.ny[b], nxg=sd.nxg[b], nyg=sd.nyg[b])

# ----------------------------------------------------------------------------
# Dense output:  value(sd, x, y, t, variable)
# ----------------------------------------------------------------------------
# 1. TIME : bracket the stored times, linear weight.
# 2. SPACE: find the block whose interior bbox contains (x,y), then bilinear
#           interpolation on that block's cell-centre grid. Ghost-cell centres
#           ARE included: at the physical boundary the value is the average of
#           the ghost centre (at -h/2 outside) and the first interior centre
#           (at +h/2), exactly as required for boundary-accurate output.
# 3. blend the two time slices.
# ----------------------------------------------------------------------------

# 1-D index + weight: returns (i, w) with interpolant (1-w)*f[i] + w*f[i+1].
function _locate(c::AbstractVector{<:Real}, q::Real)
    n = length(c)
    n == 1 && return (1, 0.0)
    q <= c[1]   && return (1, 0.0)
    q >= c[end] && return (n - 1, 1.0)
    i = clamp(searchsortedlast(c, q), 1, n - 1)
    w = (q - c[i]) / (c[i+1] - c[i])
    return (i, w)
end

# Block whose TRUE (face) bounding box contains (x,y).
#
# The true boxes tile the domain: neighbouring blocks share a face with no
# overlap and no gap, so (x,y) lands in a UNIQUE block. A point exactly on a
# shared face is caught by the FIRST matching block in the scan -- harmless,
# since the interpolant is continuous across the face. No tolerance is used:
# adjacent boxes derive the shared face from the same uniform grid, so the
# coordinates coincide and the closed inequalities leave no gap. Keeping the
# boxes un-inflated preserves the exact tiling invariant. No ghost
# extrapolation into neighbours, hence no value seam. Interpolation inside the
# block still uses xc/yc with ghost centres, so the half-cell margin near a
# face interpolates ghost <-> first interior value.
#
# `strict=true`  -> returns 0 only when (x,y) is outside EVERY true box, i.e. a
#                   genuine exterior gap (e.g. the L-shape notch).
# `strict=false` -> additionally falls back to the nearest block centre.
function _find_block(sd::SolutionData, x::Real, y::Real; strict::Bool=false)
    for k in 1:sd.nblk
        if sd.bbox[k,1] <= x <= sd.bbox[k,2] &&
           sd.bbox[k,3] <= y <= sd.bbox[k,4]
            return k                  # first matching block wins
        end
    end
    strict && return 0
    best = 1; bestd = Inf
    for k in 1:sd.nblk
        cx = 0.5*(sd.bbox[k,1]+sd.bbox[k,2]); cy = 0.5*(sd.bbox[k,3]+sd.bbox[k,4])
        d = (x-cx)^2 + (y-cy)^2
        d < bestd && (bestd = d; best = k)
    end
    return best
end

"""
    find_block(sd, x, y; strict=false) -> Int

Index of the block whose interior bounding box contains `(x, y)`. With
`strict=true`, returns `0` when the point lies in a domain gap (outside every
block); otherwise returns the nearest block.
"""
find_block(sd::SolutionData, x::Real, y::Real; strict::Bool=false) =
    _find_block(sd, x, y; strict=strict)

"""
    domain_extent(sd) -> (xmin, xmax, ymin, ymax)

Bounding box of the whole multi-block domain (union of interior block boxes).
"""
function domain_extent(sd::SolutionData)
    xmin = minimum(@view sd.bbox[:,1]); xmax = maximum(@view sd.bbox[:,2])
    ymin = minimum(@view sd.bbox[:,3]); ymax = maximum(@view sd.bbox[:,4])
    return (xmin, xmax, ymin, ymax)
end

# Bilinear interpolation of variable v on block k at time-slice it.
function _interp_space(sd::SolutionData, it::Int, k::Int, x::Real, y::Real, v::Int)
    xc = sd.xc[k]; yc = sd.yc[k]
    (ix, wx) = _locate(xc, x)
    (iy, wy) = _locate(yc, y)
    a = sd.blocks[it][k]
    f00 = a[ix,   iy,   v]; f10 = a[ix+1, iy,   v]
    f01 = a[ix,   iy+1, v]; f11 = a[ix+1, iy+1, v]
    return (1-wx)*(1-wy)*f00 + wx*(1-wy)*f10 + (1-wx)*wy*f01 + wx*wy*f11
end

"""
    value(sd, x, y, t, variable; strict=false) -> Float64

Dense-output evaluation: the value of `variable` (a name like `"rho"` or a
column index) at physical point `(x, y)` and absolute time `t`. The block and
grid cell are located automatically; interpolation is bilinear in space (ghost
cells included) and linear in time. `t` outside the stored range is clamped.

With `strict=true`, points lying in a domain gap (outside every block) return
`NaN` instead of being filled from the nearest block -- useful for masking
non-convex domains (e.g. an L-shape notch) when plotting.
"""
function value(sd::SolutionData, x::Real, y::Real, t::Real, variable; strict::Bool=false)
    v  = variable_index(sd, variable)
    k  = _find_block(sd, x, y; strict=strict)
    k == 0 && return NaN
    (it, wt) = _locate(sd.times, t)
    flo = _interp_space(sd, it, k, x, y, v)
    wt == 0.0 && return flo
    fhi = _interp_space(sd, it+1, k, x, y, v)
    return (1-wt)*flo + wt*fhi
end

"""
    point_values(sd, x, y, t; strict=false) -> NamedTuple

Dense output of ALL variables at `(x, y, t)`, returned as a NamedTuple keyed by
the variable names (e.g. `nt.rho`, `nt.vx`, `nt.Pxy`). The block and time
interval are located once and every column is interpolated the same way as
`value`. Outside the domain (`strict=true`) all entries are `NaN`.

This is the basis for plotting an arbitrary FUNCTION of the variables: pass the
returned NamedTuple to your expression, e.g.
`f = nt -> sqrt(nt.vx^2 + nt.vy^2)`.
"""
function point_values(sd::SolutionData, x::Real, y::Real, t::Real; strict::Bool=false)
    syms = sd.varsyms
    k = _find_block(sd, x, y; strict=strict)
    if k == 0
        return NamedTuple{syms}(ntuple(_ -> NaN, length(syms)))
    end
    (it, wt) = _locate(sd.times, t)
    vals = ntuple(length(syms)) do v
        flo = _interp_space(sd, it, k, x, y, v)
        wt == 0.0 ? flo : (1-wt)*flo + wt*_interp_space(sd, it+1, k, x, y, v)
    end
    return NamedTuple{syms}(vals)
end

end # module
