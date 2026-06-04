# load_solution.jl
# ----------------------------------------------------------------------------
# Loader and dense-output for the multi-block Tecplot ASCII (.dat) output of
# MultiBlockMomentsFV.
#
# Milestone 1: load the entire result (all blocks, all time instances) into a
#              single Julia data structure -- no plotting.
# Milestone 2: a "dense-output" function  value(x, y, t, variable)  that locates
#              the relevant block and grid cell for the point (x,y), and
#              interpolates both in time and in space.
#
# ---- Output file format (self-contained) ----------------------------------
#   TITLE     = "..."
#   VARIABLES = "x","y","rho",...            <- column names (x,y first)
#   ZONE T="t = <t>, Block <n>", I=<ni>, J=<nj>, F=POINT, SOLUTIONTIME=<t>, STRANDID=<id>
#   <ni*nj rows of POINT data, one column per variable>
#   ZONE ...
#
#   * Each ZONE is one (time, block) pair.  F=POINT => one row per cell centre.
#   * STRANDID groups all blocks belonging to the *same* output time instance.
#   * Each block carries ONE layer of ghost cells on every side, so for an
#     interior grid of nstx x nsty cells the zone has  I = nsty+2,  J = nstx+2.
#   * POINT ordering: the I-index (y-direction) runs fastest, the J-index
#     (x-direction) runs slowest.  Coordinates are cell *centres*, so the first
#     centre sits at xmin - h/2 (the ghost) and the last at xmax + h/2.
#
# Because the blocks have different grid sizes, the result is stored ragged
# (one array per block) rather than as one padded dense array; the conceptual
# 5-D access  data[time, block, x, y, variable]  is provided via getindex.
# ----------------------------------------------------------------------------

module LoadSolution

export SolutionData, load_solution, value, variable_index,
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
    data_start::Int   # line index of first data row
end

function _parse_header(lines)
    title = "Unknown"
    m = match(r"TITLE\s*=\s*\"([^\"]*)\"", lines[1])
    m !== nothing && (title = m.captures[1])
    vnames = String[mm.captures[1] for mm in eachmatch(r"\"([^\"]+)\"", lines[2])]
    return title, vnames
end

function _parse_zones(lines)
    zones = Zone[]
    for (i, line) in enumerate(lines)
        startswith(line, "ZONE") || continue
        bm = match(r"Block\s+(\d+)", line)
        Im = match(r"\bI\s*=\s*(\d+)", line)
        Jm = match(r"\bJ\s*=\s*(\d+)", line)
        tm = match(r"SOLUTIONTIME\s*=\s*([0-9Ee+\-.]+)", line)
        sm = match(r"STRANDID\s*=\s*(\d+)", line)
        push!(zones, Zone(
            bm === nothing ? 0   : parse(Int, bm.captures[1]),
            tm === nothing ? 0.0 : parse(Float64, tm.captures[1]),
            sm === nothing ? 0   : parse(Int, sm.captures[1]),
            Im === nothing ? 0   : parse(Int, Im.captures[1]),
            Jm === nothing ? 0   : parse(Int, Jm.captures[1]),
            i + 1))
    end
    return zones
end

# ----------------------------------------------------------------------------
# load_solution
# ----------------------------------------------------------------------------
"""
    load_solution(path) -> SolutionData

Read a complete `*.dat` output file (all blocks, all times) into a
`SolutionData`. The file is self-contained, so no input `.txt` is needed.
"""
function load_solution(path::AbstractString)
    lines = readlines(path)
    title, vnames = _parse_header(lines)
    nvar = length(vnames)

    zones = _parse_zones(lines)
    isempty(zones) && error("No ZONE records found in $path")

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

    for s in strand_order, b in blockids
        z  = zmap[(s, b)]
        it = sidx[s]; k = bidx[b]
        ni = z.I; nj = z.J            # ni = y-pts (fast), nj = x-pts (slow)
        arr = blocks[it][k]
        row = z.data_start
        @inbounds for jx in 1:nj
            for iy in 1:ni
                vals = parse.(Float64, split(lines[row]))
                for v in 1:nvar
                    arr[jx, iy, v] = vals[v]
                end
                if it == 1
                    iy == 1 && (xc[k][jx] = vals[ix_col])
                    jx == 1 && (yc[k][iy] = vals[iy_col])
                end
                row += 1
            end
        end
    end

    # interior bounding boxes (drop the half-cell ghost overhang)
    for k in 1:nblk
        dx = nxg[k] > 1 ? xc[k][2] - xc[k][1] : 0.0
        dy = nyg[k] > 1 ? yc[k][2] - yc[k][1] : 0.0
        bbox[k, 1] = xc[k][1]   + 0.5dx
        bbox[k, 2] = xc[k][end] - 0.5dx
        bbox[k, 3] = yc[k][1]   + 0.5dy
        bbox[k, 4] = yc[k][end] - 0.5dy
    end

    return SolutionData(title, vnames, times_vec, strand_order, blockids,
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

# Block whose interior bbox contains (x,y).
# `strict=true`  -> returns 0 if the point lies in no interior bbox (domain gap).
# `strict=false` -> falls back to the nearest block centre.
function _find_block(sd::SolutionData, x::Real, y::Real; strict::Bool=false)
    tol = 1e-9
    for k in 1:sd.nblk
        if x >= sd.bbox[k,1]-tol && x <= sd.bbox[k,2]+tol &&
           y >= sd.bbox[k,3]-tol && y <= sd.bbox[k,4]+tol
            return k
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

end # module
