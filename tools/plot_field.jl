# plot_field.jl
# ----------------------------------------------------------------------------
# Whole-domain plotting for MultiBlockMomentsFV results, built on the
# dense-output `value(x,y,t,variable)` from load_solution.jl.
#
# Rationale: plotting blocks individually leads to inconsistent colour scales,
# unsynchronised contour levels, and visible seams at block interfaces (cell
# centres of neighbouring blocks are half a cell apart and share no nodes).
# Instead we sample the dense output onto ONE uniform grid spanning the entire
# domain, so a single heatmap/contour call gives a consistent, seam-free image.
# Points in domain gaps (e.g. an L-shape notch) are masked to NaN and render
# blank.
#
# Backend: CairoMakie (pure-Julia, publication-quality vector output).
#   import Pkg; Pkg.add("CairoMakie")
#
# Usage:
#   include("load_solution.jl"); using .LoadSolution
#   include("plot_field.jl")
#   sd  = load_solution("../reference/square4.ref.dat")
#   fig = plot_field(sd, "rho", 0.2)                 # heatmap at t=0.2
#   save("rho.pdf", fig)
#   plot_field(sd, "T", 0.15; kind=:contourf, levels=30, outline=true)
# ----------------------------------------------------------------------------

using .LoadSolution

# ----------------------------------------------------------------------------
# Backend-agnostic sampler
# ----------------------------------------------------------------------------
"""
    sample_field(sd, variable, t; nx=400, ny=nothing, strict=true)
        -> (xs, ys, Z)

Sample the dense output of `variable` at time `t` onto a uniform grid covering
the whole domain. Returns the grid coordinate vectors `xs` (length nx), `ys`
(length ny) and `Z::Matrix` of size `(nx, ny)` with `Z[i,j] = value(xs[i],ys[j],t)`.
Domain gaps are `NaN` when `strict=true`.

`nx` defaults to 400; `ny` is chosen to preserve the domain aspect ratio if not
given. For seam-free, feature-preserving output choose a resolution at least as
fine as the finest block (see `suggest_resolution`).
"""
function sample_field(sd::SolutionData, variable, t::Real;
                      nx::Int=400, ny::Union{Int,Nothing}=nothing,
                      strict::Bool=true)
    xmin, xmax, ymin, ymax = domain_extent(sd)
    if ny === nothing
        aspect = (ymax - ymin) / (xmax - xmin)
        ny = max(2, round(Int, nx * aspect))
    end
    xs = range(xmin, xmax; length=nx)
    ys = range(ymin, ymax; length=ny)
    v  = variable_index(sd, variable)
    Z  = Matrix{Float64}(undef, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        Z[i, j] = value(sd, xs[i], ys[j], t, v; strict=strict)
    end
    return collect(xs), collect(ys), Z
end

"""
    suggest_resolution(sd) -> (nx, ny)

Resolution that matches the finest block, so the uniform sampling grid does not
under-resolve any block. Use as a starting point for `sample_field`.
"""
function suggest_resolution(sd::SolutionData)
    xmin, xmax, ymin, ymax = domain_extent(sd)
    # finest cell size over all blocks
    hx = minimum((sd.bbox[k,2]-sd.bbox[k,1]) / max(sd.nx[k],1) for k in 1:sd.nblk)
    hy = minimum((sd.bbox[k,4]-sd.bbox[k,3]) / max(sd.ny[k],1) for k in 1:sd.nblk)
    nx = max(2, ceil(Int, (xmax-xmin)/hx))
    ny = max(2, ceil(Int, (ymax-ymin)/hy))
    return (nx, ny)
end

# ----------------------------------------------------------------------------
# CairoMakie plotting
# ----------------------------------------------------------------------------
"""
    plot_field(sd, variable, t; kind=:heatmap, nx=400, ny=nothing,
               levels=20, colormap=:viridis, outline=false, strict=true,
               colorrange=nothing, title=nothing) -> Figure

Plot `variable` over the whole domain at time `t` on a single consistent grid.

Keyword arguments
- `kind`       : `:heatmap` (pixel-smooth field) or `:contourf` (filled contours)
- `nx`,`ny`    : sampling resolution (ny auto from aspect if omitted)
- `levels`     : number of filled-contour levels (only for `:contourf`)
- `colormap`   : any Makie colormap symbol
- `colorrange` : `(lo, hi)` to fix the colour scale; default = data min/max
                 (NaNs ignored) -- ONE shared scale for the whole domain
- `outline`    : draw interior block bounding boxes on top
- `mesh`       : overlay the FV mesh (cell faces, ghosts excluded). Pass `true`
                 or a NamedTuple of mesh options, e.g.
                 `mesh = (color=(:black,0.25), linewidth=0.3)`
- `strict`     : mask domain gaps to NaN (blank)
- `title`      : plot title (default: "<variable> at t = <t>")
- `rasterize`  : rasterize the field layer in vector (PDF/SVG) output to avoid
                 hairline anti-aliasing seams between cells/bands. `false` keeps
                 pure vector; `true` = 1 px/pt; an integer = density scale
                 (default 2). Does not affect PNG output.

Requires `using CairoMakie` (or GLMakie) to be loaded by the caller. Returns the
Makie `Figure`; save with `save("name.pdf", fig)`.
"""
function plot_field(sd::SolutionData, variable, t::Real;
                    kind::Symbol=:heatmap,
                    nx::Int=400, ny::Union{Int,Nothing}=nothing,
                    levels::Int=20, colormap=:viridis,
                    outline::Bool=false, mesh=false, strict::Bool=true,
                    colorrange=nothing, title=nothing, rasterize=2)
    Makie = _require_makie()

    xs, ys, Z = sample_field(sd, variable, t; nx=nx, ny=ny, strict=strict)

    vname = variable isa AbstractString ? variable : sd.varnames[variable_index(sd,variable)]
    ttl   = title === nothing ? "$(vname) at t = $(round(t; digits=4))" : title

    if colorrange === nothing
        finite = filter(isfinite, Z)
        colorrange = isempty(finite) ? (0.0, 1.0) : (minimum(finite), maximum(finite))
    end

    xmin, xmax, ymin, ymax = domain_extent(sd)
    fig = Makie.Figure(size = (760, 640))
    ax  = Makie.Axis(fig[1, 1]; title = ttl, xlabel = "x", ylabel = "y",
                     aspect = Makie.DataAspect())

    # `rasterize` renders the field as a single raster layer in vector output
    # (PDF/SVG). This suppresses the hairline white lines that CairoMakie's
    # per-polygon anti-aliasing otherwise draws between abutting heatmap cells
    # or contourf bands (each polygon anti-aliases against the white page,
    # lightening shared edges). The data itself is seam-free (blocks tile the
    # domain), so any remaining thin lines are purely this render artifact.
    # `rasterize` may be `false` (keep pure vector), `true` (1 px/pt), or an
    # integer scale (default 2 -> 2x density, crisp without seams).
    plt = if kind === :contourf
        # contourf! takes no `colorrange`; the shared scale is set by an
        # explicit level vector spanning `colorrange`.
        lev = range(colorrange[1], colorrange[2]; length = levels + 1)
        Makie.contourf!(ax, xs, ys, Z; levels = lev, colormap = colormap,
                        rasterize = rasterize)
    elseif kind === :heatmap
        Makie.heatmap!(ax, xs, ys, Z; colormap = colormap, colorrange = colorrange,
                       nan_color = :transparent, rasterize = rasterize)
    else
        error("unknown kind=$kind (use :heatmap or :contourf)")
    end

    Makie.Colorbar(fig[1, 2], plt; label = vname)
    Makie.limits!(ax, xmin, xmax, ymin, ymax)

    if mesh !== false
        mopts = mesh === true ? NamedTuple() : mesh
        draw_mesh!(ax, sd; mopts...)
    end

    if outline
        for k in 1:sd.nblk
            b = block_bounds(sd, k)
            rx = [b.xmin, b.xmax, b.xmax, b.xmin, b.xmin]
            ry = [b.ymin, b.ymin, b.ymax, b.ymax, b.ymin]
            Makie.lines!(ax, rx, ry; color = (:white, 0.6), linewidth = 0.8)
        end
    end

    return fig
end

# ----------------------------------------------------------------------------
# FV mesh (cell FACES, not centre-to-centre connections)
# ----------------------------------------------------------------------------
"""
    face_coords(sd, k) -> (xf, yf)

Cell-FACE coordinates of block `k`, ghost cells excluded. Returns `xf`
(length `nx+1`) and `yf` (length `ny+1`): the grid lines that bound the
interior finite-volume cells. The stored coordinates are cell CENTRES, so the
faces are reconstructed as the uniform interior grid from the block's interior
bounding box (`xmin..xmax` with `nx+1` lines, similarly in y).
"""
function face_coords(sd::SolutionData, k::Integer)
    b  = block_bounds(sd, k)
    xf = range(b.xmin, b.xmax; length = sd.nx[k] + 1)
    yf = range(b.ymin, b.ymax; length = sd.ny[k] + 1)
    return collect(xf), collect(yf)
end

"""
    draw_mesh!(ax, sd; blocks=1:sd.nblk, color=(:black,0.25), linewidth=0.3)

Draw the finite-volume mesh of the given `blocks` onto Makie axis `ax`. Lines
are the cell FACES (ghosts excluded), built from `face_coords`, so this is the
actual FV grid -- not a connection of cell-centre values. Block faces are drawn
as batched line segments (NaN-separated) for speed.
"""
function draw_mesh!(ax, sd::SolutionData; blocks = 1:sd.nblk,
                    color = (:black, 0.25), linewidth = 0.3)
    Makie = _require_makie()
    xpts = Float64[]; ypts = Float64[]
    for k in blocks
        xf, yf = face_coords(sd, k)
        y0, y1 = yf[1], yf[end]
        x0, x1 = xf[1], xf[end]
        for x in xf                      # vertical face lines
            push!(xpts, x, x, NaN); push!(ypts, y0, y1, NaN)
        end
        for y in yf                      # horizontal face lines
            push!(xpts, x0, x1, NaN); push!(ypts, y, y, NaN)
        end
    end
    Makie.lines!(ax, xpts, ypts; color = color, linewidth = linewidth)
    return ax
end

"""
    plot_mesh(sd; blocks=1:sd.nblk, color=:black, linewidth=0.4, title="FV mesh") -> Figure

Standalone plot of just the multi-block FV mesh (cell faces, ghosts excluded).
Useful for checking the grid layout. Returns a Makie `Figure`.
"""
function plot_mesh(sd::SolutionData; blocks = 1:sd.nblk,
                   color = :black, linewidth = 0.4, title = "FV mesh")
    Makie = _require_makie()
    xmin, xmax, ymin, ymax = domain_extent(sd)
    fig = Makie.Figure(size = (700, 640))
    ax  = Makie.Axis(fig[1, 1]; title = title, xlabel = "x", ylabel = "y",
                     aspect = Makie.DataAspect())
    draw_mesh!(ax, sd; blocks = blocks, color = color, linewidth = linewidth)
    Makie.limits!(ax, xmin, xmax, ymin, ymax)
    return fig
end

# Lazily fetch a loaded Makie-compatible module (CairoMakie or GLMakie) from Main.
function _require_makie()
    for sym in (:CairoMakie, :GLMakie, :WGLMakie, :Makie)
        if isdefined(Main, sym)
            return getfield(Main, sym)
        end
    end
    error("No Makie backend loaded. Run `using CairoMakie` (or GLMakie) first.")
end
