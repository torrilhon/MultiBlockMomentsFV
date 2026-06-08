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
# Points in domain gaps (e.g. an L-shape notch or a forward-facing step) are
# masked to NaN and can be highlighted with a solid fill (see `gap_color`).
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
    sample_field(sd, quantity, t; nx=400, ny=nothing, strict=true)
        -> (xs, ys, Z)

Sample `quantity` at time `t` onto a uniform grid covering the whole domain.
Returns the grid coordinate vectors `xs` (length nx), `ys` (length ny) and
`Z::Matrix` of size `(nx, ny)`. Domain gaps are `NaN` when `strict=true`.

`quantity` is flexible:
- a variable NAME (`"rho"`) or column index  -> that field is sampled directly;
- a FUNCTION of the variables, called per grid point with a NamedTuple of ALL
  variables at that point, so you can plot any derived quantity by name, e.g.
      sample_field(sd, nt -> sqrt(nt.vx^2 + nt.vy^2), t)     # speed |u|
      sample_field(sd, nt -> nt.p / nt.rho,           t)     # p/rho
      sample_field(sd, nt -> nt.Pxx + nt.Pyy,         t)     # trace of stress
  Outside the domain the NamedTuple is all-`NaN`, so a well-behaved `f`
  propagates `NaN` and the gap stays masked.

`nx` defaults to 400; `ny` is chosen to preserve the domain aspect ratio if not
given. For seam-free, feature-preserving output choose a resolution at least as
fine as the finest block (see `suggest_resolution`).
"""
function sample_field(sd::SolutionData, quantity, t::Real;
                      nx::Int=400, ny::Union{Int,Nothing}=nothing,
                      strict::Bool=true)
    xmin, xmax, ymin, ymax = domain_extent(sd)
    if ny === nothing
        aspect = (ymax - ymin) / (xmax - xmin)
        ny = max(2, round(Int, nx * aspect))
    end
    xs = range(xmin, xmax; length=nx)
    ys = range(ymin, ymax; length=ny)
    Z  = Matrix{Float64}(undef, nx, ny)
    if quantity isa Function
        @inbounds for j in 1:ny, i in 1:nx
            nt = point_values(sd, xs[i], ys[j], t; strict=strict)
            Z[i, j] = Float64(quantity(nt))
        end
    else
        v = variable_index(sd, quantity)
        @inbounds for j in 1:ny, i in 1:nx
            Z[i, j] = value(sd, xs[i], ys[j], t, v; strict=strict)
        end
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
# Left-out (gap) regions
# ----------------------------------------------------------------------------
"""
    gap_rectangles(sd) -> Vector{NamedTuple}

Axis-aligned rectangles of the domain that are NOT covered by any block -- the
"left-out" regions such as the forward-facing step or an L-shape notch.

The blocks tile rectangularly, so the gap is itself a union of rectangles. We
build it exactly from the block bounding boxes: lay down the grid of all block
x- and y-edges, and collect every grid cell whose centre lies inside the domain
bounding box but inside no block. Adjacent empty cells are then merged
row-by-row into maximal rectangles, so a single rectangular step comes back as
one rectangle. Each entry is `(xmin, xmax, ymin, ymax)`.
"""
function gap_rectangles(sd::SolutionData)
    xmin, xmax, ymin, ymax = domain_extent(sd)
    xe = sort!(unique(vcat(sd.bbox[:,1], sd.bbox[:,2], [xmin, xmax])))
    ye = sort!(unique(vcat(sd.bbox[:,3], sd.bbox[:,4], [ymin, ymax])))
    nxe, nye = length(xe)-1, length(ye)-1
    covered(cx, cy) = any(sd.bbox[k,1] <= cx <= sd.bbox[k,2] &&
                          sd.bbox[k,3] <= cy <= sd.bbox[k,4] for k in 1:sd.nblk)
    # empty[i,j] = grid cell (xe[i]..xe[i+1]) x (ye[j]..ye[j+1]) is a gap
    empty = falses(nxe, nye)
    for j in 1:nye, i in 1:nxe
        cx = 0.5*(xe[i]+xe[i+1]); cy = 0.5*(ye[j]+ye[j+1])
        empty[i, j] = !covered(cx, cy)
    end
    # greedy maximal-rectangle merge
    rects = NamedTuple[]
    used = falses(nxe, nye)
    for j in 1:nye, i in 1:nxe
        (empty[i,j] && !used[i,j]) || continue
        i2 = i
        while i2 < nxe && empty[i2+1, j] && !used[i2+1, j]; i2 += 1; end
        j2 = j
        while j2 < nye && all(@view empty[i:i2, j2+1]) && !any(@view used[i:i2, j2+1])
            j2 += 1
        end
        used[i:i2, j:j2] .= true
        push!(rects, (xmin=xe[i], xmax=xe[i2+1], ymin=ye[j], ymax=ye[j2+1]))
    end
    return rects
end

# ----------------------------------------------------------------------------
# CairoMakie plotting
# ----------------------------------------------------------------------------
"""
    plot_field(sd, quantity, t; label=nothing, kind=:heatmap, nx=400, ny=nothing,
               levels=20, colormap=:viridis, outline=false, strict=true,
               colorrange=nothing, title=nothing,
               gap_color=:transparent, title_box=false) -> Figure

Plot `quantity` over the whole domain at time `t` on a single consistent grid.

`quantity` is either a variable NAME/index (e.g. `"rho"`) or a FUNCTION of the
variables evaluated per grid point with a NamedTuple of all variables, e.g.
    plot_field(sd, nt -> sqrt(nt.vx^2 + nt.vy^2), t; label="|u|")
    plot_field(sd, nt -> nt.p / nt.rho,           t; label="p/rho")
When `quantity` is a function you MUST pass `label` (used for the colourbar and
the default title); for a plain variable the name is used automatically.

Keyword arguments
- `label`      : colourbar / title label; REQUIRED when `quantity` is a function
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
- `gap_color`  : fill colour for left-out regions (forward-facing step / notch).
                 Default `:transparent`; e.g. `(:red, 0.18)` for a light-red step.
- `title_box`  : if `true`, draw the title in a small white box (thin black
                 frame) anchored to the bottom-right corner of the largest
                 left-out region, instead of as an axis title on top. Has no
                 effect if there is no gap.

The colourbar is matched in height to the data-aspect plot, so the legend is
exactly as tall as the field -- not the full figure height.

Requires `using CairoMakie` (or GLMakie) to be loaded by the caller. Returns the
Makie `Figure`; save with `save("name.pdf", fig)`.
"""
function plot_field(sd::SolutionData, quantity, t::Real;
                    label=nothing,
                    kind::Symbol=:heatmap,
                    nx::Int=400, ny::Union{Int,Nothing}=nothing,
                    levels::Int=40, colormap=:deep,
                    outline::Bool=false, mesh=false, strict::Bool=true,
                    colorrange=nothing, title=nothing,
                    gap_color=:transparent, title_box::Bool=false)
    Makie = _require_makie()

    xs, ys, Z = sample_field(sd, quantity, t; nx=nx, ny=ny, strict=strict)

    # Resolve the label: a function needs an explicit `label`; a variable
    # name/index defaults to its name.
    vname = if label !== nothing
        String(label)
    elseif quantity isa Function
        error("plot_field: `quantity` is a function -- pass `label=...` for the colourbar/title")
    else
        quantity isa AbstractString ? quantity : sd.varnames[variable_index(sd, quantity)]
    end
    ttl   = title === nothing ? "$(vname) at t = $(round(t; digits=4))" : title

    if colorrange === nothing
        finite = filter(isfinite, Z)
        colorrange = isempty(finite) ? (0.0, 1.0) : (minimum(finite), maximum(finite))
    end

    xmin, xmax, ymin, ymax = domain_extent(sd)
    rects = gap_rectangles(sd)
    use_box = title_box && !isempty(rects)

    # Figure height drives the axis height; width follows from the data aspect
    # plus room for the colourbar. Picking the figure proportions from the
    # domain aspect makes the DataAspect axis fill the row vertically, so the
    # colourbar (same row) ends up the same height as the plot.
    domar = (ymax - ymin) / (xmax - xmin)
    figh  = 460
    axw   = figh / max(domar, 1e-6)
    figw  = round(Int, axw + 150)            # + colourbar & labels
    fig = Makie.Figure(size = (figw, figh))
    ax  = Makie.Axis(fig[1, 1];
                     title = use_box ? "" : ttl,
                     xlabel = "x", ylabel = "y",
                     aspect = Makie.DataAspect(),
                     xgridvisible = false, ygridvisible = false)

    plt = if kind === :contourf
        # contourf! takes no `colorrange`; the shared scale is set by an
        # explicit level vector spanning `colorrange`.
        lev = range(colorrange[1]-1e-6, colorrange[2]+1e-6; length = levels + 1)
        Makie.contourf!(ax, xs, ys, Z; levels = lev, colormap = colormap )
    elseif kind === :heatmap
        Makie.heatmap!(ax, xs, ys, Z; colormap = colormap, colorrange = colorrange,
                       nan_color = :transparent)
    else
        error("unknown kind=$kind (use :heatmap or :contourf)")
    end

    # Solid fill for the left-out region(s) -- drawn ON TOP of the field so it
    # cleanly masks any contour bleed at the gap edge (contourf has no notion
    # of NaN holes and paints a thin sliver into the step).
    if gap_color !== :transparent
        for r in rects
            Makie.poly!(ax, Makie.Point2f[(r.xmin, r.ymin), (r.xmax, r.ymin),
                                          (r.xmax, r.ymax), (r.xmin, r.ymax)];
                        color = gap_color, strokecolor = (:black, 0.35),
                        strokewidth = 0.6)
        end
    end

    # Slim colourbar matched to the plot height (Relative(1) of the axis row),
    # with compact tick/label fonts so it reads as a legend, not a second panel.
    Makie.Colorbar(fig[1, 2], plt; label = vname, width = 12,
                   ticklabelsize = 9, labelsize = 11,
                   height = Makie.Relative(1.0), tellheight = false)
    Makie.colgap!(fig.layout, 8)
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

    # Title in a small white box (thin black frame), floating INSIDE the largest
    # left-out region. The box is sized EXACTLY to the laid-out text by Makie
    # (no hand-tuned glyph widths, so any title length stays inside its frame),
    # then centred within the step's band and inset from the right edge.
    if use_box
        r = argmax(g -> (g.xmax-g.xmin)*(g.ymax-g.ymin), rects)  # largest gap
        _step_title!(fig, ax, Makie, ttl, r, (xmin, xmax, ymin, ymax), figh)
    end

    return fig
end

# Place `str` in an auto-sized white box (thin black frame) floating inside the
# gap rectangle `r`. Implemented as a nested GridLayout overlaid on the axis
# cell so Makie sizes the box to the text automatically. The box is aligned to
# the right and vertically centred within the step band, with pixel paddings
# derived from the step's fractional position so it floats comfortably inside
# the step rather than hugging the axis frame. Works for any gap that touches
# the bottom/right edges (FFS step) and degrades gracefully otherwise.
function _step_title!(fig, ax, Makie, str, r, ext, figh; fontsize = 14)
    xmin, xmax, ymin, ymax = ext
    fy = (r.ymax - r.ymin) / (ymax - ymin)          # step height fraction
    # axis pixel height ~ figh minus top/bottom protrusions; enough for layout.
    axpix   = 0.78 * figh
    bandpix = fy * axpix                            # step band height in px
    # centre the (auto-height) box in the band: bottom pad ~ half the leftover.
    botpad  = max(6.0, 0.30 * bandpix)
    gl = Makie.GridLayout(fig[1, 1]; tellwidth = false, tellheight = false,
                          halign = :right, valign = :bottom,
                          alignmode = Makie.Outside(12, 16, botpad, 12))
    Makie.Box(gl[1, 1]; color = :white, strokecolor = :black,
              strokewidth = 0.8, cornerradius = 0)
    Makie.Label(gl[1, 1], str; fontsize = fontsize, color = :black,
                padding = (10, 10, 6, 6))
    return gl
end

# ----------------------------------------------------------------------------
# FV mesh (cell FACES, not centre-to-centre connections)
# ----------------------------------------------------------------------------
"""
    face_coords(sd, k) -> (xf, yf)

Cell-FACE coordinates of block `k`, ghost cells excluded. Returns `xf`
(length `nx+1`) and `yf` (length `ny+1`): the grid lines that bound the
interior finite-volume cells. The stored coordinates are cell centres, so the
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
