# explore_ffs.jl
# ----------------------------------------------------------------------------
# Interactive entry point for the Forward-Facing Step (Mach 3) case.
#
# Mirrors explore.jl but tuned for the FFS geometry: a wide 3:1 domain with a
# rectangular step cut out of the lower-right corner. The plot below is built
# by hand (rather than via plot_field) so we can control three things the
# generic helper does not:
#   * a colourbar that is no taller than the plot axis,
#   * the left-out step region filled in light red,
#   * the title placed inside the step as a small white-boxed label.
#
# Usage in VS Code (Julia extension): same as explore.jl -- Alt+Enter on the
# first line or "Julia: Execute File in REPL". Re-running is safe (Revise).
# ----------------------------------------------------------------------------

import Pkg
Pkg.activate(@__DIR__)            # use the tools/ environment (Revise, CairoMakie)
# Pkg.instantiate()               # uncomment on first run if deps are missing

using Revise
using CairoMakie

includet(joinpath(@__DIR__, "load_solution.jl"));  using .LoadSolution
includet(joinpath(@__DIR__, "plot_field.jl"))

# --- choose the result file -------------------------------------------------
# Two reduced-output runs are provided (t = 0, 0.05, 0.1). Switch Kn here.
Kn       = "0p001"     # "0p001" or "0p005"
datafile = joinpath(@__DIR__, "..", "runs", "ffs-momentum",
                    "Grad13ForwardFacingStep_nst302400_Kn$(Kn)_tAll.dat")

# --- load the full result ---------------------------------------------------
sd = load_solution(datafile)
println(sd)
println("\nstored times    : ", sd.times)
println("variables       : ", sd.varnames)
println("suggested (nx,ny): ", suggest_resolution(sd))

# --- derived quantities (functions of the variables) ------------------------
# Each is called per grid point with a NamedTuple of ALL variables (nt.rho,
# nt.vx, nt.p, nt.Pxx, ...). Pass any of these as `var` below together with a
# matching `label`.
speed    = nt -> sqrt(nt.vx^2 + nt.vy^2)                  # |u|
mach     = nt -> sqrt(nt.vx^2 + nt.vy^2) / sqrt(nt.T)     # Mach (a = sqrt(T))

# von Mises stress from the pressure tensor. The pressure tensor is a genuine
# 3D object: Pzz = 3p - Pxx - Pyy (from p = (Pxx+Pyy+Pzz)/3), and Pxz=Pyz=0.
# The deviator is sigma = P - p*delta; sigma_vM = sqrt(3/2 * sigma:sigma):
vonMises = nt -> sqrt(1.5 * ((nt.Pxx - nt.p)^2 +
                             (nt.Pyy - nt.p)^2 +
                             (2*nt.p - nt.Pxx - nt.Pyy)^2 +   # = (Pzz - p)^2
                             2*nt.Pxy^2))

# --- plot settings ----------------------------------------------------------
# var      = vonMises      # a variable name, OR a function above (then set `label`)
# label    = "von-Mises-stress"    # required when `var` is a function (colourbar/title text)
# t        = 1.0
# kind     = :contourf
# levels   = 40
# colormap = :deep

var      = "rho"      # a variable name, OR a function above (then set `label`)
label    = "density"    # required when `var` is a function (colourbar/title text)
t        = 1.0
kind     = :contourf
levels   = 40
colormap = :viridis

# ----------------------------------------------------------------------------
# plot_ffs: contour plot specialised for the forward-facing step.
#
# The "left-out block" is the step: the part of the domain bounding box that no
# block covers. For the FFS it is the lower-right rectangle. We derive it from
# the geometry rather than hard-coding: x spans from the right edge of the
# lower strip (block 1) to the domain's right edge, and y from the domain
# bottom up to the top of the lower strip.
# ----------------------------------------------------------------------------
function plot_ffs(sd, var, t; label=nothing, kind=:contourf, levels=40,
                  colormap=:deep, nx=600, Kn=nothing, progress=:auto)
    xmin, xmax, ymin, ymax = domain_extent(sd)
    xs, ys, Z = sample_field(sd, var, t; nx=nx, strict=true, progress=progress)

    finite     = filter(isfinite, Z)
    colorrange = isempty(finite) ? (0.0, 1.0) : (minimum(finite), maximum(finite))

    # Step (left-out) rectangle: the union gap in the lower-right corner.
    # block 1 is the lower strip; its top edge is the step height and its right
    # edge is the step's left edge. Derived from geometry, not hard-coded.
    b1  = block_bounds(sd, 1)
    sx0 = 1.0 #b1.xmax            # step starts where the lower strip ends
    sx1 = 3.0 #xmax
    sy0 = 0.0 #ymin
    sy1 = 0.3333 #b1.ymax            # step top = lower-strip top

    aspect = (xmax - xmin) / (ymax - ymin)
    axh    = 300             # axis height in px; width follows the data aspect
    fig    = Figure(size = (round(Int, axh*aspect) + 140, axh + 70))

    ax = Axis(fig[1, 1]; aspect = DataAspect(),
    xgridvisible = false, ygridvisible = false)

    plt = if kind === :contourf
        lev = range(colorrange[1]-1e-6, colorrange[2]+1e-6; length = levels + 1)
        contourf!(ax, xs, ys, Z; levels = lev, colormap = colormap)
    elseif kind === :heatmap
        heatmap!(ax, xs, ys, Z; colormap = colormap, colorrange = colorrange,
                 nan_color = :transparent)
    else
        error("unknown kind=$kind")
    end

    # Fill the left-out step in light red (drawn on top so the blank NaN cells
    # in the step read as a deliberate solid region, not missing data).
    poly!(ax, Point2f[(sx0, sy0), (sx1, sy0), (sx1, sy1), (sx0, sy1)];
          color = (:red, 0.22), strokecolor = (:firebrick, 0.6), strokewidth = 0.8)

    limits!(ax, xmin, xmax, ymin, ymax)

    # Compact colourbar: a slim legend whose height is pinned to the axis. With
    # DataAspect the axis is short and wide; Aspect on the row makes the layout
    # shrink the colourbar cell to the axis height instead of stretching it.
    # Label: a function needs an explicit `label`; a variable name defaults to
    # its own name.
    vname = if label !== nothing
        String(label)
    elseif var isa Function
        error("plot_ffs: `var` is a function -- pass `label=...` for the colourbar/title")
    else
        var isa AbstractString ? var : sd.varnames[variable_index(sd, var)]
    end
    Colorbar(fig[1, 2], plt; width = 12,
             ticklabelsize = 9, labelsize = 11)
    rowsize!(fig.layout, 1, Aspect(1, 1/aspect))
    colgap!(fig.layout, 8)

    # Title moved into the step: a small white box with a thin black frame in
    # the bottom-right corner. The box is sized in data units; the text is
    # anchored to the same corner so it stays inside.
    knstr = Kn === nothing ? "" : ",  Kn = $(Kn)"
    ttl   = "G13 Moment Equations, $(vname),  t = $(round(t; digits=4))$(knstr)"
    pad   = 0.012 * (xmax - xmin)
    bw    = 0.7  * (sx1 - sx0)     # box width  (fraction of step width)
    bh    = 0.40  * (sy1 - sy0)     # box height (fraction of step height)
    bx1   = sx1 - pad;  bx0 = bx1 - bw
    by0   = sy0 + pad;  by1 = by0 + bh
    poly!(ax, Point2f[(bx0, by0), (bx1, by0), (bx1, by1), (bx0, by1)];
          color = :white, strokecolor = :black, strokewidth = 0.8)
    text!(ax, 0.5*(bx0+bx1), 0.5*(by0+by1); text = ttl,
          align = (:center, :center), fontsize = 14, color = :black,
          space = :data)

    return fig, ax, (sx0, sx1, sy0, sy1)
end

fig, ax, step = plot_ffs(sd, var, t; label=label, kind=kind, levels=levels,
                         colormap=colormap, Kn=replace(Kn, "p" => "."))
display(fig)
save(joinpath(@__DIR__, "ffs_$(var)_Kn$(Kn).png"), fig)

# ----------------------------------------------------------------------------
# Probe interactively, e.g.:
#   value(sd, 1.5, 0.5, 0.1, "T")
#   find_block(sd, 0.3, 0.1)
#   point_values(sd, 1.5, 0.5, 0.1)        # NamedTuple of all vars at a point
#
# Plot a derived quantity (function of the variables); `label` is required:
#   fig, _ = plot_ffs(sd, vonMises, t; label="σ_vM"); display(fig)
#   fig, _ = plot_ffs(sd, speed,    t; label="|u|");  display(fig)
#   fig, _ = plot_ffs(sd, mach,     t; label="Ma");   display(fig)
#
# Progress bar for fine grids (auto for large grids; force with progress=):
#   plot_ffs(sd, vonMises, t; label="σ_vM", nx=1500, progress=true)
# ----------------------------------------------------------------------------
