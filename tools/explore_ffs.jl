# explore_ffs.jl
# ----------------------------------------------------------------------------
# Interactive exploration of the FORWARD-FACING STEP (FFS) results.
#
# This is the FFS-specific sibling of explore.jl. It loads an FFS run and makes
# a polished filled-contour plot of the whole domain, with the "left-out"
# region (the step) shown as a light-red block and the title placed in a small
# white box in the bottom-right corner of that step.
#
# How to use in VS Code (Julia extension):
#   1. Open this file; make sure the Julia env is the `tools/` project
#      (status bar, or run: `import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()`).
#   2. Run the whole file (Alt+Enter on the first line / "Julia: Execute File
#      in REPL"), or step through blocks with Shift+Enter.
#   3. `sd` and the loaded modules stay live, so you can probe interactively.
#
# Re-running is safe (Revise `includet` tracks edits to load_solution.jl /
# plot_field.jl in place). Everything is path-independent via @__DIR__.
# ----------------------------------------------------------------------------

import Pkg
Pkg.activate(@__DIR__)            # use the tools/ environment (Revise, CairoMakie)
# Pkg.instantiate()               # uncomment on first run if deps are missing

using Revise                      # live-tracks the included source files
using CairoMakie                  # plotting backend (loads the FV-mesh helpers too)

includet(joinpath(@__DIR__, "load_solution.jl"));  using .LoadSolution
includet(joinpath(@__DIR__, "plot_field.jl"))

# --- choose the Knudsen number ----------------------------------------------
# Two FFS runs are available out of the box:
#   runs/ffs_Kn0p001/   (Kn = 0.001)
#   runs/ffs_Kn0p005/   (Kn = 0.005)
Kn      = 0.001                                   # 0.001 or 0.005
knstr   = replace(string(Kn), "." => "p")
datafile = joinpath(@__DIR__, "..", "runs", "ffs_Kn$(knstr)",
                    "Grad13ForwardFacingStep_nst6300_Kn$(knstr)_tAll.dat")

# --- load the full result ---------------------------------------------------
sd = load_solution(datafile)
println(sd)
println("\nstored times     : ", sd.times)
println("variables        : ", sd.varnames)
println("suggested (nx,ny) : ", suggest_resolution(sd))
println("left-out region(s): ", gap_rectangles(sd))     # the step

# --- make a polished contour plot -------------------------------------------
var = "rho"
t   = sd.times[end]
fig = plot_field(sd, var, t;
                 kind      = :contourf,
                 levels    = 30,
                 colormap  = :viridis,
                 nx        = 900,                  # fine, seam-free sampling
                 gap_color = (:red, 0.18),         # light-red step block
                 title_box = true,                 # title in a corner box
                 title     = "$var,  Kn = $Kn,  t = $t")

display(fig)                                        # inline in the VS Code plot pane
# save(joinpath(@__DIR__, "ffs_$(var)_Kn$(knstr).pdf"), fig)   # vector export
# save(joinpath(@__DIR__, "ffs_$(var)_Kn$(knstr).png"), fig; px_per_unit = 2)

# ----------------------------------------------------------------------------
# Plot-styling cheatsheet (new keywords on plot_field, all optional):
#
#   gap_color  : fill colour of the left-out region(s). `(:red, 0.18)` = light
#                red; `:transparent` (default) leaves the step blank.
#   title_box  : `true`  -> title in a small white box (thin black frame) in the
#                bottom-right corner of the step (this is the FFS look);
#                `false` -> ordinary axis title on top.
#   The colourbar is automatically matched in HEIGHT to the plotted field, so
#   the legend is exactly as tall as the plot (not the whole figure).
#
# Other ideas you might like:
#   - colormap = :balance / :vik / :turbo   for diverging or high-contrast fields
#   - mesh = (color=(:black,0.15), linewidth=0.25)   to overlay the FV grid
#   - outline = true                                 block bounding boxes
#   - swap the field: "T", "p", "vx", "vy", "qx", "qy", "Pxx", ...
# ----------------------------------------------------------------------------

# Probe interactively in the REPL after running the file, e.g.:
#   value(sd, 1.0, 0.5, t, "T")              # dense output: T at (x,y,t)
#   find_block(sd, 1.0, 0.1)                 # 0 -> point is in the step (gap)
#   domain_extent(sd)                        # whole-domain bounding box
#
#   # another field / Knudsen number:
#   fig2 = plot_field(sd, "T", t; kind=:contourf, levels=30,
#                     gap_color=(:red,0.18), title_box=true,
#                     title="T,  Kn = $Kn,  t = $t"); display(fig2)
