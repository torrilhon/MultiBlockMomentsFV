# explore.jl
# ----------------------------------------------------------------------------
# Interactive entry point for VS Code (Julia extension).
#
# How to use in VS Code:
#   1. Open this file. Make sure the Julia env is the `tools/` project:
#        - the status bar should show the env, or run in the Julia REPL:
#          `import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()`
#   2. Run the whole file:  Alt+Enter on the first line, or "Julia: Execute
#      File in REPL", or just step through blocks with Shift+Enter.
#   3. After it runs, `sd` and the loaded modules stay live in the REPL so you
#      can probe interactively (see the examples at the bottom).
#
# Re-running is safe. The source files are loaded with Revise's `includet`
# ("include + track"): Revise updates the existing definitions in place when you
# edit and save load_solution.jl / plot_field.jl -- no module is redefined, so
# there is no "two modules export SolutionData" ambiguity, and your edits take
# effect automatically without a REPL restart. (Plain `include` would redefine
# `module LoadSolution` on every re-run, shadowing the previous instance and
# producing `UndefVarError: SolutionData not defined in Main`.)
#
# Everything is path-independent via @__DIR__, so the working directory of the
# REPL does not matter.
# ----------------------------------------------------------------------------

import Pkg
Pkg.activate(@__DIR__)            # use the tools/ environment (Revise, CairoMakie)
# Pkg.instantiate()               # uncomment on first run if deps are missing

using Revise                      # live-tracks the included source files
using CairoMakie                  # plotting backend (loads the FV-mesh helpers too)

includet(joinpath(@__DIR__, "load_solution.jl"));  using .LoadSolution
includet(joinpath(@__DIR__, "plot_field.jl"))

# --- choose the result file -------------------------------------------------
#datafile = joinpath(@__DIR__, "..", "reference", "square4.ref.dat")
datafile = joinpath(@__DIR__, "..", "runs", "Lshape",
                     "Grad13Riemann2D_nst136000_Kn0p001_tAll.dat")

# --- load the full result ---------------------------------------------------
sd = load_solution(datafile)
println(sd)                       # prints times, blocks, variables, bboxes
println("\nstored times    : ", sd.times)
println("variables       : ", sd.varnames)
println("suggested (nx,ny): ", suggest_resolution(sd))

# --- make a plot ------------------------------------------------------------
fig = plot_field(sd, "qx", 0.5;
                 kind     = :contourf,
                 levels   = 40,
                 colormap = :deep,
                 outline  = false )          # add  mesh = (color=(:black,0.2), linewidth=0.25)  to overlay the FV grid

display(fig)                               # shows inline in the VS Code plot pane
# save(joinpath(@__DIR__, "rho.pdf"), fig) # uncomment to export

# ----------------------------------------------------------------------------
# Probe interactively in the REPL after running the file, e.g.:
#
#   value(sd, 0.3, -0.2, 0.15, "T")          # dense output: T at (x,y,t)
#   value(sd, 0.3, -0.2, 0.15, "rho"; strict=true)
#   find_block(sd, 0.75, 0.6)                # which block contains a point
#   block_bounds(sd, 1)                      # interior bbox of block 1
#   gridsize(sd, 3)                          # (nx,ny,nxg,nyg) of block 3
#   domain_extent(sd)                        # whole-domain bounding box
#   sd[end, 1, 10, 20, variable_index(sd,"rho")]   # raw cell value (incl. ghosts)
#
#   # another field / time / kind:
#   fig2 = plot_field(sd, "T", 0.1; kind=:heatmap); display(fig2)
#
#   # standalone FV mesh:
# figm = plot_mesh(sd, color=(:red,0.9)); display(figm)
# ----------------------------------------------------------------------------

# ---- 1D line-out: fix y and t, sweep x from -1.0 to 1.0 ----
# Assumes `sd` is loaded and `using CairoMakie` is active (as in explore.jl).

# y   = 0.0
# t   = 0.2
# var = "rho"
# Nx  = 200          # number of sample points

# xs   = range(-1.0, 1.0; length = Nx)
# vals = [value(sd, x, y, t, var; strict = true) for x in xs]   # NaN outside domain

# list of (x, val) pairs, if you want it explicitly
# xy = collect(zip(xs, vals))

# piecewise-linear 1D plot (NaNs create natural breaks where there is no data)
# fig = Figure(size = (760, 420))
# ax  = Axis(fig[1, 1]; title  = "$var  at  y = $y,  t = $t",
#            xlabel = "x", ylabel = var)
# lines!(ax, xs, vals)                 # straight segments between consecutive points
# display(fig)
# save("lineout.pdf", fig)