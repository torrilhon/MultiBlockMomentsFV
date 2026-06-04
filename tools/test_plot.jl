# End-to-end render test: actually draws with CairoMakie and saves files.
import Pkg
Pkg.add("CairoMakie")
using CairoMakie

include("load_solution.jl"); using .LoadSolution
include("plot_field.jl")

sd = load_solution(joinpath(@__DIR__, "..", "reference", "square4.ref.dat"))

fig1 = plot_field(sd, "rho", 0.2; kind=:heatmap, nx=300, outline=true)
save(joinpath(@__DIR__, "test_rho_heatmap.png"), fig1)
println("saved test_rho_heatmap.png")

fig2 = plot_field(sd, "T", 0.15; kind=:contourf, levels=30, colormap=:plasma, outline=true)
save(joinpath(@__DIR__, "test_T_contourf.png"), fig2)
println("saved test_T_contourf.png")

println("RENDER OK")
