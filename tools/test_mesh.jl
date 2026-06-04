# Render test for FV mesh plotting (faces, ghosts excluded).
using CairoMakie
include("load_solution.jl"); using .LoadSolution
include("plot_field.jl")

sd = load_solution(joinpath(@__DIR__, "..", "reference", "square4.ref.dat"))

# --- numeric checks on face reconstruction ---
xf, yf = face_coords(sd, 1)
@assert length(xf) == sd.nx[1] + 1 "block1 should have nx+1 vertical faces"
@assert length(yf) == sd.ny[1] + 1
b = block_bounds(sd, 1)
@assert isapprox(xf[1], b.xmin; atol=1e-9) && isapprox(xf[end], b.xmax; atol=1e-9)
@assert isapprox(yf[1], b.ymin; atol=1e-9) && isapprox(yf[end], b.ymax; atol=1e-9)
# faces sit half a cell beyond the first/last interior cell CENTRE
h = sd.xc[1][2] - sd.xc[1][1]
first_interior_centre = sd.xc[1][2]               # index 1 is ghost
@assert isapprox(xf[1], first_interior_centre - h/2; atol=1e-9) "face != centre - h/2"
println("face reconstruction OK  (xf: $(length(xf)) lines, yf: $(length(yf)) lines)")

# --- standalone mesh ---
fig0 = plot_mesh(sd; linewidth=0.3)
save(joinpath(@__DIR__, "test_mesh.png"), fig0); println("saved test_mesh.png")

# --- contour + mesh overlay ---
fig1 = plot_field(sd, "rho", 0.2; kind=:contourf, levels=24,
                  mesh = (color=(:black,0.25), linewidth=0.25), outline=true)
save(joinpath(@__DIR__, "test_contour_mesh.png"), fig1); println("saved test_contour_mesh.png")

println("MESH RENDER OK")
