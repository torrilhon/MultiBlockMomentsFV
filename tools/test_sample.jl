# Headless test of the sampling/masking logic in plot_field.jl
# (CairoMakie itself is not exercised here -- only the backend-agnostic parts.)
include("load_solution.jl"); using .LoadSolution

# Pull in only the sampler functions without triggering the Makie requirement.
include("plot_field.jl")

sd = load_solution(joinpath(@__DIR__, "..", "reference", "square4.ref.dat"))

println("--- domain_extent ---")
ext = domain_extent(sd)
println("  ", ext)
@assert all(isapprox.(collect(ext), [-1.0, 1.0, -1.0, 1.0]; atol=1e-6)) "extent wrong"

println("--- suggest_resolution ---")
nx, ny = suggest_resolution(sd)
println("  suggested (nx,ny) = ($nx,$ny)")
@assert nx >= 100 && ny >= 100

println("--- sample_field (coarse) ---")
xs, ys, Z = sample_field(sd, "rho", 0.2; nx=60, strict=true)
println("  size(Z) = ", size(Z), ", xs[1..3]=", xs[1:3])
@assert size(Z) == (length(xs), length(ys))
@assert all(isfinite, Z) "square4 has no gaps -> all finite expected"
# rho at t=0.2 should be a sensible positive field
@assert minimum(Z) > 0 "density should stay positive"

println("--- masking a domain gap (synthetic point outside all blocks) ---")
# square4 fully tiles [-1,1]x[-1,1]; pick a point far outside the union
vgap = value(sd, 5.0, 5.0, 0.2, "rho"; strict=true)
println("  value at (5,5) strict=true  -> ", vgap, " (expect NaN)")
@assert isnan(vgap)
vfill = value(sd, 5.0, 5.0, 0.2, "rho"; strict=false)
println("  value at (5,5) strict=false -> ", vfill, " (nearest-block fill, finite)")
@assert isfinite(vfill)

println("--- find_block / domain coverage ---")
@assert find_block(sd, -0.5, -0.5; strict=true) == 1
@assert find_block(sd,  0.75, 0.6; strict=true) == 3
@assert find_block(sd,  5.0,  5.0; strict=true) == 0
println("  block lookup OK")

println("\nSAMPLING TESTS PASSED (Makie draw not exercised here)")
