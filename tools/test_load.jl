# Quick test of load_solution.jl against the reference square4 result.
include("load_solution.jl")
using .LoadSolution

f = joinpath(@__DIR__, "..", "reference", "square4.ref.dat")
sd = load_solution(f)
println(sd)

println("\n--- metadata checks ---")
@assert sd.nt == 4            "expected 4 times, got $(sd.nt)"
@assert sd.nblk == 4          "expected 4 blocks, got $(sd.nblk)"
@assert sd.nvar == 12         "expected 12 variables, got $(sd.nvar)"
@assert sd.times == [0.0,0.1,0.2,0.3] "times wrong: $(sd.times)"
@assert sd.strandids == [1,120,240,356] "strandids wrong: $(sd.strandids)"
# block 1: interior 100x120 -> ghost 102x122
@assert (sd.nx[1],sd.ny[1]) == (100,120) "block1 interior wrong"
@assert (sd.nxg[1],sd.nyg[1]) == (102,122) "block1 ghost wrong"
# block 3: interior 50x30
@assert (sd.nx[3],sd.ny[3]) == (50,30) "block3 interior wrong"
println("metadata OK")

println("\n--- interior bounding boxes ---")
for k in 1:sd.nblk
    b = block_bounds(sd, k)
    println("  block $k: x∈[$(round(b.xmin;digits=4)),$(round(b.xmax;digits=4))]  y∈[$(round(b.ymin;digits=4)),$(round(b.ymax;digits=4))]")
end
# block1 should be x∈[-1,0.5], y∈[-1,0.2]
b1 = block_bounds(sd,1)
@assert isapprox(b1.xmin,-1.0;atol=1e-6) && isapprox(b1.xmax,0.5;atol=1e-6)
@assert isapprox(b1.ymin,-1.0;atol=1e-6) && isapprox(b1.ymax,0.2;atol=1e-6)
println("bounding boxes OK")

println("\n--- dense output checks ---")
irho = variable_index(sd,"rho")
# initial condition is uniform rho=1 everywhere (from the file head we saw)
v0 = value(sd, -0.5, -0.5, 0.0, "rho")
println("rho(-0.5,-0.5,t=0) = $v0  (expect 1.0)")
@assert isapprox(v0, 1.0; atol=1e-6)

# a point inside block 3's interior region (x∈[0.5,1], y∈[0.2,1])
v3 = value(sd, 0.75, 0.6, 0.15, "rho")
println("rho(0.75,0.6,t=0.15) = $v3  (time-interpolated between t=0.1 and 0.2)")
@assert isfinite(v3)

# time interpolation: at t halfway, value == average of the two slices at same (x,y)
xq,yq = 0.1, -0.3
fa = value(sd, xq, yq, 0.1, "T")
fb = value(sd, xq, yq, 0.2, "T")
fm = value(sd, xq, yq, 0.15, "T")
println("T at t=0.1: $fa, t=0.2: $fb, t=0.15: $fm ; avg=$(0.5*(fa+fb))")
@assert isapprox(fm, 0.5*(fa+fb); atol=1e-9) "time interpolation not linear"
println("time interpolation OK")

# spatial check: evaluating exactly at an interior cell centre reproduces the
# stored value (bilinear interpolation is exact at nodes).
k = 1
ixc = 10; iyc = 20                  # arbitrary interior indices (incl ghost offset)
xq = sd.xc[k][ixc]; yq = sd.yc[k][iyc]
stored = sd[3, k, ixc, iyc, irho]   # t-index 3 == t=0.2
interp = value(sd, xq, yq, 0.2, "rho")
println("at node ($xq,$yq): stored=$stored  interp=$interp")
@assert isapprox(stored, interp; atol=1e-8) "nodal interpolation mismatch"
println("spatial nodal interpolation OK")

println("\nALL TESTS PASSED ✅")
