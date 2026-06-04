include(joinpath(@__DIR__, "load_solution.jl"))
using .LoadSolution

sd = load_solution(joinpath(@__DIR__, "..", "reference", "square4.ref.dat"))
println("blocks=", sd.nblk, " times=", sd.times)

# square4 interfaces: vertical at x=0.5 (between B1/B2 and B4/B3),
# horizontal at y=0.2. Sample a dense line across each interface and
# count NaNs with strict=true at t=0.2.
v = variable_index(sd, "rho")
t = 0.2

# vertical interface x=0.5, scan y in [-1,1]
ys = range(-0.999, 0.999; length=500)
nanc_v = count(y -> isnan(value(sd, 0.5, y, t, v; strict=true)), ys)
println("NaNs on x=0.5 line (strict): ", nanc_v, " / ", length(ys))

# horizontal interface y=0.2, scan x in [-1,1]
xs = range(-0.999, 0.999; length=500)
nanc_h = count(x -> isnan(value(sd, x, 0.2, t, v; strict=true)), xs)
println("NaNs on y=0.2 line (strict): ", nanc_h, " / ", length(xs))

# whole grid NaN count (interior region should be fully covered)
xmin,xmax,ymin,ymax = domain_extent(sd)
gx = range(xmin, xmax; length=201); gy = range(ymin, ymax; length=200)
tot = 0; nn = 0
for x in gx, y in gy
    global tot += 1
    isnan(value(sd, x, y, t, v; strict=true)) && (global nn += 1)
end
println("Grid NaNs (strict) over square domain: ", nn, " / ", tot)

# continuity check: value just left/right of interface should be close
yp = 0.0
vl = value(sd, 0.5-1e-4, yp, t, v; strict=true)
vr = value(sd, 0.5+1e-4, yp, t, v; strict=true)
vc = value(sd, 0.5,       yp, t, v; strict=true)
println("rho at x=0.5- : ", vl)
println("rho at x=0.5  : ", vc)
println("rho at x=0.5+ : ", vr)
println("jump |L-R| = ", abs(vl-vr))

# nodal exactness still holds (interior cell centre of block 1)
xc1 = sd.xc[1]; yc1 = sd.yc[1]
ix, iy = 30, 40
xq, yq = xc1[ix], yc1[iy]
exact = sd.blocks[3][1][ix,iy,v]   # it=3 -> t=0.2
interp = value(sd, xq, yq, t, v; strict=true)
println("nodal exact=", exact, " interp=", interp, " diff=", abs(exact-interp))
