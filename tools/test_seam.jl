include(joinpath(@__DIR__, "load_solution.jl"))
using .LoadSolution

sd = load_solution(joinpath(@__DIR__, "..", "reference", "square4.ref.dat"))
println("blocks=", sd.nblk, " times=", sd.times)
xmin,xmax,ymin,ymax = domain_extent(sd)
println("domain_extent = x[$xmin,$xmax] y[$ymin,$ymax]   (expect exactly [-1,1]^2)")

v = variable_index(sd, "rho"); t = 0.2

# 1. Unique-block coverage on a fine grid + count multi-matches
function nmatch(sd,x,y)
    c=0; tol=1e-9
    for k in 1:sd.nblk
        (x>=sd.bbox[k,1]-tol && x<=sd.bbox[k,2]+tol && y>=sd.bbox[k,3]-tol && y<=sd.bbox[k,4]+tol) && (c+=1)
    end
    c
end
gx=range(xmin,xmax;length=201); gy=range(ymin,ymax;length=200)
nn=0; multi=0; zero=0
for x in gx, y in gy
    isnan(value(sd,x,y,t,v;strict=true)) && (global nn+=1)
    m=nmatch(sd,x,y); m==0 && (global zero+=1); m>=2 && (global multi+=1)
end
println("grid NaNs (strict): $nn / $(length(gx)*length(gy))")
println("points matching 0 blocks: $zero ; matching >=2 blocks (on shared faces): $multi")

# 2. continuity across x=0.5 at y=-0.3
println("\n--- rho across x=0.5 interface (y=-0.3) ---")
for x in range(0.49,0.51;length=9)
    println("x=",round(x;digits=4),"  rho=",value(sd,x,-0.3,t,v;strict=true),"  blk=",find_block(sd,x,-0.3;strict=true))
end

# 3. nodal exactness at an interior centre of block 1
xc1=sd.xc[1]; yc1=sd.yc[1]; ix,iy=30,40
exact=sd.blocks[3][1][ix,iy,v]
interp=value(sd,xc1[ix],yc1[iy],t,v;strict=true)
println("\nnodal exact=$exact interp=$interp diff=",abs(exact-interp))

# 4. physical boundary: just inside left wall x=-1, and BEYOND it -> NaN
println("\n--- physical boundary x=-1 (left wall of block 1) ---")
for x in range(-1.02,-0.98;length=9)
    println("x=",round(x;digits=4),"  rho=",value(sd,x,-0.3,t,v;strict=true))
end
