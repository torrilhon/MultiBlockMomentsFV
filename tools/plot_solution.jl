using PyPlot

function parse_solution_file(path)
    lines = readlines(path)

    title_match = match(r"TITLE\s*=\s*\"([^\"]*)\"", lines[1])
    dataset_title = title_match !== nothing ? title_match.captures[1] : "Unknown Title"

    variable_matches = collect(eachmatch(r"\"([^\"]+)\"", lines[2]))
    varnames = [m.captures[1] for m in variable_matches]
    var_dict = Dict(name => idx for (idx, name) in enumerate(varnames))

    zones = Vector{NamedTuple{(:block, :time, :I, :J, :data_start, :zone_title), Tuple{Int, Float64, Int, Int, Int, String}}}()
    for (i, line) in enumerate(lines)
        if startswith(line, "ZONE")
            title_match = match(r"T=\"([^\"]*)\"", line)
            zone_title = title_match !== nothing ? title_match.captures[1] : "Unknown Zone"

            block_match = match(r"Block\s+(\d+)", zone_title)
            block_num = block_match !== nothing ? parse(Int, block_match.captures[1]) : 0

            i_match = match(r"I=(\d+)", line)
            I = i_match !== nothing ? parse(Int, i_match.captures[1]) : 0

            j_match = match(r"J=(\d+)", line)
            J = j_match !== nothing ? parse(Int, j_match.captures[1]) : 0

            time_match = match(r"SOLUTIONTIME=([0-9Ee+\-.]+)", line)
            solution_time = time_match !== nothing ? parse(Float64, time_match.captures[1]) : 0.0

            data_start = i + 1
            push!(zones, (block=block_num, time=solution_time, I=I, J=J, data_start=data_start, zone_title=zone_title))
        end
    end

    return dataset_title, var_dict, zones, lines
end

function load_zone_data(zones, var_dict, field, time, lines)
    zones_at_time = sort([z for z in zones if z.time == time], by = z->z.block)

    field_idx = var_dict[field]

    x_all = Float64[]
    y_all = Float64[]
    field_all = Float64[]

    for zone_info in zones_at_time
        nrows = zone_info.I * zone_info.J
        data_lines = lines[zone_info.data_start : zone_info.data_start + nrows - 1]

        for line in data_lines
            values = parse.(Float64, split(strip(line)))
            push!(x_all, values[1])
            push!(y_all, values[2])
            push!(field_all, values[field_idx])
        end
    end

    return x_all, y_all, field_all
end

### define your file, field, time, and output path here

file = "Grad13Riemann2D_nst40000_Kn0p001_tAll.dat"
field = "rho"
time = 0.1
out = "" # eg "example.png"

### main code

dataset_title, var_dict, zones, lines = parse_solution_file(file)
x_all, y_all, field_all = load_zone_data(zones, var_dict, field, time, lines)

levels = range(minimum(field_all), stop=maximum(field_all), length=51)

figure(figsize=(12, 8))
contour = tricontourf(x_all, y_all, field_all, levels=levels, cmap="viridis")
xlabel("x")
ylabel("y")
title("$dataset_title: Contours of $field(x,y) at t=$(round(time, digits=3))")
axis("equal")
colorbar(contour, label=field)
tight_layout()

if !isempty(out)
    savefig(out, dpi=300, bbox_inches="tight")
    println("Saved contour image to $out")
end

matplotlib.pyplot.close()
display(gcf())
