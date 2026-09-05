function read_pos_vel_csv(file)
    lines = readlines(file)
    n = length(lines) - 1
    t, pos_z, vel_z = zeros(n), zeros(n), zeros(n)
    for (i, line) in enumerate(lines[2:end])
        parts = split(line, ',')
        t[i]      = parse(Float64, parts[1])
        pos_z[i]  = parse(Float64, parts[2])
        vel_z[i]  = parse(Float64, parts[3])
    end
    t, pos_z, vel_z
end
