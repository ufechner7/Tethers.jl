@isdefined(read_pos_vel_csv) || function read_pos_vel_csv(file)
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

# linearly interpolate y(t_src) onto t_query; t_src must be sorted ascending.
# Used to compare series that were saved on different time grids, e.g. because
# an event-triggered solver (Python/Assimulo) inserts extra points that a
# fixed-grid solver (Julia, saveat) does not.
@isdefined(interp_at) || function interp_at(t_src, y_src, t_query)
    y = zeros(length(t_query))
    j = 1
    for (i, tq) in enumerate(t_query)
        while j < length(t_src) - 1 && t_src[j+1] < tq
            j += 1
        end
        t0, t1 = t_src[j], t_src[j+1]
        y0, y1 = y_src[j], y_src[j+1]
        frac = t1 ≈ t0 ? 0.0 : (tq - t0) / (t1 - t0)
        y[i] = y0 + frac * (y1 - y0)
    end
    y
end
