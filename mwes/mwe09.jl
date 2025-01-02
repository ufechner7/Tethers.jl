using Plots, LinearAlgebra, Roots, FiniteDiff, Quaternions

# rotate a 3d vector around the y axis
function rotate_in_xz(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = cos(angle) * vec[3] + sin(angle) * vec[1]
    result
end

# rotate a 3d vector around the z axis
function rotate_in_yx(vec, angle)
    result = similar(vec)
    result[1] = cos(angle) * vec[1] + sin(angle) * vec[2]
    result[2] = cos(angle) * vec[2] - sin(angle) * vec[1]
    result[3] = vec[3]
    result
end

"""
Only on x and z axis, tether start and end point laying on x axis
"""
function generate_tether!(points, d, segments, tether_length, total_angle)
    segment_angle = total_angle / segments
    points[:, 1] .= 0
    initial_angle = -total_angle / 2 + segment_angle/2
    for i in 2:segments+1
        angle = initial_angle + (i - 2) * segment_angle
        points[1, i] = points[1, i-1] + tether_length/segments * cos(angle)
        points[3, i] = points[3, i-1] + tether_length/segments * sin(angle)
    end
    dx = points[1, end] - points[1, 1]
    dy = points[3, end] - points[3, 1]
    d[] = sqrt(dx^2 + dy^2)
    nothing
end

"""
use z axis
"""
function rotate!(points, azimuth, elevation)
    for i in eachindex(points[1, :])
        points[:, i] .= rotate_in_xz(points[:, i], elevation)
        points[:, i] .= rotate_in_yx(points[:, i], azimuth)
    end
    nothing
end

# Example usage
segments = 10 # number of segments
azimuth = 0
elevation = 0.5π * 0.9
distance = 10
start = zeros(3)
tether_length = 11 # tether length
d = Ref(0.0)
cost = Ref(0.0)
@assert tether_length > distance
tether_length = max(distance + 1e-3, tether_length)

function tether_from_distance_length!(points, distance, tether_length, segments, azimuth, elevation)
    d[] = 0.0
    cost[] = 0.0
    function f_zero!(total_angle)
        generate_tether!(points, d, segments, tether_length, total_angle)
        cost[] = d[] - distance
        return cost[]
    end
    find_zero(f_zero!, (0, 2π))
    rotate!(points, azimuth, elevation)
    return nothing
end

points = zeros(3, segments+1)
function f_jac(dx, x)
    points .= 0.0
    tether_from_distance_length!(points, x[1], x[2], segments, azimuth, elevation)
    dx .= reshape(points, 3(segments+1))
    nothing
end

x = [distance, tether_length]
J = zeros(3(segments+1), 2)
FiniteDiff.finite_difference_jacobian!(J, f_jac, x)
@show norm(J)

x_coords = points[1, :]
z_coords = points[3, :]
p = plot(x_coords, z_coords, marker=:circle, label="Tether", aspect_ratio=:equal)
scatter!([start[1], points[1, end]], [start[3], points[3, end]], label="Endpoints", color=:red)

scale_factor = 2.0
for i in 1:segments+1
    # Index into Jacobian for x and y components
    idx_x = 3i - 2
    idx_y = 3i
    
    # Velocity due to distance change (column 1)
    vx_d = J[idx_x, 1]
    vy_d = J[idx_y, 1]
    
    # Velocity due to length change (column 2)
    vx_l = J[idx_x, 2]
    vy_l = J[idx_y, 2]
    
    # Plot arrows for distance sensitivity (red)
    quiver!([x_coords[i]], [z_coords[i]], 
            quiver=([vx_d * scale_factor], [vy_d * scale_factor]), 
            color=:red, alpha=1.5, label=i==1 ? "Distance sensitivity" : "")
    
    # Plot arrows for length sensitivity (blue)
    quiver!([x_coords[i]], [z_coords[i]], 
            quiver=([vx_l * scale_factor], [vy_l * scale_factor]), 
            color=:blue, alpha=1.5, label=i==1 ? "Length sensitivity" : "")
end

xlabel!("x")
ylabel!("z")
title!("Tether with Sensitivity Vectors")

# Display the plot
display(p)
# end

# mwe()

