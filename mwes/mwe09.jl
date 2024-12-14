using Plots, LinearAlgebra, Roots, FiniteDiff


# function mwe()
function generate_tether!(points, d, n, l, start, total_angle)
    segment_angle = total_angle / n
    points[:, 1] .= start
    initial_angle = -total_angle / 2 + segment_angle/2 + π/2
    for i in 2:n+1
        angle = initial_angle + (i - 2) * segment_angle
        points[1, i] = points[1, i-1] + l/n * cos(angle)
        points[2, i] = points[2, i-1] + l/n * sin(angle)
    end
    dx = points[1, end] - points[1, 1]
    dy = points[2, end] - points[2, 1]
    d[] = sqrt(dx^2 + dy^2)
    return nothing
end

# Example usage
n = 10 # number of segments
l = 9 # segment length
start = [0.0, 0.0] # start point coordinates
desired_distance = 8.0 # desired distance
d = Ref(0.0)
cost = Ref(0.0)

function tether_from_distance_length!(points, distance, l, n, start)
    d[] = 0.0
    cost[] = 0.0
    function f_zero!(total_angle)
        generate_tether!(points, d, n, l, start, total_angle)
        cost[] = d[] - distance
        return cost[]
    end
    find_zero(f_zero!, (0, 2π))
    return nothing
end

points = zeros(2, n+1)
function f_jac(dx, x)
    points .= 0.0
    tether_from_distance_length!(points, x[1], x[2], n, start)
    dx .= reshape(points, 2(n+1))
    nothing
end

x = [desired_distance, l]
J = zeros(2(n+1), 2)
@time FiniteDiff.finite_difference_jacobian!(J, f_jac, x)
@time FiniteDiff.finite_difference_jacobian!(J, f_jac, x)
@show norm(J)

x_coords = points[1, :]
y_coords = points[2, :]
p = plot(x_coords, y_coords, marker=:circle, label="Tether", aspect_ratio=:equal)
scatter!([start[1], points[1, end]], [start[2], points[2, end]], label="Endpoints", color=:red)

scale_factor = 2.0
for i in 1:n+1
    # Index into Jacobian for x and y components
    idx_x = 2i - 1
    idx_y = 2i
    
    # Velocity due to distance change (column 1)
    vx_d = J[idx_x, 1]
    vy_d = J[idx_y, 1]
    
    # Velocity due to length change (column 2)
    vx_l = J[idx_x, 2]
    vy_l = J[idx_y, 2]
    
    # Plot arrows for distance sensitivity (red)
    quiver!([x_coords[i]], [y_coords[i]], 
            quiver=([vx_d * scale_factor], [vy_d * scale_factor]), 
            color=:red, alpha=1.5, label=i==1 ? "Distance sensitivity" : "")
    
    # Plot arrows for length sensitivity (blue)
    quiver!([x_coords[i]], [y_coords[i]], 
            quiver=([vx_l * scale_factor], [vy_l * scale_factor]), 
            color=:blue, alpha=1.5, label=i==1 ? "Length sensitivity" : "")
end

xlabel!("x")
ylabel!("y")
title!("Tether with Sensitivity Vectors")

# Display the plot
display(p)
# end

# mwe()