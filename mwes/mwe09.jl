using Plots, LinearAlgebra, Roots, ForwardDiff

function generate_tether!(points::AbstractMatrix{T}, n, s, start, total_angle) where T
    segment_angle = total_angle / n
    points[:, 1] .= start
    initial_angle = zero(T)
    for i in 2:n+1
        angle = initial_angle + (i - 2) * segment_angle
        points[1, i] = points[1, i-1] + s * cos(angle)
        points[2, i] = points[2, i-1] + s * sin(angle)
    end
    return nothing
end

function calculate_distance!(points::AbstractMatrix{T}, d::Ref{T}) where T
    dx = points[1, end] - points[1, 1]
    dy = points[2, end] - points[2, 1]
    d[] = sqrt(dx^2 + dy^2)
    return nothing
end

function distance_for_angle!(points::AbstractMatrix{T}, d::Ref{T}, n, s, start, total_angle) where T
    generate_tether!(points, n, s, start, total_angle)
    calculate_distance!(points, d)
    return nothing
end

# Example usage
n = 10 # number of segments
s = 1.0 # segment length
start = [0.0, 0.0] # start point coordinates
desired_distance = 8.0 # desired distance

# Preallocate points array
points = zeros(2, n+1)
# Preallocate distance reference
d = Ref(0.0)
cost = Ref(0.0)

# Define a function to find the root of
@inline function f!(total_angle)
    distance_for_angle!(points, d, n, s, start, total_angle)
    cost[] = d[] - desired_distance
    return cost[]
end

@time f!(0.1)
# Solve the problem using a non-allocating method
result = @time find_zero(f!, one(Float64))
println("The total_angle that results in the desired distance is: $result")

# Generate and plot the tether for the found total_angle
generate_tether!(points, n, s, start, result)
x_coords = points[1, :]
y_coords = points[2, :]
plot(x_coords, y_coords, marker=:circle, label="Tether", aspect_ratio=:equal)
scatter!([start[1], points[1, end]], [start[2], points[2, end]], label="Endpoints", color=:red)
xlabel!("x")
ylabel!("y")
title!("Tether in the xy plane")

function tether_from_distance_length!(points::AbstractMatrix{T}, distance::T, length::T, n, start) where T
    d = Ref{T}(zero(T))
    cost = Ref{T}(zero(T))
    
    function f_zero!(total_angle)
        distance_for_angle!(points, d, n, length/n, start, total_angle)
        cost[] = d[] - distance
        return cost[]
    end
    result = find_zero(f_zero!, one(T))
    return nothing
end

function f_jac!(y::AbstractVector{T}, x::AbstractVector{T}) where T
    @show x
    points = reshape(y, (2, n+1))
    tether_from_distance_length!(points, x[1], x[2], n, start)
    return nothing
end

y = reshape(points, 2*(n+1))
x = [desired_distance, result]
J = ForwardDiff.jacobian(f_jac!, y, x)
@show J[:, 1]

# Generate and plot the tether for the found total_angle
generate_tether!(points, n, s, start, result)
x_coords = points[1, :]
y_coords = points[2, :]

# Create the base plot
p = plot(x_coords, y_coords, marker=:circle, label="Tether", aspect_ratio=:equal)
scatter!([start[1], points[1, end]], [start[2], points[2, end]], label="Endpoints", color=:red)

# Extract velocities from the Jacobian
# Reshape Jacobian columns into point velocities
# First column is sensitivity to distance, second to length
scale_factor = 0.5  # Adjust this to make arrows more visible
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
            color=:red, alpha=0.5, label=i==1 ? "Distance sensitivity" : "")
    
    # # Plot arrows for length sensitivity (blue)
    # quiver!([x_coords[i]], [y_coords[i]], 
    #         quiver=([vx_l * scale_factor], [vy_l * scale_factor]), 
    #         color=:blue, alpha=1.5, label=i==1 ? "Length sensitivity" : "")
end

xlabel!("x")
ylabel!("y")
title!("Tether with Sensitivity Vectors")

# Display the plot
display(p)