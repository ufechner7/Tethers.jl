using DifferentialEquations, Plots, Plots.PlotMeasures

# Define the parameters
k = 1.0  # Spring constant
m = 1.0  # Mass of P2 and P3
g = 0.0  # No gravity in this case
F2 = 1.0  # Force applied to P2 in y direction
F3 = 1.0  # Force applied to P3 in y direction

# Initial conditions
x1, y1 = 0.0, 0.0  # P1 is fixed at origin
x2, y2 = -1.0, 2.0  # Initial position of P2
x3, y3 = 1.0, 2.0   # Initial position of P3

# Initial velocities
vx2, vy2 = 0.0, 0.0
vx3, vy3 = 0.0, 0.0

# Natural lengths of the springs (assuming rest length is the initial distance)
l1 = sqrt((x2 - x1)^2 + (y2 - y1)^2)
l2 = sqrt((x3 - x2)^2 + (y3 - y2)^2)
l3 = sqrt((x3 - x1)^2 + (y3 - y1)^2)

# Define the system of equations
function mass_spring!(du, u, t)
    x2, y2, x3, y3, vx2, vy2, vx3, vy3 = u
    
    # Compute forces due to springs
    function spring_force(xa, ya, xb, yb, l)
        dx = xb - xa
        dy = yb - ya
        length = sqrt(dx^2 + dy^2)
        force_mag = k * (length - l)
        fx = force_mag * dx / length
        fy = force_mag * dy / length
        return fx, fy
    end

    # Forces from springs
    Fx_s1, Fy_s1 = spring_force(x1, y1, x2, y2, l1)
    Fx_s2, Fy_s2 = spring_force(x2, y2, x3, y3, l2)
    Fx_s3, Fy_s3 = spring_force(x3, y3, x1, y1, l3)

    # Net forces on P2
    Fx2 = -Fx_s1 + Fx_s2
    Fy2 = -Fy_s1 + Fy_s2 + F2

    # Net forces on P3
    Fx3 = -Fx_s2 + Fx_s3
    Fy3 = -Fy_s2 + Fy_s3 + F3

    # Equations of motion
    du[1] = vx2
    du[2] = vy2
    du[3] = vx3
    du[4] = vy3
    du[5] = Fx2 / m
    du[6] = Fy2 / m
    du[7] = Fx3 / m
    du[8] = Fy3 / m
end

# Initial state
u0 = [x2, y2, x3, y3, vx2, vy2, vx3, vy3]

# Time span
tspan = (0.0, 10.0)

# Solve the differential equations
prob = ODEProblem(mass_spring!, u0, tspan)
sol = solve(prob, Tsit5(); saveat=0.05)

# Extract the solution
x2_sol = sol[1, :]
y2_sol = sol[2, :]
x3_sol = sol[3, :]
y3_sol = sol[4, :]

# Create the animation
anim = @animate for i in eachindex(x2_sol)
    # Set fixed axis limits to maintain scaling throughout the animation
    plot(xlim=(-1.5, 1.5), ylim=(-0.5, 6.0), legend=false, framestyle=:box, grid=false, size=(600, 600))

    # Scatter plot for the mass points
    scatter!([0, x2_sol[i], x3_sol[i]], [0, y2_sol[i], y3_sol[i]], markersize=5, label="")

    # Plot individual spring connections
    plot!([0, x2_sol[i]], [0, y2_sol[i]], linewidth=2, color=:blue, label="")  # P1-P2
    plot!([x2_sol[i], x3_sol[i]], [y2_sol[i], y3_sol[i]], linewidth=2, color=:red, label="")  # P2-P3
    plot!([x3_sol[i], 0], [y3_sol[i], 0], linewidth=2, color=:green, label="")  # P3-P1
end

# Save the animation
gif(anim, "mass_spring_simulation.gif", fps=30)