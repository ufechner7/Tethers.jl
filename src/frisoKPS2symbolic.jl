using Timers
tic()
using ModelingToolkit, OrdinaryDiffEq, Parameters

import Plots
toc()

# 🔹 Define struct for simulation settings
@with_kw mutable struct SimulationSettings @deftype Float64
    g_earth::Float64 = 0  # Gravity [m/s²]
    k::Float64 = 1.0  # Spring constant [N/m]
    m::Float64 = 1.0  # Mass of P2 and P3 [kg]
    F2::Float64 = 1  # External force on P2 [N]
    F3::Float64 = 1  # External force on P3 [N]
    l1::Float64 = sqrt(5)  # Rest length of spring P1-P2 [m]
    l2::Float64 = 2.0  # Rest length of spring P2-P3 [m]
    l3::Float64 = sqrt(5)  # Rest length of spring P3-P1 [m]
    duration::Float64 = 10.0  # Simulation time [s]
    save::Bool = false  # Whether to save animation frames
end

# Now, update the instance name:
se = SimulationSettings()

# 🔹 Define symbolic parameters (Using values from `Settings3`)
@parameters t k m g F2 F3 l1 l2 l3

# 🔹 Define symbolic state variables (positions & velocities)
@variables x2(t) y2(t) x3(t) y3(t) vx2(t) vy2(t) vx3(t) vy3(t)

# 🔹 Define differential operator
D = Differential(t)    

# 🔹 Define Hooke’s Law for spring force
function spring_force(xa, ya, xb, yb, l)
    dx = xb - xa
    dy = yb - ya
    length = sqrt(dx^2 + dy^2)  # Compute current length
    force_mag = k * (length - l)  # Hooke’s Law: F = -k (x - L0)
    fx = force_mag * dx / length
    fy = force_mag * dy / length
    
    return [fx, fy]
end

# 🔹 Define fixed anchor point P1 at (0,0)
x1, y1 = 0.0, 0.0

# 🔹 Compute forces symbolically
F_s1 = spring_force(x1, y1, x2, y2, l1)  # Force from P1 to P2
F_s2 = spring_force(x2, y2, x3, y3, l2)  # Force from P2 to P3
F_s3 = spring_force(x3, y3, x1, y1, l3)  # Force from P3 to P1

# 🔹 Net forces on P2
Fx2 = -F_s1[1] + F_s2[1]
Fy2 = -F_s1[2] + F_s2[2] + F2

# 🔹 Net forces on P3
Fx3 = -F_s2[1] + F_s3[1]
Fy3 = -F_s2[2] + F_s3[2] + F3

# 🔹 Equations of motion using Newton's second law: F = ma
eqs = [
    D(x2) ~ vx2,
    D(y2) ~ vy2,
    D(x3) ~ vx3,
    D(y3) ~ vy3,
    D(vx2) ~ Fx2 / m,
    D(vy2) ~ (Fy2 - m * g) / m,
    D(vx3) ~ Fx3 / m,
    D(vy3) ~ (Fy3 - m * g) / m
]

# 🔹 Create an ODE system
@named sys = ODESystem(eqs, t, [x2, y2, x3, y3, vx2, vy2, vx3, vy3], [m, k, g, F2, F3, l1, l2, l3])

# 🔹 Parameter values from `Settings3`
params = Dict(
    m => se.m,
    k => se.k,
    g => se.g_earth,
    F2 => se.F2,
    F3 => se.F3,
    l1 => se.l1,
    l2 => se.l2,
    l3 => se.l3
)

# 🔹 Initial conditions
u0 = Dict(
    x2 => -1.0, y2 => 2.0, x3 => 1.0, y3 => 2.0,
    vx2 => 0.0, vy2 => 0.0, vx3 => 0.0, vy3 => 0.0
)

# 🔹 Time span
tspan = (0.0, se.duration)

# 🔹 Convert symbolic system to an ODEProblem
prob = ODEProblem(complete(sys), u0, tspan, params)

# 🔹 Solve numerically
sol = solve(prob, Tsit5(), saveat=0.05)

# 🔹 Extract the solution
x2_sol = sol[x2, :]
y2_sol = sol[y2, :]
x3_sol = sol[x3, :]
y3_sol = sol[y3, :]
 
# 🔹 Create animation
anim = Plots.@animate for i in eachindex(x2_sol)
    Plots.plot(xlim=(-1.5, 1.5), ylim=(-0.5, 6.0), legend=false, framestyle=:box, grid=false, size=(600, 600))

    # Scatter plot for mass points
    Plots.scatter!([0, x2_sol[i], x3_sol[i]], [0, y2_sol[i], y3_sol[i]], markersize=5, label="")

    # Spring connections
    Plots.plot!([0, x2_sol[i]], [0, y2_sol[i]], linewidth=2, color=:blue, label="")  # P1-P2
    Plots.plot!([x2_sol[i], x3_sol[i]], [y2_sol[i], y3_sol[i]], linewidth=2, color=:red, label="")  # P2-P3
    Plots.plot!([x3_sol[i], 0], [y3_sol[i], 0], linewidth=2, color=:green, label="")  # P3-P1
end

# # 🔹 Save animation
Plots.gif(anim, "video/symbolic_mass_spring_simulation.gif", fps=30)
nothing