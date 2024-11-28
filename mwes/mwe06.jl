using ModelingToolkit, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D


@variables q(t)[1:4] ω(t)[1:3]
@parameters I[1:3] τ[1:3]

# Quaternion multiplication matrix
Ω = [0      -ω[1]   -ω[2]   -ω[3];
     ω[1]   0       ω[3]    -ω[2];
     ω[2]   -ω[3]   0       ω[1];
     ω[3]   ω[2]    -ω[1]   0]

# Quaternion dynamics
eqs_q = [D(q[i]) ~ 0.5 * sum(Ω[i, j] * q[j] for j in 1:4) for i in 1:4]

# Rotational dynamics (Euler's equations) https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)
eqs_ω = [
    D(ω[1]) ~ (τ[1] + (I[2] - I[3]) * ω[2] * ω[3]) / I[1],
    D(ω[2]) ~ (τ[2] + (I[3] - I[1]) * ω[3] * ω[1]) / I[2],
    D(ω[3]) ~ (τ[3] + (I[1] - I[2]) * ω[1] * ω[2]) / I[3],
]

# Combine all equations
eqs = [eqs_q; eqs_ω]

# Define the system
@mtkbuild sys = ODESystem(eqs, t)

# Initial conditions
q0 = [1.0, 0.0, 0.0, 0.0]  # Initial quaternion (identity rotation)
ω0 = [0.0, 0.0, 0.0]       # Initial angular velocity
u0 = [q0; ω0]

# Time span
tspan = (0.0, 10.0)

# Parameters: moments of inertia and applied torques
I_values = [0.1, 0.2, 0.3] # Example moments of inertia
τ_values = [0.1, 0.2, 0.3]  # Example applied torques
p = [I => I_values, τ => τ_values]

# Create and solve the ODE problem
dt = 0.02
tol = 1e-6
tspan = (0.0, 10.0)
ts    = 0:dt:10.0
prob = ODEProblem(sys, u0, tspan, p)
elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)

p = plotx(sol.t, 
    [sol[i, :] for i in 1:7]...,
    )
display(p)