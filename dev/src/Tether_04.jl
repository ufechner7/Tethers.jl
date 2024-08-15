"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0) and n tether segments. 
"""
# TODO: Distribute force correctly
# TODO: Add 2D plot

using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = 10.0                              # initial segment length            [m]
V0::Float64 = 4                                 # initial velocity of lowest mass [m/s]
segments::Int64 = 2                             # number of tether segments         [-]
POS0 = zeros(3, segments+1)
VEL0 = zeros(3, segments+1)
ACC0 = zeros(3, segments+1)
SEGMENTS0 = zeros(3, segments) 
UNIT_VECTORS0 = zeros(3, segments)
for i in 1:segments+1
    POS0[:, i] .= [0.0, 0, -(i-1)*L0]
    VEL0[:, i] .= [0.0, 0, (i-1)*V0/segments]
end
for i in 2:segments+1
    ACC0[:, i] .= G_EARTH
end
for i in 1:segments
    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
    SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
end

# defining the model, Z component upwards
@parameters mass=1.0 c_spring0=50.0 damping=0.5 l_seg=L0
@variables pos(t)[1:3, 1:segments+1]  = POS0
@variables vel(t)[1:3, 1:segments+1]  = VEL0
@variables acc(t)[1:3, 1:segments+1]  = ACC0
@variables segment(t)[1:3, 1:segments]  = SEGMENTS0
@variables unit_vector(t)[1:3, 1:segments]  = UNIT_VECTORS0
@variables norm1(t)[1:segments] = l_seg * ones(segments)
@variables rel_vel(t)[1:3, 1:segments]  = zeros(3, segments)
@variables spring_vel(t)[1:segments] = zeros(segments)
@variables c_spring(t)[1:segments] = c_spring0 * ones(segments)
@variables spring_force(t)[1:3, 1:segments] = zeros(3, segments)

eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
eqs2 = []
for i in 1:segments
    global eqs2
    eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
           norm1[i]           ~ norm(segment[:, i]),
           unit_vector[:, i]  ~ -segment[:, i]/norm1[i],
           rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
           spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
           c_spring[i]        ~ c_spring0 * (norm1[i] > l_seg),
           spring_force[:, i] ~ (c_spring[i] * (norm1[i] - l_seg) + damping * spring_vel[i]) * unit_vector[:, i],
           # TODO: the spring_force must be distributed
           acc[:, i+1]        ~ G_EARTH + spring_force[:, i] / mass]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end
eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
eqs3 = vcat(eqs1..., eqs2)
     
@named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs3))), t)
simple_sys = structural_simplify(sys)

# running the simulation
duration = 10.0
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

function plt(sol, particle)
    X = sol.t
    POS_Z = stack(sol[pos], dims=1)[:,3,particle]
    VEL_Z = stack(sol[vel], dims=1)[:,3,particle]
    C_SPRING = stack(sol[c_spring], dims=1)[:, particle-1]

    p = plot(X, [POS_Z, -L0.+0.005 .* C_SPRING], VEL_Z; xlabel="time [s]", ylabels=["pos_z [m]", "vel_z [m/s]"], 
    labels=["pos_z [m]", "c_spring", "vel_z [m/s]"], fig="segmented tether")
    display(p)
    nothing
end
# plt(sol, 2)
plt(sol, 3)
nothing

