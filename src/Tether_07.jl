"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0) and reel-out. It uses five tether segments. Added parameter NONLINEAR: If it is set to `false`, 
than a linear springs are used.
Damping now increases with the segment length. 
"""

using ModelingToolkit, OrdinaryDiffEq, PyPlot, LinearAlgebra

G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration     [m/s²]
L0::Float64 = 10.0                              # initial segment length            [m]
V0::Float64 = 4                                 # initial velocity of lowest mass [m/s]
segments::Int64 = 1                             # number of tether segments         [-]
POS0 = zeros(3, segments+1)
VEL0 = zeros(3, segments+1)
ACC0 = zeros(3, segments+1)
SEGMENTS0 = zeros(3, segments) 
UNIT_VECTORS0 = zeros(3, segments)
for i in 1:segments+1
    POS0[:, i] .= [0.0, 0, (i-1)*L0]
    VEL0[:, i] .= [0.0, 0, (i-1)*V0/segments]
end
for i in 2:segments+1
    ACC0[:, i] .= G_EARTH
end
for i in 1:segments
    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
    SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
end

# model, Z component upwards
@parameters mass=1.0 c_spring0=50.0 damping=0.5 l_seg=L0
@variables t 
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
D = Differential(t)

eqs1 = vcat(D.(pos) ~ vel,
            D.(vel) ~ acc)
eqs2 = []
for i in 1:segments
    global eqs2
    eqs2 = vcat(eqs2, segment[:, i] ~ pos[:, i+1] - pos[:, i])
    eqs2 = vcat(eqs2, norm1[i] ~ norm(segment[:, i]))
    eqs2 = vcat(eqs2, unit_vector[:, i] ~ -segment[:, i]/norm1[i])
    eqs2 = vcat(eqs2, rel_vel[:, i] ~ vel[:, i+1] - vel[:, i])
    eqs2 = vcat(eqs2, spring_vel[i] ~ -unit_vector[:, i] ⋅ vel[:, i])
    eqs2 = vcat(eqs2, c_spring[i] ~ c_spring0 * (norm1[i] > l_seg))
    eqs2 = vcat(eqs2, spring_force[:, i] ~ (c_spring[i] * (norm1[i] - l_seg) * damping * spring_vel[i]) * unit_vector[:, i])
    # TODO: the spring_force must be distributed
    eqs2 = vcat(eqs2, acc[:, i+1] ~ G_EARTH + spring_force[:, i] / mass)
end
eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
eqs = vcat(eqs1..., eqs2)
     
@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

duration = 10.0
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration
u0 = zeros(12)
u0[9] = L0
u0[12] = V0

prob = ODEProblem(simple_sys, u0, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

X = sol.t
particle = 2
POS_Z = sol(X, idxs=pos[3, particle])
VEL_Z = sol(X, idxs=vel[3, particle])
C_SPRING = sol(X, idxs=c_spring[1])

lns1 = plot(X, POS_Z, color="green", label="pos_z")
xlabel("time [s]")
ylabel("pos_z [m]")
lns2 = plot(X, L0.+0.005 .* C_SPRING, color="grey", label="c_spring")
grid(true)
legend() 
twinx()
ylabel("vel_z [m/s]") 
lns3 = plot(X, VEL_Z, color="red", label="vel_z")
lns = vcat(lns1, lns2, lns3)
labs = [l.get_label() for l in lns]
legend(lns, labs) 
PyPlot.show(block=false)
nothing

