"""
Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
for l < l_0) and reel-out. It uses five tether segments. Added parameter NONLINEAR: If it is set to `false`, 
than a linear springs are used.
Damping now increases with the segment length. 
"""

using ModelingToolkit, OrdinaryDiffEq, PyPlot, LinearAlgebra

G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration     [m/s²]
L0::Float64 = -10.0                             # initial segment length            [m]
V0::Float64 = 4                                 # initial velocity of lowest mass [m/s]
segments::Int64 = 3                             # number of tether segments         [-]
POS0 = zeros(3, segments+1)
VEL0 = zeros(3, segments+1)
ACC0 = zeros(3, segments+1)
UNIT_VECTORS0 = zeros(3, segments+1)
for i in 1:segments+1
    POS0[:, i] .= [0.0, 0, (i-1)*L0]
    VEL0[:, i] .= [0.0, 0, (i-1)*V0/segments]
    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
end
for i in 2:segments+1
    ACC0[:, i] .= G_EARTH
end

# model, Z component upwards
@parameters mass=1.0 c_spring0=50.0 damping=0.5 l_seg=L0
@variables t 
@variables pos(t)[1:3, 1:segments+1]  = POS0
@variables vel(t)[1:3, 1:segments+1]  = VEL0
@variables acc(t)[1:3, 1:segments+1]  = ACC0
# @variables unit_vector(t)[1:3, 1:segments]  = UNIT_VECTORS0
# @variables norm1(t)[1:segments] = l_seg * ones(segments)
# @variables c_spring(t) = c_spring0
# @variables spring_force(t)[1:3] = [0.0, 0.0, 0.0]
# @variables force(t) = 0.0 norm1(t) = abs(l0) spring_vel(t) = 0.0
D = Differential(t)

eqs = vcat(D.(pos) ~ vel,
           D.(vel) ~ acc,
           acc    .~ ACC0)
     
# eqs = vcat(D.(pos)      ~ vel,
#            D.(vel)      ~ acc,
#            norm1        ~ norm(pos),
#            unit_vector  ~ -pos/norm1,         # direction from point mass to origin
#            spring_vel   ~ -unit_vector ⋅ vel,
#            c_spring     ~ c_spring0 * (norm1 > abs(l0)),
#            spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
#            acc          ~ G_EARTH + spring_force/mass)

@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)
