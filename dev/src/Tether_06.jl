# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = 50.0                              # initial tether length             [m]
V_RO::Float64 = 2.0                             # reel-out speed                  [m/s]
D_TETHER::Float64 = 4                           # tether diameter                  [mm]
RHO_TETHER::Float64 = 724.0                     # density of Dyneema            [kg/m³] 
C_SPRING::Float64 = 614600.0                    # unit spring constant              [N]
DAMPING::Float64  = 473                         # unit damping constant            [Ns]
SEGMENTS::Int64 = 5                             # number of tether segments         [-]
α0 = π/10                                       # initial tether angle            [rad]
duration = 10                                   # duration of the simulation        [s]
SAVE = false                                    # save png files in folder video
mass_per_meter::Float64 = RHO_TETHER * π * (D_TETHER/2000.0)^2

# calculating consistant initial conditions
POS0 = zeros(3, SEGMENTS+1)
VEL0 = zeros(3, SEGMENTS+1)
ACC0 = zeros(3, SEGMENTS+1)
SEGMENTS0 = zeros(3, SEGMENTS) 
UNIT_VECTORS0 = zeros(3, SEGMENTS)
for i in 1:SEGMENTS+1
    l0_ = -(i-1)*L0/SEGMENTS
    POS0[:, i] .= [sin(α0) * l0_, 0, cos(α0) * l0_]
    VEL0[:, i] .= [0.0, 0, 0]
end
for i in 1:SEGMENTS
    ACC0[:, i+1] .= G_EARTH
    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
    SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
end

# defining the model, Z component upwards
@parameters c_spring0=C_SPRING/(L0/SEGMENTS) l_seg=L0/SEGMENTS
@variables pos(t)[1:3, 1:SEGMENTS+1]  = POS0
@variables vel(t)[1:3, 1:SEGMENTS+1]  = VEL0
@variables acc(t)[1:3, 1:SEGMENTS+1]  = ACC0
@variables segment(t)[1:3, 1:SEGMENTS]  = SEGMENTS0
@variables unit_vector(t)[1:3, 1:SEGMENTS]  = UNIT_VECTORS0
@variables len(t) = L0
@variables c_spring(t) = c_spring0
@variables damping(t) = DAMPING  / l_seg
@variables m_tether_particle(t) = mass_per_meter * l_seg
@variables norm1(t)[1:SEGMENTS] = l_seg * ones(SEGMENTS)
@variables rel_vel(t)[1:3, 1:SEGMENTS]  = zeros(3, SEGMENTS)
@variables spring_vel(t)[1:SEGMENTS] = zeros(SEGMENTS)
@variables c_spr(t)[1:SEGMENTS] = c_spring0 * ones(SEGMENTS)
@variables spring_force(t)[1:3, 1:SEGMENTS] = zeros(3, SEGMENTS)
@variables total_force(t)[1:3, 1:SEGMENTS] = zeros(3, SEGMENTS)

eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
eqs2 = vcat(eqs1...)

for i in SEGMENTS:-1:1
    global eqs2
    eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
           norm1[i]           ~ norm(segment[:, i]),
           unit_vector[:, i]  ~ -segment[:, i]/norm1[i],
           rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
           spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
           c_spr[i]           ~ c_spring * (norm1[i] > len/SEGMENTS),
           spring_force[:, i] ~ (c_spr[i] * (norm1[i] - (len/SEGMENTS)) + damping * spring_vel[i]) * unit_vector[:, i]]
    if i == SEGMENTS
        push!(eqs, total_force[:, i] ~ spring_force[:, i])
        push!(eqs, acc[:, i+1]       ~ G_EARTH + total_force[:, i] / 0.5*(m_tether_particle))
    else
        push!(eqs, total_force[:, i] ~ spring_force[:, i]- spring_force[:, i+1])
        push!(eqs, acc[:, i+1]       ~ G_EARTH + total_force[:, i] / m_tether_particle)
    end
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end
eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
eqs2 = vcat(eqs2, len .~ L0 + V_RO*t)
eqs2 = vcat(eqs2, c_spring .~ C_SPRING / (len/SEGMENTS))
eqs2 = vcat(eqs2, m_tether_particle .~ mass_per_meter * (len/SEGMENTS))
eqs2 = vcat(eqs2, damping  .~ DAMPING  / (len/SEGMENTS))
     
@named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs2))), t)
simple_sys = structural_simplify(sys)

# running the simulation
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

# plotting the result
function play()
    dt = 0.151
    ylim=(-1.2*(L0+V_RO*duration), 0.5)
    xlim=(-L0/2, L0/2)
    z_max = 0.0
    # text position
    xy = (L0/4.2, z_max-7)
    start = time_ns()
    i = 1
    for time in 0:dt:duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        plot2d(sol[pos][i], time; segments=SEGMENTS, xlim, ylim, xy)
        wait_until(start + 0.5*time*1e9)
    end
    nothing
end
play()
