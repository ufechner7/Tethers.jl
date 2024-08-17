# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0) and n tether segments. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers
using ModelingToolkit: t_nounits as t, D_nounits as D

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = 5.0                               # initial segment length            [m]
V0::Float64 = 2                                 # initial velocity of lowest mass [m/s]
M0::Float64 = 0.5                               # mass per particle                [kg]
C_SPRING::Float64 = 50                          # spring constant
segments::Int64 = 5                             # number of tether segments         [-]
α0 = π/10                                       # initial tether angle            [rad]
duration = 30.0                                 # duration of the simulation        [s]
POS0 = zeros(3, segments+1)
VEL0 = zeros(3, segments+1)
ACC0 = zeros(3, segments+1)
SEGMENTS0 = zeros(3, segments) 
UNIT_VECTORS0 = zeros(3, segments)
for i in 1:segments+1
    local l0
    l0 = -(i-1)*L0
    v0 = (i-1)*V0/segments
    POS0[:, i] .= [sin(α0) * l0, 0, cos(α0) * l0]
    VEL0[:, i] .= [sin(α0) * v0, 0, cos(α0) * v0]
end
for i in 2:segments+1
    ACC0[:, i] .= G_EARTH
end
for i in 1:segments
    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
    SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
end

# defining the model, Z component upwards
@parameters mass=M0 c_spring0=C_SPRING damping=0.5 l_seg=L0
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
@variables total_force(t)[1:3, 1:segments+1] = zeros(3, segments+1)

eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
eqs2 = vcat(eqs1...)
for i in segments:-1:1
    global eqs2; local eqs
    eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
           norm1[i]           ~ norm(segment[:, i]),
           unit_vector[:, i]  ~ -segment[:, i]/norm1[i],
           rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
           spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
           c_spring[i]           ~ c_spring0 * (norm1[i] > l_seg),
           spring_force[:, i] ~ (c_spring[i] * (norm1[i] - l_seg) + damping * spring_vel[i]) * unit_vector[:, i]]
    if i == segments
        push!(eqs, total_force[:, i+1] ~ spring_force[:, i])
        push!(eqs, acc[:, i+1]         ~ G_EARTH + total_force[:, i+1] / (0.5 * mass))
        push!(eqs, total_force[:, i]   ~ spring_force[:, i-1] - spring_force[:, i])
    elseif i == 1
        push!(eqs, total_force[:, i]   ~ spring_force[:, i])
        push!(eqs, acc[:, i+1]         ~ G_EARTH + total_force[:, i+1] / mass)
    else
        push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i])
        push!(eqs, acc[:, i+1]       ~ G_EARTH + total_force[:, i+1] / mass)
    end
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end
eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
     
@named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs2))), t)
simple_sys = structural_simplify(sys)

# running the simulation
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

function play()
    dt = 0.15
    ylim = (-1.2*segments*L0, 0.5)
    xlim = (-segments*L0/2, segments*L0/2)
    z_max = 0.0
    # text position
    xy = (segments*L0/4.2, z_max-3.0*segments/5)
    start = time_ns()
    i = 1
    for time in 0:dt:duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        plot2d(sol[pos][i], time; segments, xlim, ylim, xy)
        wait_until(start + 0.5*time*1e9)
    end
    nothing
end
play()
