# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0) and n tether segments. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = 5.0                               # initial segment length            [m]
V0::Float64 = 2                                 # initial velocity of lowest mass [m/s]
M0::Float64 = 0.5                               # mass per particle                [kg]
C_SPRING::Float64 = 50                          # spring constant
segments::Int64 = 5                             # number of tether segments         [-]
α0 = π/10                                       # initial tether angle            [rad]
duration = 10.0                                 # duration of the simulation        [s]
POS0 = zeros(3, segments+1)
VEL0 = zeros(3, segments+1)
for i in 1:segments+1
    local l0
    l0 = -(i-1) * L0
    v0 = (i-1) * V0/segments
    POS0[:, i] .= [sin(α0) * l0, 0, cos(α0) * l0]
    VEL0[:, i] .= [sin(α0) * v0, 0, cos(α0) * v0]
end

# defining the model, Z component upwards
@parameters mass=M0 c_spring0=C_SPRING damping=0.5 l_seg=L0
@variables pos(t)[1:3, 1:segments+1]  = POS0
@variables vel(t)[1:3, 1:segments+1]  = VEL0
@variables acc(t)[1:3, 1:segments+1]
@variables segment(t)[1:3, 1:segments]
@variables unit_vector(t)[1:3, 1:segments]
@variables norm1(t)[1:segments]
@variables rel_vel(t)[1:3, 1:segments]
@variables spring_vel(t)[1:segments]
@variables c_spring(t)[1:segments]
@variables spring_force(t)[1:3, 1:segments]
@variables total_force(t)[1:3, 1:segments+1]

# basic differential equations
eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
eqs2 = vcat(eqs1...)
# loop over all segments to calculate the spring forces
for i in segments:-1:1
    global eqs2; local eqs
    eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
           norm1[i]           ~ norm(segment[:, i]),
           unit_vector[:, i]  ~ -segment[:, i]/norm1[i],
           rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
           spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
           c_spring[i]           ~ c_spring0 * (norm1[i] > l_seg),
           spring_force[:, i] ~ (c_spring[i] * (norm1[i] - l_seg) + damping * spring_vel[i]) * unit_vector[:, i]]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end
# loop over all tether particles to apply the forces and calculate the accelerations
for i in 1:(segments+1)
    global eqs2; local eqs
    eqs = []
    if i == segments+1
        push!(eqs, total_force[:, i] ~ spring_force[:, i-1])
        push!(eqs, acc[:, i]         ~ G_EARTH + total_force[:, i] / (0.5 * mass))
    elseif i == 1
        push!(eqs, total_force[:, i] ~ spring_force[:, i])
        push!(eqs, acc[:, i]         ~ zeros(3))
    else
        push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i])
        push!(eqs, acc[:, i]         ~ G_EARTH + total_force[:, i] / mass)
    end
    eqs2 = vcat(eqs2, reduce(vcat, eqs))
end
     
@named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
simple_sys = structural_simplify(sys)

# running the simulation
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, tspan)
elapsed_time = @elapsed sol = solve(prob, KenCarp4(autodiff=false); dt, abstol=tol, reltol=tol, saveat=ts)
elapsed_time = @elapsed sol = solve(prob, KenCarp4(autodiff=false); dt, abstol=tol, reltol=tol, saveat=ts)
println("Elapsed time: $(elapsed_time) s, speed: $(round(duration/elapsed_time)) times real-time")

function play()
    dt = 0.15
    ylim = (-1.2 * segments * L0, 0.5)
    xlim = (-segments * L0/2, segments * L0/2)
    z_max = 0.0
    # text position
    xy = (segments * L0/4.2, z_max-3.0*segments/5)
    start = time_ns()
    i = 1
    for time in 0:dt:duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        plot2d(sol[pos][i], time; segments, xlim, ylim, xy)
        wait_until(start + 0.5 * time * 1e9)
    end
    nothing
end
play()
