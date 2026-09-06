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
ACC0 = zeros(3, segments+1)
for i in 1:segments+1
    local l0
    l0 = -(i-1) * L0
    v0 = (i-1) * V0/segments
    POS0[:, i] .= [sin(α0) * l0, 0, cos(α0) * l0]
    VEL0[:, i] .= [sin(α0) * v0, 0, cos(α0) * v0]
end

# defining the model, Z component upwards
@parameters mass=M0 c_spring0=C_SPRING damping=0.5 l_seg=L0
@variables pos(t)[1:3, 1:segments+1]
@variables vel(t)[1:3, 1:segments+1]  = VEL0
@variables acc(t)[1:3, 1:segments+1]  = ACC0
@variables segment(t)[1:3, 1:segments]
@variables unit_vector(t)[1:3, 1:segments]
@variables norm1(t)[1:segments]
@variables rel_vel(t)[1:3, 1:segments]
@variables spring_vel(t)[1:segments]
@variables c_spring(t)[1:segments]
@variables spring_force(t)[1:3, 1:segments]
@variables total_force(t)[1:3, 1:segments+1]

function model()
    eqs = []
    # loop over all segments to calculate the spring forces
    for i in segments:-1:1
        eqs = [
            eqs
            segment[:, i]      ~ pos[:, i+1] - pos[:, i]
            norm1[i]           ~ norm(segment[:, i])
            unit_vector[:, i]  ~ -segment[:, i]/norm1[i]
            rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i]
            spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i]
            c_spring[i]           ~ c_spring0 * (norm1[i] > l_seg)
            spring_force[:, i] ~ (c_spring[i] * (norm1[i] - l_seg) + damping * spring_vel[i]) * unit_vector[:, i]
        ]
    end
    # loop over all tether particles to apply the forces and calculate the accelerations
    for i in 1:(segments+1)
        if i == segments+1
            eqs = [
                eqs
                total_force[:, i] ~ spring_force[:, i-1]
                acc[:, i]         ~ G_EARTH + total_force[:, i] / (100 * mass)
                D(pos[:, i])      ~ vel[:, i]
                D(vel[:, i])      ~ acc[:, i]
            ]
        elseif i == 1
            eqs = [
                eqs
                total_force[:, i] ~ spring_force[:, i]
                acc[:, i]         ~ zeros(3)
                pos[:, i]         ~ zeros(3)
                vel[:, i]         ~ zeros(3)
            ]
        else
            eqs = [
                eqs
                total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i]
                acc[:, i]         ~ G_EARTH + total_force[:, i] / mass
                D(pos)[:, i]         ~ zeros(3)
                vel[:, i]         ~ zeros(3)
            ]
        end
    end
    return eqs
end

     
@named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(model())), t)
simple_sys = structural_simplify(sys)

# running the simulation
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

guesses = [
    [pos[j, i] => POS0[j, i] for j in 1:3 for i in 1:segments+1]
    [segment[j, i] => POS0[j, i+1] - POS0[j, i] for j in 1:3 for i in 1:segments]
]
prob = ODEProblem(simple_sys, nothing, tspan; guesses)
elapsed_time = @elapsed sol = solve(prob, Rodas4(autodiff=ModelingToolkit.AutoFiniteDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
elapsed_time = @elapsed sol = solve(prob, Rodas4(autodiff=ModelingToolkit.AutoFiniteDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
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
