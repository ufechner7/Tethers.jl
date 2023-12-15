# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, PyPlot, LinearAlgebra, Timers

G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration     [m/s²]
L0::Float64 = 25.0                              # initial tether length             [m]
V0::Float64 = 2                                 # initial velocity of lowest mass [m/s]
V_RO::Float64 = 2.0                             # reel-out speed                  [m/s]
M0::Float64 = 0.5                               # mass per particle                [kg]
D_TETHER::Float64 = 4                           # tether diameter                  [mm]
RHO_TETHER::Float64 = 724.0                     # densitiy of Dyneema           [kg/m³] 
C_SPRING::Float64 = 250                         # unit spring constant              [N]
DAMPING::Float64  = 2.5                         # unit damping constant            [Ns]
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
    l0 = -(i-1)*L0/segments
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

# model, Z component upwards
@parameters mass=M0 c_spring0=C_SPRING/(L0/segments) l_seg=L0/segments
@variables t 
@variables pos(t)[1:3, 1:segments+1]  = POS0
@variables vel(t)[1:3, 1:segments+1]  = VEL0
@variables acc(t)[1:3, 1:segments+1]  = ACC0
@variables segment(t)[1:3, 1:segments]  = SEGMENTS0
@variables unit_vector(t)[1:3, 1:segments]  = UNIT_VECTORS0
@variables length(t) = L0
@variables c_spring(t) = c_spring0
@variables damping(t) = DAMPING  / l_seg
@variables norm1(t)[1:segments] = l_seg * ones(segments)
@variables rel_vel(t)[1:3, 1:segments]  = zeros(3, segments)
@variables spring_vel(t)[1:segments] = zeros(segments)
@variables c_spr(t)[1:segments] = c_spring0 * ones(segments)
@variables spring_force(t)[1:3, 1:segments] = zeros(3, segments)
@variables total_force(t)[1:3, 1:segments] = zeros(3, segments)
D = Differential(t)

eqs1 = vcat(D.(pos) ~ vel,
            D.(vel) ~ acc)
eqs2 = []
for i in segments:-1:1
    global eqs2
    eqs2 = vcat(eqs2, segment[:, i] ~ pos[:, i+1] - pos[:, i])
    eqs2 = vcat(eqs2, norm1[i] ~ norm(segment[:, i]))
    eqs2 = vcat(eqs2, unit_vector[:, i] ~ -segment[:, i]/norm1[i])
    eqs2 = vcat(eqs2, rel_vel[:, i] ~ vel[:, i+1] - vel[:, i])
    eqs2 = vcat(eqs2, spring_vel[i] ~ -unit_vector[:, i] ⋅ rel_vel[:, i])
    eqs2 = vcat(eqs2, c_spr[i] ~ c_spring * (norm1[i] > length/segments))
    eqs2 = vcat(eqs2, spring_force[:, i] ~ (c_spr[i] * (norm1[i] - l_seg) + damping * spring_vel[i]) * unit_vector[:, i])
    if i == segments
        eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i])
    else
        eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i]- spring_force[:, i+1])
    end
    eqs2 = vcat(eqs2, acc[:, i+1] .~ G_EARTH + total_force[:, i] / mass)
end
eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
eqs2 = vcat(eqs2, length ~ L0 + V_RO*t)
eqs2 = vcat(eqs2, c_spring ~ C_SPRING / (length/segments))
eqs2 = vcat(eqs2, damping  ~ DAMPING  / (length/segments))
eqs = vcat(eqs1..., eqs2)
     
@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

# initial state
u0 = Dict(pos=>POS0, vel=>VEL0)

prob = ODEProblem(simple_sys, u0, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

function plot2d(sol, reltime, segments, line, sc, txt)
    index = Int64(round(reltime*50+1))
    x, z = Float64[], Float64[]
    for particle in 1:segments+1
        push!(x, (sol(sol.t, idxs=pos[1, particle]))[index])
        push!(z, (sol(sol.t, idxs=pos[3, particle]))[index])
    end
    z_max = maximum(z)
    if isnothing(line)
        line, = plot(x,z; linewidth="1")
        sc  = scatter(x, z; s=15, color="red") 
        txt = annotate("t=$(round(reltime,digits=1)) s",  
                        xy=(L0/4.2, z_max-7.0*segments/5), fontsize = 12)
        PyPlot.show(block=false)
    else
        line.set_xdata(x)
        line.set_ydata(z)
        sc.set_offsets(hcat(x,z))
        txt.set_text("t=$(round(reltime,digits=1)) s")
        gcf().canvas.draw()
    end
    line, sc, txt
end

function play()
    dt = 0.15
    ylim(-1.2*(L0+V_RO*duration), 0.5)
    xlim(-L0/2, L0/2)
    grid(true; color="grey", linestyle="dotted")
    line, sc, txt = nothing, nothing, nothing
    start = time_ns()
    for time in 0:dt:duration
        line, sc, txt = plot2d(sol, time, segments, line, sc, txt)
        wait_until(start + 0.5*time*1e9)
    end
    nothing
end
play()
