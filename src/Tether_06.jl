# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers

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
mass_per_meter::Float64 = RHO_TETHER * SEGMENTS * (D_TETHER/2000.0)^2

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
@variables t 
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
D = Differential(t)

eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
eqs2 = []
for i in SEGMENTS:-1:1
    global eqs2
    eqs2 = vcat(eqs2, segment[:, i] .~ pos[:, i+1] - pos[:, i])
    eqs2 = vcat(eqs2, norm1[i] .~ norm(segment[:, i]))
    eqs2 = vcat(eqs2, unit_vector[:, i] .~ -segment[:, i]/norm1[i])
    eqs2 = vcat(eqs2, rel_vel[:, i] .~ vel[:, i+1] - vel[:, i])
    eqs2 = vcat(eqs2, spring_vel[i] .~ -unit_vector[:, i] ⋅ rel_vel[:, i])
    eqs2 = vcat(eqs2, c_spr[i] .~ c_spring * (norm1[i] > len/SEGMENTS))
    eqs2 = vcat(eqs2, spring_force[:, i] .~ (c_spr[i] * (norm1[i] - (len/SEGMENTS)) + damping * spring_vel[i]) * unit_vector[:, i])
    if i == SEGMENTS
        eqs2 = vcat(eqs2, total_force[:, i] .~ spring_force[:, i])
        eqs2 = vcat(eqs2, acc[:, i+1] .~ G_EARTH + total_force[:, i] / 0.5*(m_tether_particle))
    else
        eqs2 = vcat(eqs2, total_force[:, i] .~ spring_force[:, i]- spring_force[:, i+1])
        eqs2 = vcat(eqs2, acc[:, i+1] .~ G_EARTH + total_force[:, i] / m_tether_particle)
    end
end
eqs2 = vcat(eqs2, acc[:, 1] .~ zeros(3))
eqs2 = vcat(eqs2, len .~ L0 + V_RO*t)
eqs2 = vcat(eqs2, c_spring .~ C_SPRING / (len/SEGMENTS))
eqs2 = vcat(eqs2, m_tether_particle .~ mass_per_meter * (len/SEGMENTS))
eqs2 = vcat(eqs2, damping  .~ DAMPING  / (len/SEGMENTS))
eqs = vcat(eqs1..., eqs2)
     
@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

# running the simulation
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration

prob = ODEProblem(simple_sys, nothing, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)

function plot2d(sol, reltime, segments, line, sc, txt, j)
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
                        xy=(L0/4.2, z_max-7), fontsize = 12)
    else
        line.set_xdata(x)
        line.set_ydata(z)
        sc.set_offsets(hcat(x,z))
        txt.set_text("t=$(round(reltime,digits=1)) s")
        gcf().canvas.draw()
    end
    if SAVE
        PyPlot.savefig("video/"*"img-"*lpad(j,4,"0"))
    end
    line, sc, txt
end

function play()
    PyPlot.close()
    dt = 0.151
    ylim(-1.2*(L0+V_RO*duration), 0.5)
    xlim(-L0/2, L0/2)
    grid(true; color="grey", linestyle="dotted")
    tight_layout(rect=(0, 0, 0.98, 0.98))
    line, sc, txt = nothing, nothing, nothing
    start = time_ns()
    mkpath("video")
    for (j, time) in pairs(0:dt:duration)
        line, sc, txt = plot2d(sol, time, SEGMENTS, line, sc, txt, j)
        wait_until(start + 0.5*time*1e9)
    end
    nothing
end
play()
