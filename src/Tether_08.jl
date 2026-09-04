# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
# New feature: A steady state solver is used to find the initial tether shape for any
# given pair of endpoints, which is then used as the initial condition for the simulation.
using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots

using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

"""
    Settings3

Simulation parameters for the 3D mass-spring tether model with a nonlinear
spring (1% stiffness for `l < l_0`), tether drag and reel-in/reel-out, used by
[`calc_initial_state`](@ref) and the main simulation loop.

# Fields
- `g_earth::Vector{Float64}`: gravitational acceleration vector [m/s²]
- `v_wind_tether::Vector{Float64}`: wind velocity acting on the tether [m/s]
- `rho`: air density [kg/m³]
- `cd_tether`: drag coefficient of the tether [-]
- `l0`: initial tether length [m]
- `v_ro`: reel-out speed [m/s]
- `d_tether`: tether diameter [mm]
- `rho_tether`: density of the tether material (Dyneema) [kg/m³]
- `c_spring`: unit spring constant [N]
- `rel_compression_stiffness`: relative compression stiffness [-]
- `damping`: unit damping constant [Ns]
- `segments::Int64`: number of tether segments [-]
- `α0`: initial tether angle [rad]
- `duration`: duration of the simulation [s]
- `save::Bool`: if `true`, save PNG files to the `video` folder
"""
@with_kw mutable struct Settings3 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [2, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 70                                      # initial tether length             [m]
    v_ro = 0.3                                   # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 10                         # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 30                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end
                              
"""
    calc_initial_state(se; p1, p2)

Compute a linearly interpolated initial position and velocity for each tether
particle between the endpoints `p1` and `p2`, for the settings `se`.

If `p2` is `nothing`, it is derived from `p1` using `se.α0` and `se.l0`; if `p1`
is `nothing`, it is derived from `p2` the same way. At least one of them must be
given.

Returns `(POS0, VEL0)`, each a `3 × (se.segments+1)` matrix; `VEL0` is all zeros.
"""
function calc_initial_state(se; p1, p2)
    isnothing(p1) && isnothing(p2) && error("at least one of p1 and p2 must be defined")
    # calculate the missing endpoint based on se.α0 and se.l0
    if isnothing(p2)
        z  = cos(se.α0) * se.l0
        y  = sin(se.α0) * se.l0
        p2 = [p1[1], p1[2] - y, p1[3] - z]
        println("p2: ", p2)
    elseif isnothing(p1)
        z  = cos(se.α0) * se.l0
        y  = sin(se.α0) * se.l0
        p1 = [p2[1], p2[2] + y, p2[3] + z]
        println("p1: ", p1)
    end
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    # use a linear interpolation between p1 and p2 for the intermediate points
    for i in 1:se.segments+1
        Δ = (p2-p1) / se.segments
        POS0[:, i] .= p1 + (i-1) * Δ
    end
    POS0, VEL0
end

"""
    model(se; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)

Build the ModelingToolkit tether model for the settings `se`, using a steady-state
solve to find a physically consistent initial tether shape between the endpoints
`p1` and `p2`.

`p1` and `p2` are the fixed or free endpoints of the tether, each a length-3 vector
or `nothing`. `fix_p1` and `fix_p2` control whether the corresponding endpoint is
held fixed (`true`) or free to accelerate under gravity and tether forces (`false`);
an endpoint that is `nothing` cannot be fixed. An endpoint that is `nothing` is
derived from the other one using `se.α0` and `se.l0`, so at least one of them must
be given.

Internally, this first builds the model with `se.v_ro` set to zero and solves for
the steady-state tether positions `POS0`, then rebuilds the model with the original
`se.v_ro` and `POS0` as initial condition.

Returns `(simple_sys, sys, pos, vel, len, c_spr)`: the compiled (`simple_sys`) and
uncompiled (`sys`) `ModelingToolkit.System`, and the symbolic variables `pos`, `vel`,
`len` and `c_spr` used to build it.
"""
function model(se; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)
    if ! isnothing(p1)
        @assert isa(p1, AbstractVector) || error("p1 must be a vector")
        @assert (length(p1) == 3)       || error("p1 must have length 3")
    else
        @assert ! fix_p1                || error("if p1 undefined it cannot be fixed")
    end

    if ! isnothing(p2)
        @assert isa(p2, AbstractVector) || error("p2 must be a vector")
        @assert (length(p2) == 3)       || error("p2 must have length 3")
    else
        @assert ! fix_p2                || error("if p2 undefined it cannot be fixed")
    end
    # straight line approximation for the tether
    POS0, VEL0 = calc_initial_state(se; p1, p2)
    # find steady state
    v_ro = se.v_ro      # save the reel-out speed
    se.v_ro = 0         # v_ro must be zero, otherwise finding the steady state is not possible
    local sol1, pos
    try
        simple_sys, sys, pos, =  model(se, p1, p2, true, true, POS0, VEL0)
        tspan = (0.0, se.duration)
        prob = ODEProblem(simple_sys, nothing, tspan)
        prob1 = SteadyStateProblem(prob)
        sol1 = solve(prob1, DynamicSS(KenCarp4(autodiff=false)))
    finally
        se.v_ro = v_ro  # restore the reel-out speed, also if the steady state solver failed
    end
    SciMLBase.successful_retcode(sol1) ||
        error("Steady state solver failed with return code $(sol1.retcode)!")
    POS0 = sol1[pos]
    # create the real model
    model(se, p1, p2, fix_p1, fix_p2, POS0, VEL0)
end

# internal helper: builds the system given precomputed initial conditions POS0, VEL0
function model(se, p1, p2, fix_p1, fix_p2, POS0, VEL0)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters rel_compression_stiffness = se.rel_compression_stiffness
    @variables begin 
        pos(t)[1:3, 1:se.segments+1]  = POS0
        vel(t)[1:3, 1:se.segments+1]  = VEL0
        acc(t)[1:3, 1:se.segments+1]
        segment(t)[1:3, 1:se.segments]
        unit_vector(t)[1:3, 1:se.segments]
        l_spring(t), c_spring(t), damping(t), m_tether_particle(t)
        len(t)[1:se.segments]
        rel_vel(t)[1:3, 1:se.segments]
        spring_vel(t)[1:se.segments]
        c_spr(t)[1:se.segments]
        spring_force(t)[1:3, 1:se.segments]
        v_apparent(t)[1:3, 1:se.segments]
        v_app_perp(t)[1:3, 1:se.segments]
        norm_v_app(t)[1:se.segments]
        half_drag_force(t)[1:3, 1:se.segments]
        total_force(t)[1:3, 1:se.segments+1]
    end

    # basic differential equations
    eqs2 = vcat([D(pos[:, i]) ~ vel[:, i] for i in axes(pos, 2)],
                [D(vel[:, i]) ~ acc[:, i] for i in axes(vel, 2)])
    # loop over all segments to calculate the spring and drag forces
    for i in 1:se.segments
        eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
               len[i]             ~ norm(segment[:, i]),
               unit_vector[:, i]  ~ -segment[:, i]/len[i],
               rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               c_spr[i]           ~ c_spring / (1+rel_compression_stiffness) 
                                     * (rel_compression_stiffness+(len[i] > l_spring)),
               spring_force[:, i] ~ (c_spr[i] * (len[i] - l_spring) 
                                     + damping * spring_vel[i]) * unit_vector[:, i],
               v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2,
               v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
               norm_v_app[i]      ~ norm(v_app_perp[:, i]),
               half_drag_force[:, i] ~ 0.25 * se.rho * se.cd_tether * norm_v_app[i] * (len[i]*se.d_tether/1000.0)
                                        * v_app_perp[:, i]]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # loop over all tether particles to apply the forces and calculate the accelerations
    for i in 1:(se.segments+1)
        eqs = []
        if i == se.segments+1
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] + half_drag_force[:, i-1])
            if isnothing(p2) || ! fix_p2
                push!(eqs, acc[:, i]         ~ se.g_earth .+ total_force[:, i] / (0.5 * m_tether_particle))
            else
                push!(eqs, acc[:, i]         ~ zeros(3))
            end
        elseif i == 1
            # spring_force[:, i] acts on particle i+1; particle i feels the reaction
            push!(eqs, total_force[:, i] ~ -spring_force[:, i] + half_drag_force[:, i])
            if isnothing(p1) || ! fix_p1
                push!(eqs, acc[:, i]     ~ se.g_earth .+ total_force[:, i] / (0.5 * m_tether_particle))
            else
                push!(eqs, acc[:, i]     ~ zeros(3))
            end
        else
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i] 
                                           + half_drag_force[:, i-1] + half_drag_force[:, i])
            push!(eqs, acc[:, i]         ~ se.g_earth .+ total_force[:, i] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # scalar equations
    eqs = [l_spring   ~ (se.l0 + se.v_ro*t) / se.segments,
    c_spring          ~ se.c_spring / l_spring,
    m_tether_particle ~ mass_per_meter * l_spring,
    damping           ~ se.damping  / l_spring]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  
        
    @named sys = System(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    tic()
    simple_sys = mtkcompile(sys) # first call: 15.0s second call: 0.02s 
    toc()
    simple_sys, sys, pos, vel, len, c_spr
end

"""
    simulate(se, simple_sys)

Simulate the tether model `simple_sys` (as returned by [`model`](@ref)) over the
duration `se.duration` using a fixed step size of 0.02s and the `FBDF` solver.

The solver is run twice; the first call triggers compilation, and only the
elapsed time of the second (timed) call is returned.

Returns a tuple `(sol, elapsed_time)` with the `ODESolution` and the elapsed
time in seconds of the second solve.
"""
function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.151
    ylim = (-1.2 * (se.l0 + se.v_ro*se.duration), 0.5)
    xlim = (-se.l0, se.l0)
    mkpath("video")
    z_max = 0.0
    # text position
    xy = (se.l0/4.2, z_max-7)
    start = time_ns()
    i = 1; j = 0
    for time in 0:dt:se.duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        plot2d(sol[pos][i], time; segments=se.segments, xlim, ylim, xy)
        if se.save
            ControlPlots.plt.savefig("video/"*"img-"*lpad(j, 4, "0"))
        end
        j += 1
        wait_until(start + 0.5 * time * 1e9)
    end
    if se.save
        println("Run the script ./bin/export_gif to create the gif file!")
    end
    nothing
end

function main(; p1=[0,0,0], p2=nothing, fix_p1=true, fix_p2=false)
    global sol, pos, vel, len, c_spr, simple_sys
    se = Settings3()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    simple_sys, sys, pos, vel, len, c_spr = model(se; p1, p2, fix_p1, fix_p2)
    sol, elapsed_time = simulate(se, simple_sys)
    if @isdefined __PC
        return sol, pos, vel, simple_sys
    end
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
    println("Equations: ", length(equations(sys)))
    println("Equations of simplified system: ", length(equations(simple_sys)))
    sol, pos, vel, simple_sys
end

main(p2=[-40,0,-47], fix_p2=false);

nothing

