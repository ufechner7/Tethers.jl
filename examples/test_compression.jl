# Step one of PlanCompression.md: investigate the equilibrium force of a tether that is
# loaded by wind, but whose end points are closer together than its unstretched length,
# i.e. that is compressed.
#
# A compressed tether buckles: it bows out in the wind direction until the aerodynamic drag
# is balanced by the (much softer, see `rel_compression_stiffness`) compression force of the
# segments. The question this script prepares is whether that force keeps its sign, i.e.
# whether the tether always stays in compression as long as the wind is significant.
#
# `main()` runs the whole investigation: vertical tethers of 6 segments with the unstretched
# lengths 1, 3, 10 and 30 m and the diameters 2, 4, 6 and 8 mm, no gravity, both end points
# fixed, at 10, 20 and 30 m/s of horizontal wind, with `l_tether_unstretched / l_tether`
# swept from 0.99 (1% extension) to 1.10 (10% compression) on the non-uniform grid of
# `ratios_around_one` — 672 operating points in about 3 minutes, written to
# `data/compression_force_vs_length.csv`.
#
# The checks and the figures are not here: `examples/plot_compression.jl` reads the CSV back
# and does them, and this script includes it at the end. That way the plots and the
# analytical formula of step two can be reworked in a second, without re-running the sweep
# and without loading ModelingToolkit.
#
# Only the unstretched length is baked into the model; the anchor positions, the wind and
# everything that follows from the diameter are parameters of `FixedEnd` and `Tether`. A
# whole strain, wind and diameter sweep therefore runs on one compiled model, so the full
# run needs only 4 `mtkcompile` calls, one per length.
using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using ADTypes: AutoFiniteDiff
using Tethers.TetherComponents: TetherSettings, set_diameter!, Tether, FixedEnd
using Tethers: TetherComponents   # for the unexported TetherComponents.mass_per_meter

"""
    OperatingPoint

One measured equilibrium of the compression experiment.

# Fields
- `v_wind`: wind speed in x direction [m/s]
- `d_tether`: tether diameter [mm]
- `l_unstretched`: unstretched tether length [m]
- `l_tether`: distance between the first and the last particle [m]
- `f_top`, `f_bot`: axial (z) force at the upper and the lower anchor [N], positive
  under tension
- `f_axial`: axial force of each segment [N], positive under tension
- `pos`: the equilibrium shape, a `3 × (segments+1)` matrix
"""
struct OperatingPoint
    v_wind::Float64
    d_tether::Float64
    l_unstretched::Float64
    l_tether::Float64
    f_top::Float64
    f_bot::Float64
    f_axial::Vector{Float64}
    pos::Matrix{Float64}
end

"""
    f_mean(op)

Mean axial force of all segments of the operating point `op` [N], positive under tension.
"""
f_mean(op::OperatingPoint) = sum(op.f_axial) / length(op.f_axial)

"""
    compression(op)

Relative compression of the operating point `op` [-]: how much longer the unstretched
tether is than the distance between its end points. Negative values mean extension.
"""
compression(op::OperatingPoint) = op.l_unstretched / op.l_tether - 1

# Magnitude of the strain steps of both sweeps, geometric and dense towards zero strain,
# because that is how the force behaves: the first run, at 6 m and 10 m/s, measured 6146 N at 1%
# extension, 3073 N at 0.5% extension, 37 N at zero strain and 4 N at 0.5% compression, but
# then only creeps from 2.8 N at 1% compression down to 0.8 N at 10% compression. A uniform
# grid spends almost all of its points on that flat tail and resolves none of the four
# decades around zero strain, so the steps shrink towards zero instead.
# The results go to `data` and not to `output`, which is in .gitignore: these CSV files are
# the input of step two, so they have to survive and be diffable.
const DATA_DIR = "data"

const EXTENSION_STEPS   = [0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002]
const COMPRESSION_STEPS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]

"""
    ratios_around_one(; below=EXTENSION_STEPS, above=COMPRESSION_STEPS)

The default grid of `l_tether_unstretched / l_tether` ratios of [`sweep_lengths`](@ref):
`1 .- below`, then `1.0`, then `1 .+ above`, i.e. from 1% extension to 10% compression.
"""
ratios_around_one(; below=EXTENSION_STEPS, above=COMPRESSION_STEPS) =
    vcat(1 .- below, 1.0, 1 .+ above)

"""
    rels_around_zero(; ext=EXTENSION_STEPS, comp=COMPRESSION_STEPS)

The default grid of relative extensions of [`sweep_distance`](@ref): `ext`, then `0.0`,
then `-comp`, i.e. the same strains as [`ratios_around_one`](@ref), but expressed as the
relative change of the distance between the end points, which is what point 1 varies.
"""
rels_around_zero(; ext=EXTENSION_STEPS, comp=COMPRESSION_STEPS) = vcat(ext, 0.0, -comp)

"""
    compression_settings(; segments=6, l0=6.0, v_wind=10.0, d_tether=4.0)

Settings of the compression experiment: `segments` segments of the total unstretched length
`l0` [m] and the diameter `d_tether` [mm], a horizontal wind of `v_wind` [m/s] in x
direction, no gravity and no reel-out.
"""
function compression_settings(; segments=6, l0=6.0, v_wind=10.0, d_tether=4.0)
    se = TetherSettings()
    se.segments       = segments
    se.l0             = l0
    se.v_ro           = 0.0                    # a steady state exists only without reel-out
    se.g_earth        = zeros(3)               # gravity off, only wind and spring forces
    se.v_wind_tether  = [v_wind, 0.0, 0.0]
    se.duration       = 10.0
    set_diameter!(se, d_tether)                # scale the spring and damping constant with d
    se
end

"""
    bow_amplitude(se, dist)

Amplitude [m] of the initial half sine wave for the end point distance `dist` [m].

For a half sine wave of the amplitude `a` over the distance `L` the arc length is
`L + π²a²/(4L)`, so the slack `se.l0 - dist` of a compressed tether is taken up by
`a = 2/π * sqrt((se.l0 - dist) * dist)`. An extended tether gets a tiny bow instead, just
to break the symmetry of the straight initial shape.
"""
function bow_amplitude(se, dist)
    slack = se.l0 - dist
    slack > 0 ? 2/π * sqrt(slack * dist) : 1e-3 * se.l0
end

"""
    initial_shape(se; p1, p2, bow)

Initial positions of all particles between the end points `p1` and `p2`: a straight line,
displaced by a half sine wave of the amplitude `bow` [m] in the wind direction.

A perfectly straight, compressed tether is an equilibrium of the spring forces, so the
steady state solver needs a shape that is already bent to find the buckled solution. The
returned tuple `(POS0, VEL0)` has the `3 × (se.segments+1)` matrices expected by
[`Tether`](@ref); `VEL0` is all zeros.
"""
function initial_shape(se; p1, p2, bow)
    n = se.segments
    POS0 = zeros(3, n+1)
    VEL0 = zeros(3, n+1)
    Δ = (p2 - p1) / n
    for i in 1:n+1
        POS0[:, i] .= p1 + (i-1) * Δ
        POS0[1, i] += bow * sin(π * (i-1)/n)   # bow out in the direction of the wind
    end
    POS0, VEL0
end

"""
    build(se; p1, p2, POS0, VEL0)

Build one [`Tether`](@ref) with both end points held by a [`FixedEnd`](@ref) and return the
simplified system.

Everything except the anchor positions and the wind is baked into the model here, so a
compiled system may only be re-used for settings whose `segments`, `l0`, `d_tether`,
`c_spring`, `damping`, `rho`, `cd_tether` and `g_earth` are unchanged. `se.v_wind_tether`
and the anchor distance may vary afterwards; [`steady_state`](@ref) passes both to the
model as parameters.
"""
function build(se; p1, p2, POS0, VEL0)
    @named tether = Tether(; se, POS0, VEL0)
    @named end1   = FixedEnd(pos0=p1)
    @named end2   = FixedEnd(pos0=p2)
    eqs = [connect(end1.flange, tether.p1),
           connect(tether.p2, end2.flange)]
    @named sys = System(eqs, t; systems=[tether, end1, end2])
    mtkcompile(sys)
end

"""
    steady_state(se, simple_sys; p2, POS0)

Solve `simple_sys` for its steady state with the lower anchor at `p2` and the initial tether
shape `POS0`, and return the resulting shape, a `3 × (se.segments+1)` matrix. Both end
points must be fixed and `se.v_ro` must be zero.

Nothing that a sweep varies is baked into the model: the anchor position is the parameter
`end2.pos_fix` of [`FixedEnd`](@ref), the wind and the tether cross section are parameters
of [`Tether`](@ref), and the shape is the initial value of the states `tether.pos_in` and
`tether.vel_in`, so one `mtkcompile` serves a whole sweep.
"""
function steady_state(se, simple_sys; p2, POS0)
    n = se.segments
    op = [simple_sys.end2.pos_fix        => collect(p2),
          simple_sys.tether.v_wind       => collect(se.v_wind_tether),
          # the four quantities that follow from the diameter, kept in step with `se`
          simple_sys.tether.d_tether     => se.d_tether,
          simple_sys.tether.c_spring_unit => se.c_spring,
          simple_sys.tether.damping_unit => se.damping,
          simple_sys.tether.mass_per_m   => TetherComponents.mass_per_meter(se),
          simple_sys.tether.pos_in       => POS0[:, 2:n],  # the inner particles are states
          simple_sys.tether.vel_in       => zeros(3, n-1)]
    prob = SteadyStateProblem(ODEProblem(simple_sys, op, (0.0, se.duration)))
    sol = solve(prob, DynamicSS(KenCarp4(autodiff=AutoFiniteDiff())))
    SciMLBase.successful_retcode(sol) ||
        error("Steady state solver failed with return code $(sol.retcode)!")
    # read the shape from `sol.original`, the ODESolution of the integration that DynamicSS
    # ran, and not from `sol` itself; see the comment in examples/Tether_10.jl
    sol.original[simple_sys.tether.pos][end]
end

"""
    segment_forces(se, POS)

Axial force of each segment and the two anchor forces of the tether shape `POS`, a
`3 × (se.segments+1)` matrix of a steady state (i.e. all velocities are zero).

This evaluates the same equations as the [`Tether`](@ref) component, but in plain Julia and
without the damping and the relative velocity terms, which vanish in the steady state.

Returns `(F_axial, F_end1, F_end2)`:
- `F_axial`: axial force of each of the `se.segments` segments [N], positive under tension
- `F_end1`, `F_end2`: force the tether exerts on the first and the last anchor [N], a
  3D vector each, including the drag of the outermost half segment
"""
function segment_forces(se, POS)
    n = se.segments
    l_seg    = se.l0 / n                       # constant, se.v_ro is zero
    c_spring = se.c_spring / l_seg
    rcs      = se.rel_compression_stiffness
    d        = se.d_tether / 1000.0
    F_axial     = zeros(n)
    F_spring    = zeros(3, n)
    F_drag_half = zeros(3, n)
    for i in 1:n
        segment = POS[:, i+1] - POS[:, i]
        len     = norm(segment)
        uv      = -segment / len               # points from particle i+1 to particle i
        # the spring is much softer under compression than under tension
        c_spr        = c_spring / (1 + rcs) * (rcs + (len > l_seg))
        F_axial[i]   = c_spr * (len - l_seg)
        F_spring[:, i] = F_axial[i] * uv
        # in the steady state the particles do not move, so v_apparent is the wind
        v_app_perp = se.v_wind_tether - (se.v_wind_tether ⋅ uv) * uv
        F_drag_half[:, i] = 0.25 * se.rho * se.cd_tether * norm(v_app_perp) * (len * d) * v_app_perp
    end
    F_end1 = -F_spring[:, 1] + F_drag_half[:, 1]
    F_end2 =  F_spring[:, n] + F_drag_half[:, n]
    F_axial, F_end1, F_end2
end

"""
    equilibrium(se, simple_sys, l_tether)

Measure one operating point of the compiled model `simple_sys`: move its lower anchor
`l_tether` [m] below the upper one, solve for the steady state and evaluate the forces.

`simple_sys` must have been built by [`build`](@ref) from settings that agree with `se` in
everything it bakes in; the anchor distance and `se.v_wind_tether` are parameters and may
differ.

Returns an [`OperatingPoint`](@ref).
"""
function equilibrium(se, simple_sys, l_tether)
    p1 = [0.0, 0.0, 0.0]
    p2 = [0.0, 0.0, -l_tether]
    POS0, _ = initial_shape(se; p1, p2, bow=bow_amplitude(se, l_tether))
    POS = steady_state(se, simple_sys; p2, POS0)
    F_axial, F_end1, F_end2 = segment_forces(se, POS)
    # tension positive: the tether pulls the upper anchor down and the lower one up
    OperatingPoint(se.v_wind_tether[1], se.d_tether, se.l0, l_tether,
                   -F_end1[3], F_end2[3], F_axial, POS)
end

"""
    report(op)

Print one line describing the operating point `op`.
"""
function report(op::OperatingPoint)
    println("v_wind: $(rpad(round(op.v_wind, digits=1), 4)) m/s, " *
            "d: $(rpad(round(op.d_tether, digits=1), 4)) mm, " *
            "l0: $(rpad(round(op.l_unstretched, digits=3), 6)) m, " *
            "l_tether: $(rpad(round(op.l_tether, digits=3), 6)) m, " *
            "compression: $(rpad(round(100*compression(op), digits=2), 6))%, " *
            "mean axial force: $(rpad(round(f_mean(op), digits=4), 12)) N, " *
            "bow: $(round(maximum(op.pos[1, :]), digits=4)) m")
end

"""
    save_csv(ops; filename)

Write the operating points `ops` to a CSV file with the columns `v_wind`, `d_tether`,
`l_tether_unstretched`, `l_tether`, the anchor forces `f_top` and `f_bot`, the mean axial
force `f_mean` and the axial force `f_seg_i` of every segment. Returns `filename`.
"""
function save_csv(ops; filename)
    isempty(ops) && error("nothing to save")
    n = length(first(ops).f_axial)
    mkpath(dirname(filename))
    open(filename, "w") do io
        println(io, "v_wind,d_tether,l_tether_unstretched,l_tether,f_top,f_bot,f_mean," *
                    join(["f_seg_$i" for i in 1:n], ","))
        for op in ops
            println(io, join([op.v_wind, op.d_tether, op.l_unstretched, op.l_tether,
                              op.f_top, op.f_bot, f_mean(op), op.f_axial...], ","))
        end
    end
    filename
end

"""
    sweep_distance(se; rels=rels_around_zero(), v_winds=[10.0], d_tethers=[4.0])

The inner sweep of one compiled model: keep the unstretched length `se.l0` and set the
distance between the end points to `(1 + rel) * se.l0`, for every diameter in `d_tethers`
[mm], every wind speed in `v_winds` [m/s] and every `rel` in `rels`; a positive `rel`
extends the tether, a negative one compresses it. See [`rels_around_zero`](@ref) for the
default, non-uniform grid.

The whole sweep runs on a single compiled model: the unstretched length is fixed by `se`,
and everything else that changes between operating points is a parameter — the position of
the lower anchor (`end2.pos_fix`), the wind (`tether.v_wind`) and the four quantities that
follow from the diameter. `se` is mutated to carry the current diameter and wind, which
after [`build`](@ref) only feed those parameters and the post-processing of
[`segment_forces`](@ref), not the model itself.

An operating point whose steady state solver fails is reported and skipped, so that one bad
point does not abort the whole sweep.

Returns the vector of [`OperatingPoint`](@ref)s.
"""
function sweep_distance(se; rels=rels_around_zero(), v_winds=[10.0], d_tethers=[4.0])
    p1 = [0.0, 0.0, 0.0]
    p2 = [0.0, 0.0, -se.l0]
    POS0, VEL0 = initial_shape(se; p1, p2, bow=1e-3*se.l0)
    simple_sys = build(se; p1, p2, POS0, VEL0)   # the only mtkcompile of this sweep
    ops = OperatingPoint[]
    for d_tether in d_tethers
        set_diameter!(se, d_tether)   # also scales se.c_spring and se.damping
        for v_wind in v_winds
            se.v_wind_tether = [v_wind, 0.0, 0.0]
            for rel in rels
                try
                    op = equilibrium(se, simple_sys, (1 + rel) * se.l0)
                    report(op)
                    push!(ops, op)
                catch e
                    @warn "no steady state for d=$d_tether mm, v_wind=$v_wind m/s, rel=$rel" exception=e
                end
            end
        end
    end
    ops
end

"""
    sweep_lengths(; l_unstretched=[1,3,10,30], ratios=ratios_around_one(),
                    v_winds=[10,20,30], d_tethers=[2,4,6,8], segments=6)

For every unstretched tether length `l0` in `l_unstretched`, every diameter in `d_tethers`
[mm] and every wind speed in `v_winds` [m/s], vary the ratio
`l_tether_unstretched / l_tether` over `ratios`, i.e. from extension (ratio < 1) to
compression (ratio > 1); see [`ratios_around_one`](@ref) for the default, non-uniform grid.
The tether always has `segments` segments.

The wind and the diameter are parameters of the model, so all of them share the compiled
system of their length; only the length itself costs an `mtkcompile`.

It is the unstretched length that is held at 1, 3, 10 and 30 m here, and the distance
`l_tether = l0 / ratio` between the end points that varies: the unstretched length is baked
into the model, the anchor position is not, so this way each length needs one `mtkcompile`
instead of one per operating point. The swept ratios are exactly the same either way, and
both lengths are written to the CSV file.

Returns the vector of [`OperatingPoint`](@ref)s of all lengths, in the order they were
measured.
"""
function sweep_lengths(; l_unstretched=[1.0, 3.0, 10.0, 30.0], ratios=ratios_around_one(),
                         v_winds=[10.0, 20.0, 30.0], d_tethers=[2.0, 4.0, 6.0, 8.0],
                         segments=6)
    ops = OperatingPoint[]
    for l0 in l_unstretched
        println("--- l_tether_unstretched = $l0 m ---")
        se = compression_settings(; segments, l0, v_wind=first(v_winds),
                                    d_tether=first(d_tethers))
        # l_tether = l0 / ratio, so the relative extension of the distance is 1/ratio - 1
        append!(ops, sweep_distance(se; rels=1 ./ ratios .- 1, v_winds, d_tethers))
    end
    ops
end

"""
    main(; l_unstretched=[1,3,10,30], ratios=ratios_around_one(), v_winds=[10,20,30],
           d_tethers=[2,4,6,8], segments=6)

Sweep the strain for each of the unstretched tether lengths `l_unstretched`, each of the
diameters `d_tethers` and each of the wind speeds `v_winds`, and write the operating points
to `data/compression_force_vs_length.csv`.

Plotting and checking that file is the job of `examples/plot_compression.jl`.

Returns the vector of [`OperatingPoint`](@ref)s.
"""
function main(; l_unstretched=[1.0, 3.0, 10.0, 30.0], ratios=ratios_around_one(),
                v_winds=[10.0, 20.0, 30.0], d_tethers=[2.0, 4.0, 6.0, 8.0], segments=6)
    t0 = time_ns()
    ops = sweep_lengths(; l_unstretched, ratios, v_winds, d_tethers, segments)
    println("Elapsed time: $(round((time_ns()-t0)/1e9, digits=1)) s for $(length(ops)) operating points")
    println("Result written to: ",
            save_csv(ops; filename=joinpath(DATA_DIR, "compression_force_vs_length.csv")))
    ops
end

ops = main();
# the checks and the figures live in their own script, which reads the CSV back, so that
# they can be reworked without re-running the sweep
include("plot_compression.jl")

nothing
