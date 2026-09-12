# How the cost of the quasi-steady model and of the dynamic mass-spring-damper model grow
# with the number of tether segments.
#
# The quasi-steady model replaces the time integration by a solve for three unknowns
# (elevation, azimuth and ground tension), whose residual integrates the tether shape
# segment by segment. The dynamic model carries 6*(segments+1) states through a stiff BDF
# integration.
#
# Both are measured per simulated second, which is the only comparison that means
# anything: the quasi-steady model has no time step of its own, so it is stepped at DT
# along a trajectory exactly as flying_circular.jl steps it, and the dynamic solve is
# divided by its simulated duration.
#
# Both models are given the same tether: same length, diameter, density and unit spring
# constant, hanging from the origin at the same elevation.
#
# Writes docs/images/qsm_vs_dynamic.png.
using Pkg
if dirname(Pkg.project().path) != normpath(joinpath(@__DIR__, ".."))
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using ModelingToolkit, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, SteadyStateDiffEq
using BenchmarkTools, LinearAlgebra, Printf, StaticArrays, Statistics
using ModelingToolkit: t_nounits as t
using ADTypes: AutoFiniteDiff
using StaticArrays: MVector
using Tethers.QuasiSteady: StaticSettings, Tether, init!, step!
using Tethers.TetherComponents: TetherSettings, set_diameter!, FixedEnd, FreeEnd, assemble_tether
import GLMakie

"number of segments the two models are compared at"
const SEGMENT_COUNTS = [4, 8, 16, 32]
"simulated duration of both runs [s]"
const DURATION = 2.0
"time step both models are sampled at [s]"
const DT = 0.02
"elevation of the tether's loose end [deg]"
const ELEVATION = 70.0
"unstretched tether length [m]"
const L_TETHER = 50.0

"""
    stats(trial)

Median, interquartile range and sample count of a `BenchmarkTools` trial, in seconds.
"""
function stats(trial)
    t = sort(trial.times) ./ 1e9
    (median = median(t), iqr = quantile(t, 0.75) - quantile(t, 0.25), n = length(t))
end

"""
    quasisteady_time(segments)

Wall clock per simulated second for the quasi-steady model with `segments` segments, as
`examples/quasisteady/flying_circular.jl` uses it: the kite walks a circular trajectory at
`DT`, and `step!` re-solves at every sample from the previous solution.

Timing a single `step!` in a loop would measure nothing: `step!` writes the converged
answer back into `te.state_vec`, so the second call onwards starts from the solution and
returns immediately. The kite has to actually move between calls.

Returns the `stats` of one revolution's solve loop, scaled to one simulated second.
"""
function quasisteady_time(segments; reps = 7)
    se = StaticSettings(; segments, elevation = ELEVATION, l_tether = L_TETHER)
    kite_distance = L_TETHER / (1 + se.slack)
    ts = 0:DT:DURATION
    # a small circle around the initial kite position, flown once per DURATION
    β, radius, ω = deg2rad(ELEVATION), 0.05 * kite_distance, 2π / DURATION
    traj = [MVector{3}(kite_distance * cos(β) + radius * sin(ω * time),
                       radius * cos(ω * time),
                       kite_distance * sin(β)) for time in ts]

    times = Float64[]
    for _ in 1:reps
        te = Tether(se)
        init!(te)
        elapsed = @elapsed for pos in traj
            step!(te, pos, te.kite_vel; tether_length = L_TETHER)
        end
        push!(times, elapsed)
    end
    sort!(times)
    (median = median(times) / DURATION,
     iqr = (quantile(times, 0.75) - quantile(times, 0.25)) / DURATION, n = length(times))
end

"""
    dynamic_time(segments)

Wall clock per simulated second for the dynamic model with `segments` segments, integrated
with `FBDF` and the analytic sparse Jacobian, plus the seconds spent building the model
(`mtkcompile` and the symbolic Jacobian). Returns `(stats, build)`, the `stats` scaled to
one simulated second.
"""
function dynamic_time(segments; seconds = 15)
    se = TetherSettings(; segments, l0 = L_TETHER, v_ro = 0.0, duration = DURATION)
    set_diameter!(se, se.d_tether)
    β = deg2rad(ELEVATION)
    p1 = zeros(3)
    p2 = [cos(β) * L_TETHER / 1.05, 0.0, sin(β) * L_TETHER / 1.05]
    POS0 = stack(p1 .+ (i - 1) / segments .* (p2 .- p1) for i in 1:segments+1)
    VEL0 = zeros(3, segments + 1)

    build = @elapsed begin
        simple_sys, = assemble_tether(se; end1 = FixedEnd(; name = :end1, pos0 = p1),
                                      end2 = FreeEnd(; name = :end2, se, pos0 = p2),
                                      POS0, VEL0)
        prob = ODEProblem(simple_sys, nothing, (0.0, se.duration); jac = true, sparse = true)
    end
    ts = 0:DT:se.duration
    trial = @benchmark solve($prob, FBDF(); dt = DT, abstol = 1e-6, reltol = 1e-6,
                             saveat = $ts) seconds=seconds
    st = stats(trial)
    (median = st.median / se.duration, iqr = st.iqr / se.duration, n = st.n), build
end

"""
    plot_scaling(segment_counts, qsm, dyn)

Plot the per-call cost of the quasi-steady model and the per-simulated-second cost of the
dynamic model against the number of segments, on logarithmic axes so that a power law in
the number of segments shows up as a straight line, and save it to
`docs/images/qsm_vs_dynamic.png`.
"""
function plot_scaling(segment_counts, qsm, dyn)
    fig = GLMakie.Figure(size = (800, 500))
    ax = GLMakie.Axis(fig[1, 1], xlabel = "tether segments", ylabel = "compute time [s]",
                      yscale = log10, xscale = log2, xticks = segment_counts,
                      title = "Quasi-steady vs. dynamic, $(Int(L_TETHER)) m tether, both per simulated second")
    for (series, label) in ((qsm, "quasi-steady, $(Int(1/DT)) step! calls"),
                            (dyn, "dynamic, FBDF with jac+sparse"))
        med = [x.median for x in series]
        GLMakie.scatterlines!(ax, segment_counts, med; label)
        # the interquartile range of the samples; usually narrower than the marker
        GLMakie.rangebars!(ax, segment_counts, [x.median - x.iqr/2 for x in series],
                           [x.median + x.iqr/2 for x in series])
    end
    GLMakie.axislegend(ax; position = :lt)
    path = joinpath(@__DIR__, "..", "..", "docs", "images", "qsm_vs_dynamic.png")
    GLMakie.save(normpath(path), fig)
    println("wrote ", normpath(path))
    fig
end

function main()
    qsm, dyn, builds = [], [], Float64[]
    @printf("%9s %8s %26s %26s %9s %10s\n", "segments", "states",
            "quasi-steady [ms/sim s]", "dynamic [ms/sim s]", "ratio", "build [s]")
    for segments in SEGMENT_COUNTS
        q = quasisteady_time(segments)
        d, b = dynamic_time(segments)
        push!(qsm, q); push!(dyn, d); push!(builds, b)
        @printf("%9d %8d %14.3f ±%-5.3f n=%-4d %14.2f ±%-5.2f n=%-4d %9.1f %10.2f\n",
                segments, 6 * (segments + 1),
                q.median * 1e3, q.iqr * 1e3, q.n,
                d.median * 1e3, d.iqr * 1e3, d.n,
                d.median / q.median, b)
    end
    plot_scaling(SEGMENT_COUNTS, qsm, dyn)
    qsm, dyn, builds
end

main()
