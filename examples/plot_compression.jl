# Plots and checks the result of the compression investigation of PlanCompression.md,
# reading it back from the CSV file that `test_compression.jl` writes.
#
# This is deliberately separate from the simulation: it needs nothing but GLMakie, so it
# loads in a second and lets the plots be reworked without re-running the 3 minute sweep.
# `test_compression.jl` includes this file at the end, so there is exactly one
# implementation of the figures. The analytical formula itself, `analytic_force`, lives in
# the package (`src/analytic_force.jl`), so that it can be reused outside the examples too.
#
#     include("examples/plot_compression.jl")    # replot data/compression_force_vs_length.csv
#
# `plot_lengths` gives the overview, one panel per wind speed and diameter, with the
# analytical formula of `analytic_force` dashed over the measurements. `plot_distance`
# drills into one operating condition and shows the individual segments and the two
# anchors. `report_sign` and `check_formula` are the two numerical checks.
# `import`, not `using`: the figures need a logarithmic force axis, which MakieControlPlots
# has not, but GLMakie and MakieControlPlots both export `plot`, and `using` both of them
# makes that name ambiguous in `Main` for every example included afterwards. Everything
# from Makie is therefore qualified below.
import GLMakie
import Tethers
using Tethers: display_if_interactive, analytic_force
using Tethers.TetherComponents: TetherSettings

"""
    Result

One row of the result file: an operating point of the compression sweep.

# Fields
- `v_wind`: wind speed in x direction [m/s]
- `d_tether`: tether diameter [mm]
- `l_unstretched`: unstretched tether length [m]
- `l_tether`: distance between the first and the last particle [m]
- `f_top`, `f_bot`: axial (z) force at the upper and the lower anchor [N], positive
  under tension
- `f_mean`: mean axial force of all segments [N]
- `f_seg`: axial force of each segment [N], positive under tension
"""
struct Result
    v_wind::Float64
    d_tether::Float64
    l_unstretched::Float64
    l_tether::Float64
    f_top::Float64
    f_bot::Float64
    f_mean::Float64
    f_seg::Vector{Float64}
end

"""
    compression(res)

Relative compression of `res` [-]: how much longer the unstretched tether is than the
distance between its end points. Negative values mean extension.
"""
compression(res::Result) = res.l_unstretched / res.l_tether - 1

"""
    extension(res)

Relative extension of `res` [-], the reciprocal view of [`compression`](@ref): how much
further apart the end points are than the unstretched tether is long.
"""
extension(res::Result) = res.l_tether / res.l_unstretched - 1

"""
    read_results(filename=joinpath("data", "compression_force_vs_length.csv"))

Read the result file written by `test_compression.jl` and return a vector of
[`Result`](@ref)s.

The columns are `v_wind`, `d_tether`, `l_tether_unstretched`, `l_tether`, `f_top`,
`f_bot`, `f_mean` and one `f_seg_i` per segment; the number of segments is taken from how
many columns there are, so a file written with a different `segments` still reads.
"""
function read_results(filename=joinpath("data", "compression_force_vs_length.csv"))
    isfile(filename) ||
        error("no result file $filename; run examples/test_compression.jl first")
    lines = filter(!isempty, strip.(readlines(filename)))
    length(lines) > 1 || error("$filename has a header but no data")
    cols = split(lines[1], ',')
    cols[1:7] == ["v_wind", "d_tether", "l_tether_unstretched", "l_tether",
                  "f_top", "f_bot", "f_mean"] ||
        error("unexpected columns in $filename: $(join(cols, ", "))")
    map(lines[2:end]) do line
        v = parse.(Float64, split(line, ','))
        length(v) == length(cols) ||
            error("$filename: expected $(length(cols)) values, got $(length(v)) in: $line")
        Result(v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8:end])
    end
end

"""
    analytic_force(se, res::Result)

[`analytic_force`](@ref) for the operating point `res`, so that it can be compared with its
measured `f_mean`. `se` only needs to hold `rho`, `cd_tether`, `d_tether` and `c_spring`; a
plain `TetherSettings()` reproduces the sweep's defaults, since `res.d_tether` is passed
through explicitly.
"""
Tethers.analytic_force(se, res::Result) =
    analytic_force(se; v_wind_perp=res.v_wind, d_segment=res.d_tether,
                     l_unstretched=res.l_unstretched, l_segment=res.l_tether,
                     segments=length(res.f_seg))

"""
    report_sign(results)

Check the claim of PlanCompression.md that the force never changes its sign, by printing
every operating point whose mean axial force is not tensile.

The drag bows a compressed tether out until its arc is longer than its unstretched length,
so the segments end up stretched even when the end points are closer together than the
tether is long; the force is expected to stay positive over the whole swept range.

Returns the vector of offending [`Result`](@ref)s, empty if the claim holds.
"""
function report_sign(results)
    bad = filter(res -> res.f_mean <= 0, results)
    if isempty(bad)
        println("The mean axial force is tensile at all $(length(results)) operating " *
                "points, minimum: $(round(minimum(res -> res.f_mean, results), digits=4)) N")
    else
        println("The mean axial force is NOT tensile at $(length(bad)) operating points:")
        for res in bad
            println("  v_wind: $(res.v_wind) m/s, d: $(res.d_tether) mm, " *
                    "l0: $(res.l_unstretched) m, l_tether: $(res.l_tether) m, " *
                    "f_mean: $(res.f_mean) N")
        end
    end
    bad
end

"""
    check_formula(results, se=TetherSettings())

Compare [`analytic_force`](@ref) with the measured mean axial force of every operating
point and print the median, 90th percentile and worst relative error.

Returns the vector of relative errors.
"""
function check_formula(results, se=TetherSettings())
    err = [analytic_force(se, res)/res.f_mean - 1 for res in results]
    a = sort(abs.(err))
    println("Analytical formula vs $(length(results)) measured points: " *
            "median $(round(100*a[cld(end,2)], digits=3))%, " *
            "p90 $(round(100*a[cld(9*end,10)], digits=3))%, " *
            "max $(round(100*a[end], digits=3))%")
    err
end

"""
    plot_lengths(results, se=TetherSettings(); min_force=1e-4)

The overview plot: the magnitude of the mean axial force over the relative compression, on
a logarithmic axis, as a grid of panels — one row per wind speed, one column per tether
diameter — with one curve per tether length, and the prediction of [`analytic_force`](@ref)
dashed in black over each of them. All panels share their axes, so the effect of the wind
is the shift down the rows and the effect of the diameter the shift across the columns.

The magnitude is plotted so that a sign change cannot break the logarithmic axis: a curve
that dives towards the clamp `min_force` [N] is a force that passes through zero, which is
exactly what the plot is meant to reveal ([`report_sign`](@ref) says whether that happened).
"""
function plot_lengths(results, se=TetherSettings(); min_force=1e-4)
    v_winds   = sort(unique(res.v_wind for res in results))
    d_tethers = sort(unique(res.d_tether for res in results))
    l0s       = sort(unique(res.l_unstretched for res in results))
    fig = GLMakie.Figure(size=(360*length(d_tethers), 260*length(v_winds)))
    axs = Matrix{GLMakie.Axis}(undef, length(v_winds), length(d_tethers))
    for (row, v_wind) in pairs(v_winds), (col, d_tether) in pairs(d_tethers)
        ax = GLMakie.Axis(fig[row, col]; yscale=log10,
                          ylabel=col == 1 ? "|mean axial force| [N]\nat $(v_wind) m/s" : "",
                          xlabel=row == length(v_winds) ? "relative compression [%]" : "",
                          title=row == 1 ? "$(d_tether) mm" : "")
        for l0 in l0s
            sel = sort(filter(res -> res.v_wind == v_wind && res.d_tether == d_tether &&
                                     res.l_unstretched == l0, results), by=compression)
            isempty(sel) && continue
            X = [100 * compression(res) for res in sel]
            # log10(0) is -Inf and would break the axis, so clamp the magnitude from below
            Y = [max(abs(res.f_mean), min_force) for res in sel]
            GLMakie.lines!(ax, X, Y; label="l0 = $l0 m")
            GLMakie.scatter!(ax, X, Y; markersize=6)
            # the analytical formula of step two, over the measured points
            GLMakie.lines!(ax, X, [max(abs(analytic_force(se, res)), min_force) for res in sel];
                           linestyle=:dash, color=:black, linewidth=1)
        end
        row == 1 && col == 1 && GLMakie.axislegend(ax; position=:rb)
        axs[row, col] = ax
    end
    GLMakie.linkaxes!(axs...)   # same scale everywhere, so the panels compare directly
    # ticks and labels only on the outer edge of the grid
    for row in eachindex(v_winds), col in eachindex(d_tethers)
        row < length(v_winds) && GLMakie.hidexdecorations!(axs[row, col]; grid=false)
        col > 1               && GLMakie.hideydecorations!(axs[row, col]; grid=false)
    end
    GLMakie.Label(fig[0, :], "Equilibrium force over tether length, wind speed and " *
                             "diameter, no gravity (dashed: analytical formula)";
                  fontsize=16)
    display_if_interactive(fig)
    fig
end

"""
    plot_distance(results; min_force=1e-4)

Drill-down into one operating condition: the axial force of every segment and the two
anchor forces over the relative extension of the tether, both on a logarithmic axis.

`results` must all share one unstretched length, diameter and wind speed, e.g. one slice of
[`read_results`](@ref); [`slice`](@ref) picks one:

    plot_distance(slice(results; l_unstretched=3.0, d_tether=4.0, v_wind=10.0))

[`plot_lengths`](@ref) shows only the mean axial force, so this is the plot to reach for
when the question is how the segments differ from each other, or how much of the anchor
force is drag rather than tension.

The force spans four decades over the swept range, so a linear axis shows nothing but the
point at 1% extension; the magnitude is plotted and clamped at `min_force` [N] for the same
reason as in [`plot_lengths`](@ref).
"""
function plot_distance(results; min_force=1e-4)
    isempty(results) && error("nothing to plot")
    results = sort(results, by=extension)
    n = length(first(results).f_seg)
    X = [100 * extension(res) for res in results]
    clamped(f) = [max(abs(f(res)), min_force) for res in results]
    fig = GLMakie.Figure()
    ax1 = GLMakie.Axis(fig[1, 1]; ylabel="|segment force| [N]", yscale=log10,
                       title="Equilibrium force of a $n-segment, " *
                             "$(first(results).l_unstretched) m x " *
                             "$(first(results).d_tether) mm tether at " *
                             "$(first(results).v_wind) m/s wind, no gravity")
    for i in 1:n
        GLMakie.lines!(ax1, X, clamped(res -> res.f_seg[i]); label="S$i")
    end
    GLMakie.axislegend(ax1; position=:rb)
    ax2 = GLMakie.Axis(fig[2, 1]; xlabel="relative extension [%], negative = compression",
                       ylabel="|anchor force| [N]", yscale=log10)
    GLMakie.lines!(ax2, X, clamped(res -> res.f_top); label="upper anchor")
    GLMakie.lines!(ax2, X, clamped(res -> res.f_bot); label="lower anchor")
    GLMakie.axislegend(ax2; position=:rb)
    GLMakie.linkxaxes!(ax1, ax2)
    GLMakie.hidexdecorations!(ax1; grid=false)
    display_if_interactive(fig)
    fig
end

"""
    slice(results; kwargs...)

The subset of `results` whose fields match every keyword, e.g.
`slice(results; l_unstretched=3.0, d_tether=4.0, v_wind=10.0)` for one curve of
[`plot_lengths`](@ref). Handy for feeding [`plot_distance`](@ref).
"""
slice(results; kwargs...) =
    filter(res -> all(getfield(res, k) == v for (k, v) in kwargs), results)

"""
    plot_results(filename=joinpath("data", "compression_force_vs_length.csv"))

Read the result file, run the two checks and show the overview plot. Returns the vector of
[`Result`](@ref)s, so that [`plot_distance`](@ref) and [`slice`](@ref) can be used on it.

Named `plot_results` and not `main` on purpose: `test_compression.jl` includes this file,
and its own `main` runs the sweep.
"""
function plot_results(filename=joinpath("data", "compression_force_vs_length.csv"))
    results = read_results(filename)
    println("Read $(length(results)) operating points from $filename")
    report_sign(results)
    check_formula(results)
    plot_lengths(results)
    results
end

results = plot_results();

nothing
