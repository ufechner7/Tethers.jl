using Pkg
if dirname(Pkg.project().path) != @__DIR__
    Pkg.activate(@__DIR__)
end

using BenchmarkTools
using Tethers: analytic_force
using Tethers.TetherComponents: TetherSettings

"""
    extension_range(filename)

`(min, max)` of the relative extension, `l_tether / l_tether_unstretched - 1`, over every
row of `filename`, the same result file `examples/plot_compression.jl` validates
[`analytic_force`](@ref) against.
"""
function extension_range(filename=joinpath(dirname(@__DIR__), "data",
                                            "compression_force_vs_length.csv"))
    lines = filter(!isempty, strip.(readlines(filename)))
    cols = split(lines[1], ',')
    i_l0 = findfirst(==("l_tether_unstretched"), cols)
    i_l  = findfirst(==("l_tether"), cols)
    ext = map(lines[2:end]) do line
        v = split(line, ',')
        parse(Float64, v[i_l]) / parse(Float64, v[i_l0]) - 1
    end
    extrema(ext)
end

const ext_min, ext_max = extension_range()

const se            = TetherSettings()
const v_wind_perp   = 20.0     # [m/s]
const d_segment     = 4.0      # [mm]
const l_unstretched = 3.0      # [m]
const segments      = 6

const n = 10_000
const extensions = ext_min .+ (ext_max - ext_min) .* rand(n)
const l_segments = l_unstretched .* (1 .+ extensions)

function run_analytic_force(se, v_wind_perp, d_segment, l_unstretched, l_segments, segments)
    for l_segment in l_segments
        analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments)
    end
end

println("Benchmarking analytic_force over $n random extensions in " *
        "[$(round(100ext_min, digits=2))%, $(round(100ext_max, digits=2))%], the range " *
        "found in data/compression_force_vs_length.csv")
b = @benchmark run_analytic_force($se, $v_wind_perp, $d_segment, $l_unstretched,
                                   $l_segments, $segments)
display(b)
println()
println("Mean time per call: $(round(BenchmarkTools.mean(b).time / n, digits=1)) ns")

nothing
