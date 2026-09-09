using Pkg
# When running this script directly (not via `examples/menu3.jl`), activate the examples
# environment (one level up) unless it is already active.
if dirname(Pkg.project().path) != normpath(joinpath(@__DIR__, ".."))
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using BenchmarkTools
using Tethers.Quasistatic: get_initial_conditions, simulate_tether

const segments = 15

# Read the initial conditions from a .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions(joinpath(@__DIR__, "..", "..", "test", "data", "input_basic_test.mat"))
simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings; prn=true)

@benchmark simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

#= On Ryzen 7950X 
julia> include("examples/quasistatic/benchmark_qsm.jl")
Iterations: 36
BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
 Range (min … max):   93.927 μs …   8.020 ms  ┊ GC (min … max): 0.00% … 98.07%
 Time  (median):     102.166 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   106.874 μs ± 117.790 μs  ┊ GC (mean ± σ):  2.18% ±  2.32%

     █▇▃▁                                                        
  ▂▃█████▇▇▆▇▆▆▆▆▆▅▅▅▅▅▅▅▅▅▅▅▅▄▅▄▄▄▄▄▄▃▃▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂ ▄
  93.9 μs          Histogram: frequency by time          133 μs <

 Memory estimate: 46.53 KiB, allocs estimate: 1118. 
 =#






