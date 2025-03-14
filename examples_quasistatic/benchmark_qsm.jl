using BenchmarkTools

include("../src/Tether_quasistatic.jl")

# Read the initial conditions from a .mat file
state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")

simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings; prn=true)
@benchmark simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

#= julia> include("examples_quasistatic/benchmark_qsm.jl")
BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
 Range (min … max):  113.274 μs …   9.631 ms  ┊ GC (min … max): 0.00% … 96.95%
 Time  (median):     116.019 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   119.865 μs ± 141.114 μs  ┊ GC (mean ± σ):  2.24% ±  2.13%

       ▄█▇▄▁                                                     
  ▁▁▂▄██████▅▄▃▃▃▃▃▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  113 μs           Histogram: frequency by time          135 μs <

 Memory estimate: 46.34 KiB, allocs estimate: 1116. 
 On Uwes Laptop on battery, fast mode.
 =#






