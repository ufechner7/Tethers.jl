using Pkg
# When running this script directly (not via `examples/menu3.jl`), activate the examples
# environment (one level up) unless it is already active.
if dirname(Pkg.project().path) != normpath(joinpath(@__DIR__, ".."))
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using BenchmarkTools, StaticArrays
using Tethers.Quasistatic: Settings, simulate_tether

const segments = 15

# Initial conditions, hardcoded. These are the values that `test/data/input_basic_test.mat`
# used to provide, so the timings below stay comparable; the .mat file itself is only a
# reference fixture for `test/test_qsm.jl` and has no business in a benchmark.
#
# state_vec is (elevation [rad], wind-frame azimuth [rad], ground tension [N]); the two
# angles are the MATLAB-convention pair (0.321750554, -0.306277369) already converted with
# `Quasistatic.matlab_to_wind`.
state_vec     = MVector{3}(1.1302856641844843, -0.7853981636973135, 1.60941384e5)
kite_pos      = MVector{3}(100.0, 100.0, 300.0)
kite_vel      = MVector{3}(0.0, 0.0, 0.0)
wind_vel      = zeros(3, segments)
tether_length = 431.66247903554

settings = Settings(rho        = 1.225,
                    g_earth    = MVector{3}(0.0, 0.0, -9.8066),
                    cd_tether  = 1.2,
                    d_tether   = 29.71973504179974,     # [mm]
                    rho_tether = 0.6729014779417218,
                    c_spring   = 8.047069220746364e7)   # E*A with E = 116 GPa, A = 693.7 mm²

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
