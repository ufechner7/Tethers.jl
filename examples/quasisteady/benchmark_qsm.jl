using Pkg
# When running this script directly (not via `examples/menu3.jl`), activate the examples
# environment (one level up) unless it is already active.
if dirname(Pkg.project().path) != normpath(joinpath(@__DIR__, ".."))
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using BenchmarkTools, StaticArrays
using StaticArrays: MVector
using Tethers.QuasiSteady: StaticSettings, Tether, step!

const segments = 15

# Initial conditions, hardcoded. These are the values that `test/data/input_basic_test.mat`
# used to provide, so the timings below stay comparable; the .mat file itself is only a
# reference fixture for `test/test_qsm.jl` and has no business in a benchmark.
kite_pos      = MVector{3}(100.0, 100.0, 300.0)
kite_vel      = MVector{3}(0.0, 0.0, 0.0)
tether_length = 431.66247903554

se = StaticSettings(segments   = segments,
                    rho        = 1.225,
                    g_earth    = MVector{3}(0.0, 0.0, -9.8066),
                    cd_tether  = 1.2,
                    d_tether   = 29.71973504179974,     # [mm]
                    rho_tether = 970.0,                 # = 0.6729 kg/m / A, Dyneema
                    c_spring   = 8.047069220746364e7)   # E*A with E = 116 GPa, A = 693.7 mm²
te = Tether(se)
# state_vec is (elevation [rad], wind-frame azimuth [rad], ground tension [N]); the two
# angles are the MATLAB-convention pair (0.321750554, -0.306277369) already converted with
# `QuasiSteady.matlab_to_wind`.
te.state_vec .= (1.1302856641844843, -0.7853981636973135, 1.60941384e5)

step!(te, kite_pos, kite_vel; tether_length, prn=true)

@benchmark step!(te, kite_pos, kite_vel; tether_length)

# `step!` reuses te.tether_pos/te.wind_vel across calls, unlike `simulate_tether` called
# directly, which allocates a fresh (3, segments) matrix every time - see how many bytes
# remain once that allocation is out of the picture.
step!(te, kite_pos, kite_vel; tether_length)   # warm up / compile
println("Bytes allocated per step!: ", @allocated(step!(te, kite_pos, kite_vel; tether_length)))

#= On Ryzen 7950X, `simulate_tether` directly (before the `Tether`/`step!` port, so this
still includes the one (3, segments) matrix allocation that `step!` now avoids by reusing
te.tether_pos - re-benchmark after the port to get numbers for `step!` itself):
julia> include("examples/quasisteady/benchmark_qsm.jl")
Iterations: 22, retcode: Success, |res|: 4.263256414560601e-14
BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
 Range (min … max):  22.793 μs …  2.311 ms  ┊ GC (min … max): 0.00% … 97.66%
 Time  (median):     23.444 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   24.000 μs ± 31.892 μs  ┊ GC (mean ± σ):  1.85% ±  1.38%

        ▁▁▄▇▅▆█▄▃▄
  ▁▁▂▄▅▆███████████▇▆▅▄▃▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  22.8 μs         Histogram: frequency by time        25.9 μs <

 Memory estimate: 6.70 KiB, allocs estimate: 63.

Was 94 μs / 46.53 KiB / 1118 allocs with 36 iterations, before

  - the residual was made allocation free (0.35 μs per evaluation, `tether_shape`),
  - the tension was moved to a logarithmic scale with a bounded trust region radius,
  - `rho_t` was read from the .mat file as a mass per unit length rather than a density,
    which had made the tether 1442x too light.

The remaining allocations are from the trust-region solve itself (NonlinearSolve/
ForwardDiff internals), not from `tether_pos` - `step!` removes the latter but does not
touch the former.
 =#
