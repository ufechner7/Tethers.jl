using Pkg
if dirname(Pkg.project().path) != @__DIR__
    Pkg.activate(@__DIR__)
end
using Test, LinearAlgebra
using Tethers: run_python
include(joinpath(@__DIR__, "test_utils.jl"))

# shared CI runners (notably macOS aarch64) are noticeably slower and noisier
# than a local dev machine, so the timing budget below needs more headroom there
const ELAPSED_TIME_LIMIT = get(ENV, "CI", "false") == "true" ? 2.0 : 1.0

@testset "Tether_06c" begin
    # without callbacks
    include("../examples/Tether_06c.jl")
    set = deepcopy(Settings2())
    set.duration = 10.0
    set.callbacks = false
    simple_sys, pos, vel = model(set)
    sol, elapsed_time = simulate(set, simple_sys)
    @test elapsed_time < ELAPSED_TIME_LIMIT
    l_tether_theoretical = set.l0 + set.v_ro * set.duration
    @test l_tether(sol, pos) ≈ l_tether_theoretical rtol=2e-3
    events = Int64(round(length(sol.t)- set.duration/set.dt)-1)
    @test events == 0

    # with callbacks
    set.callbacks = true
    simple_sys, pos, vel = model(set)
    sol, elapsed_time = simulate(set, simple_sys)
    @test elapsed_time < ELAPSED_TIME_LIMIT
    l_tether_theoretical = set.l0 + set.v_ro * set.duration
    @test l_tether(sol, pos) ≈ l_tether_theoretical rtol=2e-3
    events = Int64(round(length(sol.t)- set.duration/set.dt)-1)
    @test events >= 4 # 8 events with Rodas5, 4 events with KenCarp4

    # compare the default (with callbacks) run against the Python implementation
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Tether_06c.jl only runs main() (which writes the CSV) when __BENCH__ is
        # false; runtests.jl sets it to true for this testset, so it must be reset.
        global __BENCH__ = false
        include(joinpath(pkg_dir, "examples", "Tether_06c.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
        # Python implementation, using the implicit solver IDA
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            run_python("Tether_06c")
        end
    end

    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_06c_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_06c_python.csv"))

    # Julia's continuous_events and Python's state events fire at slightly
    # different times for this hard-switch model, so the two CSVs are not
    # guaranteed to land on identical time grids; interpolate before comparing.
    pos_z_py_interp = interp_at(t_py, pos_z_py, t_jl)
    vel_z_py_interp = interp_at(t_py, vel_z_py, t_jl)

    @test pos_z_jl ≈ pos_z_py_interp rtol=1e-2 atol=1e-2
    @test vel_z_jl ≈ vel_z_py_interp rtol=5e-2 atol=5e-2
end
nothing
