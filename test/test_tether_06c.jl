using Test, LinearAlgebra
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_06c" begin
    # without callbacks
    include("../src/Tether_06c.jl")
    set = deepcopy(Settings2())
    set.duration = 10.0
    set.callbacks = false
    simple_sys, pos, vel = model(set)
    sol, elapsed_time = simulate(set, simple_sys)
    @test elapsed_time < 1.0
    l_tether_theoretical = set.l0 + set.v_ro * set.duration
    @test l_tether(sol, pos) ≈ l_tether_theoretical rtol=2e-3
    events = Int64(round(length(sol.t)- set.duration/set.dt)-1)
    @test events == 0

    # with callbacks
    set.callbacks = true
    simple_sys, pos, vel = model(set)
    sol, elapsed_time = simulate(set, simple_sys)
    @test elapsed_time < 1.0
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
        include(joinpath(pkg_dir, "src", "Tether_06c.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
        # Python implementation, using the implicit solver IDA
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_06c.jl"))
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
