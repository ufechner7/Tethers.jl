using Test
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_03b" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation, using ModelingToolkit and the implicit solver Rodas5, with callback
        include(joinpath(pkg_dir, "src", "Tether_03b.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
        # Python implementation, using the implicit solver IDA, with callback
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_03b.jl"))
        end
    end

    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_03b_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_03b_python.csv"))

    # Python/Assimulo inserts extra rows at each event crossing, so the two CSVs
    # are not saved on the same time grid; interpolate the Python series onto
    # the Julia time grid before comparing.
    pos_z_py_interp = interp_at(t_py, pos_z_py, t_jl)
    vel_z_py_interp = interp_at(t_py, vel_z_py, t_jl)

    @test pos_z_jl ≈ pos_z_py_interp rtol=1e-3 atol=1e-3
    # vel_z is more sensitive than pos_z to how each solver steps across the
    # spring's on/off discontinuity, so it needs a looser, atol-backed tolerance
    @test vel_z_jl ≈ vel_z_py_interp rtol=1e-2 atol=5e-2
end
nothing
