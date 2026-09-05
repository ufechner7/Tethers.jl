using Test
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_03" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation, using ModelingToolkit and the implicit solver Rodas5
        include(joinpath(pkg_dir, "src", "Tether_03.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
        # Python implementation, using the implicit solver IDA
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_03.jl"))
        end
    end

    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_03_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_03_python.csv"))

    @test length(t_jl) == length(t_py)
    @test t_jl ≈ t_py atol=1e-6
    @test pos_z_jl ≈ pos_z_py rtol=1e-3
    # vel_z is more sensitive than pos_z to how each solver steps across the
    # spring's on/off discontinuity, so it needs a looser, atol-backed tolerance
    @test vel_z_jl ≈ vel_z_py rtol=1e-2 atol=1e-2
end
nothing
