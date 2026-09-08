using Test
using Tethers: run_python
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_01" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation, using ModelingToolkit and the implicit solver FBDF
        include(joinpath(pkg_dir, "examples", "Tether_01.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
        # Python implementation, using the implicit solver RADAU
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            run_python("Tether_01")
        end
    end

    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_01_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_01_python.csv"))

    @test length(t_jl) == length(t_py)
    @test t_jl ≈ t_py atol=1e-6
    @test pos_z_jl ≈ pos_z_py rtol=1e-3
    @test vel_z_jl ≈ vel_z_py rtol=1e-3
end
nothing
