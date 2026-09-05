using Test
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_02" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation, using ModelingToolkit and the implicit solver Rodas5
        include(joinpath(pkg_dir, "src", "Tether_02.jl"))
        sleep(1)
        MakieControlPlots.close("all")
        # Python implementation, using the implicit solver IDA
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_02.jl"))
        end
    end

    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_02_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_02_python.csv"))

    @test length(t_jl) == length(t_py)
    @test t_jl ≈ t_py atol=1e-6
    @test pos_z_jl ≈ pos_z_py rtol=1e-3
    # vel_z crosses zero, so a pure rtol comparison is overly strict there;
    # back it with an atol to account for small absolute solver differences
    @test vel_z_jl ≈ vel_z_py rtol=1e-2 atol=1e-2
end
nothing
