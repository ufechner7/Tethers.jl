using Test, LinearAlgebra
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_07" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation, using ModelingToolkit and the implicit solver FBDF
        include(joinpath(pkg_dir, "src", "Tether_07.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
        # Python implementation, using the implicit solver IDA
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_07.jl"))
        end
    end

    POS = stack(sol[pos], dims=1) # (time, xyz, particle)

    # p1 is fixed (acc[:, 1] ~ 0, starts at rest at the origin), so it must stay
    # fixed there for the whole simulation
    @test maximum(abs.(POS[:, :, 1])) < 1e-6

    # the simulation must produce finite results for all masses
    @test all(isfinite, POS)
    @test all(isfinite, stack(sol[vel], dims=1))

    # the wind blows in +x direction, therefore the drag must push the tether downwind:
    # the free end has to end up further downwind than it would swing without drag
    @test POS[end, 1, end] > POS[1, 1, end]

    # the two implementations must agree on the trajectory of the free end point
    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_07_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_07_python.csv"))

    @test length(t_jl) == length(t_py)
    @test t_jl ≈ t_py atol=1e-6
    @test pos_z_jl ≈ pos_z_py rtol=1e-3
    @test vel_z_jl ≈ vel_z_py rtol=1e-2 atol=1e-2
end
nothing
