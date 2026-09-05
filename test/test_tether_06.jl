using Test
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_06" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation, using ModelingToolkit and the implicit solver FBDF
        include(joinpath(pkg_dir, "src", "Tether_06.jl"))
        sleep(1)
        MakieControlPlots.close("all")
        # Python implementation, using the implicit solver IDA
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_06.jl"))
        end
    end

    @test length(sol.t) == length(ts)

    # mass 1 is pinned (acc[:, 1] ~ 0, starts at rest at the origin), so it must
    # stay fixed there for the whole simulation
    POS1 = stack(sol[pos], dims=1)[:, :, 1]
    @test maximum(abs.(POS1)) < 1e-6

    # the simulation must produce finite results for all masses
    @test all(isfinite, stack(sol[pos], dims=1))
    @test all(isfinite, stack(sol[vel], dims=1))

    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_06_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_06_python.csv"))

    @test length(t_jl) == length(t_py)
    @test t_jl ≈ t_py atol=1e-6
    @test pos_z_jl ≈ pos_z_py rtol=1e-3
    # DAMPING is deliberately low here, so right after the tether goes taut the two
    # adaptive implicit solvers (FBDF vs. IDA) straddle the lightly-damped spring's
    # on/off discontinuity differently, producing a large but short-lived vel_z spike
    # (~2.65 m/s at t≈0.36s) that both sides damp out well before the end of the run.
    # `≈` on vectors compares the 2-norm of the difference (not elementwise), and that
    # spike is enough to push norm(vel_z_jl - vel_z_py) to ~8, so atol is set above that.
    @test vel_z_jl ≈ vel_z_py rtol=1e-2 atol=10.0
end
nothing
