using Test, LinearAlgebra
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_08" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation: a SteadyStateDiffEq solve for the initial tether shape,
        # followed by a time simulation with the implicit solver FBDF
        include(joinpath(pkg_dir, "src", "Tether_08.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
        # Python implementation: the steady state is found directly with
        # scipy.optimize.least_squares, the time simulation uses SciPy's BDF
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_08.jl"))
        end
    end

    POS = stack(sol[pos], dims=1) # (time, xyz, particle)

    # p1 is fixed (fix_p1=true, acc[:, 1] ~ 0, starts at rest at the origin), so it must
    # stay fixed there for the whole simulation
    POS1 = POS[:, :, 1]
    @test maximum(abs.(POS1)) < 1e-6

    # p2 is free (fix_p2=false), so it must have moved away from its initial position
    # under gravity and tether drag
    p2_0   = POS[1, :, end]
    p2_end = POS[end, :, end]
    @test norm(p2_end - p2_0) > 1.0

    # the simulation must produce finite results for all masses
    @test all(isfinite, POS)
    @test all(isfinite, stack(sol[vel], dims=1))

    # the steady-state solver must have produced a taut initial tether shape: consecutive
    # particles must all be roughly l0/segments apart at t=0, with no folded/zero-length segments
    POS0 = POS[1, :, :]
    seg_lengths = [norm(POS0[:, i+1] - POS0[:, i]) for i in 1:size(POS0, 2)-1]
    se_default = Settings3() # main() ran with an otherwise-default Settings3
    l_seg = se_default.l0 / se_default.segments
    @test all(l -> l > 0.9 * l_seg, seg_lengths)
    @test maximum(seg_lengths) - minimum(seg_lengths) < 0.1 * l_seg

    # the two implementations must agree on the trajectory of the free end point
    t_jl, pos_z_jl, vel_z_jl = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_08_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_08_python.csv"))

    @test length(t_jl) == length(t_py)
    @test t_jl ≈ t_py atol=1e-6
    @test pos_z_jl ≈ pos_z_py rtol=1e-3
    @test vel_z_jl ≈ vel_z_py rtol=1e-2 atol=1e-2
end
nothing
