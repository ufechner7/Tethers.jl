using Test, LinearAlgebra
include(joinpath(@__DIR__, "test_utils.jl"))

@testset "Tether_10" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # the composed model: one Tether component, connected to a FixedEnd and a FreeEnd
        include(joinpath(pkg_dir, "src", "Tether_10.jl"))
        sleep(1)
        Base.invokelatest() do
            MakieControlPlots.close("all")
        end
    end

    POS = stack(sol[pos], dims=1) # (time, xyz, particle)
    se  = Base.invokelatest(TetherSettings)
    Base.invokelatest(set_diameter!, se, se.d_tether)

    # p1 is held by a FixedEnd at the origin, so it must stay there for the whole simulation
    @test maximum(abs.(POS[:, :, 1])) < 1e-6

    # p2 is held by a FreeEnd, so it must have moved away from its initial position
    @test norm(POS[end, :, end] - POS[1, :, end]) > 1.0

    # the simulation must produce finite results for all particles
    @test all(isfinite, POS)
    @test all(isfinite, stack(sol[vel], dims=1))

    # the steady-state solver must have produced a taut initial tether shape: consecutive
    # particles must all be roughly l0/segments apart at t=0
    POS0 = POS[1, :, :]
    seg_lengths = [norm(POS0[:, i+1] - POS0[:, i]) for i in 1:size(POS0, 2)-1]
    l_seg = se.l0 / se.segments
    @test all(l -> l > 0.9 * l_seg, seg_lengths)
    @test maximum(seg_lengths) - minimum(seg_lengths) < 0.1 * l_seg

    # the component must reproduce the monolithic model of example 8, which writes its
    # result to the same kind of csv file (only if that file is available)
    file_08 = joinpath(pkg_dir, "output", "Tether_08_julia.csv")
    if isfile(file_08)
        t_08, pos_z_08, vel_z_08 = read_pos_vel_csv(file_08)
        t_10, pos_z_10, vel_z_10 = read_pos_vel_csv(joinpath(pkg_dir, "output", "Tether_10_julia.csv"))
        @test length(t_08) == length(t_10)
        @test t_08 ≈ t_10 atol=1e-6
        @test pos_z_08 ≈ pos_z_10 rtol=1e-4
        @test vel_z_08 ≈ vel_z_10 rtol=1e-3 atol=1e-3
    end

    # Two tethers of half the length, joined by a point mass, are physically the same
    # system as one tether of the full length: same segment length, stiffness and damping,
    # and the knot has exactly the mass of an inner particle. So composing the same
    # component twice must give the same trajectory.
    simple_sys2, sol2 = cd(pkg_dir) do
        ssys2, _, _, _ = Base.invokelatest(model2, se; p1=[0,0,0], p2=[-40,0,-47])
        s2, _ = Base.invokelatest(simulate, se, ssys2)
        ssys2, s2
    end
    POS2 = stack([hcat(a, b[:, 2:end]) for (a, b) in
                  zip(sol2[simple_sys2.tether1.pos], sol2[simple_sys2.tether2.pos])], dims=1)
    @test size(POS2) == size(POS)
    @test maximum(abs.(POS2 - POS)) < 1e-3
end
nothing
