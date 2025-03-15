using Test
include("../src/Tether_quasistatic.jl")

if basename(pwd())=="test"
    cd("..")
end

@testset "residual_test" begin  
    state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings = get_initial_conditions("test/data/input_basic_test.mat")
    # Set up other parameters
    Ns = size(wind_vel, 2)
    buffers= [zeros(3, Ns), zeros(3, Ns), zeros(3, Ns), zeros(3, Ns), zeros(3, Ns)]
    res = zeros(3)
    # Pack in param tuple
    param = (kite_pos, kite_vel, wind_vel, tether_length, settings, buffers, true)
    # Call objective function
    Fobj, T0, pj, p0 = res!(res, state_vec, param)

    # Test type
    @test Fobj isa Vector
    @test T0 isa MVector{3,Float64}
    @test p0 isa Vector
    @test pj isa Matrix

    # Get reference values from .mat
    Fobj_ref, p0_ref, pj_ref, T0_ref = get_test_output("test/data/basic_test_results.mat")
    # Test values
    @test Fobj == Fobj_ref
    @test T0 == T0_ref
    @test pj ≈ pj_ref
    @test p0 == p0_ref
    nothing
end