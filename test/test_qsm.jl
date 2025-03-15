using Test
include("../src/Tether_quasistatic.jl")

if basename(pwd())=="test"
    cd("..")
end

"""
    get_test_output(filename)

Loads the output from the original MATLAB objective function for the tests

# Arguments
- filename: the filename of the mat file to read

# Returns
- res::MVector{3, Float64} difference between tether end and kite segment
- T0::MVector{3, Float64} force from the kite to the end of tether
- pj::(3, Ns) Matrix{Float64} x,y,z - coordinates of the Ns tether nodes
- p0::MVector{3, Float64}  x,y,z - coordinates of the kite-tether attachment
"""
function get_test_output(filename)
    vars        = matread(filename) 
    res        = MVector{3}(vec(get(vars,"Fobj", 0)))
    p0          = MVector{3}(vec(get(vars,"p0", 0)))
    pj          = get(vars,"pj", 0)
    T0          = MVector{3}(vec(get(vars,"T0", 0)))
    return res, p0, pj, T0
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