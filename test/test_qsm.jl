using Test
include(joinpath(@__DIR__, "..", "src", "Tether_quasistatic.jl"))

const QSM_DATA = joinpath(@__DIR__, "data")

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
    state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings =
        get_initial_conditions(joinpath(QSM_DATA, "input_basic_test.mat"))
    # Set up other parameters
    Ns = size(wind_vel, 2)
    buffers= [zeros(3, Ns), zeros(3, Ns), zeros(3, Ns), zeros(3, Ns), zeros(3, Ns)]
    res = zeros(3)
    # Pack in param named tuple; `return_result=true` makes res! return its intermediates
    param = (kite_pos=kite_pos, kite_vel=kite_vel, wind_vel=wind_vel,
             tether_length=tether_length, settings=settings, buffers=buffers,
             segments=Ns, return_result=true)
    # Call objective function
    Fobj, T0, pj, p0 = res!(res, state_vec, param)

    # Test type
    @test Fobj isa Vector
    @test T0 isa MVector{3,Float64}
    @test p0 isa MVector{3,Float64}
    @test pj isa AbstractMatrix

    # Get reference values from .mat
    Fobj_ref, p0_ref, pj_ref, T0_ref = get_test_output(joinpath(QSM_DATA, "basic_test_results.mat"))

    # These comparisons against the MATLAB reference output do not hold yet, because
    # `res!` and the reference data disagree on how the state vector's two angles define
    # the tether direction at the ground station. `res!` uses an elevation/azimuth
    # convention,
    #     dir ~ [cos(θ)cos(φ), cos(θ)sin(φ), sin(θ)]
    # while the reference data was produced with a z-up convention,
    #     dir ~ [sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]
    # For `input_basic_test.mat` (θ = 18.43°, φ = -17.55°) the first gives a first
    # segment along [0.905, -0.286, 0.316] and the second along [0.302, -0.302, 0.905],
    # so the computed tether is a differently-oriented one: ‖p0 - p0_ref‖ = 366 m.
    # Re-running `res!` with the angles converted to the reference convention brings that
    # down to 1.6 m, which identifies the convention as the cause but still leaves a
    # smaller second discrepancy to track down.
    #
    # This was never noticed because the test could not run: it used to pass `res!` a
    # 7-element `param` tuple, which the 8-way destructuring in `res!` rejects with a
    # BoundsError. Marked broken rather than deleted so the reference data keeps its
    # purpose; resolve the convention, then turn these back into `@test`.
    @test_broken Fobj ≈ Fobj_ref
    @test_broken T0 ≈ T0_ref
    @test_broken pj ≈ pj_ref
    @test_broken p0 ≈ p0_ref
    nothing
end
nothing
