using Test, MAT, StaticArrays, LinearAlgebra
using Tethers.QuasiSteady: get_initial_conditions
import Tethers.QuasiSteady as QSM  # for `res!`, which is internal and not exported

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
    Fobj, T0, pj, p0 = QSM.res!(res, state_vec, param)

    # Test type
    @test Fobj isa Vector
    @test T0 isa MVector{3,Float64}
    @test p0 isa MVector{3,Float64}
    @test pj isa AbstractMatrix

    # Get reference values from .mat
    Fobj_ref, p0_ref, pj_ref, T0_ref = get_test_output(joinpath(QSM_DATA, "basic_test_results.mat"))

    # `res!` and the MATLAB reference data used to disagree on how the state vector's two
    # angles define the tether direction at the ground station: `res!` uses an
    # elevation/azimuth convention, `dir ~ [cos(β)cos(φ), cos(β)sin(φ), sin(β)]`, while the
    # `.mat` fixtures store `stateVec` in a z-up convention,
    # `dir ~ [sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]`. `get_initial_conditions` now converts
    # `stateVec` from the MATLAB convention to the elevation/azimuth convention via
    # `matlab_to_wind` on load, so `state_vec` here is already directly comparable.
    # See docs/quasisteady.md, "Resolved: the angle convention", for the derivation.
    #
    # The vertical discrepancy that used to remain (~1.9 % on `T0`, ~0.4 % on `p0`) was
    # gravity after all: the `.mat` files store `T.rho_t` as a mass per unit length
    # [kg/m], while `Settings.rho_tether` is a density [kg/m^3] that the model multiplies
    # by the cross section itself, so the tether came out 1/A = 1442 times too light.
    # `get_initial_conditions` now divides by `A` on load, which yields 970.0 kg/m^3 -
    # Dyneema. With that, `T0 - Tn*dir` is `16*Ls*g*rho_t` = 2848.49 N against the
    # reference's 2848.49 N, and `T0`'s x and y already agreed to 5e-5 N out of 48525 N.
    #
    # Measured max relative error against the reference is now ~1e-16 to ~1e-15 for all
    # four quantities (double-precision roundoff), so `rtol` tightens from `1e-6` to
    # `1e-9`, leaving ~6 orders of magnitude of margin for BLAS/platform differences.
    @test Fobj ≈ Fobj_ref rtol=1e-9
    @test T0 ≈ T0_ref rtol=1e-9
    @test pj ≈ pj_ref rtol=1e-9
    @test p0 ≈ p0_ref rtol=1e-9
    nothing
end

@testset "Tether_init_step_equivalence" begin
    # `Tether`/`init!`/`step!` must reproduce `init_quasisteady`/`simulate_tether` bit for
    # bit, for a `StaticSettings` whose elevation/azimuth/l_tether/slack correspond to the
    # kite_pos/tether_length passed to the old API. `init!` derives its initial kite_pos
    # from exactly this formula, see `src/Tether_quasisteady.jl`.
    segments  = 12
    elev_deg  = 55.0
    azim_deg  = 15.0
    l_tether  = 200.0
    slack     = 0.05

    β0, φ0 = deg2rad(elev_deg), deg2rad(azim_deg)
    kite_dist = l_tether / (1 + slack)
    kite_pos0 = MVector{3}(kite_dist*cos(β0)*cos(φ0), kite_dist*cos(β0)*sin(φ0), kite_dist*sin(β0))
    kite_vel0 = MVector{3}(0.0, 0.0, 0.0)
    wind_vel0 = zeros(3, segments)

    # Old API
    state_vec_g, kite_pos_g, kite_vel_g, wind_vel_g, tether_length_g, settings_g =
        QSM.init_quasisteady(kite_pos0, l_tether; kite_vel=kite_vel0, segments, wind_vel=wind_vel0)
    state_vec_old, tether_pos_old, force_gnd_old, force_kite_old, p0_old =
        QSM.simulate_tether(state_vec_g, kite_pos_g, kite_vel_g, wind_vel_g, tether_length_g, settings_g)

    # New API
    se = QSM.StaticSettings(; segments, elevation=elev_deg, azimuth=azim_deg, l_tether, slack)
    te = QSM.Tether(se)
    QSM.init!(te)

    @test te.state_vec  ≈ state_vec_old  rtol=1e-12
    @test te.kite_pos    == kite_pos0
    @test te.tether_pos ≈ tether_pos_old rtol=1e-12
    @test te.force_gnd  ≈ force_gnd_old  rtol=1e-12
    @test te.force_kite ≈ force_kite_old rtol=1e-12
    @test te.p0         ≈ p0_old         rtol=1e-12
    @test QSM.elevation(te) == te.state_vec[1]
    @test QSM.azimuth(te)   == te.state_vec[2]
    @test QSM.tension(te)   == te.state_vec[3]

    # A subsequent step! must reproduce a subsequent simulate_tether call from the same
    # solved state.
    kite_pos_new = MVector{3}(kite_pos0 .* 1.02)
    kite_vel_new = MVector{3}(1.0, -0.5, 0.2)
    tether_length_new = 1.05 * norm(kite_pos_new)

    state_vec_old2, tether_pos_old2, force_gnd_old2, force_kite_old2, p0_old2 =
        QSM.simulate_tether(state_vec_old, kite_pos_new, kite_vel_new, wind_vel0, tether_length_new,
                            settings_g)
    QSM.step!(te, kite_pos_new, kite_vel_new; tether_length=tether_length_new)

    @test te.state_vec  ≈ state_vec_old2  rtol=1e-12
    @test te.tether_pos ≈ tether_pos_old2 rtol=1e-12
    @test te.force_gnd  ≈ force_gnd_old2  rtol=1e-12
    @test te.force_kite ≈ force_kite_old2 rtol=1e-12
    @test te.p0         ≈ p0_old2         rtol=1e-12

    # `wind_vel`/`tether_pos` are pre-allocated buffers reused across `step!` calls.
    @test size(te.wind_vel) == (3, segments)
    @test size(te.tether_pos) == (3, segments)

    # `clear!` resets the persistent state without touching `te.set`.
    QSM.clear!(te)
    @test te.state_vec == zeros(3)
    @test te.set === se
    nothing
end
nothing
