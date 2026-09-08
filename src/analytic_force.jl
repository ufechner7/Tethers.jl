"""
    analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments=Inf)

Analytical prediction of the mean axial force [N], tension positive.

`se` is a [`TetherSettings`](@ref) (or anything with the same `rho`, `cd_tether`, `d_tether`
and `c_spring` fields); `rho` and `cd_tether` are taken from it directly, and the axial
stiffness is scaled from `se.c_spring` to the given `d_segment` the same way
[`set_diameter!`](@ref) does, `EA = se.c_spring * (d_segment/se.d_tether)^2`, so `se` need
not already be set to that diameter.

A tether under a transverse load that is uniform along the chord hangs in a parabola, and a
chain of straight segments under uniform point loads is exactly that parabola sampled at
its nodes. Three relations close the system, with `L` the distance between the anchors,
`L0` the unstretched length, `EA` the axial stiffness and `w` the drag per unit length:

- sag from the force balance,        `s  = w L² / (8F)`
- arc length of the sampled parabola `ΔS = (1 - 1/n²) · 8s²/(3L)`
- and the elastic law,               `ΔS = L0 (1 + F/EA) - L`

Eliminating `s` and `ΔS` leaves a cubic in `F` with exactly one positive root, solved here
in closed form:

    (r/EA) F³ + (r-1) F² = (1 - 1/n²) w² L² / 24,   r = L0/L

The `1 - 1/n²` is the only discretisation term: an `n`-segment polyline through a parabola
is that much shorter than the smooth curve, so it needs that much more sag — and hence less
force — to take up the same slack. It is derived, not fitted.

The default `segments=Inf` gives the continuum limit, i.e. the formula for a real tether
rather than for a chain of `n` segments — `1 - 1/Inf^2` is exactly `1`, so the term simply
drops out. Pass the actual segment count to match a segmented model instead: a 6-segment
model over-predicts the continuum limit by about `1/(2n²)`, 1.4%; by 20 segments it is
0.13%.

Reproduces the measured force to 0.6% over the whole sweep; see `check_formula` in
`examples/plot_compression.jl`. The derivation is in `docs/segment_force.md`.
"""
function analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments=Inf)
    EA = se.c_spring * (d_segment/se.d_tether)^2               # axial stiffness    [N]
    w  = 0.5 * se.rho * se.cd_tether * (d_segment/1000) * v_wind_perp^2   # drag per meter [N/m]
    r  = l_unstretched / l_segment
    cn = 1 - 1/segments^2
    # F³ + a₂F² + a₀ = 0; a₁ vanishes, which is what makes the closed form short
    a2 = EA * (r - 1) / r
    a0 = -cn * EA * w^2 * l_segment^2 / (24r)
    # `a0 == 0` -- a single segment (`cn == 0`), or no transverse load -- degenerates the
    # cubic to F²(F + a2) = 0, a double root at zero. Both closed forms below break there:
    # `D` is analytically zero, so its floating-point sign is pure round-off and decides at
    # random between the two roots, and `sqrt`/`acos` at that point have an infinite
    # derivative, which fills any ForwardDiff Jacobian with NaN. Return the elastic law
    # directly instead: tension `-a2` while stretched, exactly zero once slack.
    iszero(a0) && return max(-a2, zero(a2))
    p  = -a2^2/3
    q  = 2a2^3/27 + a0
    D  = q^2/4 + p^3/27
    if D >= 0           # one real root (or a repeated root at D == 0)
        # `-q/2 ± sqrt(D)` cancels catastrophically as `p -> 0`, i.e. exactly where the
        # segment sits at its unstretched length: one of the two cube roots collapses onto
        # zero, where `cbrt` has an infinite derivative and fills a ForwardDiff Jacobian
        # with NaN. So take whichever sign *adds* rather than subtracts, and recover the
        # second root of the Cardano pair from the identity `u v = -p/3`.
        u = cbrt(-q/2 + (q <= 0 ? sqrt(D) : -sqrt(D)))
        v = iszero(u) ? zero(u) : -p/(3u)
        return u + v - a2/3
    end                 # three real roots, the physical one is the positive one
    # D < 0 forces p < 0 strictly, since p = -a2²/3 ≤ 0 always; so m below is never zero
    m  = 2sqrt(-p/3)
    th = acos(clamp(3q/(p*m), -1, 1))/3
    maximum(m*cos(th - 2π*k/3) - a2/3 for k in 0:2)
end

"""
    hooke_force(se; d_segment, l_unstretched, l_segment)

Axial force [N] of a segment under plain Hooke's law with a *constant* stiffness, tension
positive: `EA (l_segment - l_unstretched) / l_unstretched`, with `EA` scaled to `d_segment`
exactly as in [`analytic_force`](@ref).

Unlike a real tether this reference spring also pushes back when it is compressed, so it is
negative below the unstretched length. That is what makes it the yardstick
[`damping_factor`](@ref) measures the slack-capable tether against.
"""
hooke_force(se; d_segment, l_unstretched, l_segment) =
    se.c_spring * (d_segment/se.d_tether)^2 * (l_segment - l_unstretched) / l_unstretched

"""
    damping_factor(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments=Inf)

Fraction in `[0, 1]` of the nominal axial damping that a segment still carries: the quotient
of [`analytic_force`](@ref) and `abs(`[`hooke_force`](@ref)`)`, capped at one.

A slack cable does not damp axial motion, so the damper has to fade out together with the
tension rather than being switched off by hand. While the segment is stretched the sag makes
the analytical force the larger of the two -- `ΔS >= 0` in the derivation above is exactly
`F >= F_hooke` -- so the quotient saturates at `1` and the damping is left untouched. Once
the segment is shorter than its unstretched length `F_hooke` turns negative and its magnitude
grows linearly, while the wind-driven `F` stays bounded, so the quotient, and with it the
damping, decays smoothly to zero.

The cap also removes the singularity at `l_segment == l_unstretched`, where `F_hooke` is zero.

Note that this is smooth only as long as there *is* a transverse load to hold the slack
segment in tension. With `segments=1` (or `v_wind_perp=0`) `analytic_force` is exactly
`max(0, F_hooke)`, and the quotient degenerates to a hard `0`/`1` switch at the unstretched
length.
"""
function damping_factor(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments=Inf)
    f  = analytic_force(se; v_wind_perp, d_segment, l_unstretched, l_segment, segments)
    fh = abs(hooke_force(se; d_segment, l_unstretched, l_segment))
    fh <= f && return one(f)    # also catches fh == 0 at l_segment == l_unstretched
    f / fh
end
