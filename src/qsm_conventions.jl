# Conversion between the MATLAB reference angle convention and the
# elevation/azimuth convention used throughout this package. See
# `docs/quasisteady.md`, section "Resolved: the angle convention", for the
# derivation and the evidence that only these two angles need converting.
#
# Wind-frame convention (used by `res!`, `init_quasisteady`, `calc_heading`,
# `calc_clock_angle`, and the system state from KiteUtils 0.8.2 onwards):
#   - β (elevation): 0 at the horizon, 90° at zenith
#   - φ (azimuth):   positive anti-clockwise seen from above
#   dir ∝ [cos(β)cos(φ), cos(β)sin(φ), sin(β)]
#
# MATLAB reference convention (used by the `.mat` fixtures in `test/data/`):
#   - θ measured from the vertical
#   - φ as stored in `stateVec`
#   dir ∝ [sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]

"""
    matlab_to_wind(θ_m, φ_m)

Convert the tether angles at the ground station from the MATLAB reference
convention, `dir ∝ [sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]`, to the convention
used throughout this package: elevation measured up from the horizontal
plane, azimuth in the wind reference frame, positive anti-clockwise seen
from above.

# Arguments
- θ_m: MATLAB-convention angle from the vertical [rad]
- φ_m: MATLAB-convention azimuth [rad]

# Returns
- (β, φ): elevation and wind-frame azimuth [rad]
"""
function matlab_to_wind(θ_m, φ_m)
    d1 = sin(θ_m) * cos(φ_m)
    d2 = sin(φ_m)
    d3 = cos(θ_m) * cos(φ_m)
    norm_d = sqrt(d1^2 + d2^2 + d3^2)
    d1, d2, d3 = d1 / norm_d, d2 / norm_d, d3 / norm_d
    β = asin(d3)
    φ = atan(d2, d1)
    return β, φ
end

"""
    wind_to_matlab(β, φ)

Inverse of [`matlab_to_wind`](@ref): convert elevation/wind-frame-azimuth
angles back to the MATLAB reference convention, for regenerating `.mat`
reference fixtures.

# Arguments
- β: elevation, measured up from the horizontal plane [rad]
- φ: azimuth in the wind reference frame, positive anti-clockwise seen from
  above [rad]

# Returns
- (θ_m, φ_m): MATLAB-convention angle from the vertical and azimuth [rad]
"""
function wind_to_matlab(β, φ)
    d1 = cos(β) * cos(φ)
    d2 = cos(β) * sin(φ)
    d3 = sin(β)
    φ_m = asin(d2)
    θ_m = atan(d1, d3)
    return θ_m, φ_m
end
