# Create new API for quasi-steady tether model

The API should consist of:

1. a struct `Tether` with the state and the settings
2. a method `init(te::Tether, params)` that initializes the model
4. a method `step(te::Tether)` that moves the loose end point

All coordinates and angle symbols follow the KiteUtils.jl definitions:
<https://opensourceawe.github.io/KiteUtils.jl/stable/reference_frames/>

## Reference frames and symbols

The quasi-static model lives entirely in the **W (wind) reference frame**, whose
origin is the *anchor point of the tether* — which is exactly where the model
puts the ground station, at `[0, 0, 0]`. The z-axis points up, the y-axis
downwind. `kite_pos`, `kite_vel`, `wind_vel`, `tether_pos` and `p0` are all
expressed in this frame.

Angle symbols, per KiteUtils:

| Symbol | Quantity | Definition and sign |
|---|---|---|
| **β** | elevation | zero when the kite height is zero, 90° at zenith |
| **φ** | azimuth (wind frame) | positive anti-clockwise seen from above |
| **ψ** | heading / yaw | not used by the quasi-static model; reserved |

So the tether direction at the ground station is
`dir ∝ [cos(β)cos(φ), cos(β)sin(φ), sin(β)]`, which is what
`src/qsm_conventions.jl` already documents.

**θ is not a KiteUtils symbol** and must not be used for elevation. It appears
today in `tether_shape`, in the `simulate_tether` / `res!` / `init_quasistatic`
docstrings ("theta [rad]") and as `theta_init`; all of these mean β and get
renamed. The one legitimate use of θ is `θ_m` in `matlab_to_wind` /
`wind_to_matlab`, where it denotes the *MATLAB reference* angle measured from
the vertical — a genuinely different quantity, already subscripted, and kept.

## Investigate

**What are the input parameters of the quasi steady model?**

The entry point is `simulate_tether` in `src/Tether_quasistatic.jl`; everything
below is either one of its arguments or reachable from one. They fall into five
groups, and that grouping is what drives the API design below.

### 1. State / initial guess — `state_vec::MVector{3, Float64}`

`(β, φ, Tn)` at the **ground station**: elevation [rad], wind-frame azimuth
[rad], tension [N]. It is both an input (the initial guess for the nonlinear
solve) and an output, so it carries over from one step to the next.

### 2. Boundary conditions — change every step

All vectors in the W frame.

| Parameter | Type | Meaning |
|---|---|---|
| `kite_pos` | `MVector{3, Float64}` | kite position [m] |
| `kite_vel` | `MVector{3, Float64}` | kite velocity [m/s] |
| `wind_vel` | `(3, segments)` matrix | wind velocity **per segment** [m/s] |
| `tether_length` | `Float64` | unstretched tether length [m] |

`kite_vel` is not used as a velocity field directly: `tether_shape` derives
`ω = (p × v)/|p|²` and `v_parallel = v · p̂` from it and then assumes the tether
rotates rigidly with the kite. The velocity input is therefore effectively
*(radial rate, angular rate)*.

### 3. Structural — `segments`

Not a separate argument today: it is inferred as `size(wind_vel, 2)`. Segment
length is `tether_length / (segments + 1)`.

### 4. Physical properties — `settings::Settings`

| Field | Default | Unit / meaning |
|---|---|---|
| `rho` | 1.225 | kg/m³, air density |
| `g_earth` | `[0, 0, -9.81]` | m/s², W frame (only `abs(g_earth[3])` is used) |
| `cd_tether` | 0.958 | drag coefficient of the tether |
| `d_tether` | 4 | **mm**, tether diameter |
| `rho_tether` | 724 | kg/m³ — a *density*, not mass per unit length |
| `c_spring` | 614600 | N, unit spring constant (= `E*A`) |

`c_spring` is the **unit spring constant [N]** — the spring constant of a one
metre segment — which is the name the rest of the package uses
(`docs/src/theory.md`, `examples/Tether_08.jl`). It equals `E*A`, and
`get_initial_conditions` computes it that way from the MATLAB fixtures, but
`E*A` is the axial *rigidity*, not a stiffness: the spring constant of a segment
is `c_spring / l_segment` [N/m], as `Tether_08.jl` makes explicit. The
`Settings` docstring currently calls it "axial stiffness of the tether EA [N]",
which is both off-vocabulary and dimensionally loose, and is corrected below.

### 5. Solver knobs — keyword arguments

`prn = false` (print solver statistics) and `alg = DEFAULT_SOLVER` (TrustRegion
with ForwardDiff; a linear-parameterization fallback runs if it fails).

### Conclusion

Groups 3 and 4 are fixed for the lifetime of a simulation, group 1 is state that
must persist between calls, and only group 2 genuinely varies per step. The
current API forces the caller to thread all of it through by hand, which is why
`init_quasistatic` returns a six-tuple of mostly unchanged inputs.

## Design

### 1. `struct Tether`

Mutable, holding settings, structural parameters, persistent state and the
results of the last `step!`. All vectors are in the W frame:

```julia
@with_kw mutable struct Tether
    settings::Settings = Settings()
    segments::Int = 7
    # persistent state, updated by init! and step!
    state_vec::MVector{3, Float64} = zeros(MVector{3})   # (β [rad], φ [rad], Tn [N])
    # boundary conditions of the last step
    kite_pos::MVector{3, Float64} = zeros(MVector{3})
    kite_vel::MVector{3, Float64} = zeros(MVector{3})
    wind_vel::Matrix{Float64} = zeros(3, segments)
    tether_length::Float64 = 0.0
    # results of the last step
    tether_pos::Matrix{Float64} = zeros(3, segments)     # node coordinates
    force_gnd::Float64 = 0.0                             # tension at ground station [N]
    force_kite::MVector{3, Float64} = zeros(MVector{3})  # force on the tether end
    p0::MVector{3, Float64} = zeros(MVector{3})          # kite-tether attachment point
    # solver configuration
    alg = DEFAULT_SOLVER
end
```

Buffers (`tether_pos`) are allocated once at construction, so a stepping loop
stays allocation free — `simulate_tether` currently allocates a fresh
`(3, segments)` matrix on every call.

Accessors keep the KiteUtils vocabulary at the surface, so callers never index
into `state_vec` by hand:

```julia
elevation(te::Tether) = te.state_vec[1]   # β [rad]
azimuth(te::Tether)   = te.state_vec[2]   # φ [rad]
tension(te::Tether)   = te.state_vec[3]   # Tn [N]
```

### 2. `init!(te::Tether, params)`

Replaces `init_quasistatic`. Sets `te.state_vec` from the catenary solution and
stores the boundary conditions, returning `te` rather than a six-tuple:

```julia
init!(te::Tether, kite_pos, tether_length; kite_vel = zeros(MVector{3}), wind_vel = nothing)
```

- `segments` and `settings` come from `te`, not from the argument list.
- `wind_vel === nothing` keeps whatever `te` already holds (zeros by default).
- Validation moves to the constructor / a single `check_wind_vel` helper; the
  current five-branch `isnothing` cascade in `init_quasistatic` disappears.
- Internally, `theta_init` becomes `beta_init` and `phi_init` stays.
- Also runs one `step!` so that the model is in a consistent, solved state after
  `init!` — today the caller has to remember to call `simulate_tether` once
  before the loop.

### 3. `step!(te::Tether, kite_pos, kite_vel; wind_vel = te.wind_vel, tether_length = te.tether_length)`

Moves the loose end point and re-solves. `wind_vel` and `tether_length` are
keywords because they often stay constant across steps. It writes `state_vec`,
`tether_pos`, `force_gnd`, `force_kite` and `p0` into `te` and returns `te`, so
a stepping loop reads:

```julia
te = Tether(segments = 20)
init!(te, kite_pos, norm(kite_pos); kite_vel)
for ii in 1:length(gamma)
    step!(te, traj[:, ii], vel[:, ii]; tether_length = 1.05 * norm(traj[:, ii]))
    all_force_kite[:, ii] .= te.force_kite
end
```

Compare `examples/quasistatic/flying_circular.jl`, where the same loop currently
threads seven values through `simulate_tether` and manually re-assigns
`state_vec`.

### 4. Naming

`init!`/`step!` with the bang, since both mutate `te`, following the Julia
convention for functions that mutate their first argument. **Decided.**

## Decisions

- **`segments` is a `Tether` field**, not inferred from `size(wind_vel, 2)` as
  today — inferring it is fragile: a caller passing a wind field of the wrong
  width would silently change the discretization. `init!`/`step!` validate
  `wind_vel` against `te.segments` instead.
- **`simulate_tether` becomes internal.** `Tether`/`init!`/`step!` is the
  documented, exported API; `simulate_tether` (and `init_quasistatic`) stay in
  the module unexported, since tests and the MATLAB comparison scripts still
  call them directly, but `docs/quasistatic.md` and `docs/src/references.md`
  document only the new API.
- **`res!` is dropped from the API.** It already isn't exported and stays
  purely internal — the residual `tether_shape` feeds to the nonlinear solve.
  The new API gets no equivalent; `test/test_qsm.jl` keeps calling `res!`
  directly (qualified) as it does today.
- **θ → β rename applies everywhere**, not just the new API: `tether_shape`'s
  argument, `theta_init` in `init_quasistatic`, and every existing "theta
  [rad]" docstring get renamed in the same pass (implementation step 1), so
  the old and new API don't leave two conventions side by side. `θ_m` in
  `src/qsm_conventions.jl` is a different quantity and is kept.

## Implementation steps

1. Rename θ → β for elevation across the `Quasistatic` submodule: the
   `tether_shape` argument, `theta_init` in `init_quasistatic`, and the
   "theta [rad]" wording in every docstring. Leave `θ_m` in
   `src/qsm_conventions.jl` untouched. Pure rename, no behaviour change — do it
   first so the new code is written against the final vocabulary.
2. Align the `c_spring` wording with the rest of the package: the `Settings`
   docstring becomes "unit spring constant [N] (= `E*A`)" instead of "axial
   stiffness of the tether EA [N]". Optionally rename the local `EA` in
   `tether_shape` to `c_spring`, since the surrounding comment only exists to
   explain why the two are the same number. Comment-only change.
3. Add `Tether`, `init!`, `step!` and the `elevation`/`azimuth`/`tension`
   accessors to the `Quasistatic` submodule; export them. Keep `simulate_tether`
   unchanged so nothing breaks mid-way.
   Once the new API is ported and tested (steps 6-8), unexport `simulate_tether`
   and `init_quasistatic` — they remain callable (qualified) for tests and the
   MATLAB comparison scripts, but are no longer part of the documented API.
4. Re-implement `init!` on top of the existing catenary solve, dropping the
   `isnothing` cascade.
5. Make `step!` call into `tether_shape` through the existing solver setup,
   writing into the pre-allocated buffers in `te`.
6. Port `examples/quasistatic/flying_circular.jl` first — it is the one example
   with a real stepping loop, so it is the best test of whether the API is
   pleasant to use.
7. Port the remaining callers: `examples/quasistatic/run_catenary.jl`,
   `run_catenary_matlab.jl`, `force_plots.jl`, `benchmark_qsm.jl`.
8. Add tests in `test/test_qsm.jl` asserting that `init!` + `step!` reproduce the
   `init_quasistatic` + `simulate_tether` results bit for bit.
9. Confirm with `@allocated` that a `step!` in a loop does not allocate.
10. Update `docs/quasistatic.md` and `docs/src/references.md`, both of which
    document the current function-level API, and state the W frame and the β/φ
    symbols explicitly, linking the KiteUtils reference-frames page.

## Affected files

| File | Change |
|---|---|
| `src/Tether_quasistatic.jl` | θ → β rename; new struct + methods; `simulate_tether`/`init_quasistatic` unexported |
| `src/qsm_conventions.jl` | unchanged code; comment updated to name the W frame and link KiteUtils |
| `examples/quasistatic/flying_circular.jl` | port to `Tether` |
| `examples/quasistatic/run_catenary.jl` | port to `Tether` |
| `examples/quasistatic/run_catenary_matlab.jl` | port to `Tether` |
| `examples/quasistatic/force_plots.jl` | port to `Tether` |
| `examples/quasistatic/benchmark_qsm.jl` | port; add an allocation check |
| `test/test_qsm.jl` | equivalence tests for the new API |
| `docs/quasistatic.md` | document `Tether`, `init!`, `step!`, W frame, β/φ |
| `docs/src/references.md` | update the docstring listing |
