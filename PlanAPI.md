# Create new API for quasi-steady tether model

The API should consist of:

1. a struct `StaticSettings` with the model settings (do not change at runtime),
   including the initial conditions
2. a struct `Tether` with the state and the settings, a constructor te=Tether(se::StaticSettings)
3. a method `init!(te::Tether; prn=false)` that initializes the model
4. a method `step!(te::Tether, kite_pos, kite_vel; ...)` that moves the loose end point

The API follows the KiteModels.jl conventions: the settings struct carries the physical
properties *and* the initial conditions, `init!` takes no state arguments and reads
everything from `te.set`, and the per-step inputs of `step!` default to values derived
from the settings. Compare `init!(s::AKM; stiffness_factor, delta, prn, steady_state)`
and `next_step!(s, integrator; v_wind_gnd = s.set.v_wind, dt = 1/s.set.sample_freq, ...)`
in `KiteModels.jl`, and `elevation`/`azimuth`/`l_tether` in the KiteUtils `Settings`.

All coordinates and angle symbols follow the KiteUtils.jl definitions:
<https://opensourceawe.github.io/KiteUtils.jl/stable/reference_frames/>

## Reference frames and symbols

The quasi-steady model lives entirely in the **W (wind) reference frame**, whose
origin is the *anchor point of the tether* — which is exactly where the model
puts the ground station, at `[0, 0, 0]`. The z-axis points up, the y-axis
downwind. `kite_pos`, `kite_vel`, `wind_vel`, `tether_pos` and `p0` are all
expressed in this frame.

Angle symbols, per KiteUtils:

| Symbol | Quantity | Definition and sign |
|---|---|---|
| **β** | elevation | zero when the kite height is zero, 90° at zenith |
| **φ** | azimuth (wind frame) | positive anti-clockwise seen from above |
| **ψ** | heading / yaw | not used by the quasi-steady model; reserved |

So the tether direction at the ground station is
`dir ∝ [cos(β)cos(φ), cos(β)sin(φ), sin(β)]`, which is what
`src/qsm_conventions.jl` already documents.

**θ is not a KiteUtils symbol** and must not be used for elevation. It appears
today in `tether_shape`, in the `simulate_tether` / `res!` / `init_quasisteady`
docstrings ("theta [rad]") and as `theta_init`; all of these mean β and get
renamed. The one legitimate use of θ is `θ_m` in `matlab_to_wind` /
`wind_to_matlab`, where it denotes the *MATLAB reference* angle measured from
the vertical — a genuinely different quantity, already subscripted, and kept.

## Investigate

**What are the input parameters of the quasi steady model?**

The entry point is `simulate_tether` in `src/Tether_quasisteady.jl`; everything
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

### 4. Physical properties — `settings::Settings` (becomes `StaticSettings`)

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

### 5. Solver knobs

`prn = false` (print solver statistics) and `alg = DEFAULT_SOLVER` (TrustRegion
with ForwardDiff; a linear-parameterization fallback runs if it fails). `alg`
never changes at runtime and becomes a `StaticSettings` field; `prn` stays a
keyword argument of `init!`/`step!`, exactly as in `KiteModels.init!`.

### Conclusion

Groups 3, 4 and 5 are fixed for the lifetime of a simulation and become
`StaticSettings`; group 1 is state that must persist between calls; only group 2
genuinely varies per step. The current API forces the caller to thread all of it
through by hand, which is why `init_quasisteady` returns a six-tuple of mostly
unchanged inputs.

Group 2 has one exception that decides the shape of `init!`: the *initial* kite
position and tether length are needed once, to build the catenary starting guess,
and never again. KiteModels solves this by storing the initial condition in the
settings (`se.elevation`, `se.azimuth`, `se.l_tether`) so that `init!(s)` needs no
state arguments; we do the same, which is why the `StaticSettings` below carries
four initial-condition fields on top of the physical properties.

## Design

### 1. `struct StaticSettings`

Replaces `Settings`: same physical properties, plus the structural parameter, the
solver and the initial conditions. Nothing in it changes while a simulation runs.
Angles are in **degrees**, as in the KiteUtils `Settings` — radians appear only
inside `state_vec`.

```julia
@with_kw mutable struct StaticSettings @deftype Float64
    # structure
    segments::Int64 = 7
    # initial conditions (KiteUtils vocabulary)
    "initial elevation angle β                            [deg]"
    elevation = 70.0
    "initial wind-frame azimuth angle φ                   [deg]"
    azimuth = 0.0
    "initial unstretched tether length                      [m]"
    l_tether = 50.0
    "initial tether slack, l_tether = (1 + slack) * kite distance"
    slack = 0.05
    # physical properties (unchanged)
    rho = 1.225
    g_earth::MVector{3, Float64} = [0.0, 0.0, -9.81]
    cd_tether = 0.958
    d_tether = 4
    rho_tether = 724
    c_spring = 614600
    # solver
    alg = DEFAULT_SOLVER
end
```

`elevation`, `azimuth` and `l_tether` are the KiteUtils names, units and meanings
(`initial elevation angle [deg]`, …), so a caller who knows `KiteUtils.Settings`
already knows these. `slack` is the one addition, and it exists because the
catenary solve needs *both* the kite position and a tether length longer than the
straight-line distance to it: with a taut tether the coefficient goes to zero and
the solve in `init_quasisteady` degenerates. `slack = 0.05` reproduces the
`tether_length = 1.05 * norm(kite_pos)` idiom that
`examples/quasisteady/flying_circular.jl` uses today.

So the initial kite position is

```julia
kite_distance = se.l_tether / (1 + se.slack)
β, φ = deg2rad(se.elevation), deg2rad(se.azimuth)
kite_pos = kite_distance * [cos(β)cos(φ), cos(β)sin(φ), sin(β)]
```

which is the `dir` formula of `src/qsm_conventions.jl`.

`Settings` was kept briefly as a deprecated alias (`const Settings = StaticSettings`)
while `test/test_qsm.jl` and the examples were ported, since they constructed it
directly; once every caller in this repository used `StaticSettings`, the alias was
removed rather than carried forward indefinitely. **Decided.**

### 2. `struct Tether`

Mutable, holding the settings, the persistent state and the results of the last
`step!`. All vectors are in the W frame:

```julia
@with_kw mutable struct Tether
    set::StaticSettings = StaticSettings()
    # persistent state, updated by init! and step!
    state_vec::MVector{3, Float64} = zeros(MVector{3})   # (β [rad], φ [rad], Tn [N])
    # boundary conditions of the last step
    kite_pos::MVector{3, Float64} = zeros(MVector{3})
    kite_vel::MVector{3, Float64} = zeros(MVector{3})
    wind_vel::Matrix{Float64} = zeros(3, set.segments)
    tether_length::Float64 = 0.0
    # results of the last step
    tether_pos::Matrix{Float64} = zeros(3, set.segments)  # node coordinates
    force_gnd::Float64 = 0.0                             # tension at ground station [N]
    force_kite::MVector{3, Float64} = zeros(MVector{3})  # force on the tether end
    p0::MVector{3, Float64} = zeros(MVector{3})          # kite-tether attachment point
end
```

The constructor is `Tether(se::StaticSettings)`, mirroring `KPS4(kcu)`: it sizes
the buffers from `se.segments` and leaves the state at zero until `init!` runs.
The field is named `set`, as in `s.set` throughout KiteModels.

Buffers (`tether_pos`, `wind_vel`) are allocated once at construction, so a
stepping loop stays allocation free — `simulate_tether` currently allocates a
fresh `(3, segments)` matrix on every call.

Accessors keep the KiteUtils vocabulary at the surface, so callers never index
into `state_vec` by hand. They return **radians**, as `KiteModels.calc_elevation`
does:

```julia
elevation(te::Tether) = te.state_vec[1]   # β [rad]
azimuth(te::Tether)   = te.state_vec[2]   # φ [rad]
tension(te::Tether)   = te.state_vec[3]   # Tn [N]
```

### 3. `init!(te::Tether; prn = false)`

Replaces `init_quasisteady`. Takes no state arguments — everything comes from
`te.set` — and returns `te` rather than a six-tuple:

```julia
function init!(te::Tether; prn = false)
    clear!(te)                       # zero the state and the result buffers
    # kite_pos and tether_length from te.set, as above
    # catenary solve -> state_vec = (β, φ, Tn)
    step!(te, te.kite_pos, te.kite_vel; prn)   # leave the model in a solved state
    te
end
```

- `segments`, physical properties and `alg` come from `te.set`, so the
  five-branch `isnothing` cascade in `init_quasisteady` disappears entirely.
- `wind_vel` keeps whatever `te` already holds (zeros after construction); a
  caller who wants a wind field writes it into `te.wind_vel` before `init!`, or
  passes it to `step!`.
- `clear!(te::Tether)` resets the state from the settings, mirroring
  `KiteModels.clear!(s)`; `init!` calls it first, and it is also useful on its
  own to restart a simulation without rebuilding the buffers.
- Internally, `theta_init` becomes `beta_init` and `phi_init` stays.
- `init!` ends with one `step!`, so the model is in a consistent, solved state
  when it returns — today the caller has to remember to call `simulate_tether`
  once before the loop.

### 4. `step!(te::Tether, kite_pos, kite_vel; tether_length = nothing, wind_vel = nothing, prn = false)`

Moves the loose end point and re-solves. The keyword defaults follow
`KiteModels.next_step!`, which defaults its per-step inputs from the settings:

- `tether_length === nothing` → `(1 + te.set.slack) * norm(kite_pos)`, so the
  common reel-out-free case needs no argument at all.
- `wind_vel === nothing` → keep `te.wind_vel`.
- both are validated against `te.set.segments` by a single `check_wind_vel`
  helper.

It writes `state_vec`, `tether_pos`, `force_gnd`, `force_kite` and `p0` into `te`
and returns `te`, so a stepping loop reads:

```julia
se = StaticSettings(segments = 20, elevation = 60, l_tether = 52.5)
te = Tether(se)
init!(te)
for ii in 1:length(gamma)
    step!(te, traj[:, ii], vel[:, ii])
    all_force_kite[:, ii] .= te.force_kite
end
```

Compare `examples/quasisteady/flying_circular.jl`, where the same loop currently
threads seven values through `simulate_tether` and manually re-assigns
`state_vec`.

### 5. Naming

`init!`/`step!` with the bang, since both mutate `te`, following the Julia
convention for functions that mutate their first argument. **Decided.**
KiteModels calls the stepping function `next_step!`; we keep `step!`, since the
quasi-steady model has no integrator and no time step — `step!` moves the
boundary condition, it does not advance time.

## Decisions

- **The initial conditions live in `StaticSettings`**, not in the `init!`
  argument list, following `KiteModels.init!(s::AKM; …)` and the
  `elevation`/`azimuth`/`l_tether` fields of `KiteUtils.Settings`. The cost is
  that the initial kite position exists both as settings (`se.elevation`,
  `se.azimuth`, `se.l_tether`) and as state (`te.kite_pos`) — they deliberately
  diverge after the first `step!`, and `clear!` is what brings the state back to
  the settings. The benefit is that a whole simulation is described by one
  serializable settings struct, which is what makes KiteModels' settings.yaml
  workflow possible and is the reason for choosing this shape.
- **Angles in `StaticSettings` are degrees, `state_vec` is radians** — the
  KiteUtils split between what a user writes down and what the solver sees.
- **`slack` is a `StaticSettings` field.** The catenary initial guess needs a
  tether longer than the straight-line distance to the kite, so one length alone
  cannot describe the initial condition; `slack` also gives `step!` its
  `tether_length` default.
- **`segments` is a `StaticSettings` field**, not inferred from
  `size(wind_vel, 2)` as today — inferring it is fragile: a caller passing a wind
  field of the wrong width would silently change the discretization.
  `init!`/`step!` validate `wind_vel` against `te.set.segments` instead.
- **`simulate_tether` becomes internal.** `StaticSettings`/`Tether`/`init!`/
  `step!` is the documented, exported API; `simulate_tether` (and
  `init_quasisteady`) stay in the module unexported, since tests and the MATLAB
  comparison scripts still call them directly, but `docs/quasisteady.md` and
  `docs/src/references.md` document only the new API.
- **`res!` is dropped from the API.** It already isn't exported and stays
  purely internal — the residual `tether_shape` feeds to the nonlinear solve.
  The new API gets no equivalent; `test/test_qsm.jl` keeps calling `res!`
  directly (qualified) as it does today.
- **θ → β rename applies everywhere**, not just the new API: `tether_shape`'s
  argument, `theta_init` in `init_quasisteady`, and every existing "theta
  [rad]" docstring get renamed in the same pass (implementation step 1), so
  the old and new API don't leave two conventions side by side. `θ_m` in
  `src/qsm_conventions.jl` is a different quantity and is kept.

## Implementation steps

1. Rename θ → β for elevation across the `QuasiSteady` submodule: the
   `tether_shape` argument, `theta_init` in `init_quasisteady`, and the
   "theta [rad]" wording in every docstring. Leave `θ_m` in
   `src/qsm_conventions.jl` untouched. Pure rename, no behaviour change — do it
   first so the new code is written against the final vocabulary.
2. Align the `c_spring` wording with the rest of the package: the docstring
   becomes "unit spring constant [N] (= `E*A`)" instead of "axial stiffness of
   the tether EA [N]". Optionally rename the local `EA` in `tether_shape` to
   `c_spring`, since the surrounding comment only exists to explain why the two
   are the same number. Comment-only change.
3. Rename `Settings` → `StaticSettings` and add the `segments`, `elevation`,
   `azimuth`, `l_tether`, `slack` and `alg` fields; keep
   `const Settings = StaticSettings` as a deprecated alias while porting so the
   existing tests and MATLAB scripts keep working, then remove the alias once
   every caller in this repository uses `StaticSettings` directly. Document every
   field with the KiteUtils wording and units.
4. Add `Tether`, `Tether(se::StaticSettings)`, `clear!` and the
   `elevation`/`azimuth`/`tension` accessors; export them together with
   `StaticSettings`. Keep `simulate_tether` unchanged so nothing breaks mid-way.
   Once the new API is ported and tested (steps 7-9), unexport `simulate_tether`
   and `init_quasisteady` — they remain callable (qualified) for tests and the
   MATLAB comparison scripts, but are no longer part of the documented API.
5. Implement `init!(te; prn)` on top of the existing catenary solve: derive
   `kite_pos` and `tether_length` from `te.set`, drop the `isnothing` cascade,
   finish with one `step!`.
6. Implement `step!` on top of the existing solver setup, writing into the
   pre-allocated buffers in `te`, with the `tether_length`/`wind_vel` defaults
   above.
7. Port `examples/quasisteady/flying_circular.jl` first — it is the one example
   with a real stepping loop, so it is the best test of whether the API is
   pleasant to use, and its `1.05 * norm(kite_pos)` is exactly the `slack`
   default.
8. Port the remaining callers: `examples/quasisteady/run_catenary.jl`,
   `run_catenary_matlab.jl`, `force_plots.jl`, `benchmark_qsm.jl`.
9. Add tests in `test/test_qsm.jl` asserting that `init!` + `step!` reproduce the
   `init_quasisteady` + `simulate_tether` results bit for bit, for a
   `StaticSettings` whose `elevation`/`azimuth`/`l_tether`/`slack` correspond to
   the `kite_pos`/`tether_length` passed to the old API.
10. Confirm with `@allocated` that a `step!` in a loop does not allocate.
11. Update `docs/quasisteady.md` and `docs/src/references.md`, both of which
    document the current function-level API, and state the W frame and the β/φ
    symbols explicitly, linking the KiteUtils reference-frames page.

## Affected files

| File | Change |
|---|---|
| `src/Tether_quasisteady.jl` | θ → β rename; `Settings` → `StaticSettings` + new fields; new struct + methods; `simulate_tether`/`init_quasisteady` unexported |
| `src/qsm_conventions.jl` | unchanged code; comment updated to name the W frame and link KiteUtils |
| `examples/quasisteady/flying_circular.jl` | port to `Tether` |
| `examples/quasisteady/run_catenary.jl` | port to `Tether` |
| `examples/quasisteady/run_catenary_matlab.jl` | port to `Tether` |
| `examples/quasisteady/force_plots.jl` | port to `Tether` |
| `examples/quasisteady/benchmark_qsm.jl` | port; add an allocation check |
| `test/test_qsm.jl` | equivalence tests for the new API |
| `docs/quasisteady.md` | document `StaticSettings`, `Tether`, `init!`, `step!`, W frame, β/φ |
| `docs/src/references.md` | update the docstring listing |
