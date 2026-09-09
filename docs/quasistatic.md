# Quasi-static tether model: state of the port

Notes on porting the quasi-static tether model from the `andrea_quasistatic`
branch onto `main`, and on the open questions that port uncovered.

## Resolved: the angle convention

**The Julia model keeps its own convention; the MATLAB data is converted on
load.** The remaining work is to implement that — see *What to change* below.

### The convention the Julia code uses

The position of the kite is described by two angles, the azimuth angle `φ` and
the elevation angle `β`. The elevation angle is zero when the height of the
kite is zero and 90° when it is at zenith. The azimuth angle is the one in the
wind reference frame, defined **positive anti-clockwise when seen from above**
— the same convention `calc_heading()` and `calc_clock_angle()` use, and the
one written to the log file and the system state from KiteUtils 0.8.2 onwards.
(KiteUtils also knows `azimuth_north`, positive anti-clockwise, and
`azimuth_east`, positive clockwise; neither is used here.)

In this convention the tether direction at the ground station is

```julia
dir ∝ [cos(β)cos(φ), cos(β)sin(φ), sin(β)]
```

which is exactly what `res!` in `src/Tether_quasistatic.jl` already computes,
and what `init_quasistatic` already produces: `phi_init = atan(kite_pos[2],
kite_pos[1])` is anti-clockwise from above, and `theta_init = atan(z, hypot(x,
y))` is an elevation. So the state vector's `θ` **is** the elevation `β`.

### The convention the MATLAB data uses

The reference data in `test/data/` was produced with

```julia
dir ∝ [sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]
```

where `θ` is measured from the vertical. For
`test/data/input_basic_test.mat` (`θ = 18.435°`, `φ = -17.548°`, `Tn = 160941
N`, kite at `[100, 100, 300]`, 15 segments) the two readings give first
segments pointing in quite different directions:

| | first segment direction | `p0` |
| --- | --- | --- |
| `res!` as written | `[0.905, -0.286, 0.316]` | `[391.18, -123.70, 136.76]` |
| reference data | `[0.302, -0.302, 0.905]` | `[129.36, -129.36, 391.88]` |

The computed tether is a differently-oriented one, not a slightly inaccurate
one: `‖p0 - p0_ref‖ = 365.6 m` on a 431 m tether. Re-running `res!` with the
angles converted into the reference convention drops that to **1.61 m**, which
is what identifies the convention as the cause. A third candidate, the standard
spherical form `[sinθcosφ, sinθsinφ, cosθ]`, is ruled out at 90.0 m. Note that
`[sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]` and `[sin(θ), tan(φ), cos(θ)]` are the
same vector once normalised, and give identical results.

### It is a parametrisation difference, not a frame rotation

This is what makes the fix cheap. Only the two angles in `stateVec` are
affected; every *vector* quantity in the `.mat` files — `kitePos`, `kiteVel`,
`windVel`, and the reference outputs `p0`, `pj`, `T0` — is already in the same
right-handed, z-up frame the Julia code uses. The reference residual is a plain
componentwise difference of the two:

```
kitePos  [100, 100, 300]           p0_ref  [129.36, -129.36, 391.88]
Fobj_ref [-29.36, 229.36, -91.88]  ==  kitePos - p0_ref     ✓
```

If `p0_ref` were expressed in a y-flipped frame, that identity could not hold.
So there is no handedness flip to undo: MATLAB's `+φ` maps to `+y` just as the
wind-frame azimuth does. Nothing but `stateVec` needs converting, on the way in
or on the way out.

### Tolerance will not paper over this

| | max relative error | rtol needed to pass |
| --- | --- | --- |
| as written | 202 % | `2.02` |
| angles converted | 0.8 % | `0.008` |

At `rtol = 2.02` the assertion would no longer constrain anything. Only after
the convention is settled does a tolerance become meaningful — and then a
modest one does the job, though the residual 1.6 m is a second, smaller
discrepancy that still wants explaining. Land the conversion first, then chase
that 1.6 m against a tight tolerance.

### Why the loader and not `res!`

The alternative — changing `res!` to the MATLAB parametrisation — would put a
MATLAB-ism into the public state vector, alter what `state_vec` means for
`init_quasistatic`, `simulate_tether` and all six examples, and force a
conversion at every user-facing edge instead (`calc_heading`,
`calc_clock_angle`, `SysState`). `get_initial_conditions` is the only point at
which MATLAB angle data enters the package, so converting there leaves exactly
one convention in play everywhere else.

### What to change

1. Add the conversion and its inverse — the inverse is wanted as soon as the
   MATLAB reference is re-run to regenerate fixtures. `get_initial_conditions`
   is already duplicated verbatim between `src/Tether_quasistatic.jl` and
   `src/Tether_qsm_dual.jl`, so a small `src/qsm_conventions.jl` that both
   `include` beats a third copy.

   ```julia
   """
       matlab_to_wind(θ_m, φ_m)

   Convert the tether angles at the ground station from the MATLAB reference
   convention, `dir ∝ [sin(θ)cos(φ), sin(φ), cos(θ)cos(φ)]`, to the convention
   used throughout this package: elevation measured up from the horizontal
   plane, azimuth in the wind reference frame, positive anti-clockwise seen
   from above.
   """
   function matlab_to_wind(θ_m, φ_m)
       d = SVec3(sin(θ_m)*cos(φ_m), sin(φ_m), cos(θ_m)*cos(φ_m))
       d /= norm(d)
       return asin(d[3]), atan(d[2], d[1])      # elevation, azimuth
   end
   ```

2. Call it in both copies of `get_initial_conditions`:

   ```julia
   sv = vec(get(vars, "stateVec", 0))
   state_vec = MVector{3}(matlab_to_wind(sv[1], sv[2])..., sv[3])
   ```

3. **Bring `src/Tether_qsm_dual.jl` in line.** Its `res!` still computes
   `FT[1] = Tn·sinθ·cosφ`, `FT[2] = Tn·sinφ`, `FT[3] = Tn·cosθ·cosφ` — the
   MATLAB parametrisation, unconverted. The two implementations of the same
   model therefore disagree with each other today, and both are handed
   `state_vec` from the same loader (`examples/quasistatic/benchmark_qsm_dual.jl`).
   Converting in the loader *requires* this file to move to the
   elevation/azimuth form used by `Tether_quasistatic.jl`. That is an argument
   for the loader rather than against it: it collapses the two onto one
   convention instead of letting them drift further apart.

4. Turn the four `@test_broken` in `test/test_qsm.jl` back into `@test` with an
   `rtol`. The reference outputs need no conversion, per the section above.

### Caveat: the fixture's stored azimuth looks mirrored

Converted, the guess in `input_basic_test.mat` is elevation 64.76°, azimuth
**−45°**. The kite at `[100, 100, 300]` sits at elevation **64.76°**, azimuth
**+45°** — exact elevation match, mirrored azimuth. The MATLAB guess was
evidently generated from `kitePos` using a *clockwise* (`azimuth_east`)
convention, while MATLAB's own residual function treats `+φ` as `+y`, which is
why `p0_ref` has a negative y component.

For a residual unit test this is harmless: both sides evaluate the same input,
and a deliberately off-nominal guess is a reasonable thing to test a residual
at. **Do not flip the sign to make the guess point at the kite** — the stored
`p0_ref` pins the interpretation, and flipping it breaks the comparison. It
does suggest the two conventions are mixed on the MATLAB side too, which is
worth raising with whoever owns that code.

### Related: the unused frame transforms

`transformFromOtoW` / `transformFromWtoO` at the end of
`src/Tether_quasistatic.jl` are currently dead code, and their matrix carries
both a y-flip and a z-flip. If a wind-direction rotation is ever needed, that is
the second conversion site and it will have to be reconciled with the one above.

### Reproducing

```julia
include("src/Tether_quasistatic.jl")
sv, kp, kv, wv, tl, se = get_initial_conditions("test/data/input_basic_test.mat")
ref = matread("test/data/basic_test_results.mat")
Ns = size(wv, 2)
buffers() = [zeros(3, Ns) for _ in 1:5]

_, _, _, p0 = res!(zeros(3), sv, (kp, kv, wv, tl, se, buffers(), Ns, true))
norm(p0 .- vec(ref["p0"]))          # 365.6

# the same call with the angles converted to the reference convention
d = [sin(sv[1])cos(sv[2]), sin(sv[2]), cos(sv[1])cos(sv[2])]; d ./= norm(d)
sv2 = MVector(asin(d[3]), atan(d[2], d[1]), sv[3])
_, _, _, p0c = res!(zeros(3), sv2, (kp, kv, wv, tl, se, buffers(), Ns, true))
norm(p0c .- vec(ref["p0"]))         # 1.61
```

## Bugs found and fixed

**The test had never run.** `test/test_qsm.jl` passed `res!` a seven-element
`param` tuple, but `res!` destructures eight values from it, so the call died
with a `BoundsError` before reaching any assertion. This is why the convention
mismatch above went unnoticed for so long — it is a pre-existing discrepancy,
not a regression from the port. The test now builds the eight-field named
tuple, and locates its `.mat` files relative to `@__DIR__` rather than the
current working directory.

It also asserted `p0 isa Vector`, which is false: `res!` returns an
`MVector{3,Float64}`. That assertion was corrected rather than deleted.

**`GLMakie.lines!`/`scatterlines!` choked on raw `simulate_tether` output.**
`simulate_tether` returns `tether_pos` as an `MMatrix{3,segments}`, so a row
slice like `tether_pos[1,:]` is itself a `StaticArray`, not a `Base.Vector`.
GLMakie treats a fixed-size `StaticArray` passed as plot data as a single GPU
uniform rather than a per-vertex buffer, which fails shader compilation with
`Object GeometryBasics.Point{2, Float32} is not a supported uniform element
type`. `run_catenary.jl` and `flying_circular.jl` were unaffected because they
already run their `tether_pos` through `hcat` (which promotes to a plain
`Matrix`) before plotting; `run_catenary_matlab.jl` and `force_plots.jl` did
not, and crashed on their first `display_if_interactive` call. Fixed by
materializing `tether_pos` (and the `x_qs`/`y_qs` derived from it) with
`Matrix`/`collect` right after `simulate_tether` returns.

**`Tether_11.jl` was entirely dead code.** Its `main()` and the call that drove
it sat inside a `"""..."""` string literal, so including the file defined
`model` and `simulate` and then did nothing at all. `main()` was lifted out and
the plotting block ported to GLMakie.

**`Tether_11.jl` did not build under MTK 11.** Writing the differential
equations as

```julia
eqs1 = vcat(D.(pos) .~ vel, D.(vel) .~ acc)
eqs2 = vcat(eqs1...)
```

leaves a `Vector{Any}`, which the ModelingToolkit 11 `System` constructor
rejects. It now builds them per column, the way `examples/Tether_08.jl` does.

## Open question: the `maxiters` warning

`examples/Tether_11.jl` emits this during its steady-state solve, while
`examples/Tether_08.jl` does not:

```
Interrupted. Larger maxiters is needed. ...
```

The example completes and plots correctly regardless; the warning comes from
the `DynamicSS` call that produces the initial tether shape.

The obvious hypothesis is that `model` passes the prescribed `acc_p2` straight
into the steady-state solve, so the second end point accelerates forever and no
steady state exists — the same reason `se.v_ro` is zeroed on the line above.
**This was tested and is wrong**: zeroing `acc_p2` for the steady-state solve
does not remove the warning. The change was reverted, and the real cause is
still unknown.

## Verified

- The workspace root and the `test` environment both resolve against MTK 11.
- `test/test_qsm.jl`: 4 pass, 4 broken.
- All six `examples/quasistatic/` scripts run.
- `examples/Tether_11.jl` runs end to end (~67 s).
