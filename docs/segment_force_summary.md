# Axial force of a wind-loaded tether — analytical solution

A tether is held between two fixed points a distance `L` apart. Its unstretched length is
`L₀`, which may be *longer* than `L`. A steady cross wind blows perpendicular to the line
between the anchors; there is no gravity and no reel-out. What is the axial force in the
tether?

This is the continuum solution, for a tether treated as a smooth curve rather than as a
chain of discrete segments.

## Symbols

| symbol | meaning | unit |
|:------|:--------------------------------------------|:-----|
| `L`  | distance between the two anchor points | m |
| `L₀` | unstretched tether length | m |
| `r`  | `L₀/L`, the length ratio | – |
| `δ`  | `r - 1`; positive when the anchors are closer together than the tether is long | – |
| `EA` | axial stiffness of the tether | N |
| `w`  | aerodynamic drag per unit length, `½ ρ c_d d v²` | N/m |
| `s`  | sag: lateral deflection at mid-span | m |
| `F`  | axial force, positive under tension | N |

## The tether is never in compression

The result that makes a closed-form solution possible.

A slack cable carries no transverse load at all. So when the anchors are brought closer
together than the unstretched length — `δ > 0`, nominally compression — the tether does not
sit there compressed. It bows out under the wind until its arc is long enough to be taut
again, and then keeps being stretched until the tension balances the drag. Its arc is
therefore always *longer* than `L₀`, and the force is always tensile.

The whole solution below can consequently assume a single, taut stiffness `EA`; no
slack-cable or compression branch ever enters.

## Derivation

Three relations close the problem.

**Sag.** A cable under a transverse load uniform along its chord takes a parabolic shape,
with mid-span deflection set by the transverse force balance:

```math
s = \frac{w L^2}{8 F}
```

**Arc length.** For a parabola of sag `s` over a span `L`, to leading order in `s/L` the
curve is longer than its chord by

```math
\Delta S = \frac{8 s^2}{3 L}
```

**Elasticity.** That arc must be exactly as long as the tether stretched under its own
tension:

```math
\Delta S = L_0 \left(1 + \frac{F}{EA}\right) - L
```

Eliminating `s` and `ΔS` and dividing through by `L` gives a cubic in `F`:

```math
\boxed{\;\frac{r}{EA} F^3 + (r-1)\,F^2 = \frac{w^2 L^2}{24}\;}
```

Its linear term vanishes. The left-hand side is negative at `F = 0` and increases without
bound, so there is exactly one positive root — the physical solution.

## Dimensionless form

Scaling the force by the stiffness, `f = F/EA`, and defining a dimensionless load

```math
\Lambda = \frac{w L}{\sqrt{24}\,EA}
```

reduces the whole problem to one equation in one parameter pair `(δ, Λ)`:

```math
r f^3 + \delta f^2 = \Lambda^2, \qquad r = 1 + \delta
```

and, since `|δ| ≪ 1` in any realistic case, simply

```math
f^3 + \delta f^2 \approx \Lambda^2
```

The entire behaviour of the tether is the competition between two terms: `f³`, the drag
holding a bowed cable out, and `δf²`, the slack that has to be taken up.

## Closed form

With `a = EA·δ/r` and `b = -EA·w²L²/(24r)`, the cubic is `F³ + aF² + b = 0`. Substituting
`F = y - a/3` removes the quadratic term, giving `y³ + py + q = 0` with

```math
p = -\frac{a^2}{3}, \qquad q = \frac{2a^3}{27} + b, \qquad
D = \frac{q^2}{4} + \frac{p^3}{27}
```

For `D > 0` there is a single real root,

```math
F = \sqrt[3]{-\tfrac{q}{2} + \sqrt{D}}
  + \sqrt[3]{-\tfrac{q}{2} - \sqrt{D}}
  - \frac{a}{3}
```

and for `D ≤ 0` all three roots are real,

```math
F_k = 2\sqrt{-\frac{p}{3}}\;
      \cos\!\left(\frac{1}{3}\arccos\!\left(\frac{3q}{2p}\sqrt{\frac{-3}{p}}\right)
                  - \frac{2\pi k}{3}\right) - \frac{a}{3}, \quad k = 0,1,2
```

of which exactly one is positive.

## The three regimes

Each limit of the cubic is a regime with its own scaling. Note that the tether diameter
enters twice and differently — `EA ∝ d²` through the cross section, `w ∝ d` through the
frontal area — so each regime carries a different power of `d`.

| regime | condition | force | dimensionless |
|:-------------|:-------------------|:---------------------|:-------------|
| **extension** | `δ < 0`, `-δ ≫ Λ^(2/3)` | `F = EA·(1-r)/r` | `f → -δ/r` |
| **crossover** | `δ = 0` | `F = ∛(EA·w²L²/24)` | `f → Λ^(2/3)` |
| **compression** | `δ > 0`, `δ ≫ Λ^(2/3)` | `F = wL/√(24δ)` | `f → Λ/√δ` |

| regime | `d` | `v` | `L` | `δ` |
|---|---|---|---|---|
| extension | `d²` | — | — | `-δ` |
| crossover | `d⁴ᐟ³` | `v⁴ᐟ³` | `L²ᐟ³` | — |
| compression | `d¹` | `v²` | `L¹` | `δ^(-1/2)` |

**Extension** is pure Hooke's law: the tether is pulled taut, the drag is negligible against
the spring, and the force follows the cross section, `d²`. **Compression** is pure drag: the
stiffness has cancelled out entirely and the force is set by how hard the wind pushes on the
frontal area, `d¹`. The **crossover** at `δ = 0` is the geometric blend of the two — `EA ∝ d²`
and `w² ∝ d²` under a cube root give `d⁴ᐟ³`.

## Why the force never changes sign

In the compression regime,

```math
F \;\longrightarrow\; \frac{wL}{\sqrt{24\,\delta}}
```

which is strictly positive for every `δ > 0` and every non-zero wind. There is no
combination of length, diameter, wind speed or compression at which the force vanishes or
reverses; it only becomes small, as `δ^(-1/2)` flattens and `w` shrinks. The force reaches
zero solely in the limit `w → 0`, i.e. when the cross wind disappears — at which point the
tether goes genuinely slack and carries nothing.

## Range of validity

The parabolic shape and the `8s²/3L` arc length are both leading order in `s/L`, so the
solution degrades as the tether bows further. At 10% compression the sag is already about
20% of the span, which costs well under a percent; far beyond that the assumption that the
drag is uniform along the chord also weakens, because the wind meets the strongly tilted
end sections at an angle.

`F` is the **mean** tension along the arc. In a deeply bowed tether the tension is not
uniform: it is lowest at mid-span and highest at the anchors, where the transverse load has
accumulated. The anchor force can exceed the mean substantially and is not what this
solution returns.
