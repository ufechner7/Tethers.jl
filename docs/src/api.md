# API

## Tether Component
The composable components used by [Tether_10.jl](https://github.com/ufechner7/Tethers.jl/blob/main/examples/Tether_10.jl),
explained in [Re-usable tether component](@ref). They live in the submodule
`Tethers.TetherComponents`.

```@meta
CurrentModule = Tethers.TetherComponents
```

```@docs
TetherSettings
set_diameter!
mass_per_meter
l_spring
m_end
Point3D
Tether
FixedEnd
FreeEnd
```

```@meta
CurrentModule = Tethers
```

## Quasi-steady model
A quasi-steady tether model: solves for tether shape and forces given the
ground-station orientation/tension and the kite's position and velocity. See
[`docs/quasisteady.md`](https://github.com/ufechner7/Tethers.jl/blob/main/docs/quasisteady.md)
for implementation notes. It lives in the submodule `Tethers.QuasiSteady`.

```@meta
CurrentModule = Tethers.QuasiSteady
```

### Public API
`StaticSettings`, `Tether`, `init!` and `step!` are the entry point: build a
`StaticSettings`, wrap it in a `Tether`, call `init!` once, then `step!` in a loop.

```@docs
StaticSettings
Tether
init!
step!
clear!
elevation
azimuth
tension
get_initial_conditions
get_analytic_catenary
```

### Private API
Unexported internals that the public API is built on, documented here only because
this project's `checkdocs = :all` requires every docstring in the module to appear
in the manual - not part of the recommended interface.

```@docs
check_wind_vel
simulate_tether
init_quasisteady
tether_shape
res!
scaled_res
lin_res
converged
node_kinematics
segment_drag
matlab_to_wind
wind_to_matlab
```

```@meta
CurrentModule = Tethers
```

## Utilities

```@docs
display_if_interactive
run_python
Tethers.copy_bin
Tethers.copy_file
Tethers.copy_examples
Tethers.example_packages
install_examples
```
