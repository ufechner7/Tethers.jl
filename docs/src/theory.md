```@meta
CurrentModule = Tethers.TetherComponents
```

# Theory

## Inputs and outputs
A general tether model should be able to simulate a tether that connects two arbitrary
points in space. The first tutorial examples assume one (fixed) point to be at the
coordinate [0,0,0]. From [Tether with arbitrary endpoints](@ref) onwards both end points
can be chosen freely, and [Re-usable tether component](@ref) exposes them as connectors,
so that a tether can be attached to anything that defines the motion of a point.

We assume that one end of the tether is either fixed or attached to a winch, and the other
end is fixed or attached to a load that applies a force on the tether.

### Inputs
- Either or
  - p1: [x1,y1,z1] vector [m] and vel1: speed vector of point one [m/s]
  - force1: [fx1, fy1, fz1] vector of the force applied to point one [N]
- Either or
  - p2: [x2,y2,z2] vector [m] and vel2: speed vector of point two [m/s]
  - force2: [fx2, fy2, fz2] vector of the force applied to point two [N]
- v_ro: reel-out speed at point one, scalar [m/s]
- v_wind: vector of the wind speed at reference height [m/s]

So far the examples implement the first alternative of each pair: an end point is either
held at a given position (`fix_p1`/`fix_p2`), or it is free and carries a point mass that
moves under gravity and the tether forces. Driving an end point with an externally applied
force is not implemented yet.

### Outputs
- If `fix_p1`:
  - force1: [fx1, fy1, fz1] force vector, felt at point one [N]
- If `fix_p2`:
  - force2: [fx2, fy2, fz2] force vector, felt at point two [N]
- pos: vector of the position vectors of the tether particles [m]
- vel: vector of the velocity vectors of the tether particles [m/s]
- forces: vector of the scalar forces per tether segment [N]

### Configuration
The settings below are implemented as fields of the struct `Settings3` in `Tether_08.jl`,
and of [`TetherSettings`](@ref) for the re-usable component:

- segments: number of tether segments [-]
- d_tether: tether diameter [mm]
- rho_tether: tether density [kg/m³]
- c_spring: unit spring constant [N]
- damping: unit damping constant [Ns]
- rel_compression_stiffness: stiffness of a slack segment, relative to a taut one [-]
- l0: initial unstretched tether length [m]
- v_ro: reel-out speed [m/s]
- α0: initial tether angle, used to derive the second end point if only one is given [rad]
- g_earth: gravitational acceleration vector [m/s²]
- v_wind_tether: wind velocity acting on the tether [m/s]
- rho: density of the fluid at position zero and 15 °C (water, air) [kg/m³]
- cd_tether: drag coefficient of the tether [-]
- duration: duration of the simulation [s]

The initial positions and velocities of the two end points (`p1_0`, `p2_0`, `vel1_0`,
`vel2_0`) are not settings; they are passed to `model()` and `main()` as the arguments `p1`
and `p2`.

The following settings describe a height dependent wind profile. They are **not implemented
yet**; the examples use the constant wind vector `v_wind_tether` for the whole tether:

- h_ref: reference height for the wind speed [m]
- alpha: exponent of the wind profile law [-]
- z0: surface roughness [m]
- `profile_law`: integer, 1=EXP, 2=LOG, 3=EXPLOG

## Math
The tether is modelled as ``n`` segments connecting ``n+1`` point masses. The following
equations are the ones implemented by the examples with drag and reel-out
([Segmented tether with aerodynamic drag](@ref) and later).

Segment ``i`` connects the particles ``i`` and ``i+1``:

```math
\mathbf{s}_i = \mathbf{p}_{i+1} - \mathbf{p}_i, \qquad
l_i = \| \mathbf{s}_i \|, \qquad
\hat{\mathbf{u}}_i = -\frac{\mathbf{s}_i}{l_i}
```

Reeling out changes the unstretched length of every segment, and with it the segment
stiffness, the segment damping and the particle mass. With the unit spring constant ``c``,
the unit damping constant ``d`` and the tether density ``\rho_t``:

```math
l_{seg}(t) = \frac{l_0 + v_{ro}\, t}{n}, \qquad
c_{seg} = \frac{c}{l_{seg}}, \qquad
d_{seg} = \frac{d}{l_{seg}}, \qquad
m_p = \rho_t\, \pi \left(\frac{d_{tether}}{2}\right)^2 l_{seg}
```

A tether can pull, but it can hardly push. This is modelled with a much smaller stiffness
for a slack segment, with the relative compression stiffness ``\epsilon``:

```math
c_{spr,i} = \frac{c_{seg}}{1+\epsilon}
            \left( \epsilon + \begin{cases} 1 & l_i > l_{seg} \\ 0 & \text{otherwise} \end{cases} \right)
```

The spring-damper force of segment ``i``, with the speed at which the segment is stretched
``v_{s,i}``:

```math
v_{s,i} = -\hat{\mathbf{u}}_i \cdot (\mathbf{v}_{i+1} - \mathbf{v}_i), \qquad
\mathbf{F}_{s,i} = \left( c_{spr,i} \left( l_i - l_{seg} \right) + d_{seg}\, v_{s,i} \right) \hat{\mathbf{u}}_i
```

Only the component of the apparent wind perpendicular to the segment creates drag. Half of
the drag force of a segment is applied to each of its two particles, which is why the usual
factor ``\tfrac{1}{2} \rho c_d A \| \mathbf{v} \| \mathbf{v}`` appears as ``\tfrac{1}{4}``
here, with the segment area ``A = l_i d_{tether}``:

```math
\mathbf{v}_{a,i} = \mathbf{v}_w - \frac{\mathbf{v}_i + \mathbf{v}_{i+1}}{2}, \qquad
\mathbf{v}_{\perp,i} = \mathbf{v}_{a,i} - \left( \mathbf{v}_{a,i} \cdot \hat{\mathbf{u}}_i \right) \hat{\mathbf{u}}_i
```

```math
\mathbf{F}_{d,i} = \frac{1}{4} \rho\, c_d\, \| \mathbf{v}_{\perp,i} \|\, l_i\, d_{tether}\, \mathbf{v}_{\perp,i}
```

Each inner particle feels the forces of the two segments it belongs to, and gravity:

```math
m_p\, \ddot{\mathbf{p}}_i = m_p\, \mathbf{g}
    + \mathbf{F}_{s,i-1} - \mathbf{F}_{s,i}
    + \mathbf{F}_{d,i-1} + \mathbf{F}_{d,i}
    \qquad \text{for } i = 2 \ldots n
```

The two end particles belong to one segment only and therefore carry half a particle mass.
A fixed end point has zero acceleration instead; the force it would need to stay in place is
the force felt at that end of the tether.

## Model export as functional mockup unit
Functional mockup units (FMUs) are a standard to exchange models between different
simulation environments, heavily used by the car and the aerospace industries.

Simulink, Python, Modelica etc can import FMU models. They are distributed as a
zip file that contains a shared library and an XML description.

A good and detailed introduction can be found [here](https://www.iea-annex60.org/finalReport/activity_1_2.html).

For the export of Julia models as FMU the package [FMIExport](https://github.com/ThummeTo/FMIExport.jl) can be used. For importing FMU models in Python the software [PyFMI](https://jmodelica.org/pyfmi/index.html#) can be used. FMU import with Simulink is documented [here](https://nl.mathworks.com/help/simulink/ug/work-with-fmi-in-simulink.html).

Because not everybody is using Julia as main development and simulation environment we plan
to provide a tether model as FMU. The issue that used to block this,
[FMIExport#10](https://github.com/ThummeTo/FMIExport.jl/issues/10), was closed in July 2026:
FMIExport v0.6.0 contains a first prototype of co-simulation export, and the remaining work
is tracked in [FMIExport#89](https://github.com/ThummeTo/FMIExport.jl/issues/89).

**Nomenclature:**
- FMU: Functional mockup unit
- FMI: Functional mockup interface
- FMI for model exchange: A model without a solver
- [FMI for co-simulation](https://fmi-standard.org/docs/3.0.1/#_fmi_for_co_simulation_cs): A model that includes its own solver
