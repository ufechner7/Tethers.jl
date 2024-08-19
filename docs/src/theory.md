# Theory

# Math

## Inputs and outputs
A general tether model should be able to simulate a tether that connects two arbitrary points in space. The first tutorial examples assume one (fixed) point to be at the coordinate [0,0,0].

We assume that one end of the tether is either fixed or attached to a winch, and the other end is fixed or attached to a load that applies a force on the tether.

### Inputs
- Either or
  - p1: [x1,y1,z1] vector [m] and vel1: speed vector of point one [m/s]
  - force1: vector of force applied to point one
- Either or
  - p2: [x2,y2,z2] vector [m] and vel2: speed vector of point two [m/s]
  - force2: [fx2, fy2, fz2] vector of the force applied to point two [N]
- v_ro: reel-out speed at point one, scalar [m/s]
- v_wind: vector of the wind speed at reference height [m/s]

### Outputs
- If `fix_p1`:
  - force1: [fx1, fy1, fz1] force vector, felt at point one [N]
- if `fix_p2`:
  - force2: [fx1, fy1, fz1] force vector, felt at point one [N]
- pos: vector of the position vectors of the tether particles [m]
- vel: vector of the velocity vectors of the tether particles [m/s]
- forces: vector of the scalar forces per tether segment [N]

### Configuration
- segments: number of tether segments [-]
- d_tether: tether diameter [mm]
- rho_tether: tether density [kg/m³]
- c_spring: unit spring constant [N]
- damping: unit damping constant [Ns]
- l0: initial unstretched tether length [m]
- v_ro0: initial reel-out speed [m/s]
- p1_0: initial position of point one [m]
- p2_0: initial position of point two [m]
- vel1_0: initial speed vector of point one [m/s]
- vel2_0: initial speed vector of point two [m/s]
- rho: density of the fluid at position zero and 15 °C (water, air) [kg/m³]
- h_ref: reference height for the wind speed [m]
- alpha: exponent of the wind profile law [-]
- z0: surface roughness [m]
- `profile_law`: integer, 1=EXP, 2=LOG, 3=EXPLOG, 4=FAST_EXP, 5=FAST_LOG, 6=FAST_EXPLOG  

## Model export as functional mockup unit
Functional mockup units (FMUs) are a standard to exchange models between different
simulation environments, heavily used by the car and the aerospace industries.

Simulink, Python, Modelica etc can import FMU models. They are distributed as a
zip file that contains a shared library and an XML description.

A good and detailed introduction can be found [here](https://www.iea-annex60.org/finalReport/activity_1_2.html).

For the export of Julia models as FMU the package [FMIExport](https://github.com/ThummeTo/FMIExport.jl) can be used. For importing FMU models in Python the software [PyFMI](https://jmodelica.org/pyfmi/index.html#) can be used. FMU import with Simulink is documented [here](https://nl.mathworks.com/help/simulink/ug/work-with-fmi-in-simulink.html).

Because not everybody is using Julia as main development and simulation environment we plan
to provide a tether model as FMU, but this will not be possible before this [issue](https://github.com/ThummeTo/FMIExport.jl/issues/10) is resolved.

**Nomenclature:**
- FMU: Functional mockup unit
- FMI: Functional mockup interface
- FMI for model exchange: A model without a solver
- [FMI for co-simulation](https://fmi-standard.org/docs/3.0.1/#_fmi_for_co_simulation_cs): A model that includes its own solver
