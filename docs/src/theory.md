# Theory

# Math

## Inputs and outputs
A general tether model should be able to simulate a tether that connects to arbitrary points in space. The first tutorial examples assume one (fixed) point to be at the coordinate [0,0,0].

We assume that one end of the tether is either fixed or attached to a winch, and the other end
is attached to a load that applies a force on the tether.

### Inputs
- pos1: [x1,y1,z1] vector
- pos2: [x2,y2,z2] vector
- v_ro: reel-out speed at point one, scalar
- force2: [fx2, fy2, fz2] vector of the force applied to point two

### Outputs
- force1: [fx1, fy1, fz1] force vector, felt at point one
- pos: vector of the position vectors of the tether particles
- vel: vector of the velocity vectors of the tether particles
- forces: vector of the scalar forces per tether segment

### Configuration
- d_tether: tether diameter [mm]
- rho_tether: tether density [kg/m³]
- c_spring: unit spring constant [N]
- damping: unit damping constant [Ns]
- segments: number of tether segments [-]
- l0: initial unstretched tether length [m]
- v_ro0: initial

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
- FMI for model exchange: A model without solver
- FMI for co-simulation: A model that includes its own solver