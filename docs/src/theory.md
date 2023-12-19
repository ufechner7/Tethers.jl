# Theory

# Math

## Inputs and outputs

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