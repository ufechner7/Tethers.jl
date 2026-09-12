## Introduction
This package provides a tether component, and many tutorial scripts on tether modelling in Julia and Python.

A few examples where tether models can be useful:

- cranes
- undersea cables
- mooring lines of floating wind turbines
- airborne wind energy systems
- launching of sailplanes

![Tether](docs/images/Tether.gif)

Modeling of tethers and cables is difficult for several reasons. One of them is the high stiffness of the equation systems that need to be solved. I tried to implement these models with Simulink and Modelica and failed. It is possible to implement these models with Julia or Python. How to do this is explained in this tutorial. Tethers that are reeled in and out from
a winch are even more challenging to model than constant-length tethers.

A series of examples, from a simple falling mass towards a tether model, consisting of point masses connected by spring damper elements with the support of reel-out and reel-in and aerodynamic drag attached is presented.

#### Status

- all examples are of good quality and well-documented and tested
- the exported tether component is beta quality and needs further testing

## Installation

Make sure you are running a `bash` terminal (shell) and you have at least 16GB RAM (MAC, Linux and Windows supported).  
   On Windows, you can use [git for windows](https://gitforwindows.org/) which provides git AND a bash shell, but for using Julia from a `bash` terminal you must also install either  [Installation and usage of VSCode](https://ufechner7.github.io/Tethers.jl/dev/vscode/) or [Windows Terminal](https://learn.microsoft.com/en-us/windows/terminal/install). `Windows Terminal` is the simple and clean solution, `VSCode` the comfortable, fancy solution.

Check out from git:
```bash
cd repos # any folder of your choice, but without spaces in the folder name
git clone https://github.com/ufechner7/Tethers.jl
```

Build the system image:

```bash
cd repos/Tethers.jl
cd bin
./install
./create_sys_image
```

### Alternative: install as a package

If you don't want to clone the full repository, you can instead add `Tethers` to your own project and copy the examples into it. Create a new project:

```bash
mkdir test
cd test
julia --project="."
```

Then add `Tethers`.

```julia
using Pkg
pkg"add Tethers"
```

Copy the example scripts and the `bin` helper scripts to your project with:

```julia
using Tethers
install_examples()
```

This also adds the extra packages needed by the example scripts. You can now run the examples with:

```julia
include("examples/menu.jl")
```

For a faster start of the examples, build a system image (this takes a while and needs a lot of memory):

```bash
cd bin
./create_sys_image
```

`./bin/run_julia` picks it up automatically. Use `./bin/create_sys_image --update` to update all packages before building it.

## Basic example

Use the provided script to start Julia from the `Tethers.jl` folder:

```bash
cd repos/Tethers.jl
./bin/run_julia
```

From the Julia prompt, run the simulation:

```julia
include("examples/Tether_01.jl")
```

You should see a plot similar to:

![Falling mass](docs/images/FallingMass.png)

This example shows a mass that is thrown upwards, slows down and then falls.

**HINT**  
You get a menu from which you can run any of the examples by typing

```julia
menu()
```

at the Julia prompt.

**Julia code:** [Tether_01.jl](https://github.com/ufechner7/Tethers.jl/blob/main/examples/Tether_01.jl)

## Python version as comparison

From the Julia prompt execute:

```julia
run_python("Tether_01")
```

This will install Python, Matplotlib, NumPy, SciPy and CasADi and execute the script `Tether_01.py`.

**Python code:** [Tether_01.py](https://github.com/ufechner7/Tethers.jl/blob/main/examples/python/Tether_01.py)

**HINT**  
You get a menu from which you can run any of the Python examples by typing

```julia
menu2()
```

at the Julia prompt.

If you compare the Python and the Julia scripts you can see that:

- the Julia script is shorter and easier to read
- both reach the same speed, but only once each is given an analytic, sparse Jacobian

For a stiff, segmented tether the two are within about 20% of each other; see
[docs/julia_vs_python.md](https://github.com/ufechner7/Tethers.jl/blob/main/docs/julia_vs_python.md)
for the measurements and for what each ecosystem needs to get there.

Have a look at the [Examples](https://ufechner7.github.io/Tethers.jl/dev/examples/) that teach you how to construct a full tether model step by step.

## Overall comparison

Lines of code, excluding blank lines and comment lines. Execution times, which depend on
the Jacobian each side is given, are measured in
[docs/julia_vs_python.md](https://github.com/ufechner7/Tethers.jl/blob/main/docs/julia_vs_python.md)
rather than duplicated here.

| Test-case                          | LOC Julia | LOC Python |
|:-----------------------------------|:---------:|:----------:|
|Falling mass (1)                    |     36    |     48     |
|Non-linear Spring damper (3)        |     45    |     68     |
|ditto with callbacks (3b)           |     53    |    108     |
|swinging tether, 5 segments (5)     |    102    |    130     |
|Dyneema tether, reeling out (6)     |    115    |    154     |
|ditto with callbacks       (6c)     |    201    |    180     |
|Dyneema, reeling out with drag (7)  |    168    |    173     |

**Tradeoff Julia vs Python:** In Julia, the code is compiled before it is executed, which can cause about 5 to 30 seconds delay when running a simulation the first time, but speeds up the execution a lot afterward. In addition,
the Julia code is much more compact and better readable due to the use of symbolic differential equations.

Both sides need an exact, sparse Jacobian to handle the very stiff Dyneema tether at all.
Julia gets one from ModelingToolkit with `jac=true, sparse=true`; Python gets one from CasADi,
which differentiates the model expression and finds its sparsity by itself. With both compiled
and sparse the two land within about 20% of each other - see
[docs/julia_vs_python.md](https://github.com/ufechner7/Tethers.jl/blob/main/docs/julia_vs_python.md).
The remaining trade-off is setup: installing the required packages in Python takes less than a
minute, while installing and compiling the Julia software might take half an hour, because Julia
packages are distributed as source code and have to be compiled locally.

See also: [Why Julia?](https://ufechner7.github.io/2022/08/13/why-julia.html) and read the [documentation](https://ufechner7.github.io/Tethers.jl/dev/) or go straight to the [examples](https://ufechner7.github.io/Tethers.jl/dev/examples/).

## Citing

If you use Tethers.jl in your research, please cite it using the metadata in [CITATION.cff](https://github.com/ufechner7/Tethers.jl/blob/main/CITATION.cff) or the following BibTeX entry:

```bibtex
@software{fechner_tethers_jl,
  author  = {Fechner, Uwe and Bertozzi, Andrea},
  title   = {{Tethers.jl}},
  url     = {https://github.com/ufechner7/Tethers.jl},
  version = {2.0.0},
  date    = {2026-09-11},
  doi     = {10.5281/zenodo.19220562},
}
```
