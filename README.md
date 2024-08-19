# Tethers.jl
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ufechner7.github.io/Tethers.jl/dev)
[![Build Status](https://github.com/ufechner7/Tethers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ufechner7/Tethers.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Introduction

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
- the Julia examples are of good quality and well-documented and tested
- the Python examples need some more work

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
./create_sys_image
```

## Basic example
Use the provided script to start Julia from the `Tethers.jl` folder:
```bash
cd repos/Tethers.jl
./bin/run_julia
```
From the Julia prompt, run the simulation:
```julia
include("src/Tether_01.jl")
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

**Julia code:** [Tether_01.jl](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_01.jl)

## Python version as comparison
From the Julia prompt execute:
```
include("src/RunTether.jl")
```
This will install Python, Matplotlib and Assimulo and execute the script `Tether_01.py`.

**Python code:** [Tether_01.py](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_01.py)

If you compare the Python and the Julia scripts you can see that:
- the Julia script is shorter and easier to read
- Julia is about 16 times faster when running the simulation  
For a stiff, segmented tether (example 6 and 7) the Julia solvers are more than 2000 times faster than Python.

Have a look at the [Examples](https://ufechner7.github.io/Tethers.jl/dev/examples/) that teach you how to construct a full tether model step by step.

## Overall comparison
Execution time for a simulation of 10s duration with logging the state every 20ms.
Relative and absolute tolerance: $1.0^{-6}$. CPU: Ryzen 9 7950X.

| Test-case                   | Lines of code (LOC) Julia | LOC Python | Time Julia [ms] | Time Python [ms] |
|:----------------------------------|:-------------------:|:----------:|:---------------:|:---:|
|Falling mass (1)                   |     35              | 56         | 0.17            | 2.6 |
|Non-linear Spring damper (3)       |     49              | 83         | 0.64            | 20  |
|dito with callbacks (3b, 3c)       |     57              | 103        | 0.8             | 31  |
|swinging tether, 5 segments (5)    |    109              | 150        | 2.9             | 47  |
|Dyneema tether, reeling out (6)    |    125              | 160        | 4.3             | 9300 |
|dito with callbacks       (6c)     |    156              |            | 4.3             |      |
|Dyneema, reeling out with drag (7) |    175              |            | 3.3             |      |  

**Tradeoff Julia vs Python:** In Julia, the code is compiled before it is executed, which can cause about one to 10 seconds delay when running a simulation the first time, but speeds up the execution a lot afterward. In addition, Julia can run fully multithreaded, Python cannot make use of multiple CPU cores with multithreading because of the global interpreter lock. 

Furthermore, the IDA solver is hardly capable of handling a simulation with the very stiff
Dyneema tether. The Julia solvers achieve more than 2000 times the performance.

See also: [Why Julia?](https://ufechner7.github.io/2022/08/13/why-julia.html) and read the [documentation](https://ufechner7.github.io/Tethers.jl/dev/) or go straight to the [examples](https://ufechner7.github.io/Tethers.jl/dev/examples/).


