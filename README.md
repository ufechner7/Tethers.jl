# Tethers.jl
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ufechner7.github.io/Tethers.jl/dev)

## Introduction

A series of examples, from a simple falling mass towards a tether model, consisting of point masses connected by spring damper elements with support of reel-out and reel-in and aerodynamic drag attached.

-- WORK IN PROGRESS --

## Installation

1. make sure you are running a `bash` terminal (shell) and you have at least 16GB RAM (MAC, Linux and Windows supported).  
   On Windows, you can use [git for windows](https://gitforwindows.org/) which provides git AND a bash shell, but for using Julia from a `bash` terminal you must also install either  [Installation and usage of VSCode](docs/images/vscode.md) or [Windows Terminal](https://learn.microsoft.com/en-us/windows/terminal/install). `Windows Terminal` is the simple and clean solution, `VSCode` the comfortable, fancy solution.
2. install Julia 1.10 using `juliaup`, see [https://github.com/JuliaLang/juliaup](https://github.com/JuliaLang/juliaup). If `juliaup` is already installed, the following commands will do:
```
juliaup add 1.10
juliaup default 1.10 
```

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

**Julia code:** [Tether_01.jl](src/Tether_01.jl)

These differential equations define the model:
```Julia
D = Differential(t)

eqs = vcat(D.(pos) ~ vel,
           D.(vel) ~ acc,
           acc    .~ G_EARTH)
```
The term `D.(pos)` means "Apply the differential D(t) to all elements of the vector `pos`".
The second term defines that the differential of the velocity vector must be equal to the 
acceleration. For equality in symbolic equations the character `~` has to be used, because
the character `=` has the meaning "assign a value to a variable" which is not what we are doing here. The third equation means that all elements of the acceleration vector must be equal
to the elements of the gravity vector. We end up with an array of `3x3`` equations.

The next lines are:
```julia
@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)
```
This means, we create a named ordinary equation system, depending on `t`. Then we simplify
the system symbolically (order reduction). If you type `sys` in the Julia REPL (command line)
you can see that the original system had 9 equations, the second line above created a system
with ony six equations. This step helps to speed up the simulation and often also removes
algebraic loops which makes the ODE a lot simple to solve numerically later on.

## Python version as comparison
From the Julia prompt execute:
```
include("src/RunTether.jl")
```
This will install Python, Matplotlib and Assimulo and execute the script `Tether_01.py`.

**Python code:** [Tether_01.py](src/Tether_01.py)

If you compare the Python and the Julia scripts you can see that:
- the Julia script is shorter and easier to read
- Julia is about 16 times faster when running the simulation

Have a look at the [Examples] that teach you how to construct a full tether model step by step.

## Overall comparison
Execution time for a simulation of 10s duration with logging the state every 20ms.

| Testcase                    | Lines of code (LOC) Julia | LOC Python  | Time Julia [ms] | Time Python [ms] |
|:----------------------------|:-------------------:|:---:|:-----:|:---:|
|Falling mass                 |     42              | 56  | 0.17  | 2.6 |
|Non-linear Spring damper     |     61              | 83  | 0.61  | 20  |
|dito with callbacks          |     68              | 103 | 0.74  | 31  |
|swinging tether, 5 segments  |    117              | 190 | 2.90  |     |

**Tradeoff Julia vs Python:** In Julia, the code is compiled before it is executed, which can cause about one to 10 seconds delay when running a simulation the first time, but speeds up the execution a lot afterward. In addition, Julia can run fully multithreaded, Python cannot make use of multiple CPU cores with multithreading because of the global interpreter lock. 

See also: [Why Julia?](https://ufechner7.github.io/2022/08/13/why-julia.html) and read the [documentation](https://ufechner7.github.io/Tethers.jl/dev/).


