# Tethers.jl
A series of example, from a simple falling mass towards a tether model, consisting of point masses connected by spring damper elements with support of reel-out and reel-in and aerodynamic
drag attached.

-- WORK IN PROGRESS --

## Installation

1. make sure you are running a `bash` terminal (shell) and you have at least 16GB RAM (MAC, Linux and Windows supported).  
   On Windows, you can use [git for windows](https://gitforwindows.org/) which provides git AND a bash shell, but for using Julia from a `bash` terminal you must also install either  [Installation and usage of VSCode](docs/vscode.md) or [Windows Terminal](https://learn.microsoft.com/en-us/windows/terminal/install). `Windows Terminal` is the simple and clean solution, `VSCode` the comfortable, fancy solution.
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

## Running the simulation
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

![Falling mass](docs/FallingMass.png)

This example shows a mass that is thrown upwards, slows down and then falls.

[Code](src/Tether_01.jl)

## Running the Python version as comparison
From the Julia prompt execute:
```
include("src/RunTether.jl")
```
This will install Python and Matplotlib and Assimulo and execute the script `Tether_01.py`.

[Code](src/Tether_01.py)

If you compare the Python and the Julia script you can see that:
- the Julia script is shorter and easier to read
- Julia is about 16 times faster when running the simulation

## More examples
### Mass, attached to a spring damper element
From the Julia prompt, run the simulation:
```julia
include("src/Tether_02.jl")
```
![Spring damper](docs/SpringDamper.png)

Mass, attached to a spring damper element. One end of the spring at the origin, the second end
attached to the mass. Mass initially below the origin, spring un-stretched. Z-axis pointing
upwards.

[Code](src/Tether_02.jl)

### Mass, attached to a non-linear spring damper element
```julia
include("src/Tether_03.jl")
```
![Non-linear Spring damper](docs/Nonlinear.png)

Mass, attached to a spring damper element. One end of the spring at the origin, the second end
attached to the mass. Mass initially below the origin, spring un-stretched. Z-axis pointing
upwards. 

Initial velocity $4 m/s$ upwards. The compression stiffness is zero. The grey line shows that
the stiffness is zero at the beginning, and has the nominal value at the end.

[Code](src/Tether_03.jl)
