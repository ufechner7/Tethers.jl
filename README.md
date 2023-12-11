# Tethers.jl
A series of examples, from a simple falling mass towards a tether model, consisting of point masses connected by spring damper elements with support of reel-out and reel-in and aerodynamic
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

**Julia code:** [Tether_01.jl](src/Tether_01.jl)

## Running the Python version as comparison
From the Julia prompt execute:
```
include("src/RunTether.jl")
```
This will install Python and Matplotlib and Assimulo and execute the script `Tether_01.py`.

**Python code:** [Tether_01.py](src/Tether_01.py)

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

**Julia code:** [Tether_02.jl](src/Tether_02.jl)

### Mass, attached to a non-linear spring damper element
```julia
include("src/Tether_03.jl")
```
![Non-linear Spring damper](docs/Nonlinear.png)

Mass, attached to a spring damper element. One end of the spring at the origin, the second end
attached to the mass. Mass initially below the origin, spring un-stretched. Z-axis pointing
upwards. 

Initial velocity $4 m/s$ upwards. The compression stiffness is zero. The grey line shows that
the stiffness is zero at the beginning, and has the nominal value at the end. **Example:** [Tether_03.jl](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_03.jl).

The same as Python version: **Python code:** [Tether_03.py](src/Tether_03.py). 

#### Using a callback
By using a callback to detect exactly when the transition from a stiff tether segment to a loose
tether segment happens we can increase the accuracy of the simulation. **Julia code:** [Tether_03b.jl](src/Tether_03b.jl).

We only have to add the following lines of code:
```julia
function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    norm(u[1:3]) - abs(L0)
end
function affect!(integrator)
    println(integrator.t)            # Not needed, just to show that the callback works
end
cb = ContinuousCallback(condition, affect!)
```
and add the parameter `callback = cb` to the line that calls the solver:
```julia
sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts, callback = cb)
```

#### Using a callback with Python
In Python you would have to add the following attribute:
```Python
    sw0 = [vel_1[2] > 0] # array of booleans; true means the tether segment is loose (l < l_0)
```
and the following methods:
```Python
    def state_events(self, t, y, yd, sw):
        """
        This is our function that keeps track of our events. When the sign
        of any of the events has changed, we have an event.
        """
        # calculate the norm of the vector from mass1 to mass0 minus the initial segment length
        event_0 = np.linalg.norm(y[3:6]) - L_0
        return np.array([event_0])
    
    def handle_event(self, solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        state_info = event_info[0] # We are only interested in state events
        if state_info[0] != 0:     # Check if the first event function has been triggered
            if solver.sw[0]:       # If the switch is True the pendulum bounces
                print(solver.t)
```
**Example:** [Tether_03b.py](src/Tether_03b.py).  
As you can see, logging of calculated variables is not
possible with Assimulo (easy with ModelingToolkit in Julia). You need to re-calculate them
after the simulation.

#### Benchmarking
Using a callback slows the simulation down, but not much. Try it out:
```julia
include("src/Tether_03c.jl")
```
Output on a fast PC:
```
Solving the system without callback...
  0.000606 seconds (8.06 k allocations: 257.672 KiB)
Press any key...

Solving the system with callback...
  0.000741 seconds (9.93 k allocations: 365.812 KiB)
If you zoom in to the points in time where pos_z crosses -10m
you should see a difference...
```
In this example the gain of accuracy is very small, but that can be different
in other simulations. For benchmarking we call solve twice: The first call ensures that the
code is compiled, the second call measures the execution time of the code.

**Python**
The script, which executes the Python code with callbacks:
```
include("src/RunTether_03b.jl")
```
reports 31 ms for solving the problem (without printing).
Without callbacks:
```
include("src/RunTether_03.jl")
```
still 20 ms are needed.

## Comparism
| Testcase                    | Lines of code (LOC) Julia | LOC Python  | Time Julia [ms] | Time Python [ms] |
|:----------------------------|:-------------------:|:---:|:-:|:---:|
|Falling mass                 |     42              | 56  | 0.17  | 2.6  |
|Non-linear Spring damper     |     61              | 83  | 0.61  | 20  |
|dito with callbacks          |     68              | 103 | 0.74 | 31  |

**Tradeoff Julia vs Python:** In Julia the code is compiled before it is executed, that can cause about 1 to 10 seconds delay when running a simulation the first time, but speeds up the execution a lot afterwards. In addition Julia can run fully multithreaded, Python cannot really use threads because of the global interpreter lock. 