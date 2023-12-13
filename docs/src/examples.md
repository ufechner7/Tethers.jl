# Examples
A sequence of examples, from a simple mass attached to a spring-damper element to
a full segmented tether model with real-out and aerodynamic drag attached.

## Mass, attached to a spring-damper
From the Julia prompt, run the simulation:
```julia
include("src/Tether_02.jl")
```
![Spring damper](docs/images/SpringDamper.png)

Mass, attached to a spring-damper element. One end of the spring is attached at the origin, the second end is attached to the mass. Mass initially below the origin, spring un-stretched. Z-axis pointing
upwards.

**Julia code:** [Tether_02.jl](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_02.jl)

## Mass, with non-linear spring damper
```julia
include("src/Tether_03.jl")
```
![Non-linear Spring damper](docs/images/Nonlinear.png)

Mass, attached to a spring-damper element. One end of the spring is attached at the origin, and the second end is
attached to the mass. Mass initially below the origin, spring un-stretched. Z-axis pointing
upwards. 

Initial velocity $4 m/s$ upwards. The compression stiffness is zero. The grey line shows that
the stiffness is zero at the beginning, and has the nominal value at the end. **Example:** [Tether_03.jl](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_03.jl).

Thanks to the package [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) the system description is very compact and readable:
```Julia
D = Differential(t)
eqs = vcat(D.(pos)      ~ vel,
           D.(vel)      ~ acc,
           norm1        ~ norm(pos),
           unit_vector  ~ -pos/norm1,         # direction from point mass to origin
           spring_vel   ~ -unit_vector ⋅ vel,
           c_spring     ~ c_spring0 * (norm1 > abs(l0)),
           spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
           acc          ~ G_EARTH + spring_force/mass)
```

The same in Python: **Python code:** [Tether_03.py](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_03.py). 

### Using a callback
By using a callback to detect exactly when the transition from a stiff tether segment to a loose
tether segment happens we can increase the accuracy of the simulation. **Julia code:** [Tether_03b.jl](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_03b.jl).

We only have to add the following lines of code:
```julia
function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    norm(u[1:3]) - abs(L0)
end
`function affect!(integrator)
    println(integrator.t)            # Not needed, just to show that the callback works
end
cb = ContinuousCallback(condition, affect!)
```
and add the parameter `callback = cb` to the line that calls the solver:
```julia
sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts, callback = cb)
```

### Using a callback with Python
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
**Example:** [Tether_03b.py](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_03b.py).  
As you can see, logging of calculated variables is not
possible with Assimulo (easy with ModelingToolkit in Julia). You need to re-calculate them
after the simulation.

## Multi-segment tether
Using 2D arrays of variables allows to simulate a multi-segment tether:
```julia
@variables pos(t)[1:3, 1:segments+1]  = POS0
@variables vel(t)[1:3, 1:segments+1]  = VEL0
@variables acc(t)[1:3, 1:segments+1]  = ACC0
```
In this case, it is important to calculate the initial conditions of each particle such that they are physically feasible:
```julia
G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration     [m/s²]
L0::Float64 = 10.0                              # initial segment length            [m]
V0::Float64 = 4                                 # initial velocity of lowest mass [m/s]
segments::Int64 = 2                             # number of tether segments         [-]
POS0 = zeros(3, segments+1)
VEL0 = zeros(3, segments+1)
ACC0 = zeros(3, segments+1)
SEGMENTS0 = zeros(3, segments) 
UNIT_VECTORS0 = zeros(3, segments)
for i in 1:segments+1
    POS0[:, i] .= [0.0, 0, -(i-1)*L0]
    VEL0[:, i] .= [0.0, 0, (i-1)*V0/segments]
end
for i in 2:segments+1
    ACC0[:, i] .= G_EARTH
end
for i in 1:segments
    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]
    SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]
end
```
The first example of such a model is the script [Tether_04.jl](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_04.jl) which is derived from the last example.

In the script [Tether_05.jl](https://github.com/ufechner7/Tethers.jl/blob/main/src/Tether_05.jl), the spring force is distributed correctly on the two masses attached to the spring as shown here:
```julia
if i == segments
    eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i])
else
    eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i]- spring_force[:, i+1])
end
```
We loop backward over the particles, starting with the last particle, because on the last particle, only one force is acting. On particle $n-1$ two spring forces are acting in the opposite direction.

Finally, in this example, we plot the result dynamically as 2D video. Screenshot:

![Tether 2D](docs/images/Tether2d.png)


## Benchmarking
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
In this example, the gain of accuracy is very small, but that can be different
in other simulations. For benchmarking we call solve twice: The first call ensures that the
code is compiled, and the second call measures the execution time of the code.

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
still, 20 ms are needed.