var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/#Python","page":"References","title":"Python","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Assimulo offers 14 solvers with good documentation for explicit and implicit problems.","category":"page"},{"location":"references/#Julia","page":"References","title":"Julia","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"DifferentialEquations.jl offers a unified interface to about 300 different solvers from about a dozen different categories for a large range of problems. It wraps many existing open-source and commercial solvers, that have been implemented in C++ or Fortran and adds a growing number of native Julia solvers, many of them state-of-the-art.\nModelingToolkit.jl is an acausal modeling framework for automatically parallelized scientific machine learning (SciML) in Julia. A computer algebra system for integrated symbolics for physics-informed machine learning and automated transformations of differential equations. \nKiteModels.jl implements kite models, connected to a tether for airborne wind energy applications. It uses the same algorithms as this tutorial, but it is not (yet) using ModelingToolkit. \nWorking with Julia projects A must-read before creating your first project.","category":"page"},{"location":"references/#Scientific-papers","page":"References","title":"Scientific papers","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Uwe Fechner, Rolf van der Vlugt, Edwin Schreuder, Roland Schmehl. (2015). Dynamic Model of a Pumping Kite Power System  describes the tether model used in this tutorial, but also a model of a complete kite power system with experimental validation. Renewable Energy. Preprint.\nYingbo Ma, Shashi Gowda, Ranjan Anantharaman, Chris Laughman, Viral Shah, and Chris Rackauckas. (2021). ModelingToolkit: A Composable Graph Transformation System For Equation-Based Modeling.\nRackauckas, Christopher and Nie, Qing (2017). DifferentialEquations.jl–a performant and feature-rich ecosystem for solving differential equations in Julia} Journal of Open Research Software.\nD.F. Duda1, H. Fuest, T. Islam, T. Ostermann, D. Moormann1. (2022). Hybrid modeling approach for the tether of an airborne wind energy system CEAS Aeronautical Journal.","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"A sequence of examples, from a simple mass attached to a spring-damper element to a full segmented tether model with real-out and aerodynamic drag attached.","category":"page"},{"location":"examples/#Mass,-attached-to-a-spring-damper","page":"Examples","title":"Mass, attached to a spring-damper","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"From the Julia prompt, run the simulation:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"include(\"src/Tether_02.jl\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: Spring damper)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Mass, attached to a spring-damper element. One end of the spring is attached at the origin, the second end is attached to the mass. Mass initially below the origin, spring un-stretched. Z-axis pointing upwards.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Julia code: Tether_02.jl","category":"page"},{"location":"examples/#Mass,-with-non-linear-spring-damper","page":"Examples","title":"Mass, with non-linear spring damper","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"include(\"src/Tether_03.jl\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: Non-linear Spring damper)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Mass, attached to a non-linear spring-damper element. One end of the spring is attached at the origin, and the second end is attached to the mass. Mass initially below the origin, spring un-stretched. Z-axis pointing upwards. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Initial velocity 4 ms upwards. The compression stiffness is zero. The grey line shows that the stiffness is zero at the beginning, and has the nominal value at the end. Example: Tether_03.jl.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Thanks to the package ModelingToolkit.jl the system description is very compact and readable:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"D = Differential(t)\neqs = vcat(D.(pos)      ~ vel,\n           D.(vel)      ~ acc,\n           norm1        ~ norm(pos),\n           unit_vector  ~ -pos/norm1,         # direction from point mass to origin\n           spring_vel   ~ -unit_vector ⋅ vel,\n           c_spring     ~ c_spring0 * (norm1 > abs(l0)),\n           spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,\n           acc          ~ G_EARTH + spring_force/mass)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The same in Python: Python code: Tether_03.py. ","category":"page"},{"location":"examples/#Using-a-callback","page":"Examples","title":"Using a callback","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"By using a callback to detect exactly when the transition from a stiff tether segment to a loose tether segment happens we can increase the accuracy of the simulation. Julia code: Tether_03b.jl.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We only have to add the following lines of code:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0\n    norm(u[1:3]) - abs(L0)\nend\n`function affect!(integrator)\n    println(integrator.t)            # Not needed, just to show that the callback works\nend\ncb = ContinuousCallback(condition, affect!)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"and add the parameter callback = cb to the line that calls the solver:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts, callback = cb)","category":"page"},{"location":"examples/#Using-a-callback-with-Python","page":"Examples","title":"Using a callback with Python","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"In Python you would have to add the following attribute:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    sw0 = [vel_1[2] > 0] # array of booleans; true means the tether segment is loose (l < l_0)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"and the following methods:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    def state_events(self, t, y, yd, sw):\n        \"\"\"\n        This is our function that keeps track of our events. When the sign\n        of any of the events has changed, we have an event.\n        \"\"\"\n        # calculate the norm of the vector from mass1 to mass0 minus the initial segment length\n        event_0 = np.linalg.norm(y[3:6]) - L_0\n        return np.array([event_0])\n    \n    def handle_event(self, solver, event_info):\n        \"\"\"\n        Event handling. This functions is called when Assimulo finds an event as\n        specified by the event functions.\n        \"\"\"\n        state_info = event_info[0] # We are only interested in state events\n        if state_info[0] != 0:     # Check if the first event function has been triggered\n            if solver.sw[0]:       # If the switch is True the pendulum bounces\n                print(solver.t)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Example: Tether_03b.py.   As you can see, logging of calculated variables is not possible with Assimulo (easy with ModelingToolkit in Julia). You need to re-calculate them after the simulation.","category":"page"},{"location":"examples/#Multi-segment-tether","page":"Examples","title":"Multi-segment tether","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Using 2D arrays of variables allows to simulate a multi-segment tether:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@variables pos(t)[1:3, 1:segments+1]  = POS0\n@variables vel(t)[1:3, 1:segments+1]  = VEL0\n@variables acc(t)[1:3, 1:segments+1]  = ACC0","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this case, it is important to calculate the initial conditions of each particle such that they are physically feasible:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"G_EARTH     = Float64[0.0, 0.0, -9.81]          # gravitational acceleration     [m/s²]\nL0::Float64 = 10.0                              # initial segment length            [m]\nV0::Float64 = 4                                 # initial velocity of lowest mass [m/s]\nsegments::Int64 = 2                             # number of tether segments         [-]\nPOS0 = zeros(3, segments+1)\nVEL0 = zeros(3, segments+1)\nACC0 = zeros(3, segments+1)\nSEGMENTS0 = zeros(3, segments) \nUNIT_VECTORS0 = zeros(3, segments)\nfor i in 1:segments+1\n    POS0[:, i] .= [0.0, 0, -(i-1)*L0]\n    VEL0[:, i] .= [0.0, 0, (i-1)*V0/segments]\nend\nfor i in 2:segments+1\n    ACC0[:, i] .= G_EARTH\nend\nfor i in 1:segments\n    UNIT_VECTORS0[:, i] .= [0, 0, 1.0]\n    SEGMENTS0[:, i] .= POS0[:, i+1] - POS0[:, i]\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The first example of such a model is the script Tether_04.jl which is derived from the last example.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In the script Tether_05.jl, the spring force is distributed correctly on the two masses attached to the spring as shown here:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"if i == segments\n    eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i])\nelse\n    eqs2 = vcat(eqs2, total_force[:, i] ~ spring_force[:, i]- spring_force[:, i+1])\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We loop backward over the particles, starting with the last particle, because on the last particle, only one force is acting. On particle n-1 two spring forces are acting in the opposite direction.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Julia code: Tether_05.jl","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Python code: Tether_05.py","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Finally, in this example, we plot the result dynamically as 2D video. Screenshot:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: Tether 2D)","category":"page"},{"location":"examples/#Benchmarking","page":"Examples","title":"Benchmarking","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Using a callback slows the simulation down, but not much. Try it out:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"include(\"src/Tether_03c.jl\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Output on a fast PC:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Solving the system without callback...\n  0.000606 seconds (8.06 k allocations: 257.672 KiB)\nPress any key...\n\nSolving the system with callback...\n  0.000741 seconds (9.93 k allocations: 365.812 KiB)\nIf you zoom in to the points in time where pos_z crosses -10m\nyou should see a difference...","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this example, the gain of accuracy is very small, but that can be different in other simulations. For benchmarking we call solve twice: The first call ensures that the code is compiled, and the second call measures the execution time of the code.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Python The script, which executes the Python code with callbacks:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"include(\"src/RunTether_03b.jl\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"reports 31 ms for solving the problem (without printing). Without callbacks:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"include(\"src/RunTether_03.jl\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"still, 20 ms are needed.","category":"page"},{"location":"python/#Python-and-Julia-in-harmony","page":"Python and Julia","title":"Python and Julia in harmony","text":"","category":"section"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Python and Julia play very well together. You can see in the examples above that I am using Matplotlib for plotting, both in Python and in Julia. Julia has a built-in package manager. You can use it to install and remove Julia packages, but also to install or remove Python packages. That works like this:","category":"page"},{"location":"python/#Using-Python-packages-from-Julia","page":"Python and Julia","title":"Using Python packages from Julia","text":"","category":"section"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"There are three options:","category":"page"},{"location":"python/#Option-one:","page":"Python and Julia","title":"Option one:","text":"","category":"section"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Use wrapper libraries that exist for a few, very popular Python packages, e.g. PyPlot.jl for Matplotlib or SymPy.jl for SymPy. You can install them like any other Julia package, e.g.","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"]     # enter package manager mode\nadd SymPy\n<DEL> # leave the package manager","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"and on the Julia prompt:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"using SymPy","category":"page"},{"location":"python/#Option-two:","page":"Python and Julia","title":"Option two:","text":"","category":"section"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Use PyCall to use Python packages for Julia. This works for all Python packages, but  is a little bit less comfortable than option one. Example:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"using PyCall\nnp = pyimport(\"numpy\")","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Now you can use NumPy from Julia:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"julia> np.zeros(3)\n3-element Vector{Float64}:\n 0.0\n 0.0\n 0.0","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"If the package is not yet installed, you can use the notation:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"as = pyimport_conda(\"assimulo\", \"assimulo\")","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"If the command using PyCall should fail, you can execute:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"cd bin\n./build_pycall","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"to re-build it.","category":"page"},{"location":"python/#Option-three:","page":"Python and Julia","title":"Option three:","text":"","category":"section"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Using the package PythonCall. You cannot use it together with PyCall, it is the newer successor of PyCall, and it is symmetric, you can use it to call Julia from Python or Python from Julia.","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"We need to create a new project to try it out:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"cd repos # or any other folder that you use for your projects\nmkdir PythonDemo\ncd PythonDemo\njulia --project=\".\" # this creates a new, empty project","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Then enter at the Julia prompt: Example of using Python from Julia:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"re = pyimport(\"re\")   # import the re module\nwords = re.findall(\"[a-zA-Z]+\", \"PythonCall.jl is very useful!\")","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Output:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Python: ['PythonCall', 'jl', 'is', 'very', 'useful']","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Type:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"sentence = Py(\" \").join(words)","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Output:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Python: 'PythonCall jl is very useful'","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"And finally, convert this Python object to a Julia string:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"pyconvert(String, sentence)  # convert the Python string to a Julia string","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"For adding Python packages that you want to use with PythonCall use CondaPkg as explained in the next section.","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Sometimes needed: Install CondaPkg","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"] # by pressing the closing square bracket you enter the package manager mode of Julia\nadd CondaPkg # add the Python package manger\n","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Now press the \\<DEL\\> key to leave the package manager. In the Julia REPL, type:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"using CondaPkg\n]                 # enter the package manager mode\nhelp              # will show you all available commands; try for example\nconda add ipython # this will add ipython","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Press \\<DEL\\> to leave the package manager mode. In the Julia REPL, type:","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"CondaPkg.withenv() do\n    run(`ipython`)\nend","category":"page"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"You should now get an interactive Python prompt and can program in Python.","category":"page"},{"location":"python/#Further-reading","page":"Python and Julia","title":"Further reading","text":"","category":"section"},{"location":"python/","page":"Python and Julia","title":"Python and Julia","text":"Noteworthy differences Julia/Python Good to know.\nPythonCall.jl New library to call Python from Julia or Julia from Python.\nPyCall Old library to call Python from Julia.","category":"page"},{"location":"#Tethers.jl","page":"Readme","title":"Tethers.jl","text":"","category":"section"},{"location":"","page":"Readme","title":"Readme","text":"(Image: Dev)","category":"page"},{"location":"#Introduction","page":"Readme","title":"Introduction","text":"","category":"section"},{"location":"","page":"Readme","title":"Readme","text":"A series of examples, from a simple falling mass towards a tether model, consisting of point masses connected by spring damper elements with the support of reel-out and reel-in and aerodynamic drag attached.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"– WORK IN PROGRESS –","category":"page"},{"location":"#Installation","page":"Readme","title":"Installation","text":"","category":"section"},{"location":"","page":"Readme","title":"Readme","text":"make sure you are running a bash terminal (shell) and you have at least 16GB RAM (MAC, Linux and Windows supported).   On Windows, you can use git for windows which provides git AND a bash shell, but for using Julia from a bash terminal you must also install either  Installation and usage of VSCode or Windows Terminal. Windows Terminal is the simple and clean solution, VSCode the comfortable, fancy solution.\ninstall Julia 1.10 using juliaup, see https://github.com/JuliaLang/juliaup. If juliaup is already installed, the following commands will do:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"juliaup add 1.10\njuliaup default 1.10 ","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Check out from git:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"cd repos # any folder of your choice, but without spaces in the folder name\ngit clone https://github.com/ufechner7/Tethers.jl","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Build the system image:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"cd repos/Tethers.jl\ncd bin\n./create_sys_image","category":"page"},{"location":"#Basic-example","page":"Readme","title":"Basic example","text":"","category":"section"},{"location":"","page":"Readme","title":"Readme","text":"Use the provided script to start Julia from the Tethers.jl folder:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"cd repos/Tethers.jl\n./bin/run_julia","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"From the Julia prompt, run the simulation:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"include(\"src/Tether_01.jl\")","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"You should see a plot similar to:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"(Image: Falling mass)","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"This example shows a mass that is thrown upwards, slows down and then falls.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Julia code: Tether_01.jl","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"These differential equations define the model:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"D = Differential(t)\n\neqs = vcat(D.(pos) ~ vel,\n           D.(vel) ~ acc,\n           acc    .~ G_EARTH)","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"The term D.(pos) means \"Apply the differential D(t) to all elements of the vector pos\". The second term defines that the differential of the velocity vector must be equal to the  acceleration. For equality in symbolic equations the character ~ has to be used, because the character = has the meaning \"assign a value to a variable\" which is not what we are doing here. The third equation means that all elements of the acceleration vector must be equal to the elements of the gravity vector. We end up with an array of 3x3` equations.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"The next lines are:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"@named sys = ODESystem(eqs, t)\nsimple_sys = structural_simplify(sys)","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"This means we create a named ordinary equation system, depending on t. Then we simplify the system symbolically (order reduction). If you type sys in the Julia REPL (command line) you can see that the original system had 9 equations, the second line above created a system with only six equations. This step helps to speed up the simulation and often also removes algebraic loops which makes the ODE a lot simpler to solve numerically later on.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Now the parameters of the integrator are defined:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"duration = 10.0\ndt = 0.02\ntol = 1e-6\ntspan = (0.0, duration)\nts    = 0:dt:duration","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"The time step dt is the interval in which the solution shall be stored, NOT the time step of the integrator. The integrator uses a variable time step which can be much smaller or much larger as determined by the required tolerance, in this example set to tol=10^-6. The variable ts is a range object defining the sampling times for the result.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Then we define the initial condition u0. We use a dictionary of variable => value pairs to do this. In the next line, we define the ODE problem and finally, we solve it using the Rodas5 solver with the given parameters.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"u0 = Dict(vel=>[0, 0, 50.0])\n\nprob = ODEProblem(simple_sys, u0, tspan)\n@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"The macro @time measures the compilation and execution time of calling the function solve(). It is compiled only when called the first time. ","category":"page"},{"location":"#Python-version-as-comparison","page":"Readme","title":"Python version as comparison","text":"","category":"section"},{"location":"","page":"Readme","title":"Readme","text":"From the Julia prompt execute:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"include(\"src/RunTether.jl\")","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"This will install Python, Matplotlib and Assimulo and execute the script Tether_01.py.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Python code: Tether_01.py","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"If you compare the Python and the Julia scripts you can see that:","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"the Julia script is shorter and easier to read\nJulia is about 16 times faster when running the simulation","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Have a look at the Examples that teach you how to construct a full tether model step by step.","category":"page"},{"location":"#Overall-comparison","page":"Readme","title":"Overall comparison","text":"","category":"section"},{"location":"","page":"Readme","title":"Readme","text":"Execution time for a simulation of 10s duration with logging the state every 20ms. Relative and absolute tolerance: 10^-6. CPU: Ryzen 9 7950X.","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Test-case Lines of code (LOC) Julia LOC Python Time Julia [ms] Time Python [ms]\nFalling mass 42 56 0.17 2.6\nNon-linear Spring damper 61 83 0.61 20\ndito with callbacks 68 103 0.74 31\nswinging tether, 5 segments 121 148 3.5 47","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"Tradeoff Julia vs Python: In Julia, the code is compiled before it is executed, which can cause about one to 10 seconds delay when running a simulation the first time, but speeds up the execution a lot afterward. In addition, Julia can run fully multithreaded, Python cannot make use of multiple CPU cores with multithreading because of the global interpreter lock. ","category":"page"},{"location":"","page":"Readme","title":"Readme","text":"See also: Why Julia? and read the documentation.","category":"page"},{"location":"vscode/#Installation-and-usage-of-VSCode","page":"VSCode IDE","title":"Installation and usage of VSCode","text":"","category":"section"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"It is useful to install the integrated development environment VSCode, even though it is not required. You can also use any editor of your choice.","category":"page"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"VSCode provides syntax highlighting, but also the feature \"goto definition\" which can help to understand and explore the code. (Image: VSCode)","category":"page"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"You can download and install VSCode for all operating systems here.","category":"page"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"Julia development with VSCode is well documented here.","category":"page"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"I suggest to install the following Extensions:","category":"page"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"Julia\nProject Manager \nYAML\nEven Better TOML","category":"page"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"VSCode has good, integrated GIT support by default, no extension is needed for that.","category":"page"},{"location":"vscode/","page":"VSCode IDE","title":"VSCode IDE","text":"I would NOT use all the advanced features of julia-vscode, I prefer to just use the vscode terminal and launch Julia from the terminal. This makes it easy to launch Julia with any command line options and also to start and restart Julia quickly.","category":"page"}]
}
