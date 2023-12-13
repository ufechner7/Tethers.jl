## Python and Julia in harmony

Python and Julia play very well together. You could see in the examples above that I am using Matplotlib for plotting, both in Python and in Julia. Julia has a build-in package manager. You can use it install and remove Julia packages, but also to install or remove Python packages. That works like this:

### Using Python packages from Julia
There are three options:

#### Option one: 
Use wrapper libraries which exist for a few, very popular Python packages, e.g.
[PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) for Matplotlib or [SymPy.jl](https://github.com/JuliaPy/SymPy.jl) for SymPy. You can install them like any other Julia package, e.g.
```
]     # enter package manager mode
add SymPy
<DEL> # leave the package manager
```
and on the Julia prompt:
```
using SymPy
```

#### Option two:

Use PyCall to use Python packages for Julia. This works for all Python packages, but 
is a little bit less comfortable than option one. Example:

```
using PyCall
np = pyimport("numpy")
```
Now you can use NumPy from Julia:
```julia
julia> np.zeros(3)
3-element Vector{Float64}:
 0.0
 0.0
 0.0
```
If the package is not yet installed, you can use the notation:
```
as = pyimport_conda("assimulo", "assimulo")
```
If the command `using PyCall` should fail, you can execute:
```
cd bin
./build_pycall
```
to re-build it.

#### Option three:

Using the package [PythonCall](https://github.com/JuliaPy/PythonCall.jl).
You cannot use it together with `PyCall`, it is the newer successor of `PyCall`, and it is
symmetric, you can use it to call Julia from Python or Python from Julia.

We need to create a new project to try it out:
```
cd repos # or any other folder that you use for your projects
mkdir PythonDemo
cd PythonDemo
julia --project="." # this creates a new, empty project
```
Then enter at the Julia prompt:
```
]               # enter the package manger mode
add PythonCall
<BACK>          # leave the package manager mode
```
Example for using Python from Julia:
```julia
re = pyimport("re")   # import the re module
words = re.findall("[a-zA-Z]+", "PythonCall.jl is very useful!")
```
Output:
```
Python: ['PythonCall', 'jl', 'is', 'very', 'useful']
```
Type:
```julia
sentence = Py(" ").join(words)
```
Output:
```
Python: 'PythonCall jl is very useful'
```
And finally convert this Python object to a Julia string:
```julia
pyconvert(String, sentence)  # convert the Python string to a Julia string
```
For adding Python packages that you want to use with PythonCall use CondaPkg as explained
in the next section.

**Sometimes needed:** Install CondaPkg
```
] # by pressing the closing square bracket you enter the package manager mode of Julia
add CondaPkg # add the Python package manger

```
Now press the \<DEL\> key to leave the package manager.
In the Julia REPL, type:
```
using CondaPkg
]                 # enter the package manager mode
help              # will show you all available commands; try for example
conda add ipython # this will add ipython
``` 
Press \<DEL\> to leave the package manager mode.
In the Julia REPL, type:
```
CondaPkg.withenv() do
    run(`ipython`)
end
```
You should now get an interactive Python prompt and can program in Python.