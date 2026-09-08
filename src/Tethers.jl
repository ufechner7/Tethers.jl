module Tethers

using CondaPkg, Pkg

export docu, display_if_interactive, install_examples, run_python, analytic_force,
       hooke_force, damping_factor

# the re-usable, composable tether component of example 10
include("TetherComponent.jl")

# the analytical force formula of PlanCompression.md, shared by plot_compression.jl and
# test_compression.jl
include("analytic_force.jl")

LAUNCH_BROWSER = true

"""
    display_if_interactive(x)
    display_if_interactive(f, args...; kwargs...)

In the first form, display `x` (e.g. a `PlotX` figure) only when running in
an interactive session and not on CI.

In the second form, call `f(args...; kwargs...)` only when running in an
interactive session and not on CI. Use this for functions like `plot2d`
that open a window themselves and return nothing.
"""
function display_if_interactive(x)
    if isinteractive() && get(ENV, "CI", "false") == "false"
        display(x)
    end
    nothing
end

function display_if_interactive(f::Base.Callable, args...; kwargs...)
    if isinteractive() && get(ENV, "CI", "false") == "false"
        f(args...; kwargs...)
    end
    nothing
end

"""
    run_python(name)

Run the Python version of the example `name`, e.g. `run_python("Tether_01")` for the
script `examples/python/Tether_01.py`, in the Python environment managed by CondaPkg.

The script is located relative to this package, but writes its results to `output/`
relative to the current working directory, so run it from the package directory.
"""
function run_python(name)
    script = joinpath(dirname(@__DIR__), "examples", "python", "$name.py")
    isfile(script) || throw(ArgumentError("no such Python example: $script"))
    CondaPkg.withenv() do
        Base.run(`python $script`)
    end
    nothing
end

"""
    copy_examples(; overwrite=true)

Copy all example scripts (Julia and Python) to the folder `examples`
(it will be created if it doesn't exist). The `examples/Project.toml` of this
package is not copied, since it points back at this package via a relative
`[sources]` path that would not resolve in the destination. Any local
`Manifest.toml` / `Manifest-v*.toml` left over from instantiating that
Project.toml is skipped as well.
"""
function copy_examples(; overwrite=true)
    PATH = "examples"
    if ! isdir(PATH)
        mkdir(PATH)
    end
    src_path = joinpath(dirname(@__DIR__), PATH)
    files = filter(readdir(src_path)) do file
        file != "Project.toml" && !startswith(file, "Manifest")
    end
    copy_files(PATH, files; overwrite)
end

"""
    copy_bin(; overwrite=true)

Copy the helper scripts from the folder `bin` (e.g. `run_julia`, `install` and
`create_sys_image`) to the folder `bin` in the current working directory
(it will be created if it doesn't exist). Pre-built system images (`*.so`) are
not copied.
"""
function copy_bin(; overwrite=true)
    PATH = "bin"
    if ! isdir(PATH)
        mkdir(PATH)
    end
    src_path = joinpath(dirname(@__DIR__), PATH)
    files = filter(file -> !endswith(file, ".so"), readdir(src_path))
    copy_files(PATH, files; overwrite)
end

function copy_files(relpath, files; overwrite=true)
    if ! isdir(relpath)
        mkdir(relpath)
    end
    src_path = joinpath(dirname(@__DIR__), relpath)
    for file in files
        src = joinpath(src_path, file)
        dst = joinpath(relpath, file)
        if overwrite || !isfile(dst)
            cp(src, dst, force=true)
            chmod(dst, 0o774)
        end
    end
    files
end

"""
    example_packages()

Return the names of the registered packages required to run the example scripts, read
from `examples/Project.toml`. `Tethers` itself and standard library packages (currently
only `LinearAlgebra`) are excluded, since they do not need to be added.
"""
function example_packages()
    project = Pkg.TOML.parsefile(joinpath(dirname(@__DIR__), "examples", "Project.toml"))
    filter(!in(("Tethers", "LinearAlgebra")), collect(keys(project["deps"])))
end

"""
    install_examples(add_packages=true)

Install the example scripts into the current working directory.

This copies the `examples` and `bin` folders. If `add_packages` is `true`, it also
installs the packages used by the example scripts. Run the examples with:
```julia
include("examples/menu.jl")
```
"""
function install_examples(add_packages=true)
    copy_examples()
    copy_bin()
    if add_packages
        Pkg.add(example_packages())
    end
end

function docu(build=true)
    @eval using LiveServer
    # if build
    #     include("docs/make.jl")
    # end
    if Sys.islinux() && ! LAUNCH_BROWSER
        Base.run(`xdg-open "docs/build/index.html"`; wait=false)
    else
        Base.invokelatest(LiveServer.servedocs; skip_dir="docs", launch_browser=true)
    end
end

end