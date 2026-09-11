module Tethers

using CondaPkg, Pkg

export docu, display_if_interactive, install_examples, run_python

# the re-usable, composable tether component of example 10
include("TetherComponent.jl")

# the quasi-steady tether model
include("Tether_quasisteady.jl")

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
    copy_file(relpath, src_file, dst_file; overwrite=true)

Copy `src_file` from the folder `relpath` of this package to `dst_file` in the folder
`relpath` of the current working directory (it will be created if it doesn't exist).
Used for the scripts that exist in a variant for a clone of this repository and in a
variant for an installed package.
"""
function copy_file(relpath, src_file, dst_file; overwrite=true)
    if ! isdir(relpath)
        mkdir(relpath)
    end
    dst = joinpath(relpath, dst_file)
    if overwrite || !isfile(dst)
        cp(joinpath(dirname(@__DIR__), relpath, src_file), dst, force=true)
        chmod(dst, 0o774)
    end
    dst
end

"""
    copy_bin(; overwrite=true)

Copy the scripts needed to run the examples of an installed `Tethers` package to the
folders `bin` and `test` in the current working directory (they will be created if they
don't exist):

- `bin/run_julia`            starts Julia with the settings used for the examples
- `bin/create_sys_image`     builds a system image for a faster start of the examples
- `test/create_sys_image.jl` and `test/test_for_precompile.jl`, which it uses

The two scripts that build the system image exist in a variant for a clone of this
repository and in a variant for an installed package; the latter have a `2` in their name
and are copied without it. The remaining scripts in `bin` are specific to a clone of this
repository and are not copied. Pre-built system images (`*.so`) are not copied either.
"""
function copy_bin(; overwrite=true)
    copy_files("bin", ["run_julia"]; overwrite)
    copy_file("bin", "create_sys_image2", "create_sys_image"; overwrite)
    copy_files("test", ["test_for_precompile.jl"]; overwrite)
    copy_file("test", "create_sys_image2.jl", "create_sys_image.jl"; overwrite)
    [joinpath("bin", "run_julia"), joinpath("bin", "create_sys_image"),
     joinpath("test", "test_for_precompile.jl"), joinpath("test", "create_sys_image.jl")]
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