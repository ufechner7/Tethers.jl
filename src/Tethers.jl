module Tethers

using CondaPkg

export docu, display_if_interactive, run_python

# the re-usable, composable tether component of example 10
include("TetherComponent.jl")

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