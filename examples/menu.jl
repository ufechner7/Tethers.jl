using Pkg
# When running from a clone of this repository, activate the examples environment unless it
# is already active. If the examples were copied into another project with
# `install_examples()`, there is no `Project.toml` next to this script; in that case the
# environment of that project is already active and is left unchanged.
if isfile(joinpath(@__DIR__, "Project.toml")) && dirname(Pkg.project().path) != @__DIR__
    Pkg.activate(@__DIR__)
end
using REPL.TerminalMenus

examples = [("Tether_01",  "include(\"Tether_01.jl\")",  "Falling mass thrown upwards"),
            ("Tether_02",  "include(\"Tether_02.jl\")",  "Mass on a linear spring-damper"),
            ("Tether_03",  "include(\"Tether_03.jl\")",  "Mass on a non-linear spring-damper"),
            ("Tether_03b", "include(\"Tether_03b.jl\")", "Non-linear spring with callback"),
            ("Tether_03c", "include(\"Tether_03c.jl\")", "Benchmark: callback vs no callback"),
            ("Tether_04",  "include(\"Tether_04.jl\")",  "Multi-segment tether (2D arrays)"),
            ("Tether_05",  "include(\"Tether_05.jl\")",  "Segmented tether, correct force split"),
            ("Tether_06",  "include(\"Tether_06.jl\")",  "Segmented tether, reeling out"),
            ("Tether_06b", "include(\"Tether_06b.jl\")", "Reel-out tether, refactored with Settings"),
            ("Tether_06c", "include(\"Tether_06c.jl\")", "Reel-out tether with continuous callback"),
            ("Tether_07",  "include(\"Tether_07.jl\")",  "Segmented tether with aerodynamic drag"),
            ("Tether_08",  "include(\"Tether_08.jl\")",  "Tether with arbitrary/free endpoints"),
            ("Tether_09",  "include(\"Tether_09.jl\")",  "Labeled tether shape diagram for docs"),
            ("Tether_10",  "include(\"Tether_10.jl\")",  "Re-usable tether component with two end points"),
            ("Tether_11",  "include(\"Tether_11.jl\")",  "Kite flying a circular trajectory, dynamic tether model")]

name_width = maximum(length(name) for (name, _, _) in examples)
options = [rpad(name, name_width) * "  " * descr for (name, _, descr) in examples]
push!(options, "quit()")

function run_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(examples[choice][2]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end
run_menu()
