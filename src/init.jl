println("Running init...")

using MakieControlPlots
using REPL.TerminalMenus
using Tethers: run_python

examples = [("Tether_01",  "include(\"examples/Tether_01.jl\")",  "Falling mass thrown upwards"),
            ("Tether_02",  "include(\"examples/Tether_02.jl\")",  "Mass on a linear spring-damper"),
            ("Tether_03",  "include(\"examples/Tether_03.jl\")",  "Mass on a non-linear spring-damper"),
            ("Tether_03b", "include(\"examples/Tether_03b.jl\")", "Non-linear spring with callback"),
            ("Tether_03c", "include(\"examples/Tether_03c.jl\")", "Benchmark: callback vs no callback"),
            ("Tether_04",  "include(\"examples/Tether_04.jl\")",  "Multi-segment tether (2D arrays)"),
            ("Tether_05",  "include(\"examples/Tether_05.jl\")",  "Segmented tether, correct force split"),
            ("Tether_06",  "include(\"examples/Tether_06.jl\")",  "Segmented tether, reeling out"),
            ("Tether_06b", "include(\"examples/Tether_06b.jl\")", "Reel-out tether, refactored with Settings"),
            ("Tether_06c", "include(\"examples/Tether_06c.jl\")", "Reel-out tether with continuous callback"),
            ("Tether_07",  "include(\"examples/Tether_07.jl\")",  "Segmented tether with aerodynamic drag"),
            ("Tether_08",  "include(\"examples/Tether_08.jl\")",  "Tether with arbitrary/free endpoints"),
            ("Tether_09",  "include(\"examples/Tether_09.jl\")",  "Labeled tether shape diagram for docs"),
            ("Tether_10",  "include(\"examples/Tether_10.jl\")",  "Re-usable tether component with two end points")]

name_width = maximum(length(name) for (name, _, _) in examples)
options = [rpad(name, name_width) * "  " * descr for (name, _, descr) in examples]
push!(options, "quit()")

# the Python examples are all run the same way, with `run_python`, so only the name of the
# script and its description are needed here
python_examples = [("Tether_01",  "Falling mass thrown upwards"),
                    ("Tether_02",  "Mass on a linear spring-damper"),
                    ("Tether_03",  "Mass on a non-linear spring-damper"),
                    ("Tether_03b", "Non-linear spring with callback"),
                    ("Tether_04",  "Multi-segment tether (2D arrays)"),
                    ("Tether_05",  "Segmented tether, correct force split"),
                    ("Tether_06",  "Segmented tether, reeling out"),
                    ("Tether_06c", "Reel-out tether with continuous callback"),
                    ("Tether_07",  "Segmented tether with aerodynamic drag"),
                    ("Tether_08",  "Tether with arbitrary/free endpoints")]

python_name_width = maximum(length(name) for (name, _) in python_examples)
python_options = [rpad(name, python_name_width) * "  " * descr for (name, descr) in python_examples]
push!(python_options, "quit()")

function menu()
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

function menu2()
    active = true
    while active
        menu = RadioMenu(python_options, pagesize=8)
        choice = request("\nChoose Python example to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(python_options)
            run_python(python_examples[choice][1])
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end