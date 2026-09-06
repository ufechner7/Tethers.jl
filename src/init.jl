println("Running init...")

using MakieControlPlots
using REPL.TerminalMenus

examples = [("Tether_01",  "include(\"src/Tether_01.jl\")",  "Falling mass thrown upwards"),
            ("Tether_02",  "include(\"src/Tether_02.jl\")",  "Mass on a linear spring-damper"),
            ("Tether_03",  "include(\"src/Tether_03.jl\")",  "Mass on a non-linear spring-damper"),
            ("Tether_03b", "include(\"src/Tether_03b.jl\")", "Non-linear spring with callback"),
            ("Tether_03c", "include(\"src/Tether_03c.jl\")", "Benchmark: callback vs no callback"),
            ("Tether_04",  "include(\"src/Tether_04.jl\")",  "Multi-segment tether (2D arrays)"),
            ("Tether_05",  "include(\"src/Tether_05.jl\")",  "Segmented tether, correct force split"),
            ("Tether_06",  "include(\"src/Tether_06.jl\")",  "Segmented tether, reeling out"),
            ("Tether_06b", "include(\"src/Tether_06b.jl\")", "Reel-out tether, refactored with Settings"),
            ("Tether_06c", "include(\"src/Tether_06c.jl\")", "Reel-out tether with continuous callback"),
            ("Tether_07",  "include(\"src/Tether_07.jl\")",  "Segmented tether with aerodynamic drag"),
            ("Tether_08",  "include(\"src/Tether_08.jl\")",  "Tether with arbitrary/free endpoints"),
            ("Tether_09",  "include(\"src/Tether_09.jl\")",  "Labeled tether shape diagram for docs")]

name_width = maximum(length(name) for (name, _, _) in examples)
options = [rpad(name, name_width) * "  " * descr for (name, _, descr) in examples]
push!(options, "quit()")

python_examples = [("Tether_01",  "include(\"src/RunTether_01.jl\")",  "Falling mass thrown upwards"),
                    ("Tether_02",  "include(\"src/RunTether_02.jl\")",  "Mass on a linear spring-damper"),
                    ("Tether_03",  "include(\"src/RunTether_03.jl\")",  "Mass on a non-linear spring-damper"),
                    ("Tether_03b", "include(\"src/RunTether_03b.jl\")", "Non-linear spring with callback"),
                    ("Tether_04",  "include(\"src/RunTether_04.jl\")",  "Multi-segment tether (2D arrays)"),
                    ("Tether_05",  "include(\"src/RunTether_05.jl\")",  "Segmented tether, correct force split"),
                    ("Tether_06",  "include(\"src/RunTether_06.jl\")",  "Segmented tether, reeling out"),
                    ("Tether_06c", "include(\"src/RunTether_06c.jl\")", "Reel-out tether with continuous callback"),
                    ("Tether_07",  "include(\"src/RunTether_07.jl\")",  "Segmented tether with aerodynamic drag"),
                    ("Tether_08",  "include(\"src/RunTether_08.jl\")",  "Tether with arbitrary/free endpoints")]

python_name_width = maximum(length(name) for (name, _, _) in python_examples)
python_options = [rpad(name, python_name_width) * "  " * descr for (name, _, descr) in python_examples]
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
            eval(Meta.parse(python_examples[choice][2]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end