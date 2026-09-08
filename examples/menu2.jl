using Pkg
# see the comment in `menu.jl`
if isfile(joinpath(@__DIR__, "Project.toml")) && dirname(Pkg.project().path) != @__DIR__
    Pkg.activate(@__DIR__)
end
using REPL.TerminalMenus
using Tethers: run_python

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

function run_menu2()
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
run_menu2()
