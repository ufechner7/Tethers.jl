using Pkg
# see the comment in `menu.jl`
if isfile(joinpath(@__DIR__, "Project.toml")) && dirname(Pkg.project().path) != @__DIR__
    Pkg.activate(@__DIR__)
end
using REPL.TerminalMenus

quasisteady_examples = [("run_catenary",        "include(\"quasisteady/run_catenary.jl\")",             "Quasi-steady tether shape for a kite at a fixed position"),
                         ("run_catenary_matlab", "include(\"quasisteady/run_catenary_matlab.jl\")",      "Quasi-steady shape for the MATLAB reference case vs. analytic catenary"),
                         ("flying_circular",     "include(\"quasisteady/flying_circular.jl\"); main()",  "Tether shape and force while the kite flies a circular trajectory"),
                         ("force_plots",         "include(\"quasisteady/force_plots.jl\"); main()",      "Tether shape and force as a function of kite distance"),
                         ("benchmark_qsm",       "include(\"quasisteady/benchmark_qsm.jl\")",            "Benchmark: quasi-steady model, elevation/azimuth angles")]

quasisteady_name_width = maximum(length(name) for (name, _, _) in quasisteady_examples)
quasisteady_options = [rpad(name, quasisteady_name_width) * "  " * descr for (name, _, descr) in quasisteady_examples]
push!(quasisteady_options, "quit()")

function run_menu3()
    active = true
    while active
        menu = RadioMenu(quasisteady_options, pagesize=8)
        choice = request("\nChoose quasi-steady example to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(quasisteady_options)
            eval(Meta.parse(quasisteady_examples[choice][2]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end
run_menu3()
