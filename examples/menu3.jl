using Pkg
# see the comment in `menu.jl`
if isfile(joinpath(@__DIR__, "Project.toml")) && dirname(Pkg.project().path) != @__DIR__
    Pkg.activate(@__DIR__)
end
using REPL.TerminalMenus

quasistatic_examples = [("run_catenary",        "include(\"quasistatic/run_catenary.jl\")",             "Quasi-static tether shape for a kite at a fixed position"),
                         ("run_catenary_matlab", "include(\"quasistatic/run_catenary_matlab.jl\")",      "Quasi-static shape for the MATLAB reference case vs. analytic catenary"),
                         ("flying_circular",     "include(\"quasistatic/flying_circular.jl\"); main()",  "Tether shape and force while the kite flies a circular trajectory"),
                         ("force_plots",         "include(\"quasistatic/force_plots.jl\"); main()",      "Tether shape and force as a function of kite distance"),
                         ("benchmark_qsm",       "include(\"quasistatic/benchmark_qsm.jl\")",            "Benchmark: quasi-static model, elevation/azimuth angles"),
                         ("benchmark_qsm_dual",  "include(\"quasistatic/benchmark_qsm_dual.jl\")",       "Benchmark: quasi-static model, dual-number formulation")]

quasistatic_name_width = maximum(length(name) for (name, _, _) in quasistatic_examples)
quasistatic_options = [rpad(name, quasistatic_name_width) * "  " * descr for (name, _, descr) in quasistatic_examples]
push!(quasistatic_options, "quit()")

function run_menu3()
    active = true
    while active
        menu = RadioMenu(quasistatic_options, pagesize=8)
        choice = request("\nChoose quasi-static example to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(quasistatic_options)
            eval(Meta.parse(quasistatic_examples[choice][2]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end
run_menu3()
