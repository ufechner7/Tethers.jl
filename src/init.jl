println("Running init...")

import PyPlot
import PyPlot: plot as pl
import PyPlot: xlabel, ylabel, xlim, ylim, twinx, legend, figure, scatter, annotate, gcf,
               tight_layout
using REPL.TerminalMenus

function plot(x...; y...)
    res = pl(x...; y...)
    PyPlot.show(block=false)
    res
end

function grid(x...; y...)
    res = PyPlot.grid(x; y...)
    PyPlot.pause(0.01)
    PyPlot.show(block=false)
    res
end

options = ["include(\"src/Tether_01.jl\")",
           "include(\"src/Tether_02.jl\")",
           "include(\"src/Tether_03.jl\")",
           "include(\"src/Tether_04.jl\")",
           "include(\"src/Tether_05.jl\")",
           "include(\"src/Tether_06.jl\")"]

function menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end