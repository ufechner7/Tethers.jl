println("Running init...")

import PyPlot: plot as pl

function plot(x...; y...)
    pl(x...; y...)
    PyPlot.show(block=false)
end