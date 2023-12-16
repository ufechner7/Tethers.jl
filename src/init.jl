println("Running init...")

import PyPlot: plot as pl
import PyPlot: xlabel, ylabel, grid, twinx, legend, figure

function plot(x...; y...)
    res = pl(x...; y...)
    PyPlot.show(block=false)
    res
end