println("Running init...")

import PyPlot: plot as pl
import PyPlot: xlabel, ylabel, xlim, ylim, grid, twinx, legend, figure, scatter, annotate, gcf

function plot(x...; y...)
    res = pl(x...; y...)
    PyPlot.show(block=false)
    res
end