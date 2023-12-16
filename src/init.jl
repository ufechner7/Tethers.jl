println("Running init...")

import PyPlot: plot as pl
import PyPlot: xlabel, ylabel, xlim, ylim, twinx, legend, figure, scatter, annotate, gcf

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