# is included from Tether_04.jl.jl
function plot2d(x, z, reltime=0.0)
    x_max = maximum(x)
    z_max = maximum(z)
    plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    annotate!(x_max-10.0, z_max-3.0, "t=$(round(reltime,digits=1)) s")
    plot!(x, z, seriestype = :scatter) 
end