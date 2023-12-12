# is included from Tether_05.jl.jl
function plot2d(sol, reltime, segments)
    index = Int64(round(reltime*50+1))
    x = Float64[]
    z = Float64[]
    for particle in 1:segments+1
        push!(x, (sol(sol.t, idxs=pos[1, particle]))[index])
        push!(z, (sol(sol.t, idxs=pos[3, particle]))[index])
    end
    x_max = maximum(x)
    z_max = maximum(z)
    plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    annotate!(x_max+0.5, z_max-3.0, "t=$(round(reltime,digits=1)) s")
    plot!(x, z, seriestype = :scatter) 
    ylims!((-segments*10-5, 0.5))
end