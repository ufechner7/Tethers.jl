# is included from simulate.jl
function plot2d(pos, reltime=0.0; zoom=true, front=false, segments=6)
    x = Float64[] 
    z = Float64[]
    for i in eachindex(pos)
        if front
            push!(x, pos[i][2])
        else
            push!(x, pos[i][1])
        end
        push!(z, pos[i][3])
    end
    x_max = maximum(x)
    z_max = maximum(z)
    plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    annotate!(x_max-10.0, z_max-3.0, "t=$(round(reltime,digits=1)) s")
    plot!(x, z, seriestype = :scatter) 
end