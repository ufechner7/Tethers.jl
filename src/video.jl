function plot_kite(X, Y, xlim_, ylim_, lines = nothing; fig="")

    if fig != ""
        plt.figure(fig, figsize, dpi)
    end
    if isnothing(lines)
        lines=[]
        for i in 1:length(X)-1
            x1=X[i:i+1]
            y1=Y[i:i+1]
            line, = plt.plot(x1,y1; linewidth="1")
            push!(lines, line)
        end

        #sc  = plt.scatter(x, y; s=25, color="red")
    else
        i=1
        for line in lines
            println(i)
            x1=X[i:i+1]
            y1=Y[i:i+1]
            line.set_xdata(x1)
            line.set_ydata(y1)
            i+=1
        end
    end
    plt.ylim(ylim_)
    plt.xlim(xlim_)
    return lines
end



# """
#     plot2d_(pos, reltime; zoom=true, front=false, segments=6, fig="", figsize=(6.4, 4.8), dz_zoom=1.5, dz=-5.0, 
#             dx=-16.0, xlim=nothing, ylim=nothing, xy=nothing, lines, sc, txt)

# Display a video-like 2D particle system by calling `plot2d` in a loop.

# # Arguments
# - `pos`: a vector of 3D positions.
# - `reltime`: The relative time. When called the first time, set to `0.0`.
# - `zoom`: Whether to enable zooming (default: `true`).
# - `front`: Whether using a front view (default: `false`).
# - `segments`: The number of tether segments (default: `6`).
# - `fig`: The name of the figure to display (default: `""`).
# - `figsize: The size of the figure in inch (default: `(6.4, 4.8)`).
# - `dz_zoom`: The z-axis offset in zoom view (default: `1.5`).
# - `dz`: The z-axis offset in normal view (default: `-5.0`).
# - `dx`: The x-axis offset (default: `-16.0`).
# - `xlim`: The x-axis limits (default: `nothing`).
# - `ylim`: The y-axis limits (default: `nothing`) (in reality the z-axis).
# - `xy`: The x-y coordinates of the text (default: `nothing`).
# - `lines`: The lines to plot.
# - `sc`: The dots (scatter plot).
# - `txt`: The text to display.

# """
# function plot2d_(pos, reltime; zoom=true, front=false, segments=6, fig="", figsize=(6.4, 4.8), dz_zoom= 1.5, dpi=100,
#                  dz=-5.0, dx=-16.0, xlim=nothing, ylim=nothing, xy=nothing, lines, sc, txt)
#     x = Float64[] 
#     z = Float64[]
#     for i in eachindex(pos)
#         if front
#             push!(x, pos[i][2])
#         else
#             push!(x, pos[i][1])
#         end
#         push!(z, pos[i][3])
#     end
#     x_max = maximum(x)
#     x_min = minimum(x)
#     z_max = maximum(z)
#     xlabel = "x [m]"
#     if front xlabel = "y [m]" end
#     if isnothing(lines)
#         if fig != ""
#             plt.figure(fig, figsize, dpi)
#         end
#         lines=[]
#         line, = plt.plot(x,z; linewidth="1")
#         push!(lines, line)
#         sc  = plt.scatter(x, z; s=25, color="red") 
#         if zoom
#             if isnothing(xy)
#                 xy=(x_max, z_max+dz_zoom)
#             end
#             txt = plt.annotate("t=$(round(reltime,digits=1)) s",  
#                 xy, fontsize = 14)
#             if isnothing(xlim)
#                 plt.xlim(x_max-15.0, x_max+5)
#             else
#                 plt.xlim(xlim)
#             end
#             if isnothing(ylim)
#                 plt.ylim(z_max-15.0, z_max+5)
#             else
#                 plt.ylim(ylim)
#             end
#         else
#             if isnothing(xy)
#                 xy=(x_max+dx, z_max+dz)
#             end
#             txt = plt.annotate("t=$(round(reltime,digits=1)) s",  xy, fontsize = 14)
#             if isnothing(xlim)
#                 plt.xlim(0, x_max+5)
#             else
#                 plt.xlim(xlim)
#             end
#             if isnothing(ylim)
#                 plt.ylim(0, z_max+5)
#             else
#                 plt.ylim(ylim)
#             end
#         end
#         if length(pos) > segments+1
#             s=segments
#             line, = plt.plot([x[s+1],x[s+4]],[z[s+1],z[s+4]], linewidth="1"); push!(lines, line) # S6
#             line, = plt.plot([x[s+2],x[s+5]],[z[s+2],z[s+5]], linewidth="1"); push!(lines, line) # S8
#             line, = plt.plot([x[s+3],x[s+5]],[z[s+3],z[s+5]], linewidth="1"); push!(lines, line) # S7
#             line, = plt.plot([x[s+2],x[s+4]],[z[s+2],z[s+4]], linewidth="1"); push!(lines, line) # S2
#             line, = plt.plot([x[s+1],x[s+5]],[z[s+1],z[s+5]], linewidth="1"); push!(lines, line) # S5
#         end
#         plt.xlabel(xlabel, fontsize=14)
#         plt.ylabel("z [m]", fontsize=14)
#         plt.grid(true)
#         plt.grid(which="major", color="#DDDDDD")
#     else
#         lines[1].set_xdata(x)
#         lines[1].set_ydata(z)
#         if length(pos) > segments+1
#             s=segments
#             lines[2].set_xdata([x[s+1],x[s+4]]) # S6
#             lines[2].set_ydata([z[s+1],z[s+4]]) # S6
#             lines[3].set_xdata([x[s+2],x[s+5]]) # S8
#             lines[3].set_ydata([z[s+2],z[s+5]]) # S8
#             lines[4].set_xdata([x[s+3],x[s+5]]) # S7
#             lines[4].set_ydata([z[s+3],z[s+5]]) # S7
#             lines[5].set_xdata([x[s+2],x[s+4]]) # S2
#             lines[5].set_ydata([z[s+2],z[s+4]]) # S2
#             lines[6].set_xdata([x[s+1],x[s+5]]) # S5
#             lines[6].set_ydata([z[s+1],z[s+5]]) # S5
#         end
#         sc.set_offsets(hcat(x,z))
#         # xy=(x_max, z_max+8.0)
#         txt.set_text("t=$(round(reltime, RoundDown; digits=1)) s")
#         if zoom
#             txt.set_x(x_max)
#             txt.set_y(z_max+dz_zoom)
#             if isnothing(xlim)
#                 plt.xlim(x_min-5.0, x_max+5)
#             else
#                 plt.xlim(xlim)
#             end
#             if isnothing(ylim)
#                 plt.ylim(z_max-15.0, z_max+5)
#             else
#                 plt.ylim(ylim)
#             end
#         else
#             if isnothing(xy)
#                 xy=(x_max+dx, z_max+dz)
#             end
#             txt.set_x(xy[1])
#             txt.set_y(xy[2])
#             if isnothing(xlim)
#                 plt.xlim(0, x_max+5)
#             else
#                 plt.xlim(xlim)
#             end
#             if isnothing(ylim)
#                 plt.ylim(0, z_max+5)
#             else
#                 plt.ylim(ylim)
#             end
#         end
#     end
#     if front
#         plt.gca().invert_xaxis()
#     end
#     plt.tight_layout()
#     plt.pause(0.01)
#     plt.show(block=false)
#     lines, sc, txt  
# end

# plot2d__ = let lines = nothing, sc = nothing, txt = nothing  # Note: Must all be on same line as let!
#     function(pos::AbstractVector, reltime=0.0; fig="", figsize=(6.4, 4.8), dpi=100, kwargs...)
#         if reltime == 0.0
#             lines, sc, txt = nothing, nothing, nothing
#         end
#         lines, sc, txt = plot2d_(pos, reltime; lines, sc, txt, fig, figsize, dpi, kwargs...)
#     end
# end

# """
#     plot2d(pos::AbstractVector, reltime; zoom=true, front=false, segments=6, fig="", dz_zoom=1.5, 
#            dz=-5.0, dx=-16.0, xlim=nothing, ylim=nothing, xy=nothing)

# Display a video-like 2D particle system by calling `plot2d` in a loop.

# # Arguments
# - `pos`: a vector of 3D positions.
# - `reltime`: The relative time. When called the first time, set to `0.0`.
# - `zoom`: Whether to enable zooming (default: `true`).
# - `front`: Whether using a front view (default: `false`, which means side view).
# - `segments`: The number of tether segments (default: `6`).
# - `fig`: The name of the figure to display (default: `""`).
# - `dz_zoom`: The z-axis offset in zoom view (default: `1.5`).
# - `dz`: The z-axis offset in normal view (default: `-5.0`).
# - `dx`: The x-axis offset (default: `-16.0`).
# - `xlim`: The x-axis limits (default: `nothing`).
# - `ylim`: The y-axis limits (default: `nothing`) (for side view the z-axis).
# - `xy`: The x-y coordinates of the text (default: `nothing`) (for side view the z-axis).

# """
# function plot2d(pos::AbstractVector, reltime=0.0; fig="", figsize=(6.4, 4.8), dpi=100, kwargs...)
#     plot2d__(pos, reltime; fig, figsize, dpi, kwargs...)
# end
