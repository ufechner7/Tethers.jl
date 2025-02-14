using ControlPlots
include("../src/video.jl")
X=[[1,2,3], [1.1,2.2,3.3]]
Y=[[2,3,5],[2.2,3.3,4.4]]
lines,sc = nothing, nothing
xlim=(0,10)
ylim=(0,10)
for (x,y) in (X,Y)
    global lines, sc
    println(x)
    println(y)
    lines,sc = plot_kite(x,y,xlim,ylim,lines,sc)
    #plt.tight_layout()
    plt.pause(0.01)
    plt.show(block=false)
    sleep(2)
 end
 # lines=plot_kite(x,y)
 # plt.show(block=false)