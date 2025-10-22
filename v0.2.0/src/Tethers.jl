module Tethers

export docu

LAUNCH_BROWSER = true

function docu(build=true)
    @eval using LiveServer
    # if build
    #     include("docs/make.jl")
    # end
    if Sys.islinux() && ! LAUNCH_BROWSER
        Base.run(`xdg-open "docs/build/index.html"`; wait=false)
    else
        Base.invokelatest(LiveServer.servedocs; skip_dir="docs", launch_browser=true)
    end
end

end