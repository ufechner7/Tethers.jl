module Tethers

export docu

function docu(build=true)
    @eval using LiveServer
    # if build
    #     include("docs/make.jl")
    # end
    Base.invokelatest(LiveServer.servedocs; skip_dir="docs", launch_browser=true)
end

end