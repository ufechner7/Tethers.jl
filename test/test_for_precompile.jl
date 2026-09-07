__PC=true
let
    include("../examples/Tether_01.jl")
    include("../examples/Tether_02.jl")
    include("../examples/Tether_03.jl")
    include("../examples/Tether_08.jl")
    include("../examples/Tether_10.jl")

    GC.gc(true)
    let mem = Sys.free_memory() / 1024^2
        @info "Free memory: $(round(mem; digits=1)) MB"
        if haskey(ENV, "JULIA_IMAGE_THREADS")
            @info "JULIA_IMAGE_THREADS: $(ENV["JULIA_IMAGE_THREADS"])"
        else
            @info "JULIA_IMAGE_THREADS not defined!"
        end
    end
    nothing   
end

@info "Precompile script has completed execution."