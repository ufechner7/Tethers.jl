__PC=true
let 
    include("../src/init.jl")
    include("../src/Tether_01.jl")
    include("../src/Tether_02.jl")
    include("../src/Tether_03.jl")
    include("../src/Tether_08.jl")
    
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