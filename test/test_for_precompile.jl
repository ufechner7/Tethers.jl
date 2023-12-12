__PC=true
let 
    include("../src/Tether_01.jl")
    include("../src/Tether_02.jl")
    include("../src/Tether_03.jl")
    include("../src/Tether_03c.jl")
    nothing   
end

@info "Precompile script has completed execution."