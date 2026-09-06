using Test, LinearAlgebra

__BENCH__ = true;

@testset verbose=true "Tethers.jl" begin
    include("test_tether_01.jl")
    include("test_tether_02.jl")
    include("test_tether_03.jl")
    include("test_tether_03b.jl")
    include("test_tether_04.jl")
    include("test_tether_05.jl")
    include("test_tether_06.jl")
    include("test_tether_07.jl")
    include("test_tether_08.jl")
    include("test_tether_10.jl")
    include("test_tether_06c.jl")
end
nothing

