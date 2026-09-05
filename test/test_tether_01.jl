function read_tether_01_csv(file)
    lines = readlines(file)
    n = length(lines) - 1
    t, pos_z, vel_z = zeros(n), zeros(n), zeros(n)
    for (i, line) in enumerate(lines[2:end])
        parts = split(line, ',')
        t[i]      = parse(Float64, parts[1])
        pos_z[i]  = parse(Float64, parts[2])
        vel_z[i]  = parse(Float64, parts[3])
    end
    t, pos_z, vel_z
end

@testset "Tether_01" begin
    pkg_dir = dirname(@__DIR__)
    cd(pkg_dir) do
        # Julia implementation, using ModelingToolkit and the implicit solver Rodas5
        include(joinpath(pkg_dir, "src", "Tether_01.jl"))
        sleep(1)
        MakieControlPlots.close("all")
        # Python implementation, using the implicit solver RADAU
        withenv("TETHERS_BRIEF_PLOT" => "1") do
            include(joinpath(pkg_dir, "src", "RunTether_01.jl"))
        end
    end

    t_jl, pos_z_jl, vel_z_jl = read_tether_01_csv(joinpath(pkg_dir, "output", "Tether_01_julia.csv"))
    t_py, pos_z_py, vel_z_py = read_tether_01_csv(joinpath(pkg_dir, "output", "Tether_01_python.csv"))

    @test length(t_jl) == length(t_py)
    @test t_jl ≈ t_py atol=1e-6
    @test pos_z_jl ≈ pos_z_py rtol=1e-3
    @test vel_z_jl ≈ vel_z_py rtol=1e-3
end
nothing
