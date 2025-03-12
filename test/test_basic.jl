using Test, LinearAlgebra, MAT
include("../src/Tether_quasistatic.jl")

function get_init_cond(filename)
    vars = matread(filename); 
    stateVec = vec(get(vars,"stateVec", 0));
    kitePos = vec(get(vars,"kitePos", 0));
    kiteVel = vec(get(vars,"kiteVel", 0));
    windVel = get(vars,"windVel", 0);
    tetherLength = get(vars,"tetherLength", 0);


    ENVMT = get(vars,"ENVMT", 0) 
    rho_air = get(ENVMT, "rhos", 0) 
    g_earth = [0; 0; -abs(get(ENVMT, "g", 0))]      # in this way g_earth is a vector [0; 0; -9.81]

        T = get(vars,"T", 0);
    cd_tether = get(T, "CD_tether", 0) 
    d_tether = get(T, "d_tether", 0)*1000           # tether diameter                  [mm]
    rho_tether = get(T, "rho_t", 0) 
    E = get(T, "E", 0) 
    A = get(T, "A", 0)
    c_spring = E*A 

    settings = Settings(rho_air, g_earth, cd_tether, d_tether, rho_tether, c_spring)

    return stateVec, kitePos, kiteVel, windVel, tetherLength, settings
end

function get_output_objfun(filename)
    vars        = matread(filename); 
    Fobj        = vec(get(vars,"Fobj", 0))
    p0          = vec(get(vars,"p0", 0))
    pj          = get(vars,"pj", 0)
    T0          = vec(get(vars,"T0", 0)) 
    return Fobj, p0, pj, T0
end

@testset "ObjFun_test" begin  
    stateVec, kitePos, kiteVel, windVel, tetherLength, settings = get_init_cond("data/input_basic_test.mat")
    Fobj_ref, p0_ref, pj_ref, T0_ref = get_output_objfun("data/basic_test_results.mat")
    # Call objective function
    Fobj, T0, pj, p0 = objFun(stateVec, kitePos, kiteVel, windVel, tetherLength, settings)

    # Test basic types
    @test Fobj isa Vector
    @test T0 isa Vector
    @test p0 isa Vector
    @test pj isa Matrix

    # Test values
    @test Fobj == Fobj_ref
    @test T0 == T0_ref
    @test pj ≈ pj_ref
    @test p0 == p0_ref
    nothing
end