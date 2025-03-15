using LinearAlgebra, StaticArrays, BenchmarkTools

#= BenchmarkTools.Trial: 10000 samples with 186 evaluations per sample.
 Range (min … max):  552.016 ns …  33.646 μs  ┊ GC (min … max): 0.00% … 97.83%
 Time  (median):     556.575 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   570.646 ns ± 463.443 ns  ┊ GC (mean ± σ):  1.27% ±  1.83%

  ▁ ▁██▆▅▃▂▁▁▁▁▁ ▂▄▃▂▂▂▂▁▁▁▁▄▄▃▂▁▁▁▁▁ ▁▁     ▁                  ▂
  █████████████████████████████████████████████▇█▇▆▆▇▆▆▆▆▆▆▆▆▅▆ █
  552 ns        Histogram: log(frequency) by time        606 ns <

 Memory estimate: 80 bytes, allocs estimate: 2. =#

struct Settings 
    rho::Float64
    g_earth::SVector{3, Float64}
    cd_tether::Float64
    d_tether::Float64 
    rho_tether::Float64
    c_spring::Float64   
end

function objFun!(res, stateVec, kitePos, kiteVel, windVel, tetherLength, settings)

    # settings follows same syntax as Settings3 in Tether_09.jl
    g  = abs(settings.g_earth[3])  # in this function, g is considered a scalar 

    Ns = length(windVel)           # number of masses - windVel is a 3xNs matrix
    Ls = tetherLength / (Ns+1)        # segment length
    mj = settings.rho_tether*Ls     # segment mass

    p = kitePos                     # input kite position
    p_unit = normalize(p)
    v = kiteVel                     # input kite velocity

    # Extract states
    theta   = stateVec[1]
    phi     = stateVec[2]
    Tn      = stateVec[3]

    # Compute omega_t
    omega_t = cross(p / sum(abs2, p), v);

    # Tether cross section
    A = pi/4*(settings.d_tether/1000)^2 # [m2] 
    # Compute Young's modulus    
    E = settings.c_spring/A

    # # First element from ground station (Ns)
    # TODO: a lot of repeated computation below
    FT = Tn .* SA[sin(theta)* cos(phi), sin(phi), cos(theta)*cos(phi)]
    pj = Ls .* SA[sin(theta)* cos(phi), sin(phi), cos(theta)*cos(phi)]
    pj_unit = normalize(pj)
    vj = dot(v, p_unit) * p_unit + cross(omega_t, pj)
    aj = cross(omega_t, cross(omega_t, pj))

    # Drag calculation first element
    v_a_p = vj - windVel[end]
    if all(abs.(v_a_p) .< 1e-3)
        Fd = SA[0.0, 0.0, 0.0]
    else
        v_a_p_t = dot(pj_unit, v_a_p) * pj_unit
        v_a_p_n = v_a_p - v_a_p_t
        Fd = (-0.5 * settings.rho * Ls * settings.d_tether * settings.cd_tether * 
            norm(v_a_p_n) * v_a_p_n )# particle drag
    end

    # All other segments and masses except for segment connected to the kite
    for ii = Ns:-1:2
        if ii == Ns
            FT = (mj+0.5*mj) * aj + FT - Fd + SA[0, 0, (mj+0.5*mj)*g]
        else
            FT = mj*aj + FT - Fd + SA[0, 0, mj*g]
        end

        FT_norm = norm(FT)
        l_i_1 = (FT_norm / (E*A) + 1) * Ls
        
        # Drag calculation
        v_a_p = vj - windVel[ii] # moving this up here

        Δpj = (l_i_1/FT_norm) * FT
        pj += Δpj
        vj = dot(v, p_unit) * p_unit + cross(omega_t, pj)
        aj = cross(omega_t, cross(omega_t, pj))
        
        if all(abs.(v_a_p) .< 1e-3)
            Fd = SA[0.0, 0.0, 0.0]
        else
            v_a_p_t = (dot(Δpj, v_a_p) / sum(abs2, Δpj)) * Δpj
            v_a_p_n = v_a_p - v_a_p_t;
            Fd = -0.5 * settings.rho * Ls * settings.d_tether * settings.cd_tether * norm(v_a_p_n) * v_a_p_n # particle drag
        end
    end

    T0 = (mj+0.5*mj) .* aj + FT - Fd + SA[0, 0, (mj+0.5*mj)*g]
    l_i_1 = (norm(T0)/(E*A) + 1) * Ls
    p0 = pj + l_i_1 .* (T0 / norm(T0))
    
    res .= p - p0
    nothing
end

stateVec = rand(SVector{3, Float64})
kitePos = SA_F64[100, 100, 300]
kiteVel = SA_F64[0, 0, 0]
windVel = rand(SVector{3, Float64}, 15)
tetherLength = 500.0
settings = Settings(1.225, SA_F64[0, 0, -9.806], 0.9, 4, 0.85, 500000)
res = zeros(3)

objFun!(res, stateVec, kitePos, kiteVel, windVel, tetherLength, settings)
@benchmark objFun!(res, stateVec, kitePos, kiteVel, windVel, tetherLength, settings)
