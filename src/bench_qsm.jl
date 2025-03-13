using LinearAlgebra, StaticArrays, BenchmarkTools

struct Settings 
    rho::Float64
    g_earth::MVector{3, Float64}
    cd_tether::Float64
    d_tether::Float64                               # tether diameter                  [mm]
    rho_tether::Float64                             # density of Dyneema            [kg/m³]
    c_spring::Float64   
end


function objFun!(res, stateVec, kitePos, kiteVel, windVel, tetherLength, settings, buffers)
    g = abs(settings.g_earth[3])
    Ns = size(windVel, 2)
    Ls = tetherLength / (Ns + 1)
    mj = settings.rho_tether * Ls
    drag_coeff = -0.5 * settings.rho * Ls * settings.d_tether * settings.cd_tether
    A = π/4 * (settings.d_tether/1000)^2
    E = settings.c_spring / A

    # Preallocate arrays
    FT = buffers[1]
    Fd = buffers[2]
    pj = buffers[3]
    vj = buffers[4]
    aj = buffers[5]

    # Unpack state variables
    θ, φ, Tn = stateVec[1], stateVec[2], stateVec[3]

    # Precompute common values
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    cosφ = cos(φ)
    norm_p = norm(kitePos)
    p_unit = kitePos ./ norm_p
    v_parallel = dot(kiteVel, p_unit)
    
    # First element calculations
    FT[1, Ns] = Tn * sinθ * cosφ
    FT[2, Ns] = Tn * sinφ
    FT[3, Ns] = Tn * cosθ * cosφ

    pj[1, Ns] = Ls * sinθ * cosφ
    pj[2, Ns] = Ls * sinφ
    pj[3, Ns] = Ls * cosθ * cosφ

    # Velocity and acceleration calculations
    ω = cross(kitePos / norm_p^2, kiteVel) # 3 alloc
    a = cross(ω, @view(pj[:, Ns]))         # 3 alloc
    b = cross(ω, cross(ω, @view(pj[:, Ns])))
    vj[:, Ns] .= v_parallel * p_unit + a
    aj[:, Ns] .= b

    # Drag calculation for first element
    v_a_p1 = vj[1, Ns] - windVel[1, Ns]
    v_a_p2 = vj[2, Ns] - windVel[2, Ns]
    v_a_p3 = vj[3, Ns] - windVel[3, Ns]

    if all(x -> abs(x) < 1e-3, (v_a_p1, v_a_p2, v_a_p3))
        Fd[:, Ns] .= 0.0
    else
        dir1, dir2, dir3 = pj[1, Ns]/Ls, pj[2, Ns]/Ls, pj[3, Ns]/Ls
        v_dot_dir = v_a_p1*dir1 + v_a_p2*dir2 + v_a_p3*dir3
        v_a_p_t1 = v_dot_dir * dir1
        v_a_p_t2 = v_dot_dir * dir2
        v_a_p_t3 = v_dot_dir * dir3

        v_a_p_n1 = v_a_p1 - v_a_p_t1
        v_a_p_n2 = v_a_p2 - v_a_p_t2
        v_a_p_n3 = v_a_p3 - v_a_p_t3

        norm_v_a_p_n = sqrt(v_a_p_n1^2 + v_a_p_n2^2 + v_a_p_n3^2)
        coeff = drag_coeff * norm_v_a_p_n

        Fd[1, Ns] = coeff * v_a_p_n1
        Fd[2, Ns] = coeff * v_a_p_n2
        Fd[3, Ns] = coeff * v_a_p_n3
    end

    # Process other segments
    @inbounds for ii in Ns:-1:2
        # Tension force calculations
        if ii == Ns
            mj_total = 1.5mj
            g_term = mj_total * g
        else
            mj_total = mj
            g_term = mj * g
        end

        @views for k in 1:3
            FT[k, ii-1] = mj_total * aj[k, ii] + FT[k, ii] - Fd[k, ii]
        end
        FT[3, ii-1] += g_term

        # Position calculations
        ft_norm = sqrt(FT[1, ii-1]^2 + FT[2, ii-1]^2 + FT[3, ii-1]^2)
        l_i_1 = (ft_norm/(E*A) + 1) * Ls
        ft_dir = FT[1, ii-1]/ft_norm, FT[2, ii-1]/ft_norm, FT[3, ii-1]/ft_norm

        pj[1, ii-1] = pj[1, ii] + l_i_1 * ft_dir[1]
        pj[2, ii-1] = pj[2, ii] + l_i_1 * ft_dir[2]
        pj[3, ii-1] = pj[3, ii] + l_i_1 * ft_dir[3]

        # Velocity and acceleration
        a = cross(ω, @view(pj[:, ii-1]))           # 28 allocations
        b = cross(ω, cross(ω, @view(pj[:, ii-1]))) # 28 allocations
        vj[:, ii-1] .= v_parallel * p_unit + a
        aj[:, ii-1] .= b

        # Drag calculations
        v_a_p1 = vj[1, ii] - windVel[1, ii]
        v_a_p2 = vj[2, ii] - windVel[2, ii]
        v_a_p3 = vj[3, ii] - windVel[3, ii]

        if all(x -> abs(x) < 1e-3, (v_a_p1, v_a_p2, v_a_p3))
            Fd[:, ii-1] .= 0.0
        else
            dx = pj[1, ii-1] - pj[1, ii]
            dy = pj[2, ii-1] - pj[2, ii]
            dz = pj[3, ii-1] - pj[3, ii]
            segment_norm = sqrt(dx^2 + dy^2 + dz^2)
            dir1 = dx/segment_norm
            dir2 = dy/segment_norm
            dir3 = dz/segment_norm

            v_dot_dir = v_a_p1*dir1 + v_a_p2*dir2 + v_a_p3*dir3
            v_a_p_t1 = v_dot_dir * dir1
            v_a_p_t2 = v_dot_dir * dir2
            v_a_p_t3 = v_dot_dir * dir3

            v_a_p_n1 = v_a_p1 - v_a_p_t1
            v_a_p_n2 = v_a_p2 - v_a_p_t2
            v_a_p_n3 = v_a_p3 - v_a_p_t3

            norm_v_a_p_n = sqrt(v_a_p_n1^2 + v_a_p_n2^2 + v_a_p_n3^2)
            coeff = drag_coeff * norm_v_a_p_n

            Fd[1, ii-1] = coeff * v_a_p_n1
            Fd[2, ii-1] = coeff * v_a_p_n2
            Fd[3, ii-1] = coeff * v_a_p_n3
        end
    end

    # Final ground connection calculations
    T0_1 = 1.5mj*aj[1,1] + FT[1,1] - Fd[1,1]
    T0_2 = 1.5mj*aj[2,1] + FT[2,1] - Fd[2,1]
    T0_3 = 1.5mj*aj[3,1] + FT[3,1] - Fd[3,1] + 1.5mj*g
    T0_norm = sqrt(T0_1^2 + T0_2^2 + T0_3^2)
    
    l_i_1 = (T0_norm/(E*A) + 1) * Ls
    T0_dir1 = T0_1/T0_norm
    T0_dir2 = T0_2/T0_norm
    T0_dir3 = T0_3/T0_norm

    p0 = [pj[1,1] + l_i_1*T0_dir1, 
          pj[2,1] + l_i_1*T0_dir2,
          pj[3,1] + l_i_1*T0_dir3]

    res .= kitePos - p0
    nothing
end

stateVec::MVector{3} = MVector{3}(rand(3,))
kitePos::SVector{3} = SVector{3}([100, 100, 300])
kiteVel::SVector{3} = SVector{3}([0, 0, 0])
windVel::SMatrix{3, 15} = SMatrix{3, 15}(rand(3,15))
tetherLength = 500
settings::Settings = Settings(1.225, [0, 0, -9.806], 0.9, 4, 0.85, 500000)
Ns::Int64 = size(windVel, 2)
buffers::Vector{Matrix{Float64}} = [zeros(3, Ns), zeros(3, Ns), zeros(3, Ns), zeros(3, Ns), zeros(3, Ns)]
res::Vector{Float64} = zeros(3)

objFun!(res, stateVec, kitePos, kiteVel, windVel, tetherLength, settings, buffers)
@benchmark objFun!(res, stateVec, kitePos, kiteVel, windVel, tetherLength, settings, buffers)



