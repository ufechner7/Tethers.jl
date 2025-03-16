using LinearAlgebra, StaticArrays, ADTypes, NonlinearSolve, MAT, ForwardDiff, PreallocationTools
# experimental version for using automated differentiation

const MVec3 = MVector{3, Float64}
const Dual = ForwardDiff.Dual

"""
    Settings

Contains the environmental and tether properties

# Fields
  - rho::Float64: density of air [kg/m³] 
  - g_earth::MVector{Float64}: gravitational acceleration [m/s]
  - cd_tether::Float64: drag ceofficient of the tether
  - d_tether::Float64: diameter of the tether [mm]
  - rho_tether::Float64: density of the tether (Dyneema) [kg/m³]
  - c_spring::Float64: axial stiffness of the tether EA [N] 
""" 

struct Settings 
    rho::Float64
    g_earth::MVector{3, Float64}
    cd_tether::Float64
    d_tether::Float64                              
    rho_tether::Float64                             
    c_spring::Float64   
end

"""
    simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

Function to determine the tether shape and forces, based on a quasi-static model.

# Arguments
- state_vec::MVector{3, Float64}: state vector (theta [rad], phi [rad], Tn [N]);  
  tether orientation and tension at ground station
- kite_pos::MVector{3, Float64}: kite position vector in wind reference frame
- kite_vel::MVector{3, Float64}: kite velocity vector in wind reference frame
- wind_vel:: (3, Ns) MMatrix{Float64} wind velocity vector in wind reference frame for each Ns node of the tether
- tether_length: tether length
- settings:: Settings struct containing environmental and tether parameters: see [Settings](@ref)

# Returns
- state_vec::MVector{3, Float64}: state vector (theta [rad], phi [rad], Tn [N]);  
  tether orientation and tension at ground station 
- tether_pos::Matrix{Float64}: x,y,z - coordinates of the tether nodes
- force_gnd::Float64: Line tension at the ground station
- force_kite::Vector{Float64}: force from the kite to the end of tether
- p0::Vector{Float64}:  x,y,z - coordinates of the kite-tether attachment
"""
function simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings; prn=false)
    Ns = size(wind_vel, 2)
    buffers = [DiffCache(zeros(3, Ns)), DiffCache(zeros(3, Ns)), DiffCache(zeros(3, Ns)), DiffCache(zeros(3, Ns)), DiffCache(zeros(3, Ns)), ]
    
    # Pack parameters in param named tuple - false sets res! for in-place solution
    param = (kite_pos=kite_pos, kite_vel=kite_vel, wind_vel=wind_vel, 
         tether_length=tether_length, settings=settings, buffers=buffers, 
         returnFlag=false)
    # Define the nonlinear problem
    dual_vec = Dual.(state_vec, 0)
    prob = NonlinearProblem(res!, dual_vec, param)
    # Solve the problem with TrustRegion method
    sol = solve(prob, TrustRegion(); show_trace = Val(false)) 

    iterations = sol.stats.nsteps  # Field name may vary; verify with `propertynames(sol)`
    state_vec = sol.u
    if prn
        println("Iterations: ", iterations)
    end

    # # Set the returnFlag to true so that res! returns outputs
    # param = (; param..., returnFlag=true)
    # res = zeros(Dual, 3)
    # res, force_kite, tether_pos, p0 = res!(res, state_vec, param)

    # force_gnd = state_vec[3]
    # state_vec, tether_pos, force_gnd, force_kite, p0
end


"""
    res!(res, state_vec, param)

Calculates difference between tether end and kite given tether ground segment orientation 
and magnitude.

# Arguments
- res: 3 element vector, difference between tether end and kite segment
- state_vec::state vector (theta [rad], phi [rad], Tn [N]) - tether orientation and tension at ground station
- par:: 7-elements tuple:
    - kite_pos::MVector{3, Float64} kite position vector in wind reference frame
    - kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
    - wind_vel::MMatrix{Float64} wind velocity vector in wind reference frame for each Ns node of the tether
    - tether_length: tether length
    - settings:: Settings struct containing enviromental and tether parameters: see [Settings](@ref)
    - buffers:: (5, ) Vector{Matrix{Float64}}  Vector of (3, Ns) Matrix{Float64} empty matrices for preallocation
    - returnFlag:: Boolean to determine use for in-place optimization or for calculating returns

# Returns (if returnFlag==true)
- res::residual of length 3, difference between tether end and kite segment
- T0::Vector{Float64} force from the kite to the end of tether
- pj:: (3, Ns) Matrix{Float64} x,y,z - coordinates of the Ns tether nodes
- p0::Vector{Float64}  x,y,z - coordinates of the kite-tether attachment

# Example usage
state_vec = rand(3,)
kite_pos = [100, 100, 300] 
kite_vel = [0, 0, 0]
wind_vel = rand(3,15)
tether_length = 500
settings = Settings(1.225, [0, 0, -9.806], 0.9, 4, 0.85, 500000)
res!(res, state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)
"""
function res!(res, state_vec, param)
    kite_pos, kite_vel, wind_vel, tether_length, settings, buffers, returnFlag = param
    g = abs(settings.g_earth[3])
    Ns = size(wind_vel, 2)
    Ls = tether_length / (Ns + 1)
    mj = settings.rho_tether * Ls
    drag_coeff = -0.5 * settings.rho * Ls * settings.d_tether * settings.cd_tether
    A = π/4 * (settings.d_tether/1000)^2
    E = settings.c_spring / A

    # Preallocate arrays
    FT = get_tmp(buffers[1], state_vec)
    Fd = get_tmp(buffers[2], state_vec)
    pj = get_tmp(buffers[3], state_vec)
    vj = get_tmp(buffers[4], state_vec)
    aj = get_tmp(buffers[5], state_vec)

    # Unpack state variables
    θ, φ, Tn = state_vec[1], state_vec[2], state_vec[3]

    # Precompute common values
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    cosφ = cos(φ)
    norm_p = norm(kite_pos)
    p_unit = kite_pos ./ norm_p
    v_parallel = dot(kite_vel, p_unit)
    
    # First element calculations
    FT[1, Ns] = Tn * sinθ * cosφ
    FT[2, Ns] = Tn * sinφ
    FT[3, Ns] = Tn * cosθ * cosφ

    pj[1, Ns] = Ls * sinθ * cosφ
    pj[2, Ns] = Ls * sinφ
    pj[3, Ns] = Ls * cosθ * cosφ

    # Velocity and acceleration calculations
    ω = cross(kite_pos / norm_p^2, kite_vel) # 3 alloc
    a = cross(ω, (@view(pj[:, Ns])))         # 3 alloc
    b = cross(ω, cross(ω, (@view(pj[:, Ns]))))
    vj[:, Ns] .= v_parallel * p_unit + a
    aj[:, Ns] .= b

    # Drag calculation for first element
    v_a_p1 = vj[1, Ns] - wind_vel[1, Ns]
    v_a_p2 = vj[2, Ns] - wind_vel[2, Ns]
    v_a_p3 = vj[3, Ns] - wind_vel[3, Ns]

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

        for k in 1:3
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
        a = cross(ω, (@view(pj[:, ii-1])))           # 28 allocations
        b = cross(ω, cross(ω, (@view(pj[:, ii-1])))) # 28 allocations
        vj[:, ii-1] .= v_parallel * p_unit + a
        aj[:, ii-1] .= b

        # Drag calculations
        v_a_p1 = vj[1, ii] - wind_vel[1, ii]
        v_a_p2 = vj[2, ii] - wind_vel[2, ii]
        v_a_p3 = vj[3, ii] - wind_vel[3, ii]

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

    res .= kite_pos - p0
    if returnFlag
        return res, MVector(T0_1, T0_2, T0_3), pj, p0
    else
        nothing
    end
end

"""
    get_initial_conditions(filename)

Loads the initialization data for the basic examples and tests

# Arguments
- filename: the filename of the mat file to read

# Returns
- state_vec::MVector{3, Float64} state vector (theta [rad], phi [rad], Tn [N]) - tether orientation and tension at ground station
- kite_pos::MVector{3, Float64} kite position vector in wind reference frame
- kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
- wind_vel::MMatrix{3, Ns, Float64} wind velocity vector in wind reference frame for each Ns node of the tether
- tether_length: Float64 tether length
- settings::Settings struct containing enviromental and tether parameters: see [Settings](@ref)
"""
function get_initial_conditions(filename)
    vars = matread(filename) 
    state_vec = MVector{3}(vec(get(vars,"stateVec", 0)))
    kite_pos = MVector{3}(vec(get(vars,"kitePos", 0)))
    kite_vel = MVector{3}(vec(get(vars,"kiteVel", 0)))
    wind_vel = get(vars,"windVel", 0)
    tether_length = get(vars,"tetherLength", 0)

    ENVMT = get(vars,"ENVMT", 0) 
    rho_air = get(ENVMT, "rhos", 0) 
    g_earth = [0; 0; -abs(get(ENVMT, "g", 0))]      # in this way g_earth is a vector [0; 0; -9.81]

    T = get(vars,"T", 0)
    cd_tether = get(T, "CD_tether", 0) 
    d_tether = get(T, "d_tether", 0)*1000           # tether diameter                  [mm]
    rho_tether = get(T, "rho_t", 0) 
    E = get(T, "E", 0) 
    A = get(T, "A", 0)
    c_spring = E*A 

    settings = Settings(rho_air, g_earth, cd_tether, d_tether, rho_tether, c_spring)

    return state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings
end
