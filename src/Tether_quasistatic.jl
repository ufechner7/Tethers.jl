using LinearAlgebra, StaticArrays, ADTypes, NonlinearSolve, MAT, Parameters#, QuadGK

const MVec3 = MVector{3, Float64}
const SVec3 = SVector{3, Float64}
# const segments = 15

# Iterations: 36
# BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
#  Range (min … max):  84.169 μs …  2.388 ms  ┊ GC (min … max): 0.00% … 93.18%
#  Time  (median):     89.110 μs              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   92.513 μs ± 47.339 μs  ┊ GC (mean ± σ):  1.10% ±  2.09%

#     ▂█▇▆▂▁                                                     
#   ▁▃██████▆▆▆▅▅▅▄▅▄▄▄▃▃▄▃▃▃▃▃▃▄▃▃▃▃▃▃▃▃▃▃▃▃▂▂▂▃▂▂▂▂▂▂▂▂▂▂▂▁▁▁ ▃
#   84.2 μs         Histogram: frequency by time         108 μs <

#  Memory estimate: 41.20 KiB, allocs estimate: 971.

"""
    Settings

Contains the environmental and tether properties

# Fields
  - rho::Float64: density of air [kg/m³] 
  - g_earth::MVector{Float64}: gravitational acceleration [m/s]
  - cd_tether::Float64: drag coefficient of the tether
  - d_tether::Float64: diameter of the tether [mm]
  - rho_tether::Float64: density of the tether (Dyneema) [kg/m³]
  - c_spring::Float64: axial stiffness of the tether EA [N] 
""" 

@with_kw mutable struct Settings @deftype Float64
    rho = 1.225
    g_earth::MVector{3, Float64} = [0.0, 0.0, -9.81]
    cd_tether = 0.958                            
    d_tether = 4                                 
    rho_tether = 724                             
    c_spring = 614600                            
end

"""
    simulate_tether(state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings)

Function to determine the tether shape and forces, based on a quasi-static model.

# Arguments
- state_vec::MVector{3, Float64}: state vector (theta [rad], phi [rad], Tn [N]);  
  tether orientation and tension at ground station
- kite_pos::MVector{3, Float64}: kite position vector in wind reference frame
- kite_vel::MVector{3, Float64}: kite velocity vector in wind reference frame
- wind_vel:: (3, segments) MMatrix{Float64} wind velocity vector in wind reference frame for each segment of the tether
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
    segments = size(wind_vel)[2]
    buffers= [MMatrix{3, segments}(zeros(3, segments)), MMatrix{3, segments}(zeros(3, segments)), MMatrix{3, segments}(zeros(3, segments)), 
              MMatrix{3, segments}(zeros(3, segments)), MMatrix{3, segments}(zeros(3, segments))]
    
    # Pack parameters in param named tuple - false sets res! for in-place solution
    param = (kite_pos=kite_pos, kite_vel=kite_vel, wind_vel=wind_vel, 
         tether_length=tether_length, settings=settings, buffers=buffers, segments = segments, 
         return_result=false)
    # Define the nonlinear problem
    prob = NonlinearProblem(res!, state_vec, param)
    # Solve the problem with TrustRegion method
    sol = solve(prob, TrustRegion(autodiff=AutoFiniteDiff()), show_trace=Val(false)) 

    iterations = sol.stats.nsteps  # Field name may vary; verify with `propertynames(sol)`
    state_vec = sol.u
    if prn
        println("Iterations: ", iterations)
    end
    # Set the return_result to true so that res! returns outputs
    param = (; param..., return_result=true)
    res = MVector(0.0, 0, 0)
    res, force_kite, tether_pos, p0 = res!(res, state_vec, param)

    force_gnd = state_vec[3]
    state_vec, tether_pos, force_gnd, force_kite, p0
end


"""
    res!(res, state_vec, param)

Calculates difference between tether end and kite given tether ground segment orientation 
and magnitude.

# Arguments
- res::Vector{Float64} difference between tether end and kite segment
- state_vec::MVector{3, Float64} state vector (theta [rad], phi [rad], Tn [N]);
  tether orientation and tension at ground station
- par:: 7-elements tuple:
    - kite_pos::MVector{3, Float64} kite position vector in wind reference frame
    - kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
    - wind_vel::MMatrix{Float64} wind velocity vector in wind reference frame for each segment of the tether
    - tether_length: tether length
    - settings:: Settings struct containing environmental and tether parameters: see [Settings](@ref)
    - buffers:: (5, ) Vector{Matrix{Float64}}  Vector of (3, segments) Matrix{Float64} empty matrices for preallocation
    - segments:: number of tether segments
    - return_result:: Boolean to determine use for in-place optimization or for calculating returns

# Returns (if return_result==true)
- res::Vector{Float64} difference between tether end and kite segment
- T0::Vector{Float64} force from the kite to the end of tether
- pj:: (3, segments) Matrix{Float64} x,y,z - coordinates of the tether nodes
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
    kite_pos, kite_vel, wind_vel, tether_length, settings, buffers, segments, return_result = param
    g = abs(settings.g_earth[3])
    Ls = tether_length / (segments + 1)
    drag_coeff = -0.5 * settings.rho * Ls * settings.d_tether * settings.cd_tether
    A = π/4 * (settings.d_tether/1000)^2
    mj = settings.rho_tether * Ls * A
    E = settings.c_spring / A

    # Preallocate arrays
    FT = buffers[1]
    Fd = buffers[2]
    pj = buffers[3]
    vj = buffers[4]
    aj = buffers[5]

    # Unpack state variables (elevation, azimuth)
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
    FT[1, segments] = Tn * cosθ * cosφ # cos(elevation)cos(azimuth)
    FT[2, segments] = Tn * cosθ * sinφ # cos(elevation)sin(azimuth)
    FT[3, segments] = Tn * sinθ        # sin(elevation)

    pj[1, segments] = Ls * cosθ * cosφ
    pj[2, segments] = Ls * cosθ * sinφ
    pj[3, segments] = Ls * sinθ


    # Velocity and acceleration calculations
    ω = cross(kite_pos / norm_p^2, kite_vel)
    a = cross(ω, SVec3(pj[:, segments]))         
    b = cross(ω, cross(ω, SVec3(pj[:, segments])))
    vj[:, segments] .= v_parallel * p_unit + a
    aj[:, segments] .= b

    # Drag calculation for first element
    v_a_p1 = vj[1, segments] - wind_vel[1, segments]
    v_a_p2 = vj[2, segments] - wind_vel[2, segments]
    v_a_p3 = vj[3, segments] - wind_vel[3, segments]

    if all(x -> abs(x) < 1e-3, (v_a_p1, v_a_p2, v_a_p3))
        Fd[:, segments] .= 0.0
    else
        dir1, dir2, dir3 = pj[1, segments]/Ls, pj[2, segments]/Ls, pj[3, segments]/Ls
        v_dot_dir = v_a_p1*dir1 + v_a_p2*dir2 + v_a_p3*dir3
        v_a_p_t1 = v_dot_dir * dir1
        v_a_p_t2 = v_dot_dir * dir2
        v_a_p_t3 = v_dot_dir * dir3

        v_a_p_n1 = v_a_p1 - v_a_p_t1
        v_a_p_n2 = v_a_p2 - v_a_p_t2
        v_a_p_n3 = v_a_p3 - v_a_p_t3

        norm_v_a_p_n = sqrt(v_a_p_n1^2 + v_a_p_n2^2 + v_a_p_n3^2)
        coeff = drag_coeff * norm_v_a_p_n

        Fd[1, segments] = coeff * v_a_p_n1
        Fd[2, segments] = coeff * v_a_p_n2
        Fd[3, segments] = coeff * v_a_p_n3
    end

    # Process other segments
    @inbounds for ii in segments:-1:2
        # Tension force calculations
        if ii == segments
            mj_total = 1.5mj
            g_term = mj_total * g
        else
            mj_total = mj
            g_term = mj * g
        end

        FT[:, ii-1] .= mj_total * aj[:, ii] + FT[:, ii] - Fd[:, ii]
        FT[3, ii-1] += g_term

        # Position calculations
        ft_norm = sqrt(FT[1, ii-1]^2 + FT[2, ii-1]^2 + FT[3, ii-1]^2)
        l_i_1 = (ft_norm/(E*A) + 1) * Ls
        ft_dir = FT[1, ii-1]/ft_norm, FT[2, ii-1]/ft_norm, FT[3, ii-1]/ft_norm

        pj[1, ii-1] = pj[1, ii] + l_i_1 * ft_dir[1]
        pj[2, ii-1] = pj[2, ii] + l_i_1 * ft_dir[2]
        pj[3, ii-1] = pj[3, ii] + l_i_1 * ft_dir[3]

        # Velocity and acceleration
        a = cross(ω, SVec3(pj[:, ii-1]))           
        b = cross(ω, cross(ω, SVec3(pj[:, ii-1])))
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

    p0 = MVector(pj[1,1] + l_i_1*T0_dir1, 
                 pj[2,1] + l_i_1*T0_dir2,
                 pj[3,1] + l_i_1*T0_dir3)

    res .= kite_pos - p0
    if return_result
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
- state_vec::MVector{3, Float64} state vector (theta [rad], phi [rad], Tn [N])  
  tether orientation and tension at ground station
- kite_pos::MVector{3, Float64} kite position vector in wind reference frame
- kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
- wind_vel::MMatrix{3, segments, Float64} wind velocity vector in wind reference frame for each segment of the tether
- tether_length: Float64 tether length
- settings::Settings struct containing environmental and tether parameters: see [Settings](@ref)
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

"""
    init_quasistatic(kite_pos, tether_length; kite_vel = nothing, segments = nothing, wind_vel = nothing, settings = nothing)

Initialize the quasi-static tether model providing an initial guess for the state vector based on the numerical solution of the catenary equation

# Arguments
- kite_pos::MVector{3, Float64} kite position vector in wind reference frame
- tether_length: Float64 tether length
- kite_vel::MVector{3, Float64} kite velocity vector in wind reference frame
- segments::Int number of tether segments
- wind_vel::MMatrix{3, segments, Float64} wind velocity vector in wind reference frame for each segment of the tether
- settings::Settings struct containing environmental and tether parameters: see [Settings](@ref)

# Returns
- state_vec::MVector{3, Float64} state vector (theta [rad], phi [rad], Tn [N])  
  tether orientation and tension at ground station
"""
function init_quasistatic(kite_pos, tether_length; kite_vel = nothing, segments = nothing, wind_vel = nothing, settings = nothing)
    # Some basic checks
    @assert isa(kite_pos, MVector{3}) || error("kite_pos must be a MVector of size (3,1)")
    if isnothing(kite_vel) 
        kite_vel = MVector{3}([0.0, 0.0, 0.0])
    end

    if isnothing(segments) && isnothing(wind_vel)
        segments = 7
        wind_vel = zeros(3, segments)
    elseif isnothing(segments) && !isnothing(wind_vel)
        @assert size(wind_vel)[1] == 3 || error("wind_vel should have 3 rows!")
        segments = size(wind_vel)[2]
    elseif !isnothing(segments) && isnothing(wind_vel)
        wind_vel = zeros(3, segments)
    elseif !isnothing(segments) || !isnothing(wind_vel)
        @assert size(wind_vel)[1] == 3 || error("wind_vel should have 3 rows!")
        @assert size(wind_vel)[2] == segments || error("wind_vel should have the same number of columns as segments!")
    end

    if isnothing(settings)
        settings = Settings()
    else
        @assert typeof(settings) == Settings || error("settings should be of type Settings!")
    end

    kite_dist = norm(kite_pos)

    # azimuth angle calculation
    phi_init = atan(kite_pos[2], kite_pos[1])        
    
    # tension definition
    tension = 0.0002*settings.c_spring
    function solve_catenary(kite_pos, tether_length, segments)  
        hvec = kite_pos[1:2]    
        h = norm(hvec)
        v = kite_pos[3]
    
        # Function for nonlinear solver
        function f!(res, coeff, param)    
            tether_length, v, h = param  
            res[] = sqrt(tether_length^2 - v^2) - (2 * sinh(coeff[] * h / 2) / coeff[])       
        end
    
        u0 = [0.1]  # Initial guess as a scalar in an array
    
        # Define and solve nonlinear problem
        param = (tether_length, v, h)
        prob = NonlinearProblem(f!, u0, param)
        coeff = solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff()), show_trace=Val(false)) 
        coeff_val = coeff[]  # Extract scalar value
        
        # Adjust catenary solution to specific case
        X = LinRange(0, h, segments)
        angle1 = atan(hvec[1], hvec[2])
        XY = [sin(angle1) * X'; cos(angle1) * X']
        x_min = -(1 / 2) * (log((tether_length + v) / (tether_length - v)) / coeff_val - h)
        bias = -cosh(-x_min * coeff_val) / coeff_val    
        # Compute z-coordinates of catenary
        z_catenary = cosh.((X .- x_min) .* coeff_val) ./ coeff_val .+ bias
        x_catenary = XY[1, :]
        y_catenary = XY[2, :]    
        return x_catenary, y_catenary, z_catenary
    end

    # Solve the catenary equation
    x_catenary, y_catenary, z_catenary = solve_catenary(kite_pos, tether_length, segments)  
    # Calculate the elevation angle
    theta_init = atan(z_catenary[2], sqrt(x_catenary[2]^2 + y_catenary[2]^2))    

    # Assemble state vector
    state_vec = MVector{3}([theta_init, phi_init, tension])        
    
    return state_vec, kite_pos, kite_vel, wind_vel, tether_length, settings
end

"""
    get_analytic_catenary(filename)

Loads the analytic catenary curve for the 2D catenary example

# Arguments
- filename: the filename of the mat file to read

# Returns
- x_cat: x coordinates of the catenary curve
- y_cat: x coordinates of the catenary curve
"""
function get_analytic_catenary(filename)
    vars        = matread(filename)
    vars        = get(vars, "analytic_catenary", 0)
    x_cat       = vec(get(vars, "x", 0))
    y_cat       = vec(get(vars, "y", 0))
    return x_cat, y_cat
end

function transformFromOtoW(windDirection_rad,vec_O)
    M_WO = [cos(windDirection_rad) sin(windDirection_rad) 0; sin(windDirection_rad) -cos(windDirection_rad) 0;  0 0 -1]
    vec_W = M_WO*vec_O
    return vec_W
end

function transformFromWtoO(windDirection_rad,vec_W)
    M_OW = [cos(windDirection_rad) sin(windDirection_rad) 0; sin(windDirection_rad) -cos(windDirection_rad) 0; 0 0 -1]
    vec_O = M_OW*vec_W
    return vec_O
end