using LinearAlgebra

struct Settings 
    rho::Float64
    g_earth::Vector{Float64}
    cd_tether::Float64
    d_tether::Float64                                 # tether diameter                  [mm]
    rho_tether::Float64                             # density of Dyneema            [kg/m³]
    c_spring::Float64   
end

stts = Settings(1.225, [0; 0; -9.81], 0.9, 4.0, 750.0, 614600.0)
stateVec = [0; 0; 0]
kitePos = [10; 10; 40]
kiteVel = [1; 1; 1]
windVel =  [10 10 10 10 10;
            10 10 10 10 10;
            10 10 10 10 10]

tetherLength = 100

function objFun(stateVec, kitePos, kiteVel, windVel, tetherLength, settings)

    # settings follows same syntax as Settings3 in Tether_09.jl
    g  = settings.g_earth[3]   # in this function, g is considered a scalar 

    Ns = size(windVel, 2)           # number of masses - windVel is a 3xNs matrix
    Ls = tetherLength/(Ns+1)        # segment length
    mj = settings.rho_tether*Ls     # segment mass

    p = kitePos                     # input kite position
    v = kiteVel                     # input kite velocity

    # Extract states
    theta   = stateVec[1]
    phi     = stateVec[2]
    Tn      = stateVec[3]

    # Compute omega_t
    omega_t = cross(p/(norm(p)^2),v);

    # Tether cross section
    A = pi/4*(settings.d_tether/1000)^2 # [m2] 
    # Compute Young's modulus    
    E = settings.c_spring/A

    FT = zeros(3,Ns) # Tension forces
    Fd = zeros(3,Ns) # Drag forces
    pj = zeros(3,Ns) # Mass positions
    vj = zeros(3,Ns) # Mass velocities
    aj = zeros(3,Ns) # Mass accellerations

    # First element from ground station (Ns)
    FT[:,Ns] = Tn.*[sin(theta)* cos(phi); sin(phi); cos(theta)*cos(phi)]
    pj[:,Ns] = Ls.*[sin(theta)* cos(phi); sin(phi); cos(theta)*cos(phi)]
    vj[:,Ns] = dot(v,p/norm(p)) .* (p/norm(p)) + cross(omega_t,pj[:,Ns])
    aj[:,Ns] = cross(omega_t,cross(omega_t,pj[:,Ns]))

    # Drag calculation first element
    v_a_p = vj[:,Ns] - windVel[:,Ns]
    if all(abs.(v_a_p) .< 1e-3)
        Fd[:,Ns] = [0;0;0]
    else
        v_a_p_t = (dot((pj[:,Ns])/norm(pj[:,Ns]),v_a_p)*
            ((pj[:,Ns])/norm(pj[:,Ns])))
        v_a_p_n = v_a_p - v_a_p_t
        Fd[:,Ns] = (-0.5 * settings.rho * Ls * settings.d_tether * settings.cd_tether * 
            norm(v_a_p_n) * v_a_p_n )# particle drag
    end

    # All other segments and masses except for segment connected to the kite
    for ii = Ns:-1:2
        if ii == Ns
            FT[:,ii-1] = (mj+0.5*mj)*aj[:,ii] + FT[:,ii] - Fd[:,ii] + [0;0;(mj+0.5*mj)*g]
        else
            FT[:,ii-1] = mj*aj[:,ii] + FT[:,ii] - Fd[:,ii] + [0;0;mj*g]
        end

        l_i_1 = (norm(FT[:,ii-1])/(E*A) + 1)*Ls
        
        pj[:,ii-1] = pj[:,ii] + l_i_1.*(FT[:,ii-1]/norm(FT[:,ii-1]))       
        vj[:,ii-1] = dot(v,p./norm(p)) * (p./norm(p)) + cross(omega_t,pj[:,ii-1])
        aj[:,ii-1] = cross(omega_t,cross(omega_t,pj[:,ii-1]))
        
        # Drag calculation
        v_a_p = vj[:,ii] - windVel[:,ii]
        if all(abs.(v_a_p) .< 1e-3)
            Fd[:,ii-1] = [0;0;0]
        else
            v_a_p_t = (dot((pj[:,ii-1]-pj[:,ii])/norm(pj[:,ii-1]-pj[:,ii]),v_a_p) * 
                ((pj[:,ii-1]-pj[:,ii])/norm(pj[:,ii-1]-pj[:,ii])))
            v_a_p_n = v_a_p - v_a_p_t;
            Fd[:,ii-1] = -0.5 * settings.rho * Ls * settings.d_tether * settings.cd_tether * norm(v_a_p_n) * v_a_p_n # particle drag
        end
    end

    T0 = (mj+0.5*mj).*aj[:,1] + FT[:,1] - Fd[:,1] + [0;0;(mj+0.5*mj)*g];
    l_i_1 = (norm(T0)/(E*A) + 1)*Ls;
    p0 = pj[:,1] + l_i_1.*(T0/norm(T0));
    
    Fobj = p-p0;      

    return Fobj, T0, pj, p0
end

Fobj, T0, pj, p0 = objFun(stateVec, kitePos, kiteVel, windVel, 
                            tetherLength, stts)