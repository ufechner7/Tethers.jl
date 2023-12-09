# Example two: Falling mass, attached to non-linear spring without compression stiffness,
# initially moving upwards with 4 m/s.
using ModelingToolkit, OrdinaryDiffEq, PyPlot, LinearAlgebra

G_EARTH  = Float64[0.0, 0.0, -9.81]    # gravitational acceleration [m/s²]
L0 = -10.0                             # initial spring length      [m]
V0 = 4                                 # initial velocity           [m/s]

# model, Z component upwards
@parameters mass=1.0 c_spring=50.0 damping=0.5 l0=L0
@variables t pos(t)[1:3] = [0.0, 0.0,  L0]
@variables   vel(t)[1:3] = [0.0, 0.0,  V0] 
@variables   acc(t)[1:3] = G_EARTH
@variables unit_vector(t)[1:3]  = [0.0, 0.0, -sign(L0)]
@variables spring_force(t)[1:3] = [0.0, 0.0, 0.0]
@variables force(t) = 0.0 norm1(t) = abs(l0) spring_vel(t) = 0.0
D = Differential(t)

eqs = vcat(D.(pos)      ~ vel,
           D.(vel)      ~ acc,
           norm1        ~ norm(pos),
           unit_vector  ~ -pos/norm1,         # direction from point mass to origin
           spring_vel   ~ -unit_vector ⋅ vel,
           spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector,
           acc          ~ G_EARTH + spring_force/mass)

@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

duration = 10.0
dt = 0.02
tol = 1e-6
tspan = (0.0, duration)
ts    = 0:dt:duration
u0 = zeros(6)
u0[3] = L0
u0[6] = V0

prob = ODEProblem(simple_sys, u0, tspan)
@time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, tstops=ts, saveat=ts)

X = sol.t
POS_Z = sol(X, idxs=pos[3])
VEL_Z = sol(X, idxs=vel[3])

plot(X, POS_Z, color="green")
xlabel("time [s]")
grid(true)
twinx()
ylabel("vel_z [m/s]") 
plot(X, VEL_Z, color="red") 
PyPlot.show(block=false)
nothing

#= # Calculate the vector res1, that depends on the velocity and the acceleration.
# The drag force of each segment is distributed equaly on both particles.
function calc_res(s::KPS3, pos1, pos2, vel1, vel2, mass, veld, result, i)
    s.segment .= pos1 - pos2
    height = (pos1[3] + pos2[3]) * 0.5
    rho = calc_rho(s.am, height)               # calculate the air density
    rel_vel = vel1 - vel2                # calculate the relative velocity
    s.av_vel .= 0.5 * (vel1 + vel2)
    norm1 = norm(s.segment)
    s.unit_vector .= normalize(s.segment) # unit vector in the direction of the tether
    # # look at: http://en.wikipedia.org/wiki/Vector_projection
    # # calculate the relative velocity in the direction of the spring (=segment)
    spring_vel = s.unit_vector ⋅ rel_vel

    k2 = 0.05 * s.c_spring * s.stiffness_factor             # compression stiffness tether segments
    if norm1 - s.segment_length > 0.0
        s.spring_force .= (s.c_spring * s.stiffness_factor * (norm1 - s.segment_length) + s.damping * spring_vel) .* s.unit_vector
    else
        s.spring_force .= k2 * ((norm1 - s.segment_length) + (s.damping * spring_vel)) .* s.unit_vector
    end
    s.seg_area = norm1 * s.set.d_tether/1000.0
    s.last_v_app_norm_tether = calc_drag(s, s.av_vel, s.unit_vector, rho, s.v_app_perp, s.seg_area)
    s.force .= s.spring_force + 0.5 * s.last_tether_drag

    if i == s.set.segments+1 # add the drag of the bridle lines
        s.bridle_area =  s.set.l_bridle * s.set.d_line/1000.0
        s.last_v_app_norm_tether = calc_drag(s, s.av_vel, s.unit_vector, rho, s.v_app_perp, s.bridle_area)
        s.force .+= s.last_tether_drag  
    end
   
    s.total_forces .= s.force + s.last_force
    s.last_force .= 0.5 * s.last_tether_drag - s.spring_force
    s.forces[i] .= s.total_forces
    acc = s.total_forces ./ mass # create the vector of the spring acceleration
    # result .= veld - (s.acc + SVector(0,0, -G_EARTH)) # Python code, wrong
    result .= veld - (SVector(0, 0, -G_EARTH) - acc)
    nothing
end =#