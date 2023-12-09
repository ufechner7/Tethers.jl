# This is the old code from KiteSim.jl that simulates the tether numerically without
# using ModelingToolkit, just as reference.

# Calculate the vector res1, that depends on the velocity and the acceleration.
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
end