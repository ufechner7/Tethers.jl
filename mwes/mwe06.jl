using ModelingToolkit, ControlPlots, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

# Pre-calculate moments of inertia
function calculate_inertia(m, x, y, z)
    n = length(m)
    Ix = sum(m[i] * (y[i]^2 + z[i]^2) for i in 1:n)
    Iy = sum(m[i] * (x[i]^2 + z[i]^2) for i in 1:n)
    Iz = sum(m[i] * (x[i]^2 + y[i]^2) for i in 1:n)
    return [Ix, Iy, Iz]
end

# Pre-calculate moments of inertia
function calculate_full_inertia(m, x, y, z)
    n = length(m)
    I = zeros(3,3)
    for i in 1:n
        # Diagonal terms
        I[1,1] += m[i] * (y[i]^2 + z[i]^2)  # Ixx
        I[2,2] += m[i] * (x[i]^2 + z[i]^2)  # Iyy
        I[3,3] += m[i] * (x[i]^2 + y[i]^2)  # Izz
        # Products of inertia
        I[1,2] = I[2,1] -= m[i] * x[i] * y[i]  # Ixy
        I[1,3] = I[3,1] -= m[i] * x[i] * z[i]  # Ixz
        I[2,3] = I[3,2] -= m[i] * y[i] * z[i]  # Iyz
    end
    return I
end

# Calculate full inertia tensor and find principal axes
function get_principal_axes(m, x, y, z)
    I = calculate_full_inertia(m, x, y, z)
    eigenvals, eigenvecs = eigen(I)
    
    # Sort by magnitude
    p = sortperm(eigenvals)
    eigenvals = eigenvals[p]
    eigenvecs = eigenvecs[:, p]
    
    # Ensure right-handed coordinate system
    if det(eigenvecs) < 0
        eigenvecs[:, 3] *= -1
    end
    
    return eigenvals, eigenvecs
end

function to_principal(v_B, Q)
    return Q * v_B
end

function to_body(v_P, Q)
    return Q' * v_P
end

# Example values
n = 5  # Number of squares
m_values = [1.0, 1.0, 1.0, 1.0, 1.0]  # Example masses
x_values = [1.0, 1.0, 1.0, 1.0, 1.0]  # Example x-coordinates
y_values = [1.0, 0.0, -1.0, -2.0, -3.0]  # Example y-coordinates
z_values = [-1.0, 0.0, 0.0, 0.0, -1.0]   # Example z-coordinates

# Calculate total mass
total_mass = sum(m_values)

# Calculate center of mass
com_x = sum(m_values .* x_values) / total_mass
com_y = sum(m_values .* y_values) / total_mass
com_z = sum(m_values .* z_values) / total_mass

# Update values centered around origin
m_values = m_values  # Example masses
x_values .-= com_x  # Center x coordinates
y_values .-= com_y  # Center y coordinates
z_values .-= com_z   # Center z coordinates

# Calculate and apply transformation
principal_moments, Q = get_principal_axes(m_values, x_values, y_values, z_values)

# q in body frame, ω in principal frame
@variables x(t)[1:3] v(t)[1:3] q(t)[1:4] ω(t)[1:3] τ(t)[1:3] z(t)[1:3] ω_body(t)[1:3] r_world(t)[1:3,1:n] v_world(t)[1:3,1:n]
@parameters m F_b[1:3,1:n] F_p[1:3,1:n] r_p[1:3,1:n] I[1:3] c_t c_r

# Function to rotate vector by quaternion
function rotate_by_quaternion(vec, q)
    w, x, y, z = q[1], q[2], q[3], q[4]
    [1-2*(y^2+z^2)    2*(x*y-w*z)     2*(x*z+w*y);
     2*(x*y+w*z)      1-2*(x^2+z^2)   2*(y*z-w*x);
     2*(x*z-w*y)      2*(y*z+w*x)     1-2*(x^2+y^2)] * vec
end

function conj(q)
    w, x, y, z = q[1], q[2], q[3], q[4]
    return [w, -x, -y, -z]
end

# Quaternion multiplication matrix
Ω = [0       -ω[1]   -ω[2]   -ω[3];
     ω[1]    0       ω[3]    -ω[2];
     ω[2]    -ω[3]   0       ω[1];
     ω[3]    ω[2]    -ω[1]   0]


eqs_x = [D(x[i]) ~ v[i] for i in 1:3]
F_total = sum(F_b[:, i] for i in 1:n)
eqs_v = [D(v[i]) ~ F_total[i]/m - c_t*v[i] for i in 1:3]

τ = sum(cross(rotate_by_quaternion(r_p[:, i], q), F_p[:, i]) for i in 1:n)

# error_constant = 1e3
# qnorm_error = 1 - sum(q[i]^2 for i in 1:4)
eqs_q = [D(q[i]) ~ 0.5 * sum(Ω[i, j] * q[j] for j in 1:4) for i in 1:4]
# eqs_q = [D(q[i]) ~ 0.5 * sum(Ω[i, j] * q[j] for j in 1:4) + 0.5*qnorm_error*q[i]/error_constant for i in 1:4]

# Rotational dynamics (Euler's equations)
@show I[2]
eqs_ω = [
    D(ω[1]) ~ (τ[1] + (I[2] - I[3]) * ω[2] * ω[3]) / I[1] - c_r*ω[1]
    D(ω[2]) ~ (τ[2] + (I[3] - I[1]) * ω[3] * ω[1]) / I[2] - c_r*ω[2]
    D(ω[3]) ~ (τ[3] + (I[1] - I[2]) * ω[1] * ω[2]) / I[3] - c_r*ω[3]
    z ~ rotate_by_quaternion([0, 0, 1], q)
    ω_body ~ Q' * ω
]

ω_world = rotate_by_quaternion(ω, q)
v_eqs = [
    [r_world[:,i] ~ x + rotate_by_quaternion(r_p[:,i], q) for i in 1:n]
    [v_world[:,i] ~ v + rotate_by_quaternion(cross(ω, r_p[:,i]), q) for i in 1:n]]

# Combine all equations
eqs = [v_eqs; eqs_x; eqs_v; eqs_q; eqs_ω]

@mtkbuild sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs)), t) # ; continuous_events = [0 ~ 0] => [q ~ q / norm(q)]

# Initial conditions
x0 = zeros(3)
v0 = zeros(3)
q0 = [1.0, 0.0, 0.0, 0.0]  # Initial quaternion (identity rotation)
ω0 = [0.0, 0.0, 0.0]       # Initial angular velocity
# z0 = zeros(3)
u0 = [x0; v0; q0; ω0]

# Parameters: pre-calculated moments of inertia and applied torques
r_body = [x_values y_values z_values]'
r_principal = [to_principal(r_body[:,i], Q) for i in 1:n]
r_principal = hcat(r_principal...)


F_values = zeros(3, n)
# F_values[3, :] .= -9.81  # Gravity on each point
F_values[1, end] += 1
F_principal = [to_principal(F_values[:,i], Q) for i in 1:n]
F_principal = hcat(F_principal...)

c_critical = 2.0 * sqrt(principal_moments[1])  # Using largest moment
println("Critical damping coefficient: ", c_critical)

p = [m => sum(m_values), 
     F_p => F_principal, 
     F_b => F_values,
     r_p => r_principal,  # transformed coordinates
     I => principal_moments,  # diagonal moments only
     c_t => 0.0,
     c_r => c_critical]


# Create and solve the ODE problem
dt = 0.02
tol = 1e-8
tspan = (0.0, 10)
ts    = 0:dt:tspan[2]
prob = ODEProblem(sys, u0, tspan, p)

sol = solve(prob, QBDF(autodiff=false),
    dt=dt, 
    abstol=tol, 
    reltol=tol,
    saveat=ts)
@time sol = solve(prob, QBDF(autodiff=false), 
    dt=dt, 
    abstol=tol, 
    reltol=tol,
    saveat=ts)
    
p = plotx(sol.t, 
    [sol[i, :] for i in 7:13]...;
    ylabels=string.(unknowns(sys))[7:13]
    )
display(p)

@show sol[ω_body][50]
@show principal_moments Q
nothing