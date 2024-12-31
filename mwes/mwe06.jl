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
    # Sort by magnitude of moments
    p = sortperm(eigenvals)
    return eigenvals[p], eigenvecs[:, p]
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
y_values = [-1.0, 0.0, 1.0, 2.0, 3.0]  # Example y-coordinates
z_values = [1.0, 0.0, 0.0, 0.0, 1.0]   # Example z-coordinates

# Calculate center of mass
com_x = sum(m_values .* x_values) / total_mass
com_y = sum(m_values .* y_values) / total_mass
com_z = sum(m_values .* z_values) / total_mass

# Update values centered around origin
m_values = m_values  # Example masses
x_values .-= com_x  # Center x coordinates
y_values .-= com_y  # Center y coordinates
z_values .-= com_z   # Center z coordinates

# Calculate total mass
total_mass = sum(m_values)

# Calculate and apply transformation
principal_moments, Q = get_principal_axes(m_values, x_values, y_values, z_values)

# # Pre-calculate inertia
# I_values = calculate_inertia(m_values, x_values, y_values, z_values)
# @show I_values
# @show calculate_full_inertia(m_values, x_values, y_values, z_values)

@variables x(t)[1:3] v(t)[1:3] q(t)[1:4] ω(t)[1:3] τ(t)[1:3] z(t)[1:3]
@parameters m F[1:3,1:n] r[1:3,1:n] I[1:3] c_t c_r

# Function to rotate vector by quaternion
function rotate_by_quaternion(vec, q)
    w, x, y, z = q[1], q[2], q[3], q[4]
    [1-2*(y^2+z^2)    2*(x*y-w*z)     2*(x*z+w*y);
     2*(x*y+w*z)      1-2*(x^2+z^2)   2*(y*z-w*x);
     2*(x*z-w*y)      2*(y*z+w*x)     1-2*(x^2+y^2)] * vec
end

# Quaternion multiplication matrix
Ω = [0       -ω[1]   -ω[2]   -ω[3];
     ω[1]    0       ω[3]    -ω[2];
     ω[2]    -ω[3]   0       ω[1];
     ω[3]    ω[2]    -ω[1]   0]


eqs_x = [D(x[i]) ~ v[i] for i in 1:3]
F_total = sum(F[:, i] for i in 1:n)
eqs_v = [D(v[i]) ~ F_total[i]/m - c_t*v[i] for i in 1:3]

τ = sum(cross(rotate_by_quaternion(r[:, i], q), F[:, i]) for i in 1:n)

# Normalize quaternion in dynamics
function normalize_quaternion(q)
    norm = sqrt(sum(q[i]^2 for i in 1:4))
    return [q[i]/norm for i in 1:4]
end

error_constant = 0.1
qnorm_error = 1 - sum(q[i]^2 for i in 1:4)
eqs_q = [D(q[i]) ~ 0.5 * sum(Ω[i, j] * q[j] for j in 1:4) + 0.5*qnorm_error*q[i]/error_constant for i in 1:4]

# Rotational dynamics (Euler's equations)
@show I[2]
eqs_ω = [
    D(ω[1]) ~ (τ[1] + (I[2] - I[3]) * ω[2] * ω[3]) / I[1] - c_r*ω[1]
    D(ω[2]) ~ (τ[2] + (I[3] - I[1]) * ω[3] * ω[1]) / I[2] - c_r*ω[2]
    D(ω[3]) ~ (τ[3] + (I[1] - I[2]) * ω[1] * ω[2]) / I[3] - c_r*ω[3]
    z ~ rotate_by_quaternion([0, 0, 1], q)
]

# Combine all equations
eqs = [eqs_x; eqs_v; eqs_q; eqs_ω]

@mtkbuild sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs)), t)

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
F_values[3, :] .= -9.81  # Gravity on each point
F_values[3, end] -= 1
F_principal = [to_principal(F_values[:,i], Q) for i in 1:n]
F_principal = hcat(F_principal...)

p = [m => sum(m_values), 
     F => F_principal, 
     r => r_principal,  # transformed coordinates
     I => principal_moments,  # diagonal moments only
     c_t => 0.0,
     c_r => 0.0]


# Create and solve the ODE problem
dt = 0.02
tol = 1e-6
tspan = (0.0, 100.0)
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
