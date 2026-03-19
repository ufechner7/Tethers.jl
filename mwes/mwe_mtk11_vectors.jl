# Diagnostic script to understand how MTK 11 handles vector equations.
# Run with: include("mwes/mwe_mtk11_vectors.jl")
using ModelingToolkit, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

println("=== MTK version: ", pkgversion(ModelingToolkit))

# ---- helpers ----------------------------------------------------------------
function show_eqs(label, eqs)
    println("\n--- $label ---")
    println("  typeof: ", typeof(eqs))
    println("  length: ", length(eqs))
    for (i, eq) in enumerate(eqs)
        println("  [$i] ", typeof(eq), "  =>  ", eq)
    end
end

# ---- small 2-particle, 1-segment model --------------------------------------
@parameters mass=1.0 c_spring=50.0 damping=0.5 l_seg=5.0
@variables pos(t)[1:3, 1:2] = zeros(3, 2)
@variables vel(t)[1:3, 1:2] = zeros(3, 2)
@variables acc(t)[1:3, 1:2]
@variables segment(t)[1:3]
@variables unit_vector(t)[1:3]
@variables norm1(t)
@variables rel_vel(t)[1:3]
@variables spring_vel(t)
@variables c_spr(t)
@variables spring_force(t)[1:3]

G_EARTH = [0.0, 0.0, -9.81]

# Step 1: differential equations (array equations via broadcasting)
eqs1 = vcat(D.(pos) .~ vel,
            D.(vel) .~ acc)
show_eqs("eqs1  (D.(pos) .~ vel  etc.)", eqs1)

# Step 2: vcat(eqs1...)
eqs2 = vcat(eqs1...)
show_eqs("eqs2 = vcat(eqs1...)", eqs2)

# Step 3: a mixed list of scalar + array equations (like the inner loop)
eqs_loop = [segment       ~ pos[:, 2] - pos[:, 1],   # array equation
            norm1          ~ norm(segment),             # scalar equation
            unit_vector    ~ -segment / norm1,          # array equation
            rel_vel        ~ vel[:, 2] - vel[:, 1],    # array equation
            spring_vel     ~ -unit_vector ⋅ rel_vel,   # scalar equation
            c_spr          ~ c_spring * (norm1 > l_seg), # scalar equation
            spring_force   ~ (c_spr * (norm1 - l_seg) + damping * spring_vel) * unit_vector, # array
            acc[:, 2]      ~ G_EARTH .+ spring_force / mass,  # array equation
            acc[:, 1]      ~ zeros(3)]                         # array equation

show_eqs("eqs_loop  (mixed scalar+array)", eqs_loop)

# Step 4: reduce(vcat, eqs_loop) -- what type do we get?
eqs_reduced = reduce(vcat, eqs_loop)
show_eqs("reduce(vcat, eqs_loop)", eqs_reduced)

# Step 5: Symbolics.scalarize.() applied element-wise
eqs_scalarized = Symbolics.scalarize.(eqs_loop)
println("\n--- Symbolics.scalarize.(eqs_loop) ---")
println("  typeof: ", typeof(eqs_scalarized))
for (i, e) in enumerate(eqs_scalarized)
    println("  [$i] typeof=", typeof(e), "  eltype=", eltype(e))
end

# Step 6: reduce(vcat, Symbolics.scalarize.(eqs_loop)) -- the final pattern
eqs_final = reduce(vcat, Symbolics.scalarize.(eqs_loop))
show_eqs("reduce(vcat, Symbolics.scalarize.(eqs_loop))", eqs_final)

# Step 7: root cause - eqs2 is Vector{BasicSymbolicImpl}, loop eqs are Vector{Equation}
# When vcat-ed, Julia upcasts to Vector{Any} → final reduce(vcat, scalarize.(...)) breaks.
println("\n=== ROOT CAUSE ANALYSIS ===")
println("eqs2 eltype:              ", eltype(eqs2))
println("eqs_loop eltype:          ", eltype(eqs_loop))
eqs_mixed = vcat(eqs2, eqs_loop)
println("vcat(eqs2, eqs_loop) type: ", typeof(eqs_mixed))

# Fix A: scalarize_to_vec always returning a vector, then explicitly convert to Vector{Equation}
scalarize_to_vec(e) = (r = Symbolics.scalarize(e); isa(r, AbstractVector) ? r : [r])

println("\n--- Fix A2: explicit convert(Vector{Equation}, ...) ---")
try
    eqs_fixA_raw = reduce(vcat, scalarize_to_vec.(vcat(eqs2, eqs_loop)))
    println("  raw type: ", typeof(eqs_fixA_raw))
    eqs_fixA = convert(Vector{Equation}, eqs_fixA_raw)
    println("  converted type: ", typeof(eqs_fixA), "  length: ", length(eqs_fixA))
    @named sys = ODESystem(eqs_fixA, t)
    println("  ODESystem SUCCESS, equations: ", length(equations(sys)))
catch e; println("  FAILED: ", e); end

# Fix B2: column-wise differential equations (avoids 2D D(pos) issue)
println("\n--- Fix B2: column-wise D(pos[:,i]) ~ vel[:,i] ---")
try
    eqs_diff = vcat([D(pos[:, i]) ~ vel[:, i] for i in axes(pos, 2)],
                    [D(vel[:, i]) ~ acc[:, i] for i in axes(vel, 2)])
    println("  typeof(eqs_diff): ", typeof(eqs_diff), "  length: ", length(eqs_diff))
    println("  eltype: ", eltype(eqs_diff))
    eqs_fixB2 = vcat(eqs_diff, eqs_loop)
    println("  vcat type: ", typeof(eqs_fixB2))
    eqs_final_B2 = reduce(vcat, Symbolics.scalarize.(eqs_fixB2))
    println("  after scalarize/reduce type: ", typeof(eqs_final_B2), "  length: ", length(eqs_final_B2))
    @named sys = ODESystem(eqs_final_B2, t)
    println("  ODESystem SUCCESS, equations: ", length(equations(sys)))
catch e; println("  FAILED: ", e); end

# Fix C: build eqs2 typed as Vector{Equation} by splatting into typed comprehension
println("\n--- Fix C: Equation[eqs1[:]...] to retype the differential equations ---")
try
    eqs_diff_typed = Equation[eqs1[:]...]
    println("  typeof: ", typeof(eqs_diff_typed), "  length: ", length(eqs_diff_typed))
    eqs_fixC = vcat(eqs_diff_typed, eqs_loop)
    println("  vcat type: ", typeof(eqs_fixC))
    eqs_final_C = reduce(vcat, Symbolics.scalarize.(eqs_fixC))
    println("  after scalarize/reduce type: ", typeof(eqs_final_C), "  length: ", length(eqs_final_C))
    @named sys = ODESystem(eqs_final_C, t)
    println("  ODESystem SUCCESS, equations: ", length(equations(sys)))
catch e; println("  FAILED: ", e); end


