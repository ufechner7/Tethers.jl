# Example one: Falling mass.
using Groebner, ModelingToolkit, Optimization, OptimizationOptimJL, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

G_EARTH::Vector{Float64} = [0.0, 0.0, -9.81]    # gravitational acceleration     [m/s²]
L0::Float64 = -10.0                             # initial spring length      [m]

# defing the model, Z component upwards
@parameters mass=1.0 c_spring=50.0 damping=0.5 l0=L0
@variables pos(t)[1:3] = [0.0, 0.0,  L0]
@variables vel(t)[1:3] = [0.0, 0.0,  0.0] 
@variables acc(t)[1:3]
@variables unit_vector(t)[1:3]
@variables spring_force(t)[1:3]
@variables norm1(t) spring_vel(t)
@variables total_acc(t)
@variables distance(t)

eqs = [
    pos          ~ distance * [0.0, 0.0, 1.0]
    vel          ~ [0.0, 0.0, 0.0]
    norm1        ~ norm(pos)
    unit_vector  ~ -pos / norm1         # direction from point mass to origin
    spring_vel   ~ -unit_vector ⋅ vel
    spring_force ~ (c_spring * (norm1 - abs(l0)) + damping * spring_vel) * unit_vector
    acc          ~ G_EARTH + spring_force/mass
    total_acc    ~ norm(acc)
]
eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

function merge_eqs(eqs, opt_var)
    found_vars = Dict{Num, Int}()
    subs = Dict{Num, Num}()
    eqs = Vector{Any}(eqs)
    allunknowns = Symbolics.OrderedSet()
    ps = Symbolics.OrderedSet()

    opt_i = 1
    # for i in eachindex(eqs)
    #     # eqs[i] = eqs[i].lhs - eqs[i].rhs
    #     ModelingToolkit.collect_vars!(allunknowns, ps, eqs[i], t)
    #     # for var in get_variables(eqs[i])
    #     #     if !(haskey(found_vars, var))
    #     #         found_vars[var] = i
    #     #         @show eqs[i] var typeof(var)
    #     #         @show solve_for(0 ~ eqs[i], var)
    #     #         eqs[i] = var ~ solve_for(0 ~ eqs[i], var)
    #     #         break
    #     #     end
    #     # end
    # end
    # println(allunknowns)
    # println(ps)
    # for i in eachindex(eqs)
    #     if eqs[i].lhs 
    # end
    for i in eachindex(eqs)
        @assert length(Symbolics.get_variables(eqs[i].lhs)) == 1
        if occursin(eqs[i].lhs, opt_var)
            opt_i = i
        end
        subs[eqs[i].lhs] = eqs[i].rhs
        eqs[i] = eqs[i].rhs
    end
    eqs[opt_i] = Symbolics.fixpoint_sub(eqs[opt_i], subs)
    return eqs[opt_i]
end

opt_eq = merge_eqs(eqs, total_acc)

@mtkbuild sys = OptimizationSystem(opt_eq, [distance], [mass, c_spring, l0])
prob = OptimizationProblem(sys, [L0], [mass => 1.0, c_spring => 50.0, l0 => L0], grad = true, hess = true)
@time sol = solve(prob, Optim.NelderMead(); maxiters=10_000)


