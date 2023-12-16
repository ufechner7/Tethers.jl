using ModelingToolkit, LinearAlgebra

@variables t pos(t)[1:3] = [0.0, 0.0, 10.0]

# working
# pos/norm(pos)

# not working
@register_symbolic LinearAlgebra.normalize(pos)

normalize(pos)