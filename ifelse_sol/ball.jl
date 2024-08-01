using ModelingToolkit, OrdinaryDiffEq, Plots

## Test bouncing ball with equation affect
@variables t x(t)=1 v(t)=0 d_x(t)=0
@variables switch(t)=1
D = Differential(t)

# root_eqs = [x ~ 0]
# affect   = [v ~ -v, switch ~ 1-switch]
continuous_events = [
    [x ~ 0] => [v ~ -v]
    # [x ~ 0.5] => [v ~ 2*v]
]

d_x=v+x
# v -= x

eqs = [
    D(x) ~ d_x
    # v~v+0.1
    D(v) ~ -9.8
    D(switch) ~ 0
]

@named ball = ODESystem(eqs, t; continuous_events = continuous_events)

ball = structural_simplify(ball)

tspan = (0.0,5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob,TRBDF2())
print(0 <= minimum(sol[x]) <= 1e-10) # the ball never went through the floor but got very close
plot(sol)