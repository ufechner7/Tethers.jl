using ModelingToolkit, OrdinaryDiffEq, Plots

## Test bouncing ball with equation affect
@variables t x(t)=1 v(t)=0
@variables switch(t)=1
D = Differential(t)

root_eqs = [x ~ 0]
affect   = [v ~ -v, switch ~ 1-switch]
# root_eqs2 = [x ~ 0.5]
# affect2 = [switch ~ 1 - switch]

@named ball = ODESystem([
    D(x) ~ v
    D(v) ~ -9.8
    D(switch) ~ 0
], t; continuous_events = [root_eqs => affect])

ball = structural_simplify(ball)

tspan = (0.0,5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob,Tsit5())
print(0 <= minimum(sol[x]) <= 1e-10) # the ball never went through the floor but got very close
plot(sol)