using OrdinaryDiffEq
using ModelingToolkit
# using Plots
import IfElse

@variables t
vars = @variables y(t)
pars = @parameters uMax uMin u
@variables ifCond1(t) = false
@variables ifCond2(t) = false
D = Differential(t)

continuous_events = [
    (uMax - t ~ 0) => [ifCond1 ~ true, ifCond2 ~ false]
    (uMin - t ~ 0) => [ifCond1 ~ false, ifCond2 ~ true]
]


eqs = [
  D(y) ~ IfElse.ifelse(ifCond1 == true, uMax, IfElse.ifelse(ifCond2 == true, uMin, u)),
  D(ifCond1) ~ 0,
  D(ifCond2) ~ 0,
]

@named ifequationDer = ODESystem(eqs, t; continuous_events)
ifequationDer = structural_simplify(ifequationDer)
tspan = (0.0,12.0)
prob = ODEProblem(ifequationDer,
                  [y => 0.0, ifCond1=>false, ifCond2=>false],
                  tspan,
                  [uMax => 10.0,
                   uMin => 2.0,
                   u => 4.0])
sol = solve(prob,Tsit5())

# plot(sol)