using MakieControlPlots, GLMakie

x = 1:5
y = 1:5
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, x, y; color=:black)
scatter!(ax, x, y)
hidedecorations!(ax)
hidespines!(ax)
display(fig)
