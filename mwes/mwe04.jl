using ControlPlots

x = 1:5
y = 1:5
plt.plot(x,y, color="black")
plt.scatter(x,y)
ax = plt.gca()
ax.set_axis_off()
plt.show()

