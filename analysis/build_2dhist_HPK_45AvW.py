import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.array([0.1, 0.25, 0.5, 0.75, 1])
Y = np.array([1, 2, 3, 4, 5])
X, Y = np.meshgrid(X, Y)
Z = np.array([ 0.07361370709130977,
0.0692044036772698,
0.06515310064525154,
0.060581932293235295,
0.056172580395264245,
0.027647347915331393,
0.024102657473828053,
0.02013471169093374,
0.018547968827322767,
0.015125969599530083,
0.006006230735636041,
0.004396620876297627,
0.003248692614151628,
0.0020022031985840742,
0.0017156403191523801,
0.0021149398389965703,
0.00024201944751032228,
0.00020768157450246015,
0.0,
0.0001344203842663807,
0.002029742917420471,
6.746027531154551e-05,
3.359432520247127e-05,
3.354896152252795e-05,
3.372335486551842e-05,]).reshape(5, 5) * 100
# print(Z)

# Plot the surface.
# ax.scatter3D(X, Y, Z)
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)

# Customize the z axis.
# ax.set_zlim(0, 7)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
cax = fig.add_axes([0.4, .87, 0.35, 0.03])

fig.colorbar(surf, cax=cax, shrink=0.5,orientation='horizontal', aspect=5)
ax.set_xlabel('Thickness of the trenches [mkm]')
ax.set_ylabel('Trenches Depth [mkm]')
ax.set_zlabel('Probability of crosstalks at 7V (%)')
plt.title('HPK Crosstalk Probability 45mkm Avalanche Width')
plt.show()



