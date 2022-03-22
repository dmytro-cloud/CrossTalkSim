import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.array([0.1, 0.25, 0.5, 0.75, 1])
Y = np.array([1, 2, 3, 4, 5])
X, Y = np.meshgrid(X, Y)
Z = np.array([ 0.06164432840644401,
0.058714426519476975,
0.05679063202022547,
0.05400740154178384,
0.0505432447671994,
0.023008353838290585,
0.02201921752979639,
0.018763863578835808,
0.016887661877577637,
0.01611547076168048,
0.005414628555371925,
0.003739976224924085,
0.0027809198903344726,
0.002379524744731568,
0.0019169693817192215,
0.0016182134207806324,
0.0,
0.00041494117492872157,
6.754858477090933e-05,
3.3501043051346725e-05,
0.002179832894256533,
0.0004504200647729461,
0.0,
3.3619917286944506e-05,
0.0 ]).reshape(5, 5) * 100
# print(Z)

# Plot the surface.
# ax.scatter3D(X, Y, Z)
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)

# Customize the z axis.
ax.set_zlim(0, 7)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
cax = fig.add_axes([0.4, .87, 0.35, 0.03])

fig.colorbar(surf, cax=cax, shrink=0.5,orientation='horizontal', aspect=5)
ax.set_xlabel('Thickness of the trenches [mkm]')
ax.set_ylabel('Trenches Depth [mkm]')
ax.set_zlabel('Probability of crosstalks at 7V (%)')
plt.title('HPK Crosstalk Probability')
plt.show()



