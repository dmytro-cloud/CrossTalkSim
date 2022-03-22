import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.array([0.1, 0.25, 0.5, 0.75, 1])
Y = np.array([1, 2, 3, 4, 5])
X, Y = np.meshgrid(X, Y)
Z = np.array([ 0.05425565958136484,
0.05329130811592519,
0.04810514640041961,
0.05053598849964893,
0.04595749718845093,
0.023302588204896477,
0.020312873035867134,
0.01669402427915583,
0.014047389938882042,
0.013602995159919442,
0.00490615016954555,
0.0030807482817091537,
0.00261875736379765,
0.0018197986588925166,
0.0012812058224497221,
0.002730737083766591,
6.716846757213841e-05,
3.3440920339535496e-05,
6.708785735002964e-05,
0.00017546414688991095,
0.0020709533795054236,
3.3423258644305804e-05,
0.0002085695024941351,
3.351041490872946e-05,
0.0,]).reshape(5, 5) * 100
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
plt.title('HPK Crosstalk Probability 0.1mkm AvD')
plt.show()



