#!/usr/bin/python3

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

# Transform Cartesian to spherical
def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    return np.array([r, theta, phi])

# Transform spherical to Cartesian coordinates
def sph2cart(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x, y, z])

filename = '../output/HPK_full_1.00_d1_AvW45.root'

file = uproot.open(filename) #"../output/test.root")

data = file["ntp"].arrays(["photon_initial_px", "photon_initial_py", "photon_initial_pz"], library='pd')
data['r'], data["theta"], data["phi"] = cart2sph(data.photon_initial_px, data.photon_initial_py, data.photon_initial_pz)
print(data)
plt.hist(data.phi, label=r'$\phi$')
plt.xlabel(r'$\varphi$')
plt.title(r'Initial momentum $\varphi$ of the generated photons')
plt.show()
plt.hist(np.cos(data.theta), label=r'$\cos \theta$')
plt.xlabel(r'$\cos \theta$')
plt.title(r'Initial momentum $\cos \theta$ of the generated photons')
plt.show()