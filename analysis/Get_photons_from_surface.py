#!/usr/bin/python

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

filename = sys.argv[1]
filename_params = sys.argv[2]
file_pararms = open(filename_params, 'r')
params = []
for line in file_pararms:
	params.append(float(line))
file_pararms.close()

PitchWidth = params[0]
AvWidth = params[1]
trench_width = params[2]
trench_depth = params[3]
si_thickness = params[4]
d_star = 0.000005 # params[5]
electron_width = 0.000145 #0.001 - HPK param # params[6]
w_star = 0.0022 # 0.004 - HPK # params[7]
holes_width = w_star - electron_width

file = uproot.open(filename) #"../output/test.root")

data = file["ntp"].arrays(["second_last_x", "second_last_y", "second_last_z",
 "photon_last_x", "photon_last_y", "photon_last_z", "photon_initial_wl"], library='pd')

cut = data["photon_last_z"] < - si_thickness / 2 - 0.01
data_cut = data[cut]

# Find the direction of the selected events and make a cut on theta angle

data_cut['r_x'] = data_cut.photon_last_x - data_cut.second_last_x
data_cut['r_y'] = data_cut.photon_last_y - data_cut.second_last_y
data_cut['r_z'] = data_cut.photon_last_z - data_cut.second_last_z

# Surface n is (0, 0, -1) since the front side is at the negative end of coordinte system
data_cut['cos_theta'] = - data_cut.r_z / np.sqrt(data_cut.r_x**2 + data_cut.r_y**2 + data_cut.r_z**2)

# Assume 45 degrees apperture
plt.hist(data_cut.cos_theta, bins=100)
plt.show()

data_cut = data_cut[data_cut['cos_theta'] > np.sqrt(2) / 2]

plt.hist(data_cut.photon_initial_wl, bins=100)
plt.xlabel('Wavelength [nm]')
plt.show()

print(len(data))
print(len(data_cut))
plt.hist(data_cut.second_last_z, bins=100, range=(- si_thickness / 2 - 0.01, 0))
plt.xlabel('Z, [mm]')
plt.title('Depth distribution')
plt.show()
plt.hist2d(data_cut.second_last_x, data_cut.second_last_y, range=( (-PitchWidth, PitchWidth),
(-PitchWidth, PitchWidth)), norm=mpl.colors.LogNorm(), bins=50)
plt.xlabel('X, [mm]')
plt.ylabel('Y, [mm]')
plt.title('One cell')
plt.colorbar()
plt.show()

plt.hist2d(data_cut.second_last_x, data_cut.second_last_y, norm=mpl.colors.LogNorm(), bins=100)
plt.colorbar()
plt.title('The whole board')
plt.xlabel('X, [mm]')
plt.ylabel('Y, [mm]')
plt.show()


# Convolution with nearby cells probability


print(np.full(100, 5))

