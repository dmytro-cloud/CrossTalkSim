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

data = file["ntp"].arrays(["photon_initial_x", "photon_initial_y", "photon_initial_z",
 "photon_last_x", "photon_last_y", "photon_last_z"], library='pd')
# print(data)
plt.hist2d(data.photon_last_x, data.photon_last_y, bins=1000, norm=mpl.colors.LogNorm(), range=((-0.15, 0.15), (-0.15, 0.15)))
# plt.zscale('log')
plt.colorbar()
plt.xlabel('X, [mm]')
plt.ylabel('Y, [mm]')
plt.show()

# plt.hist2d(data.photon_last_z, data.photon_last_y, bins=1000, norm=mpl.colors.LogNorm())
# plt.colorbar()
# plt.xlabel('Z, [mm]')
# plt.ylabel('Y, [mm]')
# plt.show()

# 3d Plot
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# ax.scatter(data.photon_last_x[:10000], data.photon_last_y[:10000], data.photon_last_z[:10000], marker='.')

# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# ax.set_zlim(-0.15, -0.13)
# ax.set_xlim(-0.04, 0.04)
# ax.set_ylim(-0.05, -0.1)
# plt.show()

# FBK 
# PitchWidth = 0.035
# AvWidth = 0.032
# trench_width = 0.001
# trench_depth = 0.003
# si_thickness = 0.3

# z_cut_si = (data['photon_last_z'] < si_thickness / 2) & (data['photon_last_z'] > -si_thickness / 2) 
z_cut_tr_electron = (data['photon_last_z'] < -si_thickness / 2 + d_star + electron_width) & (data['photon_last_z'] > -si_thickness / 2 + d_star) 

x_shift = y_shift = PitchWidth - AvWidth / 2
height_cut = (  AvWidth/2 > data['photon_last_y'] )  &  ( data['photon_last_y'] > -AvWidth/2 )
right_cut = (data['photon_last_x'] > x_shift) & (data['photon_last_x'] < x_shift + AvWidth)
left_cut = (data['photon_last_x'] < -x_shift) & (data['photon_last_x'] > -x_shift - AvWidth)

width_cut = (  AvWidth/2 > data['photon_last_x'] )  &  ( data['photon_last_x'] > -AvWidth/2 )
top_cut = ((data['photon_last_y'] > y_shift) & (data['photon_last_y'] < y_shift + AvWidth))
bot_cut = ((data['photon_last_y'] < -y_shift) & (data['photon_last_y'] > -y_shift - AvWidth))

right_sens_electrons = data[right_cut & height_cut & z_cut_tr_electron]
left_sens_electrons = data[left_cut & height_cut & z_cut_tr_electron]
top_sens_electrons = data[top_cut & width_cut & z_cut_tr_electron]
bot_sens_electrons = data[bot_cut & width_cut & z_cut_tr_electron]

right_top_sens_electrons = data[top_cut & right_cut & z_cut_tr_electron]
left_top_sens_electrons = data[top_cut & left_cut & z_cut_tr_electron]
right_bot_sens_electrons = data[bot_cut & right_cut & z_cut_tr_electron]
left_bot_sens_electrons = data[bot_cut & left_cut & z_cut_tr_electron]

z_cut_tr_holes = (data['photon_last_z'] < -si_thickness / 2 + d_star + electron_width + holes_width) & (data['photon_last_z'] > -si_thickness / 2 + d_star + electron_width) 

right_sens_holes = data[right_cut & height_cut & z_cut_tr_holes]
left_sens_holes = data[left_cut & height_cut & z_cut_tr_holes]
top_sens_holes = data[top_cut & width_cut & z_cut_tr_holes]
bot_sens_holes = data[bot_cut & width_cut & z_cut_tr_holes]

right_top_sens_holes = data[top_cut & right_cut & z_cut_tr_holes]
left_top_sens_holes = data[top_cut & left_cut & z_cut_tr_holes]
right_bot_sens_holes = data[bot_cut & right_cut & z_cut_tr_holes]
left_bot_sens_holes = data[bot_cut & left_cut & z_cut_tr_holes]

l_limit = -0.12
r_limit = 0.12
bins = 1000

# plt.hist2d(right_sens_holes.photon_last_x, right_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_sens_holes.photon_last_x, left_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens_holes.photon_last_x, top_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens_holes.photon_last_x, bot_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())

# plt.hist2d(right_top_sens_holes.photon_last_x, right_top_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_top_sens_holes.photon_last_x, left_top_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(right_bot_sens_holes.photon_last_x, right_bot_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_bot_sens_holes.photon_last_x, left_bot_sens_holes.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.show()

# plt.hist2d(right_sens_holes.photon_last_x, right_sens_holes.photon_last_z, bins=bins,  norm=mpl.colors.LogNorm())
# plt.hist2d(left_sens_holes.photon_last_x, left_sens_holes.photon_last_z, bins=bins, norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens_holes.photon_last_x, top_sens_holes.photon_last_z, bins=bins, norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens_holes.photon_last_x, bot_sens_holes.photon_last_z, bins=bins, norm=mpl.colors.LogNorm())

# plt.hist2d(right_top_sens_holes.photon_last_x, right_top_sens_holes.photon_last_z, bins=bins, norm=mpl.colors.LogNorm())
# plt.hist2d(left_top_sens_holes.photon_last_x, left_top_sens_holes.photon_last_z, bins=bins, norm=mpl.colors.LogNorm())
# plt.hist2d(right_bot_sens_holes.photon_last_x, right_bot_sens_holes.photon_last_z, bins=bins, norm=mpl.colors.LogNorm())
# plt.hist2d(left_bot_sens_holes.photon_last_x, left_bot_sens_holes.photon_last_z, bins=bins, norm=mpl.colors.LogNorm())
# plt.show()


# print('Right prob: ', len(right_sens) / len(data))
# print('Left prob: ', len(left_sens) / len(data))
# print('Top prob: ', len(top_sens) / len(data))
# print('Bot prob: ', len(bot_sens) / len(data))

# print('Right top prob: ', len(right_top_sens) / len(data))
# print('Left top prob: ', len(left_top_sens) / len(data))
# print('Right bot prob: ', len(right_bot_sens) / len(data))
# print('Left bot prob: ', len(left_bot_sens) / len(data))

# Calculate all pitch excluding trenches

x_shift = y_shift = PitchWidth/2 + trench_width
height_cut = (  PitchWidth/2 - trench_width > data['photon_last_y'] )  &  ( data['photon_last_y'] > -PitchWidth/2 + trench_width )
right_cut = (data['photon_last_x'] > x_shift) & (data['photon_last_x'] < x_shift + (PitchWidth - 2 * trench_width)) 
left_cut = (data['photon_last_x'] < -x_shift) & (data['photon_last_x'] > -x_shift - (PitchWidth - 2 * trench_width))

width_cut = (  PitchWidth/2 - trench_width > data['photon_last_x'] )  &  ( data['photon_last_x'] > -PitchWidth/2 + trench_width )
top_cut = ((data['photon_last_y'] > y_shift) & (data['photon_last_y'] < y_shift + (PitchWidth - 2 * trench_width)))
bot_cut = ((data['photon_last_y'] < -y_shift) & (data['photon_last_y'] > -y_shift - (PitchWidth - 2 * trench_width)))

# Electron avalanches
right_sens_electrons = data[right_cut & height_cut & z_cut_tr_electron]
left_sens_electrons = data[left_cut & height_cut & z_cut_tr_electron]
top_sens_electrons = data[top_cut & width_cut & z_cut_tr_electron]
bot_sens_electrons = data[bot_cut & width_cut & z_cut_tr_electron]

right_top_sens_electrons = data[top_cut & right_cut & z_cut_tr_electron]
left_top_sens_electrons = data[top_cut & left_cut & z_cut_tr_electron]
right_bot_sens_electrons = data[bot_cut & right_cut & z_cut_tr_electron]
left_bot_sens_electrons = data[bot_cut & left_cut & z_cut_tr_electron]

# Holes avalanches
right_sens_holes = data[right_cut & height_cut & z_cut_tr_holes]
left_sens_holes = data[left_cut & height_cut & z_cut_tr_holes]
top_sens_holes = data[top_cut & width_cut & z_cut_tr_holes]
bot_sens_holes = data[bot_cut & width_cut & z_cut_tr_holes]

right_top_sens_holes = data[top_cut & right_cut & z_cut_tr_holes]
left_top_sens_holes = data[top_cut & left_cut & z_cut_tr_holes]
right_bot_sens_holes = data[bot_cut & right_cut & z_cut_tr_holes]
left_bot_sens_holes = data[bot_cut & left_cut & z_cut_tr_holes]

# plt.hist2d(right_sens_electrons.photon_last_x, right_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_sens_electrons.photon_last_x, left_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens_electrons.photon_last_x, top_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens_electrons.photon_last_x, bot_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())

# plt.hist2d(right_top_sens_electrons.photon_last_x, right_top_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_top_sens_electrons.photon_last_x, left_top_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(right_bot_sens_electrons.photon_last_x, right_bot_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_bot_sens_electrons.photon_last_x, left_bot_sens_electrons.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.colorbar()
# plt.xlabel('X, [mm]')
# plt.ylabel('Y, [mm]')
# plt.show()

# print(right_sens.photon_last_z)
# all_sens = np.concatenate((right_sens.photon_last_z, left_sens.photon_last_z, top_sens.photon_last_z,
# 	bot_sens.photon_last_z, right_top_sens.photon_last_z, left_top_sens.photon_last_z,
# 	right_bot_sens.photon_last_z, left_bot_sens.photon_last_z))

# plt.hist(all_sens, bins=30)
# plt.xlabel('Z, [mm]')
# plt.title('Depth distribution of absorbed photons')
# plt.show()

# plt.hist2d(left_sens.photon_last_x, left_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens.photon_last_x, top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens.photon_last_x, bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())

# plt.hist2d(right_top_sens.photon_last_x, right_top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_top_sens.photon_last_x, left_top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(right_bot_sens.photon_last_x, right_bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_bot_sens.photon_last_x, left_bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())

prob_array_electrons = []
prob_array_electrons.append(len(right_sens_electrons) / len(data))
prob_array_electrons.append(len(left_sens_electrons) / len(data))
prob_array_electrons.append(len(top_sens_electrons) / len(data))
prob_array_electrons.append(len(bot_sens_electrons) / len(data))
prob_array_electrons.append(len(right_top_sens_electrons) / len(data))
prob_array_electrons.append(len(left_top_sens_electrons) / len(data))
prob_array_electrons.append(len(right_bot_sens_electrons) / len(data))
prob_array_electrons.append(len(left_bot_sens_electrons) / len(data))

prob_array_holes = []
prob_array_holes.append(len(right_sens_holes) / len(data))
prob_array_holes.append(len(left_sens_holes) / len(data))
prob_array_holes.append(len(top_sens_holes) / len(data))
prob_array_holes.append(len(bot_sens_holes) / len(data))
prob_array_holes.append(len(right_top_sens_holes) / len(data))
prob_array_holes.append(len(left_top_sens_holes) / len(data))
prob_array_holes.append(len(right_bot_sens_holes) / len(data))
prob_array_holes.append(len(left_bot_sens_holes) / len(data))

print('Everything excluding trenches X and Y with thickness cut. Electrons avalanches')
print('Right prob: ', len(right_sens_electrons) / len(data))
print('Left prob: ', len(left_sens_electrons) / len(data))
print('Top prob: ', len(top_sens_electrons) / len(data))
print('Bot prob: ', len(bot_sens_electrons) / len(data))

print('Right top prob: ', len(right_top_sens_electrons) / len(data))
print('Left top prob: ', len(left_top_sens_electrons) / len(data))
print('Right bot prob: ', len(right_bot_sens_electrons) / len(data))
print('Left bot prob: ', len(left_bot_sens_electrons) / len(data))

print('Everything excluding trenches X and Y with thickness cut. Holes avalanches')
print('Right prob: ', len(right_sens_holes) / len(data))
print('Left prob: ', len(left_sens_holes) / len(data))
print('Top prob: ', len(top_sens_holes) / len(data))
print('Bot prob: ', len(bot_sens_holes) / len(data))

print('Right top prob: ', len(right_top_sens_holes) / len(data))
print('Left top prob: ', len(left_top_sens_holes) / len(data))
print('Right bot prob: ', len(right_bot_sens_holes) / len(data))
print('Left bot prob: ', len(left_bot_sens_holes) / len(data))

P_1_electron = sum(prob_array_electrons)
P_1_holes = sum(prob_array_holes)
# print('Photon in 1 closest cell (electron): ', sum(prob_array_electrons) * 100, '%')
# print('Photon in 1 closest cell (holes): ', P_1_holes * 100, '%')

P_2_electron  = 0
for i in range(8):
	for j in range(i):
		P_2_electron += prob_array_electrons[i] * prob_array_electrons[j]
# print('Photon in 2 closest cells (electron) ', P_2_electron * 100 * 2, '%') # Delete 2 factor?

P_2_hole  = 0
for i in range(8):
	for j in range(i):
		P_2_hole += prob_array_holes[i] * prob_array_holes[j]
# print('Photon in 2 closest cells (electron) ', P_2_hole * 100 * 2, '%') # Delete 2 factor?

# Calculate probability of crosstalk

voltages = [2, 3, 4, 5, 6, 7]

# Probabilities for electron and holes generating avalanche
electron_prob_FBK = np.array([0.75, 0.87, 0.93, 0.96, 0.98, 0.99])
holes_prob_FBK = np.array([0.07, 0.1, 0.12, 0.15, 0.175, 0.19])

electron_prob_HPK = np.array([0.7, 0.82, 0.89, 0.93, 0.96, 0.98])
holes_prob_HPK = np.array([0.07, 0.11, 0.14, 0.17, 0.2, 0.24])

# Number of photons
gain_FBK = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7]) * 1e-12 / 1.602 / 1e-19 # Read from plot Giacomo Gallina. Development of a single vacuum ultra-violet photon-sensing solution for nEXO. PhD thesis, University of British Columbia, 2021.
gain_HPK = np.array([1.6, 2.2, 3, 3.6, 4.4, 5.1]) * 1e6 # https://www.osti.gov/servlets/purl/1557100
numer_of_photons_per_carrier_FBK = 4.04e-6 # photon/e
numer_of_photons_per_carrier_HPK = 8.71e-6 # photon/e
photon_num_FBK = gain_FBK * numer_of_photons_per_carrier_FBK
photon_num_HPK = gain_HPK * numer_of_photons_per_carrier_HPK

total_prob_FBK = photon_num_FBK * (electron_prob_FBK * P_1_electron + holes_prob_FBK * P_1_holes)
total_prob_HPK = photon_num_HPK * (electron_prob_HPK * P_1_electron + holes_prob_HPK * P_1_holes)

# print(total_prob_FBK * 100)
print(total_prob_FBK[-1], ',', sep='')
plt.errorbar(voltages, total_prob_HPK * 100, yerr=0.02e-6 * gain_FBK, linestyle='None', marker='o')
# plt.errorbar(voltages, photon_num_HPK * 0.00002 * 100, yerr=0.02e-6 * gain_HPK, linestyle='None', marker='o')
plt.xlabel('Voltage [V]')
plt.ylabel('Probability of crosstalk [%]')
plt.show()











