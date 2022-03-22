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

file = uproot.open(filename) #"../output/test.root")

data = file["ntp"].arrays(["photon_initial_x", "photon_initial_y", "photon_initial_z",
 "photon_last_x", "photon_last_y", "photon_last_z"], library='pd')
print(data)
# plt.hist2d(data.photon_last_x, data.photon_last_y, bins=1000, norm=mpl.colors.LogNorm())#, range=((-0.3, 0.3), (-0.3, 0.3)))
# # plt.zscale('log')
# plt.colorbar()
# plt.xlabel('X, [mm]')
# plt.ylabel('Y, [mm]')
# plt.show()

# plt.hist2d(data.photon_last_z, data.photon_last_y, bins=1000, norm=mpl.colors.LogNorm())
# plt.colorbar()
# plt.xlabel('Z, [mm]')
# plt.ylabel('Y, [mm]')
# plt.show()

# FBK 
# PitchWidth = 0.035
# AvWidth = 0.032
# trench_width = 0.001
# trench_depth = 0.003
# si_thickness = 0.3

# z_cut_si = (data['photon_last_z'] < si_thickness / 2) & (data['photon_last_z'] > -si_thickness / 2) 
z_cut_tr = (data['photon_last_z'] < -si_thickness / 2 + trench_depth) & (data['photon_last_z'] > -si_thickness / 2) 

x_shift = y_shift = PitchWidth - AvWidth / 2
height_cut = (  AvWidth/2 > data['photon_last_y'] )  &  ( data['photon_last_y'] > -AvWidth/2 )
height_cut = (  AvWidth/2 > data['photon_last_y'] )  &  ( data['photon_last_y'] > -AvWidth/2 ) 
right_cut = (data['photon_last_x'] > x_shift) & (data['photon_last_x'] < x_shift + AvWidth)
left_cut = (data['photon_last_x'] < -x_shift) & (data['photon_last_x'] > -x_shift - AvWidth)

width_cut = (  AvWidth/2 > data['photon_last_x'] )  &  ( data['photon_last_x'] > -AvWidth/2 )
top_cut = ((data['photon_last_y'] > y_shift) & (data['photon_last_y'] < y_shift + AvWidth))
bot_cut = ((data['photon_last_y'] < -y_shift) & (data['photon_last_y'] > -y_shift - AvWidth))

right_sens = data[right_cut & height_cut & z_cut_tr]
left_sens = data[left_cut & height_cut & z_cut_tr]
top_sens = data[top_cut & width_cut & z_cut_tr]
bot_sens = data[bot_cut & width_cut & z_cut_tr]

right_top_sens = data[top_cut & right_cut & z_cut_tr]
left_top_sens = data[top_cut & left_cut & z_cut_tr]
right_bot_sens = data[bot_cut & right_cut & z_cut_tr]
left_bot_sens = data[bot_cut & left_cut & z_cut_tr]

l_limit = -0.12
r_limit = 0.12
bins = 1000

# plt.hist2d(right_sens.photon_last_x, right_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_sens.photon_last_x, left_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens.photon_last_x, top_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens.photon_last_x, bot_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())

# plt.hist2d(right_top_sens.photon_last_x, right_top_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_top_sens.photon_last_x, left_top_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(right_bot_sens.photon_last_x, right_bot_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_bot_sens.photon_last_x, left_bot_sens.photon_last_y, bins=bins, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.show()

print('Right prob: ', len(right_sens) / len(data))
print('Left prob: ', len(left_sens) / len(data))
print('Top prob: ', len(top_sens) / len(data))
print('Bot prob: ', len(bot_sens) / len(data))

print('Right top prob: ', len(right_top_sens) / len(data))
print('Left top prob: ', len(left_top_sens) / len(data))
print('Right bot prob: ', len(right_bot_sens) / len(data))
print('Left bot prob: ', len(left_bot_sens) / len(data))

# Calculate all pitch excluding trenches

x_shift = y_shift = PitchWidth/2 + trench_width
height_cut = (  PitchWidth/2 - trench_width > data['photon_last_y'] )  &  ( data['photon_last_y'] > -PitchWidth/2 + trench_width )
right_cut = (data['photon_last_x'] > x_shift) & (data['photon_last_x'] < x_shift + (PitchWidth - 2 * trench_width)) 
left_cut = (data['photon_last_x'] < -x_shift) & (data['photon_last_x'] > -x_shift - (PitchWidth - 2 * trench_width))

width_cut = (  PitchWidth/2 - trench_width > data['photon_last_x'] )  &  ( data['photon_last_x'] > -PitchWidth/2 + trench_width )
top_cut = ((data['photon_last_y'] > y_shift) & (data['photon_last_y'] < y_shift + (PitchWidth - 2 * trench_width)))
bot_cut = ((data['photon_last_y'] < -y_shift) & (data['photon_last_y'] > -y_shift - (PitchWidth - 2 * trench_width)))

right_sens = data[right_cut & height_cut & z_cut_tr]
left_sens = data[left_cut & height_cut & z_cut_tr]
top_sens = data[top_cut & width_cut & z_cut_tr]
bot_sens = data[bot_cut & width_cut & z_cut_tr]

right_top_sens = data[top_cut & right_cut & z_cut_tr]
left_top_sens = data[top_cut & left_cut & z_cut_tr]
right_bot_sens = data[bot_cut & right_cut & z_cut_tr]
left_bot_sens = data[bot_cut & left_cut & z_cut_tr]

# plt.hist2d(right_sens.photon_last_x, right_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_sens.photon_last_x, left_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens.photon_last_x, top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens.photon_last_x, bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())

# plt.hist2d(right_top_sens.photon_last_x, right_top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_top_sens.photon_last_x, left_top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(right_bot_sens.photon_last_x, right_bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_bot_sens.photon_last_x, left_bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
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
# plt.hist2d(left

prob_array = []
prob_array.append(len(right_sens) / len(data))
prob_array.append(len(left_sens) / len(data))
prob_array.append(len(top_sens) / len(data))
prob_array.append(len(bot_sens) / len(data))
prob_array.append(len(right_top_sens) / len(data))
prob_array.append(len(left_top_sens) / len(data))
prob_array.append(len(right_bot_sens) / len(data))
prob_array.append(len(left_bot_sens) / len(data))

print('Everything excluding trenches X and Y with thickness cut')
print('Right prob: ', len(right_sens) / len(data))
print('Left prob: ', len(left_sens) / len(data))
print('Top prob: ', len(top_sens) / len(data))
print('Bot prob: ', len(bot_sens) / len(data))

print('Right top prob: ', len(right_top_sens) / len(data))
print('Left top prob: ', len(left_top_sens) / len(data))
print('Right bot prob: ', len(right_bot_sens) / len(data))
print('Left bot prob: ', len(left_bot_sens) / len(data))

P_1 = sum(prob_array)
print('Photon in 1 closest cell: ', sum(prob_array) * 100, '%')

P_3  = 0
for i in range(8):
	for j in range(i):
		P_3 += prob_array[i] * prob_array[j]
print('Photon in 2 closest cells ', P_3 * 100 * 2, '%') # Delete 2 factor?

# Calculate probability of crosstalk

voltages = [2, 3, 4, 5, 6, 7]
gain_FBK = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7]) * 1e-12 / 1.602 / 1e-19 # Read from plot Giacomo Gallina. Development of a single vacuum ultra-violet photon-sensing solution for nEXO. PhD thesis, University of British Columbia, 2021.
gain_HPK = np.array([1.6, 2.2, 3, 3.6, 4.4, 5.1]) * 1e6 # https://www.osti.gov/servlets/purl/1557100
numer_of_photons_per_carrier_FBK = 4.04e-6 # photon/e
numer_of_photons_per_carrier_HPK = 8.71e-6 # photon/e
photon_num_FBK = gain_FBK * numer_of_photons_per_carrier_FBK
photon_num_HPK = gain_HPK * numer_of_photons_per_carrier_HPK
print(photon_num_FBK * P_1 * 100)
plt.errorbar(voltages, photon_num_FBK * P_1 * 100, yerr=0.02e-6 * gain_FBK, linestyle='None', marker='o')
# plt.errorbar(voltages, photon_num_HPK * 0.00002 * 100, yerr=0.02e-6 * gain_HPK, linestyle='None', marker='o')
plt.xlabel('Voltage [V]')
plt.ylabel('Probability of crosstalk [%]')
plt.show()












