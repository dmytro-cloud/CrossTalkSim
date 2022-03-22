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

z_cut_si = (data['photon_last_z'] < si_thickness / 2) & (data['photon_last_z'] > -si_thickness / 2) 

x_shift = y_shift = PitchWidth - AvWidth / 2
height_cut = (  AvWidth/2 > data['photon_last_y'] )  &  ( data['photon_last_y'] > -AvWidth/2 )
height_cut = (  AvWidth/2 > data['photon_last_y'] )  &  ( data['photon_last_y'] > -AvWidth/2 ) 
right_cut = (data['photon_last_x'] > x_shift) & (data['photon_last_x'] < x_shift + AvWidth)
left_cut = (data['photon_last_x'] < -x_shift) & (data['photon_last_x'] > -x_shift - AvWidth)

width_cut = (  AvWidth/2 > data['photon_last_x'] )  &  ( data['photon_last_x'] > -AvWidth/2 )
top_cut = ((data['photon_last_y'] > y_shift) & (data['photon_last_y'] < y_shift + AvWidth))
bot_cut = ((data['photon_last_y'] < -y_shift) & (data['photon_last_y'] > -y_shift - AvWidth))

right_sens = data[right_cut & height_cut & z_cut_si]
left_sens = data[left_cut & height_cut & z_cut_si]
top_sens = data[top_cut & width_cut & z_cut_si]
bot_sens = data[bot_cut & width_cut & z_cut_si]

right_top_sens = data[top_cut & right_cut & z_cut_si]
left_top_sens = data[top_cut & left_cut & z_cut_si]
right_bot_sens = data[bot_cut & right_cut & z_cut_si]
left_bot_sens = data[bot_cut & left_cut & z_cut_si]

l_limit = -0.12
r_limit = 0.12

# plt.hist2d(right_sens.photon_last_x, right_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_sens.photon_last_x, left_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens.photon_last_x, top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens.photon_last_x, bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())

# plt.hist2d(right_top_sens.photon_last_x, right_top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_top_sens.photon_last_x, left_top_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(right_bot_sens.photon_last_x, right_bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_bot_sens.photon_last_x, left_bot_sens.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
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

right_sens = data[right_cut & height_cut & z_cut_si]
left_sens = data[left_cut & height_cut & z_cut_si]
top_sens = data[top_cut & width_cut & z_cut_si]
bot_sens = data[bot_cut & width_cut & z_cut_si]

right_top_sens = data[top_cut & right_cut & z_cut_si]
left_top_sens = data[top_cut & left_cut & z_cut_si]
right_bot_sens = data[bot_cut & right_cut & z_cut_si]
left_bot_sens = data[bot_cut & left_cut & z_cut_si]

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
prob_array = []
prob_array.append(len(right_sens) / len(data))
prob_array.append(len(left_sens) / len(data))
prob_array.append(len(top_sens) / len(data))
prob_array.append(len(bot_sens) / len(data))
prob_array.append(len(right_top_sens) / len(data))
prob_array.append(len(left_top_sens) / len(data))
prob_array.append(len(right_bot_sens) / len(data))
prob_array.append(len(left_bot_sens) / len(data))

print('Everything excluding trenches X and Y')
print('Right prob: ', len(right_sens) / len(data))
print('Left prob: ', len(left_sens) / len(data))
print('Top prob: ', len(top_sens) / len(data))
print('Bot prob: ', len(bot_sens) / len(data))

print('Right top prob: ', len(right_top_sens) / len(data))
print('Left top prob: ', len(left_top_sens) / len(data))
print('Right bot prob: ', len(right_bot_sens) / len(data))
print('Left bot prob: ', len(left_bot_sens) / len(data))

print('Photon in 1 closest cell: ', sum(prob_array) * 100, '%')

P_3  = 0
for i in range(8):
	for j in range(i):
		P_3 += prob_array[i] * prob_array[j]
print('Photon in 2 closest cells ', P_3 * 100 * 2, '%') # Remove 2 factor?

# Volume under trenches
trench_x_shift = trench_y_shift = PitchWidth/2
# trench_x_limit = trench_y_limit = PitchWidth/2 + trench_width
z_cut_trench = (data['photon_last_z'] > -si_thickness / 2 + trench_depth + trench_depth/3 ) & (data['photon_last_z'] < si_thickness/2)

height_cut = (  PitchWidth/2 - trench_width > data['photon_last_y'] )  &  ( data['photon_last_y'] > -PitchWidth/2 + trench_width )
right_under_trench_cut = (data['photon_last_x'] > trench_x_shift) & (data['photon_last_x'] < trench_x_shift + trench_width) 
left_under_trench_cut = (data['photon_last_x'] < -trench_x_shift) & (data['photon_last_x'] > -trench_x_shift - trench_width) 

width_cut = (  PitchWidth/2 - trench_width > data['photon_last_x'] )  &  ( data['photon_last_x'] > -PitchWidth/2 + trench_width )
top_under_trench_cut = ((data['photon_last_y'] > trench_y_shift) & (data['photon_last_y'] < trench_y_shift + trench_width))
bot_under_trench_cut = ((data['photon_last_y'] < -trench_y_shift) & (data['photon_last_y'] > -trench_y_shift - trench_width))

right_sens_trench = data[right_under_trench_cut & height_cut & z_cut_trench]
left_sens_trench = data[left_under_trench_cut & height_cut & z_cut_trench]
top_sens_trench = data[top_under_trench_cut & width_cut & z_cut_trench]
bot_sens_trench = data[bot_under_trench_cut & width_cut & z_cut_trench]

print('Everything under trenches')
print('Right prob: ', len(right_sens_trench) / len(data))
print('Left prob: ', len(left_sens_trench) / len(data))
print('Top prob: ', len(top_sens_trench) / len(data))
print('Bot prob: ', len(bot_sens_trench) / len(data))
prob_array = []
prob_array.append( (len(right_sens) + len(right_sens_trench))/ len(data))
prob_array.append( (len(left_sens) + len(left_sens_trench)) / len(data))
prob_array.append( (len(top_sens) + len(top_sens_trench)) / len(data))
prob_array.append( (len(bot_sens) + len(bot_sens_trench)) / len(data))
prob_array.append(len(right_top_sens) / len(data))
prob_array.append(len(left_top_sens) / len(data))
prob_array.append(len(right_bot_sens) / len(data))
prob_array.append(len(left_bot_sens) / len(data))

print('Photon in 1 closest cell: ', sum(prob_array) * 100, '%')

P_3  = 0
for i in range(8):
	for j in range(i):
		P_3 += prob_array[i] * prob_array[j]
print('Photon in 2 closest cells ', P_3 * 100 * 2, '%') # delete 2 factor?

# plt.hist2d(right_sens_trench.photon_last_z, right_sens_trench.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(left_sens_trench.photon_last_z, left_sens_trench.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(top_sens_trench.photon_last_z, top_sens_trench.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.hist2d(bot_sens_trench.photon_last_z, bot_sens_trench.photon_last_y, bins=1000, range=((l_limit, r_limit), (l_limit, r_limit)), norm=mpl.colors.LogNorm())
# plt.colorbar()
# plt.show()
# data = file["NumberElectrons"].arrays(["electron_initial_x", "electron_initial_y", "electron_initial_z",
#  "electron_last_x", "electron_last_y", "electron_last_z"], library='pd')

# print(data.electron_initial_x, data.electron_initial_y)
# plt.hist2d(data.electron_initial_x, data.electron_initial_y, bins=1000, norm=mpl.colors.LogNorm())
# # plt.zscale('log')
# plt.show()