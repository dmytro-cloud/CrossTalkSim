#!/usr/bin/python

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd

filename = sys.argv[1]
filename_params = sys.argv[2]
filename_uniform = sys.argv[3]

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

# Functions for CS prob 
def Apply_cut(data, x_boundaries, y_boundaries, z_boundaries):
	cut_x = ((data['photon_last_x'] < x_boundaries[1]) & (data['photon_last_x'] > x_boundaries[0]))
	cut_y = ((data['photon_last_y'] < y_boundaries[1]) & (data['photon_last_y'] > y_boundaries[0]))
	cut_z = ((data['photon_last_z'] < z_boundaries[1]) & (data['photon_last_z'] > z_boundaries[0]))
	return data[cut_x & cut_y & cut_z]

def Get_data_after_cuts(data, holes_av_region, cells_region=3):
	av_data = {}
	# Av region count
	for x_i in range(cells_region):
		for y_j in range(cells_region):				
			left_cut = PitchWidth/2 + trench_width/2 + (x_i-1) * PitchWidth
			right_cut = PitchWidth/2 + x_i * PitchWidth - trench_width / 2
			bottom_cut = PitchWidth/2 + trench_width/2 + (y_j - 1) * PitchWidth
			top_cut = PitchWidth/2 + y_j * PitchWidth - trench_width / 2
			z_bottom_cut = -si_thickness / 2 + d_star
			z_top_cut =  -si_thickness / 2 + d_star + electron_width

			if holes_av_region:
				z_top_cut = -si_thickness / 2 + d_star + electron_width + holes_width
				z_bottom_cut = -si_thickness / 2 + d_star + electron_width

			if x_i == 0 and y_j == 0:
				left_cut = - PitchWidth/2
				bottom_cut = - PitchWidth/2
				right_cut = top_cut = PitchWidth/2

			av_data[str(x_i) + ':' + str(y_j)] = Apply_cut(data, (left_cut, right_cut), (bottom_cut, top_cut), (z_bottom_cut, z_top_cut))
			av_data[str(-x_i) + ':' + str(-y_j)] = Apply_cut(data, (-right_cut, -left_cut), (-top_cut, -bottom_cut), (z_bottom_cut, z_top_cut))
			av_data[str(-x_i) + ':' + str(y_j)] = Apply_cut(data, (-right_cut, -left_cut), (bottom_cut, top_cut), (z_bottom_cut, z_top_cut))
			av_data[str(x_i) + ':' + str(-y_j)] = Apply_cut(data, (left_cut, right_cut), (-top_cut, -bottom_cut), (z_bottom_cut, z_top_cut))
	return av_data

# def Get_small_cell_data(data, )

# Functions for mapping
def Load_and_cut_data(filename, cos_theta = np.sqrt(2) / 2):

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

	# # Assume 45 degrees apperture
	# plt.hist(data_cut.cos_theta, bins=100)
	# plt.show()

	data_cut = data_cut[data_cut['cos_theta'] > cos_theta]
	return data_cut

# Shifts the initial data. Make probability cut if needed
def Shift_xy(data, x_shift=0, y_shift=0, probability=1.):
	tmp = data.copy()
	tmp.second_last_x += x_shift * PitchWidth 
	tmp.second_last_y += y_shift * PitchWidth

	# Might be better to switch to the actual number to not get confused if the files have different number of events
	return tmp[:int(len(tmp) // (1. / probability))]

electrons_av_data = Get_data_after_cuts(data, False)
holes_av_data = Get_data_after_cuts(data, True)

# Plot to actually test it ( by eye )
res_electrons = []
for key, d in electrons_av_data.items():
	res_electrons.append(d)

res_holes = []
for key, d in holes_av_data.items():
	res_holes.append(d)

# Ttal data for holes and electrons separately
res_electrons = pd.concat(res_electrons)
res_holes = pd.concat(res_holes)
# plt.hist2d(res_electrons.photon_last_x, res_electrons.photon_last_y, bins=1000, norm=mpl.colors.LogNorm(), range=((-0.15, 0.15), (-0.15, 0.15)))
# plt.show()
print('Right top:', len(electrons_av_data['1:1']) / len(data))
print('Right:', len(electrons_av_data['1:0']) / len(data))
print('Left:', len(electrons_av_data['-1:0']) / len(data))

# Excluding the initial cell
P_1_electron = len(res_electrons) / len(data) - len(electrons_av_data['0:0']) / len(data)
P_1_holes = (len(res_holes) - len(holes_av_data['0:0']))/len(data)

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

electron_total_prob_FBK = 1 - np.power(1 - electron_prob_FBK * P_1_electron, photon_num_FBK)
holes_total_prob_FBK = 1 - np.power(1 - holes_prob_FBK * P_1_holes, photon_num_FBK)
total_prob_FBK = electron_total_prob_FBK + holes_total_prob_FBK
total_prob_FBK_correct = 1 - np.power(1 - electron_prob_FBK * P_1_electron - holes_prob_FBK * P_1_holes, photon_num_FBK)
total_prob_HPK = photon_num_HPK * (electron_prob_HPK * P_1_electron + holes_prob_HPK * P_1_holes)

# print(total_prob_FBK * 100)
print(total_prob_FBK[-1], ',', sep='')
# plt.errorbar(voltages, total_prob_FBK * 100, yerr=0.02e-6 * gain_FBK, linestyle='None', marker='o', label='Full Optical CS probability')
plt.errorbar(voltages, total_prob_FBK_correct * 100, yerr=0.02e-6 * gain_FBK, linestyle='None', marker='o', label='Full Optical CS probability')
plt.errorbar(voltages, electron_total_prob_FBK * 100, yerr=0.02e-6 * gain_FBK, linestyle='None', marker='x', label='Electron driven CS probability')
# plt.errorbar(voltages, photon_num_HPK * 0.00002 * 100, yerr=0.02e-6 * gain_HPK, linestyle='None', marker='o')
plt.legend()
plt.xlabel('Voltage [V]')
plt.ylabel('Probability of crosstalk [%]')
plt.show()



#### IMAGE CALCULATIONS ###	
print(photon_num_FBK[-1])
# As an example for 7V overvoltage
full_av_data_prob = {}
for key, _ in electrons_av_data.items():
	full_av_data_prob[key] = 1 - np.power(1 - (len(electrons_av_data[key]) * 0.99 + len(holes_av_data[key]) * 0.19) / len(data), photon_num_FBK[-1])

data_cut = Load_and_cut_data(filename)
data_cut_uniform = Load_and_cut_data(filename_uniform)

# Use not only electrons AV data but also holes 
print( full_av_data_prob['0:' + str(1)]) 
# Convolution with nearby cells probability
res = [data_cut]

for i in [-2, -1, 0, 1, 2]:
	for j in [-2, -1, 0, 1, 2]:
		res.append(Shift_xy(data_cut_uniform, i, j, full_av_data_prob[str(i) + ':' + str(j)]))

print(full_av_data_prob['0:0']  )

# Self cross-talk in the cell (?)
# res.append(Shift_xy(data_cut_uniform, 0, 0, 0.01))

res = pd.concat(res)
plt.hist2d(res.second_last_x, res.second_last_y, norm=mpl.colors.LogNorm(), bins=75, range=((-0.15, 0.15), (-0.15, 0.15)))
plt.colorbar()
plt.title('The whole board')
plt.xlabel('X, [mm]')
plt.ylabel('Y, [mm]')
plt.show()








