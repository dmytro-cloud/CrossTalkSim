#!/usr/bin/python

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd

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
			if x_i == 0 and y_j == 0:
				continue
			left_cut = PitchWidth/2 + trench_width/2 + (x_i-1) * PitchWidth
			right_cut = PitchWidth/2 + x_i * PitchWidth - trench_width / 2
			bottom_cut = PitchWidth/2 + trench_width/2 + (y_j - 1) * PitchWidth
			top_cut = PitchWidth/2 + y_j * PitchWidth - trench_width / 2
			z_bottom_cut = -si_thickness / 2 + d_star
			z_top_cut =  -si_thickness / 2 + d_star + electron_width

			if holes_av_region:
				z_top_cut = -si_thickness / 2 + d_star + electron_width + holes_width
				z_bottom_cut = -si_thickness / 2 + d_star + electron_width

			av_data[str(x_i) + ':' + str(y_j)] = Apply_cut(data, (left_cut, right_cut), (bottom_cut, top_cut), (z_bottom_cut, z_top_cut))
			av_data[str(-x_i) + ':' + str(-y_j)] = Apply_cut(data, (-right_cut, -left_cut), (-top_cut, -bottom_cut), (z_bottom_cut, z_top_cut))
			av_data[str(-x_i) + ':' + str(y_j)] = Apply_cut(data, (-right_cut, -left_cut), (bottom_cut, top_cut), (z_bottom_cut, z_top_cut))
			av_data[str(x_i) + ':' + str(-y_j)] = Apply_cut(data, (left_cut, right_cut), (-top_cut, -bottom_cut), (z_bottom_cut, z_top_cut))
	return av_data

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

P_1_electron = len(res_electrons) / len(data)
P_1_holes = len(res_holes) / len(data)
# print('Photon in 1 closest cell (electron): ', sum(prob_array_electrons) * 100, '%')
# print('Photon in 1 closest cell (holes): ', P_1_holes * 100, '%')

# P_2_electron  = 0
# for i in range(8):
# 	for j in range(i):
# 		P_2_electron += prob_array_electrons[i] * prob_array_electrons[j]
# # print('Photon in 2 closest cells (electron) ', P_2_electron * 100 * 2, '%') # Delete 2 factor?

# P_2_hole  = 0
# for i in range(8):
# 	for j in range(i):
# 		P_2_hole += prob_array_holes[i] * prob_array_holes[j]
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

electron_total_prob_FBK = photon_num_FBK * electron_prob_FBK * P_1_electron
holes_total_prob_FBK = photon_num_FBK * holes_prob_FBK * P_1_holes
total_prob_FBK = electron_total_prob_FBK + holes_total_prob_FBK
total_prob_HPK = photon_num_HPK * (electron_prob_HPK * P_1_electron + holes_prob_HPK * P_1_holes)

# print(total_prob_FBK * 100)
print(total_prob_FBK[-1], ',', sep='')
plt.errorbar(voltages, total_prob_FBK * 100, yerr=0.02e-6 * gain_FBK, linestyle='None', marker='o', label='Full Optical CS probability')
plt.errorbar(voltages, electron_total_prob_FBK * 100, yerr=0.02e-6 * gain_FBK, linestyle='None', marker='x', label='Electron driven CS probability')
# plt.errorbar(voltages, photon_num_HPK * 0.00002 * 100, yerr=0.02e-6 * gain_HPK, linestyle='None', marker='o')
plt.legend()
plt.xlabel('Voltage [V]')
plt.ylabel('Probability of crosstalk [%]')
plt.show()











