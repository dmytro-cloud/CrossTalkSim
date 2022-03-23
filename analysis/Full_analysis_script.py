#!/usr/bin/python

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import json
import time

start = time.time()

# filename = sys.argv[1] # File that you will obtaine probabilities of the crosstalk from and put as initial light (100% prob for appearing)
# filename_params = sys.argv[2]

# Load geometry parameters
filename_params = '../Data/Geometry/Test_data.json' #sys.argv[2]
file_params = open(filename_params)
sipm_params = json.load(file_params)
file_params.close()

# Load electrons map
filename_electrons = ""
if sipm_params["number_of_electrons_generated"] > 1000:
    filename_electrons = "../Data/Geant4_electrons_input/" + str( int(sipm_params["number_of_electrons_generated"] / 1000)) + "k_points.json";
else:
    filename_electrons = "../Data/Geant4_electrons_input/" + str( int(sipm_params["number_of_electrons_generated"] / 1000)) + "_points.json";

file_electrons = open(filename_electrons)
electrons_map = json.load(file_electrons)
file_electrons.close()

# To mkm from mm
for key, value in sipm_params.items():
    if key == '_typename' or key == 'n_bins' or key == 'n_cells':
        continue
    sipm_params[key] *= 1e-3
print(sipm_params)

PitchWidth = sipm_params['pitch_width']
AvWidth = sipm_params['av_width']
trench_width = sipm_params['trench_width']
trench_depth = sipm_params['trench_depth']
si_thickness = sipm_params['si_thickness']
number_of_bins = sipm_params['n_bins']
number_of_cells = sipm_params['n_cells']
d_star = sipm_params['d_star'] # params[5]
electron_width = sipm_params['electron_width'] #0.001 - HPK param # params[6]
w_star = sipm_params['w_star'] # 0.004 - HPK # params[7]
holes_width = w_star - electron_width

# Functions for CS prob 
def Apply_cut(data, x_boundaries, y_boundaries, z_boundaries):
    cut_x = ((data['photon_last_x'] < x_boundaries[1]) & (data['photon_last_x'] > x_boundaries[0]))
    cut_y = ((data['photon_last_y'] < y_boundaries[1]) & (data['photon_last_y'] > y_boundaries[0]))
    cut_z = ((data['photon_last_z'] < z_boundaries[1]) & (data['photon_last_z'] > z_boundaries[0]))
    return data[cut_x & cut_y & cut_z]

def Get_data_little_cuts(data, x_cut, y_cut, z_cut, binning=3):
    cell_data = {}
    for ii in range(binning):
        for jj in range(binning):
            small_l_cut = x_cut[0] + ii * PitchWidth / binning
            small_r_cut = x_cut[0] + (ii + 1) * PitchWidth / binning
            small_b_cut = y_cut[0] + jj * PitchWidth / binning
            small_t_cut = y_cut[0] + (jj + 1) * PitchWidth / binning
            cell_data[str(ii) + ':' + str(jj)] = Apply_cut(data, (small_l_cut, small_r_cut), (small_b_cut, small_t_cut), z_cut)
    return cell_data

# Sounds super silly but it works only for odd numbers so far
def Get_data_after_cuts(data, holes_av_region, cells=5, binning=3):
    cells_region = cells // 2 + 1
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

            av_data[str(x_i) + ':' + str(y_j)] = Get_data_little_cuts(data, (left_cut, right_cut), (bottom_cut, top_cut), (z_bottom_cut, z_top_cut), binning=binning)
            av_data[str(-x_i) + ':' + str(-y_j)] = Get_data_little_cuts(data, (-right_cut, -left_cut), (-top_cut, -bottom_cut), (z_bottom_cut, z_top_cut), binning=binning)
            av_data[str(-x_i) + ':' + str(y_j)] = Get_data_little_cuts(data, (-right_cut, -left_cut), (bottom_cut, top_cut), (z_bottom_cut, z_top_cut), binning=binning)
            av_data[str(x_i) + ':' + str(-y_j)] = Get_data_little_cuts(data, (left_cut, right_cut), (-top_cut, -bottom_cut), (z_bottom_cut, z_top_cut), binning=binning)
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
    tmp.photon_last_x += x_shift * PitchWidth 
    tmp.second_last_y += y_shift * PitchWidth
    tmp.photon_last_y += y_shift * PitchWidth

    # To avoid division by 0
    if probability == 0:
        return tmp[:0]

    # Might be better to switch to the actual number to not get confused if the files have different number of events
    return tmp[:int(len(tmp) // (1. / probability))]

# file = uproot.open(filename) #"../output/test.root")

# data = file["ntp"].arrays(["photon_initial_x", "photon_initial_y", "photon_initial_z",
#  "photon_last_x", "photon_last_y", "photon_last_z"], library='pd')

# electrons_av_data = Get_data_after_cuts(data, False, cells=number_of_cells, binning=number_of_bins)
# # print(electrons_av_data['0:0'])
# holes_av_data = Get_data_after_cuts(data, True, cells=number_of_cells, binning=number_of_bins)

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

# Load binned uniform data
full_av_data_prob = {}
for cell_i in range(-(number_of_cells//2), number_of_cells//2 + 1):
    for cell_j in range(-(number_of_cells//2), number_of_cells//2 + 1):
        for bin_i in range(number_of_bins):
            tmp = []
            for bin_j in range(number_of_bins):
                electrons_map_key_bin = str(bin_i) + ':' + str(bin_j)
                electrons_map_key_cell = str(cell_i) + ':' + str(cell_j)
                # Do not load and calculate everything if probability to appear close to 0
                if electrons_map[electrons_map_key_cell][electrons_map_key_bin] < 1e-9:
                    continue

                file = uproot.open('../Data/Geant4_output/FBK_8bins/FBK_binned_' + str(bin_i) + '_' + str(bin_j) + '.root') #"../output/test.root")
                data = file["ntp"].arrays(["photon_initial_x", "photon_initial_y", "photon_initial_z", "second_last_x", "second_last_y", "second_last_z",
                 "photon_last_x", "photon_last_y", "photon_last_z"], library='pd')

                tmp_electrons_av_data = Get_data_after_cuts(Shift_xy(data, cell_i, cell_j), False, cells=number_of_cells, binning=number_of_bins)
                tmp_holes_av_data = Get_data_after_cuts(Shift_xy(data, cell_i, cell_j), True, cells=number_of_cells, binning=number_of_bins)

                print(electrons_map_key_cell, electrons_map_key_bin)
                # Calculate probabilities
                for key, item in tmp_electrons_av_data.items():
                    if key not in full_av_data_prob:
                        full_av_data_prob[key] = {}
                    for key2, _ in item.items():
                        # low CS rate for initial cell
                        # With the shift need to adjust this coordinate to the initial cell of the electron
                        if key == electrons_map_key_cell:
                            if key2 not in full_av_data_prob[key]:
                                full_av_data_prob[key][key2] = electrons_map[electrons_map_key_cell][electrons_map_key_bin] * 0.001 * (1 - np.power(1 - (len(tmp_electrons_av_data[key][key2]) * 0.99 + len(tmp_holes_av_data[key][key2]) * 0.19) / len(data), photon_num_FBK[-1]))
                            else:
                                full_av_data_prob[key][key2] += electrons_map[electrons_map_key_cell][electrons_map_key_bin] * 0.001 * (1 - np.power(1 - (len(tmp_electrons_av_data[key][key2]) * 0.99 + len(tmp_holes_av_data[key][key2]) * 0.19) / len(data), photon_num_FBK[-1]))
                        else:
                            if key2 not in full_av_data_prob[key]:
                                full_av_data_prob[key][key2] = electrons_map[electrons_map_key_cell][electrons_map_key_bin] * ( 1 - np.power(1 - (len(tmp_electrons_av_data[key][key2]) * 0.99 + len(tmp_holes_av_data[key][key2]) * 0.19) / len(data), photon_num_FBK[-1]))
                            else:
                                full_av_data_prob[key][key2] += electrons_map[electrons_map_key_cell][electrons_map_key_bin] * ( 1 - np.power(1 - (len(tmp_electrons_av_data[key][key2]) * 0.99 + len(tmp_holes_av_data[key][key2]) * 0.19) / len(data), photon_num_FBK[-1]))


CS_probability = 0
for i in range(-(number_of_cells//2), number_of_cells//2 + 1):
    for j in range(-(number_of_cells//2), number_of_cells//2 + 1):
        for ii in range(number_of_bins):
            for jj in range(number_of_bins):
                CS_probability += full_av_data_prob[str(i) + ':' + str(j)][str(ii) + ':' + str(jj)]

weights = []
x = []
y = []
for ii in range(8):
    for jj in range(8):
        weights.append(full_av_data_prob['0:0'][str(ii) + ':' + str(jj)])
        x.append(ii)
        y.append(jj)

plt.hist2d(x, y, weights=weights, bins=8)
plt.savefig('CS_probability.png')
# plt.show()

# print(full_av_data_prob)
print("CS probability: ", CS_probability)
end = time.time()
print('Time to calculate CS prob: ', end - start)
#### IMAGE CALCULATIONS ### 
# As an example for 7V overvoltage

# Convolution with nearby cells probability
res = []
for i in range(-(number_of_cells//2), number_of_cells//2 + 1):
    for j in range(-(number_of_cells//2), number_of_cells//2 + 1):
        for ii in range(number_of_bins):
            for jj in range(number_of_bins):
                if electrons_map[str(i) + ':' + str(j)][str(ii) + ':' + str(jj)] < 1e-9:
                    res.append(pd.DataFrame())
                    continue
                data_tmp = Load_and_cut_data('../Data/Geant4_output/FBK_8bins/FBK_binned_' + str(ii) + '_' + str(jj) + '.root')
                # Taking into account electron data
                res.append(Shift_xy(data_tmp, i, j, electrons_map[str(i) + ':' + str(j)][str(ii) + ':' + str(jj)]))

# Load binned uniform data
uniform_binned_data = []
for i in range(number_of_bins):
    tmp = []
    for j in range(number_of_bins):
        tmp.append(Load_and_cut_data('../Data/Geant4_output/FBK_8bins/FBK_binned_' + str(i) + '_' + str(j) + '.root'))
    uniform_binned_data.append(tmp)

#I have to generate data for every of the n cells and use it.
for i in range(-(number_of_cells//2), number_of_cells//2 + 1):
    for j in range(-(number_of_cells//2), number_of_cells//2 + 1):
        for ii in range(number_of_bins):
            for jj in range(number_of_bins):
                res.append(Shift_xy(uniform_binned_data[ii][jj], i, j, full_av_data_prob[str(i) + ':' + str(j)][str(ii) + ':' + str(jj)]))

print(full_av_data_prob['0:0'])

# Self cross-talk in the cell (?)
# res.append(Shift_xy(data_cut_uniform, 0, 0, 0.01))

res = pd.concat(res)
plt.hist2d(res.second_last_x, res.second_last_y, norm=mpl.colors.LogNorm(), bins=200, range=((-0.15, 0.15), (-0.15, 0.15)))
plt.colorbar()
plt.title('The whole board')
plt.xlabel('X, [mm]')
plt.ylabel('Y, [mm]')
plt.savefig('100k_point_8bins.png')
plt.show()





