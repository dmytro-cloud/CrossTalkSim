#!/usr/bin/python

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import json

# Functions for prob comoputing
def Apply_cut(data, x_boundaries, y_boundaries, z_boundaries = (- np.inf, np.inf)):  # We dont use z cut for the electron.
    cut_x = ((data['finalX'] < x_boundaries[1]) & (data['finalX'] > x_boundaries[0]))
    cut_y = ((data['finalY'] < y_boundaries[1]) & (data['finalY'] > y_boundaries[0]))
    cut_z = ((data['finalZ'] < z_boundaries[1]) & (data['finalZ'] > z_boundaries[0]))
    return data[cut_x & cut_y]

# Binned data
def Get_data_little_cuts(data, x_cut, y_cut, z_cut = (- np.inf, np.inf), binning=3):
    cell_data = {}
    for ii in range(binning):
        for jj in range(binning):
            small_l_cut = x_cut[0] + ii * pitch_width / binning
            small_r_cut = x_cut[0] + (ii + 1) * pitch_width / binning
            small_b_cut = y_cut[0] + jj * pitch_width / binning
            small_t_cut = y_cut[0] + (jj + 1) * pitch_width / binning
            cell_data[str(ii) + ':' + str(jj)] = Apply_cut(data, (small_l_cut, small_r_cut), (small_b_cut, small_t_cut), z_cut)
    return cell_data

def Get_data_after_cuts(data, cells=5, nbins=3):
    cells_region = cells // 2 + 1
    av_data = {}
    # Av region count
    for x_i in range(cells_region):
        for y_j in range(cells_region):             
            left_cut = pitch_width/2 + trench_width/2 + (x_i-1) * pitch_width
            right_cut = pitch_width/2 + x_i * pitch_width - trench_width / 2
            bottom_cut = pitch_width/2 + trench_width/2 + (y_j - 1) * pitch_width
            top_cut = pitch_width/2 + y_j * pitch_width - trench_width / 2
            # z_bottom_cut = si_thickness
            # Means that electron reached a trench depth
            z_top_cut = si_thickness - trench_depth

            # if holes_av_region:
            #     z_top_cut = -si_thickness / 2 + d_star + electron_width + holes_width
            #     z_bottom_cut = -si_thickness / 2 + d_star + electron_width

            if x_i == 0 and y_j == 0:
                left_cut = - pitch_width/2
                bottom_cut = - pitch_width/2
                right_cut = top_cut = pitch_width/2

            av_data[str(x_i) + ':' + str(y_j)] = Get_data_little_cuts(data, (left_cut, right_cut), (bottom_cut, top_cut), (z_top_cut, np.inf), binning=nbins)
            av_data[str(-x_i) + ':' + str(-y_j)] = Get_data_little_cuts(data, (-right_cut, -left_cut), (-top_cut, -bottom_cut), (z_top_cut, np.inf), binning=nbins)
            av_data[str(-x_i) + ':' + str(y_j)] = Get_data_little_cuts(data, (-right_cut, -left_cut), (bottom_cut, top_cut), (z_top_cut, np.inf), binning=nbins)
            av_data[str(x_i) + ':' + str(-y_j)] = Get_data_little_cuts(data, (left_cut, right_cut), (-top_cut, -bottom_cut), (z_top_cut, np.inf), binning=nbins)
    return av_data
 

# Load geometry. Geometry is in mm -> change to mkm
# Unify units later
filename_params = '../Data/Geometry/Test_data.json' #sys.argv[2]
file_params = open(filename_params)
sipm_params = json.load(file_params)
file_params.close()
# for key, value in sipm_params.items():
#     if key == '_typename' or key == 'n_bins' or key == 'n_cells':
#         continue
#     sipm_params[key] *= 1e3
print(sipm_params)

pitch_width = sipm_params['pitch_width']
av_width = sipm_params['av_width']
trench_width = sipm_params['trench_width']
trench_depth = sipm_params['trench_depth']
si_thickness = sipm_params['si_thickness']
number_of_cells = sipm_params['n_cells']


filename = "" # root with electrons information
if sipm_params["number_of_electrons_generated"] > 1000:
    filename = "../Data/ElectronDriftOutput/" + str( int(sipm_params["number_of_electrons_generated"] / 1000)) + "k_points.root";
else:
    filename = "../Data/ElectronDriftOutput/" + str( int(sipm_params["number_of_electrons_generated"] / 1000)) + "_points.root";

file = uproot.open(filename) #"../output/test.root"

data = file["ntp"].arrays(["initialX", "initialY", "initialZ", "initialT", "finalX", "finalY", "finalZ", "finalT"], library='pd')
plt.hist2d(data.finalX, data.finalY, bins=100)
plt.show()

number_of_electrons = len(data)
cutted_data = Get_data_after_cuts(data, cells=sipm_params["n_cells"], nbins=sipm_params["n_bins"])

survived_electrons = 0
json_output = {}
for i in range(-(number_of_cells//2), number_of_cells//2 + 1):
    for j in range(-(number_of_cells//2), number_of_cells//2 + 1):
        value = cutted_data[str(i) + ":" + str(j)]
        json_output[str(i) + ":" + str(j)] = {}
        for ii in range(sipm_params["n_bins"]):
            for jj in range(sipm_params["n_bins"]):
                bins_value = value[str(ii) + ":" + str(jj)]
                if bins_value.empty:
                    json_output[str(i) + ":" + str(j)][str(ii) + ":" + str(jj)] = 0
                else:
                    survived_electrons += len(bins_value) / number_of_electrons
                    json_output[str(i) + ":" + str(j)][str(ii) + ":" + str(jj)] = len(bins_value) / number_of_electrons
 

# print(json_output)
json_string = json.dumps(json_output)
# print(json_string)
# Using a JSON string

outputFilename = ""
if sipm_params["number_of_electrons_generated"] > 1000:
    outputFilename = "../Data/Geant4_electrons_input/" + str( int(sipm_params["number_of_electrons_generated"] / 1000)) + "k_points.json";
else:
    outputFilename = "../Data/Geant4_electrons_input/" + str( int(sipm_params["number_of_electrons_generated"] / 1000)) + "_points.json";

with open(outputFilename, 'w') as outfile:
    outfile.write(json_string)

print(survived_electrons)


