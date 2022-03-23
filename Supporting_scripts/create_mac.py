#!/usr/bin/python

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import json

filename_params = '../Data/Geometry/Test_data.json' #sys.argv[2]
file_params = open(filename_params)
sipm_params = json.load(file_params)
file_params.close()
# for key, value in sipm_params.items():
#     if key == '_typename' or key == 'n_bins' or key == 'n_cells':
#         continue
#     sipm_params[key] *= 1e3
print(sipm_params)

output_string = '''
# Sets some default verbose

/control/verbose 1
/run/verbose 1
#/run/setCut 0.0001 mm
/event/verbose 0
/tracking/verbose 0

# Output file name
/SiPM/output/filename ../Data/Geant4_output/'''

output_string += 'test.root' + '\n' #sipm_params['output_filename'] + '\n'

output_string += '''
# Adding geometery control(s)
/SiPM/det/Pitch ''' + str(sipm_params['pitch_width']) + ' um\n'
output_string += '/SiPM/det/SiTh ' + str(sipm_params['si_thickness']) + ' um\n'

output_string += '''# Same for generator
/SiPM/gen/SiTh ''' +  str(sipm_params['si_thickness']) + ' um\n\n'

output_string += '/SiPM/det/TrD ' + str(sipm_params['trench_depth']) + ' um\n'  		# Trench Depth
output_string += '/SiPM/det/AvW ' + str(sipm_params['av_width']) + ' um\n' 		# Avalanche d. Width
output_string += '''# Same for generator
/SiPM/gen/AvWidthGen ''' + str(sipm_params['av_width']) + ' um\n\n'     # Avalanche Width for generator

output_string += '/SiPM/det/AvTh ' + str(sipm_params['avalanche_thickness']) + ' um\n'
output_string += '/SiPM/gen/AvTh ' + str(sipm_params['avalanche_thickness']) + ' um\n'
output_string += '/SiPM/det/AvD ' + str(sipm_params['avalanche_depth']) + ' um\n'           
output_string += '/SiPM/gen/AvD ' + str(sipm_params['avalanche_depth']) + ' um\n\n'

output_string += '/SiPM/det/SiO2FTh ' + str(sipm_params['si_o_front_thickness']) + ' um\n'
output_string += '/SiPM/det/SiO2BTh ' + str(sipm_params['si_o_back_thickness']) + ' um\n'
output_string += '/SiPM/det/TrBW ' + str(sipm_params['trench_width']) + ' um\n' 		# Trench Bottom Width
output_string += '/SiPM/det/TrTW ' + str(sipm_params['trench_width']) + ' um\n' 		# So far the same as bottom
if sipm_params['polysi_for_trenches'] == 0:
	output_string += '/SiPM/det/PolySi false\n'
else:
	output_string += '/SiPM/det/PolySi true\n'

output_string += '''
/run/initialize

# Uncomment for visualization
#/vis/open HepRepFile
#/vis/scene/create
#/vis/scene/add/volume
#/vis/sceneHandler/attach
#/vis/scene/add/trajectories

# Store trajectories for tracking
/tracking/storeTrajectory 1

/gun/energy 0.01 keV

/run/beamOn 10000 

'''

f = open('../Data/Geant4_mac_files/test.mac', 'w')
f.write(output_string)
f.close()
# print(output_string)












