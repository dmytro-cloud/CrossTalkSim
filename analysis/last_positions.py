#!/usr/bin/python

import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import json
import time

# Disable indexing warnings
pd.options.mode.chained_assignment = None  # default='warn'


start = time.time()

# filename = sys.argv[1] # File that you will obtaine probabilities of the crosstalk from and put as initial light (100% prob for appearing)
# filename_params = sys.argv[2]
# filename = '../Data/Geant4_output/FBK_8bins_W/FBK_W_binned_4_4.root'
filename = '../Data/Geant4_output/Backside_4bins_W/Backside_W_binned_2_2.root'
# filename = '../Data/Geant4_output/Backside_4bins_Si/Backside_polysi_binned_2_2.root'
file = uproot.open(filename) #"../output/test.root")
data = file["ntp"].arrays(["second_last_x", "second_last_y", "second_last_z",
     "photon_last_x", "photon_last_y", "photon_last_z", "photon_initial_wl"], library='pd')

plt.hist2d(data.photon_last_x, data.photon_last_y, bins=100, norm=mpl.colors.LogNorm(),  range=((-0.4, 0.4), (-0.4, 0.4)))
plt.show()
plt.hist2d(data.photon_last_x, data.photon_last_z, bins=100, norm=mpl.colors.LogNorm(),  range=((-0.4, 0.4), (-0.02, 0.02)))
plt.show()