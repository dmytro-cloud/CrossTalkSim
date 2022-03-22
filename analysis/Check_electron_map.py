import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import json


filename_electrons = '../Data/Geant4_electrons_input/100k_point_8bins.json' #sys.argv[2]
file_electrons = open(filename_electrons)
electrons_map = json.load(file_electrons)
file_electrons.close()

weights = []
x = []
y = []
for ii in range(8):
	for jj in range(8):
		weights.append(electrons_map['0:0'][str(ii) + ':' + str(jj)])
		x.append(ii)
		y.append(jj)

plt.hist2d(x, y, weights=weights, bins=8)
plt.show()