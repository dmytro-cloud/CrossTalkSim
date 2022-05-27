# CrossTalkSim
SiPM optical cross-talk simulation

Scripts that work with the data produced by Geant4 part (https://gitlab.triumf.ca/phaar/g4sipm/-/tree/crosstalkJulyUoA) of the code.

Files sctructure:

Electron_drift - contains code to describe electron propagation within SiPM and create a map of avalanches

Data - contains Geometry of a SiPM and also intermediate outputs of scripts like a map of avalanches

Supporting_scripts - additional scripts that help to operate. Currently contains only a script to generate input mac file to geant4 simulation

Analysis - contains scripts to analyse the data.
	Full_analysis_script.py - the main analysis script that calculates crosstalk and ammount of emitted light.