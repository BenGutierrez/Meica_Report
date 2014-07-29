#!/usr/bin/env python
"""
Gutierrez, B.  Generates three png files for each accepted components along with time series of these components
in each png file.

sample input
meica_figure.py ~/Documents/NIFTI/WE/meica.WE_TE1234Run_4/20140528_094651MPRAGE1isoG2s013a1001_ns_at.nii.gz 
~/Documents/NIFTI/WE/meica.WE_TE1234Run_4/WE_meica_Run_4_mefc.nii.gz ~/Documents/NIFTI/WE/meica.WE_TE1234Run_4/TED/feats_OC2.nii
 ~/Documents/NIFTI/WE/meica.WE_TE1234Run_4/TED/meica_mix.1D
"""
import matplotlib.pyplot as plt
import accepted_figures
import numpy as np
import sys
import os

com = accepted_figures.center_of_mass(anatomical = sys.argv[1], overlay = sys.argv[2], threshold_map = sys.argv[3])
maps = accepted_figures.collect_data(anatomical = sys.argv[1], overlay = sys.argv[2], threshold_map = sys.argv[3])
path = os.getcwd()
if len(sys.argv) == 5:
	TS = str(sys.argv[4])
else:
 	TS = ''
os.system('rm -r png_dump')
os.mkdir('png_dump')
os.chdir(path+'/png_dump')

# for i in range(maps[1].shape[3]):
for i in range(1):
	plot1 = accepted_figures.overlay_axial(maps, com, threshold = 1.96, alpha = 0.8, 
		index = i, series = TS)
  	plot2 = accepted_figures.overlay_cornial(maps, com, threshold = 1.96, alpha = 0.8, 
  		index = i, series = TS)
  	plot3 = accepted_figures.overlay_sagital(maps, com, threshold = 1.96, alpha = 0.8, 
  		index = i, series = TS)
os.chdir(path)