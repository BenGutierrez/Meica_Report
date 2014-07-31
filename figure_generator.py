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
import kappa_vs_rho
import sys
import os
import analysis
fig = plt.figure()
components = kappa_vs_rho.file_parse(sys.argv[5])
accept, reject, middle, ignore = kappa_vs_rho.split_components(sys.argv[5], components)
plot = kappa_vs_rho.kappa_vs_rho_plot(accept, reject, middle, ignore)

com = accepted_figures.center_of_mass(anatomical = sys.argv[1], overlay = sys.argv[2])
maps = accepted_figures.collect_data(anatomical = sys.argv[1], overlay = sys.argv[2], threshold_map = sys.argv[3])
# path = os.getcwd()
# if len(sys.argv) == 6:
# 	TS = str(sys.argv[4])
# else:
#  	TS = ''
# os.system('rm -r png_dump')
# os.mkdir('png_dump')
# os.chdir(path+'/png_dump')
# plt.savefig('kappa_vs_rho')
threshold = 1.96
# N = 0
# for i in range(maps[1].shape[3]):
# 	if i in accept:
# 		plot1 = accepted_figures.overlay_axial(maps, com, threshold = threshold, alpha = 0.8, 
# 			index = i, series = TS, accepted = 0, Counter = N)
#   		#plot2 = accepted_figures.overlay_coronal(maps, com, threshold = threshold, alpha = 0.8, 
#   			#index = i, series = TS, accepted = 0)
#   		plot3 = accepted_figures.overlay_sagital(maps, com, threshold = threshold, alpha = 0.8, 
#   			index = i, series = TS, accepted = 0, Counter = N)
#   		N = N + 1
# 	else: 
# 		plot1 = accepted_figures.overlay_axial(maps, com, threshold = threshold, alpha = 0.8, 
# 			index = i, series = TS, accepted = 1)
#   		#plot2 = accepted_figures.overlay_coronal(maps, com, threshold = threshold, alpha = 0.8, 
#   			#index = i, series = TS, accepted = 1)
#   		plot3 = accepted_figures.overlay_sagital(maps, com, threshold = threshold, alpha = 0.8, 
#   			index = i, series = TS, accepted = 1)
# os.chdir(path)

analysis.analysis_rst(accept, reject, middle, ignore, threshold)

