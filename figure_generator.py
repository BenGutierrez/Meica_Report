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

# startdir
# anat
# prefix
# setname
# os.system
threshold = 1.96
path = os.getcwd()
os.system('rm -r png_dump')
os.mkdir('png_dump')
os.chdir(path+'/png_dump')
components = kappa_vs_rho.file_parse(sys.argv[5])
plot = kappa_vs_rho.kr_vs_component(sys.argv[5])
plt.savefig('kappa_rho_vs_components')
plt.close()
fig = plt.figure()
accept, reject, middle, ignore = kappa_vs_rho.split_components(sys.argv[5], components)
plot = kappa_vs_rho.kappa_vs_rho_plot(accept, reject, middle, ignore)
plt.savefig('kappa_vs_rho')
plt.close()
accepted_figures.tsnr('/Users/gutierrezbe/Documents/NIFTI/WE/meica.WE_TE1234Run_4/WE_meica_Run_4_tsoc.nii.gz','/Users/gutierrezbe/Documents/NIFTI/WE/meica.WE_TE1234Run_4/WE_meica_Run_4_medn.nii.gz')
maps = accepted_figures.collect_data(anatomical = sys.argv[1], overlay = sys.argv[2], threshold_map = sys.argv[3])

if len(sys.argv) == 6:
	TS = str(sys.argv[4])
else:
 	TS = ''

plot2 = accepted_figures.montage(maps, accept, threshold = threshold, alpha = 0.8, 
	series = TS, Axial = 0, Sagital = 0, Coronal = 1)
os.chdir(path)

analysis.diagnostics_rst()
analysis.analysis_rst(accept, reject, middle, ignore, threshold, sys.argv[5])

