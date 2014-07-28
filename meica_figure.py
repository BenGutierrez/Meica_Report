#!/usr/bin/env python
"""
Gutierrez, B.  call on other functions in order to make a meica report form

"""
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib
import numpy as np
import sys
import kappa_vs_rho
import overlay_mass
from parse import parse
import os

fig = plt.figure()
gs0 = gridspec.GridSpec(2,2)
ax1 = fig.add_subplot(gs0[0,0])
components = kappa_vs_rho.file_parse(sys.argv[1])
accept, reject, middle, ignore = kappa_vs_rho.split_components(sys.argv[1], components)
plot = kappa_vs_rho.kappa_vs_rho_plot(accept, reject, middle, ignore)
N = 0
(anatomical_data, overlay_data, threshold_map, ant_hdr, hdr, com) = overlay_mass.center_of_mass(anatomical = sys.argv[2], overlay = sys.argv[3], threshold_map = sys.argv[4])

path = parse('{:S}'+'/feats_OC2.nii',sys.argv[4])[0]+'/meica_mix.1D'

gs00 = gridspec.GridSpecFromSubplotSpec(6,6, subplot_spec = gs0[1,0],wspace = -.5, hspace = 0)
for i in range(0,6,2):
	fig.add_subplot(gs00[i:i+2,0:2])
	plot1 = overlay_mass.overlay_axial(anatomical_data, overlay_data, threshold_map, com, threshold = 1.96, alpha = 0.8,
	 index = N, overlay_hdr = hdr, anatomical_hdr = ant_hdr)
	#fig.add_subplot(gs00[i+2,0:2])
	#os.system("csh -c 1dplot %(path)s[%(N)d]" % {'path': path, 'N': N})
	#print("csh -c 1dplot %(path)s[%(N)d]" % {'path': path, 'N': N})
	#subprocess.call('csh',shell=True)
	fig.add_subplot(gs00[i:i+2,2:4])
	plot2 = overlay_mass.overlay_cornial(anatomical_data, overlay_data, threshold_map, com, threshold = 1.96, alpha = 0.8,
	 index = N, overlay_hdr = hdr, anatomical_hdr = ant_hdr)
	fig.add_subplot(gs00[i:i+2,4:6])
	plot3 = overlay_mass.overlay_sagital(anatomical_data, overlay_data, threshold_map, com, threshold = 1.96, alpha = 0.8,
	 index = N, overlay_hdr = hdr, anatomical_hdr = ant_hdr)
	N = N + 1

plt.show()