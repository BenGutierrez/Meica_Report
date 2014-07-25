#!/usr/bin/env python
"""
Gutierrez, B.  call on other functions in order to make a meica report form

"""
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import sys
import kappa_vs_rho
import overlay_mass

fig = plt.figure(1)
gs = gridspec.GridSpec(6,5)
ax1 = fig.add_subplot(gs[0:3,0:3])
components = kappa_vs_rho.file_parse(sys.argv[1])
accept, reject, middle, ignore = kappa_vs_rho.split_components(sys.argv[1], components)
plot = kappa_vs_rho.kappa_vs_rho_plot(accept, reject, middle, ignore)
N = 0
anatomical_data, overlay_data, ant_hdr, hdr, com = overlay_mass.center_of_mass(anatomical = sys.argv[2], overlay = sys.argv[3])


for i in range(3,6):
	fig.add_subplot(gs[i,0])
	plot1 = overlay_mass.overlay_axial(anatomical_data, overlay_data, com, threshold=12.34, alpha =0.8,
	 index=N, overlay_hdr = hdr, anatomical_hdr = ant_hdr)
	fig.add_subplot(gs[i,1])
	plot2 = overlay_mass.overlay_cornial(anatomical_data, overlay_data, com, threshold=12.34, alpha=0.8,
	 index=N, overlay_hdr = hdr, anatomical_hdr = ant_hdr)
	fig.add_subplot(gs[i,2])
	plot3 = overlay_mass.overlay_sagital(anatomical_data, overlay_data, com, threshold=12.34, alpha=0.8,
	 index=N, overlay_hdr = hdr, anatomical_hdr = ant_hdr)
	N = N + 1

plt.show()