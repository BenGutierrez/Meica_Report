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
import parse as parse
from matplotlib.ticker import MaxNLocator




fig = plt.figure(1)
gs = gridspec.GridSpec(3,3)
ax1 = fig.add_subplot(gs[0,0])
components = kappa_vs_rho.file_parse(sys.argv[1])
accept, reject, middle, ignore = kappa_vs_rho.split_components(sys.argv[1], components)
plot = kappa_vs_rho.kappa_vs_rho_plot(accept, reject, middle, ignore)

ax2 = fig.add_subplot(gs[1,0])
plot = overlay_mass.overlay(sys.argv[2],sys.argv[3],threshold=12.34, index = 7)

plt.show()