#!/usr/bin/env python
"""
Gutierrez, B.  Displayes an anatomical image with a super imposed
functional image overlayed.  Files must be in NIFTI format.
"""

"""
The function overlay accepts an two datasets.  The first used in meica.py
as the anatomical, the second from the selection of output files from meica (i.e medn, mefc, mefl, tsoc).
The threshold can be set to any positive number that controls how much of the overlay is seen.  Alpha 
is the blending of the image, i.e the transparency of the overlay and can be anythng from 0-1.
Index is the component number of the overlay that is to be viewd.  For x,y,z, represent the slice number
to be viewed in the x,y,z dimension.  Only sone should be defined in the function at a time.

import matplotlib.pyplot as plt
import nibabel as ni
import numpy as np
import scipy

overlay_mass.overlay(anatomical= '/Users/gutierrezbe/Documents/NIFTI/WE/meica.WE_TE1234Run_5/20140528_094651MPRAGE1isoG2s013a1001_ns_at.nii.gz',
        overlay= '/Users/gutierrezbe/Documents/NIFTI/WE/meica.WE_TE1234Run_5/WE_meica_Run_5_mefl.nii.gz',threshold =12.34,
        index =7 ,alpha=0.6)

what
"""
import matplotlib.pyplot as plt
import numpy as np
import nibabel as ni
import scipy.ndimage

def overlay(anatomical, overlay, threshold, index, alpha = 0.8, axial = 0, cornial = 0, sagital = 0):

	anatomical = ni.load(anatomical)
	overlay = ni.load(overlay)
	anatomical_data = anatomical.get_data()
	overlay_data = overlay.get_data()
	com = np.empty(shape = (overlay_data.shape[3],4))
	for i in range(overlay_data.shape[3]):
		com[i,0] = i
		com[i,1:4] = np.asarray(scipy.ndimage.measurements.center_of_mass(abs(overlay_data[:,:,:,i])))

	#gs = gridspec.GridSpec(1,3)
	if axial == 0:
		#ax1 = fig.add_subplot(gs[0,0])
		overlay_axial(anatomical_data, overlay_data, com, threshold, alpha, index)
	if cornial == 0:
		#ax2 = fig.add_subplot(gs[0,1])
		overlay_cornial(anatomical_data, overlay_data, com, threshold, alpha, index)
	if sagital == 0: 
		#ax3 = fig.add_subplot(gs[0,2])
		overlay_sagital(anatomical_data, overlay_data, com, threshold, alpha, index)

def overlay_axial(anatomical_data, overlay_data, com, threshold, alpha, index):
	
	fig = plt.figure()
	n = com[index,3]/overlay_data.shape[2]
	threshold_overlay = np.empty(shape=(overlay_data.shape[0],overlay_data.shape[1]))

	for i , j in np.ndenumerate(overlay_data[:,:,round(com[index,3]),index]):
		if abs(overlay_data[i[0],i[1],round(com[index,3]),index]) <= threshold:
			threshold_overlay[i[0]][i[1]] = np.nan
		else:
			threshold_overlay[i[0]][i[1]] = abs(overlay_data[i[0],i[1],round(com[index,3]),index])

	plt.imshow(anatomical_data[:,:,round(anatomical_data.shape[2]*n)].T,cmap = 'Greys_r',origin='lower')
	plt.imshow(threshold_overlay.T,cmap = 'RdYlGn',extent=(0,anatomical_data.shape[0],0,anatomical_data.shape[1]),
		alpha=alpha,vmin=threshold,origin='lower', interpolation='bilinear')
	threshold_overlay= [0]

def overlay_cornial(anatomical_data, overlay_data, com, threshold, alpha, index):
	fig = plt.figure()
	n = com[index,2]/overlay_data.shape[1]
	threshold_overlay = np.empty(shape=(overlay_data.shape[0],overlay_data.shape[2]))

	for i , j in np.ndenumerate(overlay_data[:,round(com[index,2]),:,index]):
		if abs(overlay_data[i[0],round(com[index,2]),i[1],index]) <= threshold:
			threshold_overlay[i[0]][i[1]] = np.nan
		else:
			threshold_overlay[i[0]][i[1]] = abs(overlay_data[i[0],round(com[index,2]),i[1],index])

	plt.imshow(anatomical_data[:,round(anatomical_data.shape[1]*n),:].T,cmap = 'Greys_r',origin='lower')
	plt.imshow(threshold_overlay.T,cmap = 'RdYlGn',extent=(0,anatomical_data.shape[0],0,anatomical_data.shape[2]),
		alpha=alpha,vmin=threshold,origin='lower',interpolation='bilinear')
	threshold_overlay= [0]

def overlay_sagital(anatomical_data, overlay_data, com, threshold, alpha, index):
	fig = plt.figure()
	n = com[index,1]/overlay_data.shape[0]
	threshold_overlay = np.empty(shape=(overlay_data.shape[1],overlay_data.shape[2]))

	for i , j in np.ndenumerate(overlay_data[round(com[index,1]),:,:,index]):
		if abs(overlay_data[round(com[index,3]),i[0],i[1],index]) <= threshold:
			threshold_overlay[i[0]][i[1]] = np.nan
		else:
			threshold_overlay[i[0]][i[1]] = abs(overlay_data[round(com[index,3]),i[0],i[1],index])


	plt.imshow(anatomical_data[round(anatomical_data.shape[0]*n),:,:].T,cmap = 'Greys_r',origin='lower')
	plt.imshow(threshold_overlay.T[:,::-1],cmap = 'RdYlGn',extent=(anatomical_data.shape[1],0,0,anatomical_data.shape[2]),
		alpha=alpha,vmin=threshold,origin='lower',interpolation='bilinear')