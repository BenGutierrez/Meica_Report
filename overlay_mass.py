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
"""
import matplotlib.pyplot as plt
import numpy as np
import nibabel as ni
import scipy.ndimage
from matplotlib import gridspec

def center_of_mass(anatomical, overlay):
	anatomical = ni.load(anatomical)
	overlay = ni.load(overlay)
	anatomical_data = anatomical.get_data()
	overlay_data = overlay.get_data()
	anatomical_hdr = anatomical.get_header()
	overlay_hdr = overlay.get_header()
	com = np.empty(shape = (overlay_data.shape[3],4))
	for i in range(overlay_data.shape[3]):
		com[i,0] = i
		com[i,1:4] = np.asarray(scipy.ndimage.measurements.center_of_mass(abs(overlay_data[:,:,:,i])))
	return(anatomical_data, overlay_data, anatomical_hdr, overlay_hdr, com)

def overlay_axial(anatomical_data, overlay_data, com, threshold, alpha, index, overlay_hdr, anatomical_hdr):
	
	n = com[index,3]/overlay_data.shape[2]
	threshold_overlay = np.empty(shape = (overlay_data.shape[0],overlay_data.shape[1]))
	a_xmax = anatomical_data[:,:,round(anatomical_data.shape[2]*n)].shape[0]*anatomical_hdr['srow_x'][0] + anatomical_hdr['srow_x'][3]
	a_ymax = anatomical_data[:,:,round(anatomical_data.shape[2]*n)].shape[1]*anatomical_hdr['srow_y'][1] + anatomical_hdr['srow_y'][3]
	xmax = overlay_data.shape[0]*overlay_hdr['srow_x'][0] + overlay_hdr['srow_x'][3]
	ymax = overlay_data.shape[1]*overlay_hdr['srow_y'][1] + overlay_hdr['srow_y'][3]
	a_xmin = anatomical_hdr['srow_x'][3]
	a_ymin = anatomical_hdr['srow_y'][3]
	xmin = overlay_hdr['srow_x'][3]
	ymin = overlay_hdr['srow_y'][3]

	for i , j in np.ndenumerate(overlay_data[:,:,round(com[index,3]),index]):
	 	if abs(overlay_data[i[0],i[1],round(com[index,3]),index]) <= threshold:
	 		threshold_overlay[i[0]][i[1]] = np.nan
	 	else:
	 		threshold_overlay[i[0]][i[1]] = abs(overlay_data[i[0],i[1],round(com[index,3]),index])

	plt.imshow(anatomical_data[:,:,round(anatomical_data.shape[2]*n)].T, cmap = 'Greys_r', 
		origin = 'lower', interpolation = 'nearest', extent = [a_xmin,a_xmax,a_ymin,a_ymax])
	plt.imshow(threshold_overlay.T, cmap = 'RdYlGn', extent = [xmin,xmax,ymin,ymax],
		alpha = alpha, vmin = threshold, origin = 'lower', interpolation='nearest')
	plt.axis('off')

def overlay_cornial(anatomical_data, overlay_data, com, threshold, alpha, index, overlay_hdr, anatomical_hdr):
	n = com[index,2]/overlay_data.shape[1]
	threshold_overlay = np.empty(shape = (overlay_data.shape[0],overlay_data.shape[2]))
	a_xmax = anatomical_data[:,round(anatomical_data.shape[1]*n),:].shape[0]*anatomical_hdr['srow_x'][0] + anatomical_hdr['srow_x'][3]
	a_zmax = anatomical_data[:,round(anatomical_data.shape[1]*n),:].shape[1]*anatomical_hdr['srow_z'][2] + anatomical_hdr['srow_z'][3]
	xmax = overlay_data.shape[0]*overlay_hdr['srow_x'][0] + overlay_hdr['srow_x'][3]
	zmax = overlay_data.shape[2]*overlay_hdr['srow_z'][2] + overlay_hdr['srow_z'][3]
	a_xmin = anatomical_hdr['srow_x'][3]
	a_zmin = anatomical_hdr['srow_z'][3]
	xmin = overlay_hdr['srow_x'][3]
	zmin = overlay_hdr['srow_z'][3]

	for i , j in np.ndenumerate(overlay_data[:,round(com[index,2]),:,index]):
		if abs(overlay_data[i[0],round(com[index,2]),i[1],index]) <= threshold:
			threshold_overlay[i[0]][i[1]] = np.nan
		else:
			threshold_overlay[i[0]][i[1]] = abs(overlay_data[i[0],round(com[index,2]),i[1],index])

	plt.imshow(anatomical_data[:,round(anatomical_data.shape[1]*n),:].T,cmap = 'Greys_r', 
		origin = 'lower', interpolation = 'nearest', extent = [a_xmin,a_xmax,a_zmin,a_zmax])
	plt.imshow(threshold_overlay.T,cmap = 'RdYlGn',extent = [xmin,xmax,zmin,zmax],
		alpha = alpha, vmin = threshold, origin = 'lower', interpolation = 'nearest')
	plt.axis('off')

def overlay_sagital(anatomical_data, overlay_data, com, threshold, alpha, index, overlay_hdr, anatomical_hdr):
	n = com[index,1]/overlay_data.shape[0]
	threshold_overlay = np.empty(shape = (overlay_data.shape[1],overlay_data.shape[2]))
	a_ymax = anatomical_data[round(anatomical_data.shape[0]*n),:,:].shape[0]*anatomical_hdr['srow_y'][1] + anatomical_hdr['srow_y'][3]
	a_zmax = anatomical_data[round(anatomical_data.shape[0]*n),:,:].shape[1]*anatomical_hdr['srow_z'][2] + anatomical_hdr['srow_z'][3]
	ymax = overlay_data.shape[1]*overlay_hdr['srow_y'][1] + overlay_hdr['srow_y'][3]
	zmax = overlay_data.shape[2]*overlay_hdr['srow_z'][2] + overlay_hdr['srow_z'][3]
	a_ymin = anatomical_hdr['srow_y'][3]
	a_zmin = anatomical_hdr['srow_z'][3]
	ymin = overlay_hdr['srow_y'][3]
	zmin = overlay_hdr['srow_z'][3]

	for i , j in np.ndenumerate(overlay_data[round(com[index,1]),:,:,index]):
		if abs(overlay_data[round(com[index,3]),i[0],i[1],index]) <= threshold:
			threshold_overlay[i[0]][i[1]] = np.nan
		else:
			threshold_overlay[i[0]][i[1]] = abs(overlay_data[round(com[index,3]),i[0],i[1],index])

	plt.imshow(anatomical_data[round(anatomical_data.shape[0]*n),:,:].T[:,::-1],cmap = 'Greys_r',
		origin ='lower', interpolation = 'nearest', extent = [a_ymax,a_ymin,a_zmin,a_zmax])
	plt.imshow(threshold_overlay.T[:,::-1],cmap = 'RdYlGn', extent = [ymax,ymin,zmin,zmax],
		alpha = alpha, vmin = threshold, origin = 'lower', interpolation = 'nearest')
	plt.axis('off')




