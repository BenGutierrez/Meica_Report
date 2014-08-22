#!/usr/bin/env python
"""
Gutierrez, B.  Contains a function for finding the center of mass of 
fMRI image componants.  Also contains functions for overlaying 
a 2D fMRI image onto a 2D anatomical image in the axial, sagital, and cortial
plane
"""
import matplotlib as mpl
mpl.use('Agg')

from matplotlib.colors import LinearSegmentedColormap
from numpy.core.umath_tests import inner1d
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec
import matplotlib.pyplot as plt
from parse import parse
import nibabel as ni
import numpy as np
import os

cdict = {'red':  ((0.0,  0.0, 0.0),
						   (0.2, 0.6, 0.6),
						   (0.35, 0.9, 0.9),
						   (0.5, 1.0, 1.0),
		                   (1.0, 1.0, 1.0)),

		         'green': ((0.0, 1.0, 1.0),
		         		   (0.5, 1.0, 1.0),
		         		   (0.65, 0.9, 0.9),
		                   (0.9, 0.5, 0.5),
		                   (1.0, 0.0, 0.0)),

		         'blue':  ((0.0, 0.0, 0.0),
		                   (1.0, 0.0, 0.0))
		        }

RGB = LinearSegmentedColormap('RGB', cdict)



"""
Removes two dimensional slices from 3d matrix image that contain all non-zero elements
image: three dim array
axis: axial that the two dimmensional plane fixes, i.e. (0,1) fixes z
"""
def file_parse(file):
	txt_file = open(str(file))
	components = [None,'#REJ ','#MID ','#IGN ']
	txt_file.seek(0)
	while components[0] == None:
		line = txt_file.readline()
		if line == '':
			print('unable to parse the ctab txt_file')
			break
		else:	
			components[0] = parse('#ACC {:S} {}\n', line)
	for i in range(1,4):
		components[i] = parse(components[i]+'{:S} {}\n', txt_file.readline())
	for i in range(4):
		components[i] = map(int, components[i][0].split(','))
	txt_file.close()
	return components

"""
Take the infomation from the comp file and sperate it into 4 2D matricies
for the accepted, rejected, middle, ignore bins.
"""
def split_components(comp_table_title,components):
	a = 0
	b = 0
	c = 0
	d = 0
	comp_table = np.loadtxt(str(comp_table_title))
	accept = np.zeros(shape = (len(components[0]),5))
	reject = np.zeros(shape = (len(components[1]),5))
	middle = np.zeros(shape = (len(components[2]),5))
	ignore = np.zeros(shape = (len(components[3]),5))
	for i in comp_table[:,0]:
		if i in components[0]:
			accept[a,:] = comp_table[i,:]
			a = a + 1
		elif i in components[1]:
			reject[b,:] = comp_table[i,:]
			b = b + 1
		elif i in components[2]:
			middle[c,:] = comp_table[i,:]
			c = c + 1
		elif i in components[3]:
			ignore[d,:] = comp_table[i,:]
			d = d + 1
	return accept,reject,middle,ignore

def mask(image, axis):
	im_mask = np.zeros(image.shape)
	image[np.isnan(image)] = 0
	im_mask[image != 0] = 1
	if axis == (0,1):
		im_mask = np.sum(im_mask, axis = 0)
		im_mask = np.sum(im_mask, axis = 0)
		return image[:,:,im_mask > 0]
	elif axis == (1,2):
		im_mask = np.sum(im_mask, axis = 1)
		im_mask = np.sum(im_mask, axis = 1)
		return image[im_mask > 0,:,:]
	elif axis == (0,2):
		im_mask = np.sum(im_mask, axis = 0)
		im_mask = np.sum(im_mask, axis = 1)
		return image[:,im_mask > 0,:]

"""
set up floodfill algorithm.  acts as a clustering algorithm.  x,y,z designates where to begin algorithm in matrix
matrix: 2D array
x: integer index of matrix
y: integer index of matrix
z: integer index of matrix
"""
def flood(matrix, x, y, z):
	N = 0
	itemindex = [[],[],[]]
	original = matrix.copy()
	N = floodfill(matrix, x, y, z, N, itemindex)
	if N >= 20: # twenty is arbitrary choosen
		return original
	else:
		return matrix
"""
Rest of the flood fill algorithm.
itemindex: array containing location of non-zero elements of matrix next to (x,y,z)
"""
def floodfill(matrix, x, y, z, N, itemindex):
	itemindex[0] = itemindex[0][1:len(itemindex[0]) + 1]
	itemindex[1] = itemindex[1][1:len(itemindex[0]) + 1]
	itemindex[2] = itemindex[2][1:len(itemindex[0]) + 1]
	N += np.sum(matrix[max(x-1,0):min(x+2,matrix.shape[0]),max(y-1,0):min(y+2,matrix.shape[1]),max(z-1,0):min(z+2,matrix.shape[2])])
	if N <= 19: #20 -1
		matrix[x,y,z] = 0
		test = np.where(matrix[max(x-1,0):min(x+2,matrix.shape[0]-1),max(y-1,0):min(y+2,matrix.shape[1]-1),max(z-1,0):min(z+2,matrix.shape[2]-1)] == 1)
		itemindex[0].extend(test[0]+ x + (2 - (min(x+2,matrix.shape[0]) - max(x - 1,0)))) #(2-...) is a fudge factor for being on the edge of the matrix
		itemindex[1].extend(test[1]+ y + (2 - (min(y+2,matrix.shape[1]) - max(y - 1,0))))
		itemindex[2].extend(test[2]+ z + (2 - (min(z+2,matrix.shape[2]) - max(z - 1,0))))
		matrix[max(x-1,0):min(x+2,matrix.shape[0]),max(y-1,0):min(y+2,matrix.shape[1]),max(z-1,0):min(z+2,matrix.shape[2])] = 0
		if len(itemindex[0])>0:#ensures that if no more 1's around (x,y,z) not calling function again
			N = floodfill(matrix,itemindex[0][0],itemindex[1][0],itemindex[2][0], N, itemindex)
	return N

"""
Creates a montage of greyscale 10 images of axial, sagital and coronal views of meica components along with a time series of the component.
accepted components also get overlayed onto the anatomical image and floodfill is performed on statiscially significant voxels.
maps: array of 5 elements [anat dataset, mefl dataset, feats dataset, anat header, overlay header]
accept: array accepted components
threshold: float statiscial threshold
alpha: float alpha blending value
series: string meica_mix.1D path
axial: true or false. plot axial or not
sagital: true or false. plot sagital or not
coronal: true or false. plot coronal or not
"""
def montage(maps, accept, threshold, alpha, series ='', Axial = 0, Sagital = 0, Coronal = 0):
	anat = maps[0]
	overlay = maps[1]
	threshold_data = maps[2]
	anat_hdr = maps[3]
	overlay_hdr = maps[4]
	l = 0
	ax_extreme = np.zeros(8)
	ax = ['x','y']
	cor_extreme = np.zeros(8)
	co = ['x','z']
	sag_extreme = np.zeros(8)
	sa = ['y','z']

	for j in [0,1]:
		if Axial == 0:
			ax_extreme[j+6] = anat.shape[j]*anat_hdr['srow_%s' % ax[j]][j] + anat_hdr['srow_%s' % ax[j]][3]#setting up talairach coordinates
			ax_extreme[j+4] = overlay.shape[j]*overlay_hdr['srow_%s' % ax[j]][j] + overlay_hdr['srow_%s' % ax[j]][3]
			ax_extreme[j+2] = anat_hdr['srow_%s' % ax[j]][3]
			ax_extreme[j] = overlay_hdr['srow_%s' % ax[j]][3]
		if Sagital == 0:
			sag_extreme[j+6] = anat.shape[j+1]*anat_hdr['srow_%s' % sa[j]][j+1] + anat_hdr['srow_%s' % sa[j]][3]
			sag_extreme[j+4] = overlay.shape[j+1]*overlay_hdr['srow_%s' % sa[j]][j+1] + overlay_hdr['srow_%s' % sa[j]][3]
			sag_extreme[j+2] = anat_hdr['srow_%s' % sa[j]][3]
			sag_extreme[j] = overlay_hdr['srow_%s' % sa[j]][3]
		if Coronal == 0:
			cor_extreme[j+6] = anat.shape[j*2]*anat_hdr['srow_%s' % co[j]][j*2] + anat_hdr['srow_%s' % co[j]][3]
			cor_extreme[j+4] = overlay.shape[j*2]*overlay_hdr['srow_%s' % co[j]][j*2] + overlay_hdr['srow_%s' % co[j]][3]
			cor_extreme[j+2] = anat_hdr['srow_%s' % co[j]][3]
			cor_extreme[j] = overlay_hdr['srow_%s' % co[j]][3]
	for i in range(overlay.shape[3]):
		fig = plt.figure(figsize = (3.2*5,9 - (Axial + Sagital + Coronal)*2))
		gs0 = gridspec.GridSpec(3 - (Axial + Sagital + Coronal),10)
		if i in accept:
			overlay_acc = threshold_data[:,:,:,l].copy()
			overlay_acc = np.absolute(overlay_acc)
			overlay_acc[overlay_acc < threshold] = 0 #threshold mefl by feats dataset
			overlay_mask = np.zeros(overlay_acc.shape)
			overlay_mask[overlay_acc != 0] = 1

			itemindex = np.where(overlay_mask == 1)
			for j in range(len(itemindex[0])):
				overlay_mask = flood(overlay_mask,itemindex[0][j],itemindex[1][j],itemindex[2][j])#floodfill
			overlay_acc[overlay_mask == 0] = np.nan
			for j in range(10):#plot montage of accept component onto anatomical
				if Axial == 0:
					ax1 = fig.add_subplot(gs0[0,j])
					plt.imshow(anat[:,:,anat.shape[2]*j*.1].T, cmap = 'Greys_r', 
						origin = 'lower', interpolation = 'nearest', extent = [ax_extreme[2],ax_extreme[6],ax_extreme[3],ax_extreme[7]])
					bar = plt.imshow(overlay_acc[:,:,overlay_acc.shape[2]*j*.1].T, cmap = RGB, extent = [ax_extreme[0],ax_extreme[4],ax_extreme[1],ax_extreme[5]],
						alpha = alpha, origin = 'lower', interpolation = 'gaussian', vmin = threshold, vmax = 5)
					plt.axis('off')
				if Sagital == 0:
					ax2 = fig.add_subplot(gs0[1 - Axial,j])
					plt.imshow(anat[anat.shape[0]*j*.1,::-1,:].T, cmap = 'Greys_r', 
						origin = 'lower', interpolation = 'nearest', extent = [sag_extreme[2],sag_extreme[6],sag_extreme[3],sag_extreme[7]])
					bar = plt.imshow(overlay_acc[overlay_acc.shape[0]*j*.1,::-1,:].T, cmap = RGB, extent = [sag_extreme[0],sag_extreme[4],sag_extreme[1],sag_extreme[5]],
						alpha = alpha, origin = 'lower', interpolation = 'gaussian', vmin = threshold, vmax = 5)
					plt.axis('off')
				if Coronal == 0:
					ax3 = fig.add_subplot(gs0[2 - (Axial + Sagital),j])
					plt.imshow(anat[:,anat.shape[1]*j*.1,:].T, cmap = 'Greys_r', 
						origin = 'lower', interpolation = 'nearest', extent = [cor_extreme[2],cor_extreme[6],cor_extreme[3],cor_extreme[7]])
					bar = plt.imshow(overlay_acc[:,overlay_acc.shape[1]*j*.1,:].T, cmap = RGB, extent = [cor_extreme[0],cor_extreme[4],cor_extreme[1],cor_extreme[5]],
						alpha = alpha, origin = 'lower', interpolation = 'gaussian', vmin = threshold, vmax =5)
					plt.axis('off')
			l += 1
			
			gs1 = gridspec.GridSpec(1,1)
			ax4 = fig.add_subplot(gs1[0,0])
			if Axial + Sagital + Coronal == 0:
				fig.subplots_adjust(top = 0.3)
				gs0.tight_layout(fig, h_pad = -5, w_pad = 0.5, rect = [0,.3,.95,1])
			elif Axial + Sagital + Coronal == 1:
				fig.subplots_adjust(top = 0.21)
				gs0.tight_layout(fig, h_pad = -13, w_pad = 1.5, rect = [0,.2,.95,1])
			elif Axial + Sagital + Coronal == 2:
				fig.subplots_adjust(top = 0.21)
				gs0.tight_layout(fig, rect = [0,.2,.95,1])
			time_series = np.loadtxt(series)
		 	plt.plot(np.arange(time_series.shape[0]),time_series[:,i])
		 	plt.xlabel('Time (TR)', fontsize = 12)
		 	plt.ylabel('Arbitrary BOLD units', fontsize = 12)
			gs1.tight_layout(fig, rect = [0,0,.95,.35 + (Axial + Sagital + Coronal) * .1])
			right = max(gs0.right, gs1.right)
			left = max(gs0.left, gs1.left)
			gs0.update(left = left, right = right)
			gs1.update(left = left, right = right)
			cbar_ax = fig.add_axes([(gs0.right + ((gs0.right + 1)/2 - gs0.right)/2), gs1.bottom, .01, gs0.top - (gs1.bottom * 2)])
			fig.colorbar(bar, cax = cbar_ax)
			plt.ylabel('absoulte z score', fontsize = 12, rotation = 270)
			N = str(i)
			while len(N) < len(str(overlay.shape[3])):
				N = '0' + N
			plt.savefig('Accepted_Component_' + N)
			plt.close()
		
		fig = plt.figure(figsize = (3.2*5,9 - (Axial + Sagital + Coronal)*2))
		gs0 = gridspec.GridSpec(3 - (Axial + Sagital + Coronal),10)

		overlay_z = mask(overlay[:,:,:,i],(0,1))
		overlay_x = mask(overlay[:,:,:,i],(1,2))
		overlay_y = mask(overlay[:,:,:,i],(0,2))
		contrast = overlay_z[overlay_z != 0]
		maximum = np.percentile(contrast,98)
		minimum = np.percentile(contrast,2)

		for j in range(10):#plot greyscale component montage
			if Axial == 0:
				ax1 = fig.add_subplot(gs0[0,j])
				plt.imshow(overlay_z[:,:,overlay_z.shape[2]*j*.1].T, cmap = 'Greys_r',
					origin = 'lower', vmin = minimum, vmax = maximum)
				plt.axis('off')
			if Sagital == 0:
				ax2 = fig.add_subplot(gs0[1 - Axial,j])
				plt.imshow(overlay_y[overlay_y.shape[0]*j*.1,::-1,:].T, cmap = 'Greys_r',
					origin = 'lower', vmin = minimum, 
					vmax = maximum)
				plt.axis('off')
			if Coronal == 0:
				ax3 = fig.add_subplot(gs0[2 - (Axial + Sagital),j])
				plt.imshow(overlay_x[:,overlay_x.shape[1]*j*.1,:].T, cmap = 'Greys_r',
					origin = 'lower', vmin = minimum, 
					vmax = maximum)
				plt.axis('off')
		if Axial + Sagital + Coronal == 0:
			fig.subplots_adjust(top=0.3)
			gs0.tight_layout(fig, h_pad = -4, w_pad = 1, rect = [0,.3,1,1])
		elif Axial + Sagital + Coronal == 1:
			fig.subplots_adjust(top=0.21)
			gs0.tight_layout(fig, h_pad = -13, w_pad = 2.5, rect = [0,.2,1,1])
		elif Axial + Sagital + Coronal == 2:
			fig.subplots_adjust(top=0.15)
			gs0.tight_layout(fig, w_pad = 4, rect = [0,.2,1,1]) 

		gs1 = gridspec.GridSpec(1,1)
		ax4 = fig.add_subplot(gs1[0,0])
		time_series = np.loadtxt(series)
	 	plt.plot(np.arange(time_series.shape[0]),time_series[:,i])
	 	plt.xlabel('Time (TR)', fontsize = 12)
		plt.ylabel('Arbitrary BOLD units', fontsize = 12)
		gs1.tight_layout(fig, rect = [0,0,1,.35 + (Axial + Sagital + Coronal) * .1])
		right = max(gs0.right, gs1.right)
		left = max(gs0.left, gs1.left)
		gs0.update(left = left, right = right)
		gs1.update(left = left, right = right)

		N = str(i)
		while len(N) < len(str(overlay.shape[3])):
			N = '0' + N
		plt.savefig('Component_' + N)
		print ('++ figures created for Component %s' % N)
		plt.close()

def coreg(setname,anat):
	fig = plt.figure(figsize = (3.2*5,4))
	gs0 = gridspec.GridSpec(1,3)
	os.chdir('../%s' % setname)
	if os.path.isfile('ocv_uni_vrm_e3.nii'):
		os.system('rm -f ocv_uni_vrm*')
	if '.nii.gz' in anat[-7:]:
		anat_name = anat[:-7]
	if '.nii' in anat[-4:]:
		anat_name = anat[:-4]
	os.system('3dcalc -a ocv_uni_vr.nii.gz -b eBvrmask.nii.gz -expr "step(b)*a" -prefix ocv_uni_vrm')
	os.system('3drefit -view orig ocv_uni_vrm+tlrc.')
	os.system('@AddEdge ocv_uni_vrm+orig. %s+orig.' % anat_name)
	os.system('3dcalc -a ocv_uni_vrm_e3+orig -expr %s -prefix ocv_uni_vrm_e3.nii' % "a")
	anatomical = ni.load(anat).get_data()
	overlay = ni.load('ocv_uni_vrm_e3.nii').get_data()
	os.chdir('../png_dump')
	overlay[overlay == 0] = np.nan

	fig = plt.figure(figsize = (3.2*5,4))
	gs0 = gridspec.GridSpec(1,10)
	for i in np.arange(0,1,.1):#plot montage of accept component onto 
		ax1 = fig.add_subplot(gs0[0,int(i*10)])
		plt.imshow(anatomical[:,::-1,anatomical.shape[2]*i].T, cmap = 'Greys_r')
		plt.imshow(overlay[:,::-1,overlay.shape[2]*i].T, alpha = 0.8, cmap = RGB)
		plt.axis('off')
	gs0.tight_layout(fig, w_pad = -2)
	fig.subplots_adjust(right = 0.9)
	plt.savefig('coregistration')
	plt.close()

"""
Collects data from the given data sets.
"""
def collect_data(anatomical, overlay, threshold_map):
	anatomical = ni.load(anatomical)
	overlay = ni.load(overlay)
	threshold = ni.load(threshold_map)
	anat_data = anatomical.get_data()
	overlay_data = overlay.get_data()
	threshold_data = threshold.get_data()
	anat_hdr = anatomical.get_header()
	overlay_hdr = overlay.get_header()
	return(anat_data, overlay_data, threshold_data, anat_hdr, overlay_hdr)

"""
Makes TSNR figure of medn and medn/tsoc datasets
tsoc: string path to tsoc dataset
medn: string path to medn dataset
"""
def tsnr(tsoc,medn):

	RGB = LinearSegmentedColormap('RGB', cdict)

	medn_data = ni.load(medn).get_data()
	tsoc_data = ni.load(tsoc).get_data()
	medn_tsnr = medn_data.mean(-1)/medn_data.std(-1)
	tsoc_tsnr = tsoc_data.mean(-1)/tsoc_data.std(-1)
	frac = medn_tsnr/tsoc_tsnr

	medn_tsnr = mask(image = medn_tsnr, axis = (0,1))#remove z slices without nonzero elements
	tsoc_tsnr = mask(image = tsoc_tsnr, axis = (0,1))
	frac_tsnr = mask(image = frac, axis = (0,1))#remove z slices without nonzero elements
	medn_tsnr[medn_tsnr == 0] = np.nan
	tsoc_tsnr[tsoc_tsnr == 0] = np.nan
	frac_tsnr[frac_tsnr == 0] = np.nan
	background = np.zeros((medn_tsnr[:,:,0].T).shape)
	tsnr_options = [medn_tsnr,tsoc_tsnr,frac_tsnr]
	for j in range(3):
		tsnr = tsnr_options[j]
		fig = plt.figure(figsize = (3.2*5,4)) #medn tsnr figure
		gs0 = gridspec.GridSpec(1,10)
		#plot montage of medn TSNR
		tsnr_mask = tsnr[np.isnan(tsnr) == False]
		if j in [0,2]: 
			maximum = np.percentile(tsnr_mask,95)
			minimum = np.percentile(tsnr_mask,5)
		for i in range(0,10):
			ax1 = fig.add_subplot(gs0[0,i])
			plt.imshow(background, cmap = 'Greys_r')
			plot = plt.imshow(tsnr[:,::-1,i*.1*tsnr.shape[2]].T, vmin = minimum, vmax = maximum, cmap = RGB)
			plt.axis('off')
		gs0.tight_layout(fig, w_pad = -1, rect = [0,0,0.95,1])
		cbar = fig.add_axes([(gs0.right + ((gs0.right + 1)/2 - gs0.right)/2), gs0.bottom, .01, gs0.top - gs0.bottom])
		fig.subplots_adjust(right = 0.9)
		fig.colorbar(plot, cax = cbar)
		if j == 0:
			plt.savefig('medn_tsnr')
		elif j ==1:
			plt.savefig('tsoc_tsnr')
		else:
			plt.savefig('tsnr_ratio')
		plt.close()

	#plot histogram of the TSNR of medn
	medn_mask = medn_tsnr[np.isnan(medn_tsnr) == False]
	tsoc_mask = tsoc_tsnr[np.isnan(tsoc_tsnr) == False]
	frac_mask = frac_tsnr[np.isnan(frac_tsnr) == False]
	fig = plt.figure()
	plt.hist(medn_mask, bins = 100, range = [np.percentile(medn_mask,0.01),np.percentile(medn_mask,99.99)])
	plt.title('TSNR medn', fontsize = 15)
	plt.xlabel('TSNR', fontsize = 15)
	plt.ylabel('Frequency', fontsize = 15)
	plt.savefig('medn_tsnr_hist')
	plt.close()
	fig = plt.figure()
	plt.hist(tsoc_mask, bins = 100, range = [np.percentile(tsoc_mask,0.01),np.percentile(tsoc_mask,99.99)])
	plt.title('TSNR tsoc', fontsize = 15)
	plt.xlabel('TSNR', fontsize = 15)
	plt.ylabel('Frequency', fontsize = 15)
	plt.savefig('tsoc_tsnr_hist')
	plt.close()

	#plot histogram of the TSNR ratio of medn/tsnr
	fig = plt.figure()
	plt.hist(frac_mask, bins = 100, range = [np.percentile(frac_mask,0.01),np.percentile(frac_mask,99.99)])
	plt.title('TSNR medn / TSNR tsoc', fontsize = 15)
	plt.xlabel('TSNR ratio', fontsize = 15)
	plt.ylabel('Frequency', fontsize = 15)
	plt.savefig('tsnr_ratio_hist')
	plt.close()

"""
plot kappa vs rho and represent %varaince by the size of the markers.
accept: array of all accepted components
reject: array of all rejected components
middle: array of all middle kappa components
ignore: array of all ignore components
"""
def kappa_vs_rho_plot(accept,reject,middle,ignore):
	plt.figure(2)# this simple figure is created an removed in order to take the legend from it.  
	#plt.legend has issue where marker size in legend is propoertional to marker size in plot
	trial_1 = plt.scatter(1,1, c = 'r', marker = 'o')
	trial_2 = plt.scatter(1,1, c = 'b', marker = '^')
	trial_3 = plt.scatter(1,1, c = 'g', marker = 'v')
	trial_4 = plt.scatter(1,1, c = 'c', marker = '*')
	plt.close(2)
	fig = plt.figure()
	plt.title('ME-ICA Analysis, ' + r'$\kappa$' + ' vs ' + r'$\rho$', fontsize = 14)
	ACC = plt.scatter(accept[:,1], accept[:,2], c = 'r', marker = 'o', s = 50*accept[:,4]) 
	REJ = plt.scatter(reject[:,1], reject[:,2], c = 'b', marker = '^', s = 50*reject[:,4])
	MID = plt.scatter(middle[:,1], middle[:,2], c = 'g', marker = 'v', s = 50*middle[:,4])
	IGN = plt.scatter(ignore[:,1], ignore[:,2], c = 'c', marker = '*', s = 50*ignore[:,4])
	plt.legend((trial_1,trial_2,trial_3,trial_4),('Accepted','Rejected','Middle',
		'Ignore'), scatterpoints = 1, loc = 'upper right', markerscale = 2)
	plt.gca().xaxis.set_major_locator( MaxNLocator(nbins = 5))
	plt.tick_params(axis = 'x', which = 'both', top = 'off')
	plt.tick_params(axis = 'y', which = 'both', right = 'off')
	plt.xlabel(r'$\kappa$', fontsize = 15)
	plt.ylabel(r'$\rho$', fontsize = 15)
	plt.savefig('kappa_vs_rho')
	plt.close()

"""
plot kappa and rho vs their component number.
comptable title: path to the ctab file
"""
def kr_vs_component(comp_table_title):
	fig = plt.figure()
	components = np.loadtxt(str(comp_table_title))
	plt.figure()
	plt.title('ME-ICA Analysis, ' + r'$\kappa$' + ' and ' + r'$\rho$' + ' vs Component Rank', fontsize = 14)
	plt.ylabel(r'$\kappa$' ', ' + r'$\rho$' , fontsize = 15)
	plt.xlabel('Component Rank' , fontsize = 15)
	kappa = plt.plot(components[:,0], components[:,1])
	rho = plt.plot(components[:,0], components[:,2])
	plt.legend((r'$\kappa$', r'$\rho$'))
	plt.savefig('kappa_rho_vs_components')
	plt.close()	


def correlation(startdir,setname,nsmprage,threshold):
	cdict = {'red':  	  ((0.0, 0.0, 0.0),
						   (0.6, 0.0, 0.0),
						   (0.7, 0.6, 0.6),
						   (0.75, 0.9, 0.9),
						   (0.8, 1.0, 1.0),
		                   (1.0, 1.0, 1.0)),

		         'green': ((0.0, 0.0, 0.0),
		         		   (0.6, 1.0, 1.0),
		         		   (0.8, 1.0, 1.0),
		         		   (0.85, 0.9, 0.9),
		                   (0.95, 0.5, 0.5),
		                   (1.0, 0.0, 0.0)),

		         'blue':  ((0.0, 1.0, 1.0),
		         		   (0.5, 0.0, 0.0),
		                   (1.0, 0.0, 0.0))
		        }

	RGB = LinearSegmentedColormap('RGB', cdict)
	for i in ['pC', 'mPFC']:
		if i == 'pC':
			MNI = np.array([[0],[-53],[26]])

		if i == 'mPFC':
			MNI = np.array([[0],[52],[6]])
		ext = np.zeros(8)
		axial = ['x','y']
		beta = ni.load('%s/%s/TED/betas_hik_OC.nii' % (startdir, setname)).get_data()
		beta_header = ni.load('%s/%s/TED/betas_hik_OC.nii' % (startdir, setname)).get_header()
		anatomical = ni.load(nsmprage).get_data()
		anat_hdr = ni.load(nsmprage).get_header()
		b = beta_header['quatern_b'] 
		c = beta_header['quatern_c']
		d = beta_header['quatern_d']
		a = np.sqrt(1.0-(b*b+c*c+d*d))
		R = np.array([[a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*b*d+2*a*c], [2*b*c+2*a*d, a*a+c*c-b*b-d*d,
		 2*c*d-2*a*b], [2*b*d-2*a*c, 2*c*d+2*a*b, a*a+d*d-c*c-b*b]])
		if np.linalg.det(R) == -1:
			q_fac = -1
		else:
			q_fac = 1

		native_offset = np.dot(np.linalg.inv(R),MNI) - np.array([[beta_header['qoffset_x']] ,[beta_header['qoffset_y']], [beta_header['qoffset_z']]])
		native = np.divide(native_offset,np.array([[beta_header['pixdim'][1]], [beta_header['pixdim'][2]], [q_fac * beta_header['pixdim'][3]] ]))

		for j in [0,1]:
			ext[j+6] = anatomical.shape[j]*anat_hdr['srow_%s' % axial[j]][j] + anat_hdr['srow_%s' % axial[j]][3]#setting up talairach coordinates
			ext[j+4] = beta.shape[j]*beta_header['srow_%s' % axial[j]][j] + beta_header['srow_%s' % axial[j]][3]
			ext[j+2] = anat_hdr['srow_%s' % axial[j]][3]
			ext[j] = beta_header['srow_%s' % axial[j]][3]

		seed = beta[native[0,0],native[1,0],native[2,0],:] #Shape: 1 x n_c
		seed = np.asarray([seed])

		mask = beta.prod(-1) != 0
		beta_mask = beta[mask] #Shape: n_v x n_c
		norm = np.array([np.sqrt(inner1d(beta_mask,beta_mask))]).T
		norm_seed = np.sqrt(np.dot(seed,seed.T))
		beta_norm= np.divide(beta_mask,norm)
		R = np.array([inner1d(beta_norm,seed/norm_seed)]).T
		z = np.arctanh(R)*np.sqrt((beta.shape[3] - 3)) #z_value 1D array
		z[np.isnan(z)] = 5
		z[z > 5] = 5
		z[z < -5] = -5

		z_scores = np.zeros(mask.shape)
		z_threshold = np.zeros(mask.shape)
		z_scores[mask == True] = z[:,0]
		z_threshold[np.absolute(z_scores) > threshold] = 1

		itemindex = np.where(z_threshold == 1)
		for j in range(len(itemindex[0])):
			z_threshold = flood(z_threshold,itemindex[0][j],itemindex[1][j],itemindex[2][j])#floodfill
		z_scores[z_threshold == 0 ] = np.nan
		fig = plt.figure(figsize = (3.2*5,4))
		gs0 = gridspec.GridSpec(1,10)
		minimum = min(z_scores[np.isnan(z_scores) == False])
		maximum = max(z_scores[np.isnan(z_scores) == False])
		for j in np.arange(0,1,.1):#plot montage of accept component onto 
			ax1 = fig.add_subplot(gs0[0,int(j*10)])
			plt.imshow(anatomical[:,::-1,anatomical.shape[2]*j].T, cmap = 'Greys_r', extent = [ext[2],ext[6],ext[3],ext[7]])
			plot = plt.imshow(z_scores[:,::-1,z_scores.shape[2]*j].T, alpha = 0.8, cmap = RGB , extent = [ext[0],ext[4],ext[1],ext[5]],
			 interpolation = 'gaussian', vmin = minimum, vmax = maximum)
			plt.axis('off')
		gs0.tight_layout(fig, w_pad = -2.7, rect = [0,0,0.95,1])
		cbar = fig.add_axes([(gs0.right + ((gs0.right + 1)/2 - gs0.right)/2), gs0.bottom, .01, gs0.top - gs0.bottom])
		fig.colorbar(plot, cax = cbar)
		plt.ylabel('z score',fontsize = 12, rotation = 270)
		plt.savefig('%s_correlation' % i)
		plt.close()

