#!/usr/bin/env python
"""
Gutierrez, B.  Contains all functions for creating the figures
of the mecia report
"""

import matplotlib as mpl
mpl.use('Agg')

from matplotlib.colors       import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import ImageGrid
from multiprocessing         import Pool
from numpy.core.umath_tests  import inner1d
from matplotlib.ticker       import MaxNLocator
from matplotlib              import gridspec
from matplotlib              import pylab

import matplotlib.pyplot as plt
import seaborn as sns
import nibabel as ni
import numpy as np
import subprocess
import sys
import os


cdict = {'red':   ((0.0,  0.0,  0.0),
                   (0.2,  0.6,  0.6),
                   (0.35, 0.9,  0.9),
                   (0.5,  1.0,  1.0),
                   (1.0,  1.0,  1.0)),

         'green': ((0.0,  1.0,  1.0),
                   (0.5,  1.0,  1.0),
                   (0.65, 0.9,  0.9),
                   (0.9,  0.5,  0.5),
                   (1.0,  0.0,  0.0)),

         'blue':  ((0.0,  0.0,  0.0),
                   (1.0,  0.0,  0.0))}

GYR = LinearSegmentedColormap('GYR', cdict)

def FFT(Denoised_components_ts, TR, N, outputDir):
    t      = Denoised_components_ts[:,N]
    FFT    = abs(np.fft.fft(t))
    FFT    = np.fft.fftshift(FFT)

    freq      = np.fft.fftfreq(t.size,float(TR))
    freq      = np.fft.fftshift(freq)
    freq_axis = freq[np.where(freq == 0)[0][0]:]

    fig  = plt.figure(figsize= (8,4))

    plt.plot(freq_axis,FFT[-(freq_axis.shape)[0]:])
    plt.title('FFT of the Time Series', fontsize = 15)
    plt.xlabel('Frequency(Hz)' , fontsize = 15)
    plt.ylabel('Amplitude' , fontsize = 15)
    plt.savefig('%s/Report_Figures/FFT_Component_' % outputDir + str(N).zfill(len(str(Denoised_components_ts.shape[1])))
        , bbox_inches='tight')
    plt.close()

"""
Primary usage is for finding upper and lower bounds for non-zero 2d slices so as to 
display only useful information when displaying montages of three dimmensional images.  
Also returns an image with 2D zero-filled slices removed, but this is not always the more useful option.
image: 3D array
axis: String Dimmension to look for 2D zero-filled slices in, i.e, 'x','y','z'.
"""
def ImageDimBounds(image, axis):
    img_mask = np.zeros(image.shape)
    image[np.isnan(image)] = 0
    img_mask[image != 0] = 1

    lower_bound = 0
    upper_bound = 0

    if axis == 'x':
        first_dim  = 1
        second_dim = 1
    if axis == 'z':
        first_dim  = 0
        second_dim = 0
    if axis == 'y':
        first_dim  = 0
        second_dim = 1

    img_bounds = np.sum(img_mask, axis = first_dim)
    img_bounds = np.sum(img_bounds, axis = second_dim)

    i = img_bounds[0]
    j = img_bounds[img_bounds.shape[0]-1]

    while i == 0:
        lower_bound += 1
        if lower_bound == img_bounds.shape[0]:
            lower_bound = 0
            break
        else:
            i = img_bounds[lower_bound]

    while j == 0:
        upper_bound += 1
        if upper_bound == img_bounds.shape[0]:
            upper_bound = 0
            img_bounds = np.ones(img_bounds.shape)
            break
        else:
            j = img_bounds[img_bounds.shape[0]-1-upper_bound]

    if axis == 'x':
        cropped_img = image[img_bounds > 0,:,:]
    if axis == 'y':
        cropped_img = image[:,img_bounds > 0,:]
    if axis == 'z':
        cropped_img = image[:,:,img_bounds > 0]

    return (cropped_img,lower_bound,img_bounds.shape[0]-1-upper_bound)

def montage_control(TED, outputDir, Denoised_components, Denoised_components_ts, Ncpu):
    
    generic_image = np.ones(shape=Denoised_components[:,:,:,0].shape)
    fig           = plt.figure(figsize=(12,4))
    grid1         = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,axes_pad=0.0,share_all=True,cbar_mode=None)
    for j in range(10):
        grid1[j].imshow(generic_image[:,:,0].T, cmap = 'Greys_r', vmin = 0, vmax = 2)
        grid1[j].axes.get_xaxis().set_ticks([])
        grid1[j].axes.get_yaxis().set_ticks([])
    plt.savefig('%s/Report_Figures/Axial_GS_Component_XX' % (outputDir) , bbox_inches='tight', dpi=150)
    plt.close()

    fig = plt.figure(figsize=(12,4))
    grid1 = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,axes_pad=0.0,share_all=True,cbar_mode=None)
    for j in range(10):
        grid1[j].imshow(generic_image[0,:,:].T, cmap = 'Greys_r', vmin = 0, vmax = 2)
        grid1[j].axes.get_xaxis().set_ticks([])
        grid1[j].axes.get_yaxis().set_ticks([])
    plt.savefig('%s/Report_Figures/Sagittal_GS_Component_XX' % (outputDir), bbox_inches='tight', dpi=150)
    plt.close()
    
    fig = plt.figure(figsize=(12,4))
    grid1 = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,axes_pad=0.0,share_all=True,cbar_mode=None)
    for j in range(10):
        grid1[j].imshow(generic_image[:,0,:].T, cmap = 'Greys_r', vmin = 0, vmax = 2)
        grid1[j].axes.get_xaxis().set_ticks([])
        grid1[j].axes.get_yaxis().set_ticks([])
    plt.savefig('%s/Report_Figures/Coronal_GS_Component_XX' % (outputDir), bbox_inches='tight', dpi=150)
    plt.close()

    pool = Pool(processes=Ncpu)
    pool.map(gs_montage, [{'Denoised_components':Denoised_components,'N':str(i).zfill(len(str(Denoised_components.shape[3]))),
     'outputDir':outputDir,'Denoised_components_ts': Denoised_components_ts} for i in range(Denoised_components.shape[3])])

    p = subprocess.Popen(['3dinfo','-tr', '%s/betas_OC.nii' % TED], stdout = subprocess.PIPE, stderr = subprocess.PIPE)# retrieve TR
    TR = p.communicate()[0][:-2]

    for N in range(Denoised_components.shape[3]):
        FFT(Denoised_components_ts, TR, N, outputDir)

"""
Creates a montage of greyscale 10 images of axial, Sagittal and coronal views of meica components along with a time series of the component.
accepted components also get overlayed onto the anatomical image and floodfill is performed on statiscially significant voxels.
maps: array of 5 elements [anat dataset, mefl dataset, feats dataset, anat header, overlay header]
overlay: mefl dataset
Denoised_components_ts: string meica_mix.1D path
i: component number
N: component number as string of standardized length
"""

def gs_montage(item):
    Denoised_components    = item['Denoised_components']
    N                      = item['N']
    outputDir              = item['outputDir']
    Denoised_components_ts = item['Denoised_components_ts']
    i = int(N)

    contrast           = Denoised_components[Denoised_components[:,:,:,i] != 0,i]#fix contrast
    tmp,lower,upper     = ImageDimBounds(Denoised_components[:,:,:,i],'z')#find range of indicies for each axis that contain non-zero values for displaying to the user.
    ax_montage_spacing  = np.linspace(lower,upper,10)/Denoised_components.shape[2]
    tmp,lower,upper     = ImageDimBounds(Denoised_components[:,:,:,i],'x')
    sag_montage_spacing = np.linspace(lower,upper,10)/Denoised_components.shape[0]
    tmp,lower,upper     = ImageDimBounds(Denoised_components[:,:,:,i],'y')
    cor_montage_spacing = np.linspace(lower,upper,10)/Denoised_components.shape[1]
    maximum             = np.percentile(contrast, 98)
    minimum             = np.percentile(contrast, 2)

    #plot greyscale component montage
    fig = plt.figure(figsize = (12,4))
    grid1 = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,axes_pad=0.0,share_all=True,cbar_mode=None)
    for j in range(10):
        grid1[j].imshow(Denoised_components[:,:,int((Denoised_components.shape[2]-1)*ax_montage_spacing[j]),i].T, cmap = 'Greys_r',
         vmin = minimum, vmax = maximum)
        grid1[j].axes.get_xaxis().set_ticks([])
        grid1[j].axes.get_yaxis().set_ticks([])
    plt.savefig('%s/Report_Figures/Axial_GS_Component_' % outputDir + N, bbox_inches='tight', dpi=150)
    plt.close()
    
    fig = plt.figure(figsize = (12,4))
    grid2 = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,axes_pad=0.0,share_all=True,cbar_mode=None)
    for j in range(10):
        grid2[j].imshow(Denoised_components[int((Denoised_components.shape[0]-1)*sag_montage_spacing[j]),:,::-1,i].T, cmap = 'Greys_r',
         vmin = minimum, vmax = maximum)
        grid2[j].axes.get_xaxis().set_ticks([])
        grid2[j].axes.get_yaxis().set_ticks([])
    plt.savefig('%s/Report_Figures/Sagittal_GS_Component_' % outputDir + N, bbox_inches='tight', dpi=150)
    plt.close()
    
    fig = plt.figure(figsize = (12,4))
    grid3 = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,axes_pad=0.0,share_all=True,cbar_mode=None)
    for j in range(10):
        grid3[j].imshow(Denoised_components[:,int((Denoised_components.shape[1]-1)*cor_montage_spacing[j]),::-1,i].T, cmap = 'Greys_r',
         vmin = minimum, vmax = maximum)
        grid3[j].axes.get_xaxis().set_ticks([])
        grid3[j].axes.get_yaxis().set_ticks([])
    plt.savefig('%s/Report_Figures/Coronal_GS_Component_' % outputDir + N, bbox_inches='tight', dpi=150)
    plt.close()

    fig = plt.figure(figsize= (8,4))
    plt.plot(np.arange(Denoised_components_ts.shape[0]),Denoised_components_ts[:,i])
    plt.title('Time Series of the Component', fontsize = 15)
    plt.xlabel('Time (TR)', fontsize = 15)
    plt.ylabel('Arbitrary BOLD units', fontsize = 15)
    plt.xlim([0,Denoised_components_ts.shape[0]-1])
    plt.savefig('%s/Report_Figures/TimeSeries_' % outputDir + N, bbox_inches='tight')
    plt.close()
    
    print("++ INFO [Figures]: Component [%s] figure completed." % N)

"""
Makes TSNR figures of medn, tsoc, and medn/tsoc datasets
tsoc: string path to tsoc dataset
medn: string path to medn dataset
"""
def tsnr(tsoc_data, medn_data, outputDir):
    medn_data[medn_data == 0] = np.nan
    tsoc_data[tsoc_data == 0] = np.nan

    medn_tsnr = medn_data.mean(-1)/medn_data.std(-1)
    tsoc_tsnr = tsoc_data.mean(-1)/tsoc_data.std(-1)
    frac = medn_tsnr/tsoc_tsnr

    medn_tsnr,lower,upper = ImageDimBounds(image = medn_tsnr, axis = 'z')#remove z slices without nonzero elements
    tsoc_tsnr,lower,upper = ImageDimBounds(image = tsoc_tsnr, axis = 'z')
    frac_tsnr,lower,upper = ImageDimBounds(image = frac, axis = 'z')#remove z slices without nonzero elements
    medn_tsnr[medn_tsnr == 0] = np.nan
    tsoc_tsnr[tsoc_tsnr == 0] = np.nan
    frac_tsnr[frac_tsnr == 0] = np.nan
    background = np.zeros((medn_tsnr[:,:,0].T).shape)

    fig = plt.figure(figsize = (12,2))#plot montage of medn TSNR
    grid = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,cbar_pad=0.05,axes_pad=0.0,
        share_all=True,cbar_mode="single")

    tsnr_mask = medn_tsnr[np.isnan(medn_tsnr) == False]
    maximum = np.percentile(tsnr_mask,95)
    minimum = np.percentile(tsnr_mask,5)
    SaveTSNR=np.zeros((medn_tsnr.shape[1],medn_tsnr.shape[0],10))
    for i in range(10):
        grid[i].imshow(background, cmap = 'Greys_r')
        plot = grid[i].imshow(medn_tsnr[:,:,int(i*.1*medn_tsnr.shape[2])].T, vmin = minimum, vmax = maximum, cmap = GYR)
        SaveTSNR[:,:,i]= medn_tsnr[:,:,int(i*.1*medn_tsnr.shape[2])].T
        grid[i].axes.get_xaxis().set_ticks([])
        grid[i].axes.get_yaxis().set_ticks([])
    cb1 = grid.cbar_axes[0].colorbar(plot)
    grid.cbar_axes[0].tick_params(labelsize=9)
    grid.cbar_axes[0].get_yaxis().labelpad = 15
    cb1.set_label_text('TSNR', fontsize=9, rotation=270)
    plt.savefig('%s/Report_Figures/medn_tsnr' % outputDir, bbox_inches='tight', dpi=150)
    np.save('%s/axialized_nifti/medn_tsnr' % (outputDir),SaveTSNR)
    np.savetxt('%s/axialized_nifti/tsnr_thresholds.txt' % (outputDir),np.array([minimum,maximum]))
    plt.close()

    fig = plt.figure(figsize = (12,2))#plot montage of medn TSNR
    grid = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,cbar_pad=0.05,axes_pad=0.0,
        share_all=True,cbar_mode="single")
    SaveTSNR=np.zeros((tsoc_tsnr.shape[1],tsoc_tsnr.shape[0],10))
    for i in range(10):
        grid[i].imshow(background, cmap = 'Greys_r')
        grid[i].imshow(tsoc_tsnr[:,:,int(i*.1*tsoc_tsnr.shape[2])].T, vmin = minimum, vmax = maximum, cmap = GYR)
        SaveTSNR[:,:,i]= tsoc_tsnr[:,:,int(i*.1*tsoc_tsnr.shape[2])].T
        grid[i].axes.get_xaxis().set_ticks([])
        grid[i].axes.get_yaxis().set_ticks([])
    cb2 = grid.cbar_axes[0].colorbar(plot)
    grid.cbar_axes[0].tick_params(labelsize=9)
    grid.cbar_axes[0].get_yaxis().labelpad = 15
    cb2.set_label_text('TSNR', fontsize=9, rotation=270)
    plt.savefig('%s/Report_Figures/tsoc_tsnr' % outputDir, bbox_inches='tight', dpi=150)
    np.save('%s/axialized_nifti/tsoc_tsnr' % (outputDir),SaveTSNR)
    plt.close()

    fig = plt.figure(figsize = (12,2))#plot montage of medn TSNR
    grid = ImageGrid(fig, 111 , nrows_ncols=(1,10),cbar_location='right',add_all=True,cbar_pad=0.05,axes_pad=0.0,
        share_all=True,cbar_mode="single")
    tsnr_mask = frac_tsnr[np.isnan(frac_tsnr) == False]
    maximum = np.percentile(tsnr_mask,95)
    minimum = np.percentile(tsnr_mask,5)
    SaveTSNR=np.zeros((frac_tsnr.shape[1],frac_tsnr.shape[0],10))
    for i in range(10):
        grid[i].imshow(background, cmap = 'Greys_r')
        plot = grid[i].imshow(frac_tsnr[:,:,int(i*.1*frac_tsnr.shape[2])].T, vmin = minimum, vmax = maximum, cmap = GYR)
        SaveTSNR[:,:,i]= frac_tsnr[:,:,int(i*.1*frac_tsnr.shape[2])].T
        grid[i].axes.get_xaxis().set_ticks([])
        grid[i].axes.get_yaxis().set_ticks([])
    cb3 = grid.cbar_axes[0].colorbar(plot)
    grid.cbar_axes[0].tick_params(labelsize=9)
    grid.cbar_axes[0].get_yaxis().labelpad = 15
    cb3.set_label_text('TSNR ratio', fontsize=9, rotation=270)
    plt.savefig('%s/Report_Figures/tsnr_ratio' % outputDir, bbox_inches='tight', dpi=150)
    np.save('%s/axialized_nifti/ratio_tsnr' % (outputDir),SaveTSNR)
    np.savetxt('%s/axialized_nifti/ratio_tsnr_thresholds.txt' % (outputDir),np.array([minimum,maximum]))
    plt.close()

    #plot histogram of the TSNR of medn
    medn_mask = medn_tsnr[np.isnan(medn_tsnr) == False]
    tsoc_mask = tsoc_tsnr[np.isnan(tsoc_tsnr) == False]
    frac_mask = frac_tsnr[np.isnan(frac_tsnr) == False]
    fig = plt.figure()
    plt.hist(medn_mask, bins = 100, range = [np.percentile(medn_mask,5),np.percentile(medn_mask,95)])
    plt.title('TSNR medn', fontsize = 15)
    plt.xlabel('TSNR', fontsize = 15)
    plt.ylabel('Frequency', fontsize = 15)
    plt.savefig('%s/Report_Figures/medn_tsnr_hist' % outputDir)
    plt.close()

    fig = plt.figure()
    plt.hist(tsoc_mask, bins = 100, range = [np.percentile(tsoc_mask,5),np.percentile(tsoc_mask,95)])
    plt.title('TSNR tsoc', fontsize = 15)
    plt.xlabel('TSNR', fontsize = 15)
    plt.ylabel('Frequency', fontsize = 15)
    plt.savefig('%s/Report_Figures/tsoc_tsnr_hist' % outputDir)
    plt.close()

    #plot histogram of the TSNR ratio of medn/tsnr
    fig = plt.figure()
    plt.hist(frac_mask, bins = 100, range = [np.percentile(frac_mask,5),np.percentile(frac_mask,95)])
    plt.title('TSNR medn / TSNR tsoc', fontsize = 15)
    plt.xlabel('TSNR ratio', fontsize = 15)
    plt.ylabel('Frequency', fontsize = 15)
    plt.savefig('%s/Report_Figures/tsnr_ratio_hist' % outputDir)
    plt.close()
    print("++ INFO [Figures]: TSNR figures created.")
    return("Median meica denoised TSNR:  %s\nMedian optimally combined TSNR:   %s\nMedian denoised over optimally combined TSNR ratio:   %s"
     % (np.percentile(medn_mask,50),np.percentile(tsoc_mask,50),np.percentile(frac_mask,50)))

def kappa_vs_rho_plot(accepted,rejected,middle_kappa,ignored,ctab,outputDir,TED_dir):
    plt.figure(2)# this simple figure is created and removed in order to take the legend from it.  
    #plt.legend has issue where marker size in legend is propoertional to marker size in plot
    trial_1 = plt.scatter(1,1, c = 'b', marker = 'o')
    trial_2 = plt.scatter(1,1, c = 'r', marker = 'v')
    trial_3 = plt.scatter(1,1, c = 'g', marker = '^')
    trial_4 = plt.scatter(1,1, c = 'c', marker = '*')
    plt.close(2)
    fig = plt.figure()
    plt.title('ME-ICA Analysis, ' + r'$\kappa$' + ' vs ' + r'$\rho$', fontsize = 14)
    ACC = plt.scatter(ctab[accepted,1], ctab[accepted,2], c = 'b', marker = 'o', s = 50 * ctab[accepted,4]) 
    REJ = plt.scatter(ctab[rejected,1], ctab[rejected,2], c = 'r', marker = 'v', s = 50 * ctab[rejected,4])
    MID = plt.scatter(ctab[middle_kappa,1], ctab[middle_kappa,2], c = 'g', marker = '^', s = 50 * ctab[middle_kappa,4])
    IGN = plt.scatter(ctab[ignored,1], ctab[ignored,2], c = 'c', marker = '*', s = 50 * ctab[ignored,4])
    plt.legend((trial_1, trial_2, trial_3, trial_4),('accepted','rejected','middle_kappa',
        'ignored'), scatterpoints = 1, loc = 'upper right', markerscale = 2)
    plt.gca().xaxis.set_major_locator(MaxNLocator(nbins = 5))
    plt.tick_params(axis = 'x', which = 'both', top = 'off')
    plt.tick_params(axis = 'y', which = 'both', right = 'off')
    plt.xlabel(r'$\kappa$', fontsize = 15)
    plt.ylabel(r'$\rho$', fontsize = 15)
    plt.savefig('%s/Report_Figures/kappa_vs_rho.png'  % outputDir)
    plt.close()
    txt_file = open(TED_dir + '/' + 'comp_table.txt')
    lines = txt_file.readlines()
    for i in range(len(lines)):
        if '#Dataset variance explained by ICA' in lines[i]:
            index = i
    if len(rejected) == 0:
        reject_high_k = ''
    else:
        reject_high_k = rejected[0]
    if len(middle_kappa) == 0:
        middle_high_k = ''
    else:
        middle_high_k = middle_kappa[0]
    print("++ INFO [Figures]: kappa vs rho figures created.")
    return("Number of accepted components:   %s\nNumber of rejected components:   %s\nNumber of middle kappa components:   %s\nNumber of ignored components:   %s\naccepted variance:   %s" % 
        (len(accepted),len(rejected),len(middle_kappa),len(ignored),sum(ctab[accepted,4])) + 
        "\nRejected variance:   %s\nMiddle Kappa variance:   %s\nIgnored variance:   %s\nVariance explained by ICA:   %s\nLargest variance accepted component:   %s\nHighest kappa rejected component:   %s\nHighest kappa middle kappa component:   %s"
         % (sum(ctab[rejected,4]),sum(ctab[middle_kappa,4]),sum(ctab[ignored,4]),lines[index][41:-1],accepted[np.argmax(ctab[accepted,4])],reject_high_k,middle_high_k))
"""
plot kappa and rho vs their component number.
comptable title: path to the ctab file
"""
def kr_vs_component(ctab,outputDir):
    fig = plt.figure()
    ctab
    plt.title('ME-ICA Analysis, ' + r'$\kappa$' + ' and ' + r'$\rho$' + ' vs Component Rank', fontsize = 14)
    plt.ylabel(r'$\kappa$' ', ' + r'$\rho$', fontsize = 15)
    plt.xlabel('Component Rank', fontsize = 15)
    kappa = plt.plot(ctab[:,0], ctab[:,1])
    rho = plt.plot(ctab[:,0], ctab[:,2])
    plt.legend((r'$\kappa$', r'$\rho$'))
    plt.savefig('%s/Report_Figures/kappa_rho_vs_components' % outputDir)
    plt.close()
    print("++ INFO [Figures]: kappa and rho vs component number figures created.")

def motion(outputDir, motion_file):
    motion = np.loadtxt('%s' % motion_file)
    grad = np.zeros(motion.shape)
    for Nc in range(motion.shape[1]):
        grad[:,Nc] = np.gradient(motion[:,Nc])
    deriv = np.linalg.norm(grad,ord=2,axis = 1)

    plt.figure()
    fig, ax = plt.subplots(nrows=6, ncols=1, sharex=True, sharey = True)
    ax[0].set_title('Subject Motion', fontsize = 14)
    plt.xlabel('Time (TR)', fontsize = 14)
    plt.ylabel('distance (mm)', fontsize = 14)
    label = ['x', 'y', 'z', 'pitch', 'roll', 'yaw']
    for i in range(6):
        ax[i].plot(np.arange(motion.shape[0]),motion[:,5-i])
        ax[i].plot(np.arange(motion.shape[0]),np.zeros(shape = (motion.shape[0],1)), ls = '--', c = 'k')
        ax[i].set_ylabel(label[i])

    ax[0].set_yticks([ax[0].get_ylim()[0], 0, ax[0].get_ylim()[1]])
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.savefig('%s/Report_Figures/motion_plot' % outputDir)
    plt.close()

    plt.figure()
    plt.plot(np.arange(motion.shape[0]),deriv[:])
    plt.xlabel('Time (TR)', fontsize = 14)
    plt.ylabel('rate (mm/TR)', fontsize = 14)
    plt.title('rate of motion')
    plt.savefig('%s/Report_Figures/motion_rate' % outputDir)
    plt.close()
    print("++ INFO [Figures]: Motion figures created.")
    itemindex = np.where(np.absolute(motion)==np.max(np.absolute(motion)))
    return("Max head displacement in any one dirrection:   %s\nTime of max head displacement (TR):   %s\nMax rate of head motion:   %s\nTime of max head motion rate (TR):   %s" % (np.max(np.absolute(motion)),itemindex[0][0],
        np.max(np.absolute(deriv)),np.argmax(np.absolute(deriv))))

