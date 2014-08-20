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
import meica_figures
import numpy as np
import rst_files
import os


startdir = '/Users/gutierrezbe/Documents/NIFTI/WE'
nsmprage = '/Users/gutierrezbe/Documents/NIFTI/WE/20140528_094651MPRAGE1isoG2s013a1001_ns_at.nii.gz'
outprefix = 'prefix'
setname = 'meica.WE_TE1234label'
threshold = 1.96
MNI = True

os.chdir(startdir)
# os.system('rm -rf %s/sphinx' % startdir)
# os.system('mkdir %s/sphinx' % startdir)
# print '!!!!hit enter for sphinx-quickstart prompts unless otherwise stated below:!!!!'
# print '> Root path for the documentation [.]: "sphinx"'
# print '> Project name: (users choice)'
# print '> Author name(s): (users choice)'
# print '> Project version: (users choice)'
# print '> Project release [m]: (users choice)'
# print '> pngmath: include math, rendered as PNG images (y/n) [n]: "y"'
# os.system('sphinx-quickstart')

components = meica_figures.file_parse('%s/%s_ctab.txt' % (startdir,outprefix))
maps = meica_figures.collect_data(anatomical = nsmprage, overlay = '%s/%s_mefl.nii.gz' % (startdir,outprefix), threshold_map = '%s/%s/TED/feats_OC2.nii' % (startdir,setname))
accept, reject, middle, ignore = meica_figures.split_components('%s/%s_ctab.txt' % (startdir,outprefix), components)

os.system('rm -rf %s/png_dump' % startdir)
os.mkdir('%s/png_dump' % startdir)
os.chdir('%s/png_dump' % startdir)

#make figures
print('++ making figures')
print('++ ignore tight_layout : falling back to Agg renderer message')
meica_figures.kr_vs_component('%s/%s_ctab.txt' % (startdir,outprefix))
meica_figures.kappa_vs_rho_plot(accept, reject, middle, ignore)
meica_figures.tsnr('%s/%s_tsoc.nii.gz' % (startdir,outprefix),'%s/%s_medn.nii.gz' % (startdir,outprefix))
print('++ this set of figures may take a while')
meica_figures.montage(maps, accept, threshold, 0.8, '%s/%s/TED/meica_mix.1D' % (startdir,setname), Axial = 0, Sagital = 0, Coronal = 1)
if nsmprage != '':
	meica_figures.coreg(setname,nsmprage[len(startdir)+1:])
	if MNI == True:
		meica_figures.correlation(startdir, setname, nsmprage, threshold = 2.1)

os.chdir('%s/sphinx' % startdir)

#make .rst files for sphinx to use
print('++ occupying sphinx directory with .rst files')
os.system('rm -f index.rst')#remove sphinx version of index
rst_files.diagnostics_rst(nsmprage)
rst_files.index_rst()
rst_files.intro_rst()
rst_files.analysis_rst(accept, reject, middle, ignore, threshold, '%s/%s_ctab.txt' % (startdir,outprefix))
if nsmprage != '' and MNI == True:
	rst_files.correlation_rst()



#run sphinx
os.system('make html')
os.system('make latex')
