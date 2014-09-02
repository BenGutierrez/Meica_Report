#!/usr/bin/env python
"""
Gutierrez, B.  Generates three png files for each accepted components along with time series of these components
in each png file.

"""

import commands
import sys
import os

#Run dependency check
def dep_check():
	print '++ Checking system for dependencies...'
	fails = 0
	matplotlib_installed = 0
	nibabel_installed = 0
	numpy_installed = 0
	parse_installed = 0

	try:
		import numpy
		numpy_installed = 1
	except:
		print "*+ Can't import Numpy! Please check Numpy installation for this Python version."
		fails += 1

	try:
		import matplotlib
		matplotlib_installed = 1
	except:
		print "*+ Can't import Matplotlib! Please check Matplotlib installation."
		fails += 1

	try:
		import nibabel
		nibabel_installed = 1
	except:
		print "*+ Can't import Nibabel! Please check Nibabel installation."
		fails += 1

	try:
		import parse
		parse_installed = 1
	except:
		print "*+ Can't import Parse! Please check Parse installation."
		fails += 1

	if numpy_installed:
		print " + Numpy version: %s" % (numpy.__version__)
		if float('.'.join(numpy.__version__.split('.')[0:2]))<1.5:
			fails += 1
			print "*+ Numpy version is too old! Please upgrade to Numpy >=1.5.x!"
		import numpy.__config__ as nc
		if nc.blas_opt_info == {}:
			fails += 1
			print "*+ Numpy is not linked to BLAS! Please check Numpy installation."

	afnicheck = commands.getstatusoutput("3dinfo")
	if afnicheck[0]!=0:
		print "*+ Can't run AFNI binaries. Make sure AFNI is on the path!"
		fails += 1
	elif not afnicheck[1].__contains__('Alternate Alternative Usage'):
		print "*+ This seems like an old version of AFNI. Please upgrade to latest version of AFNI."
		fails += 1
	if fails == 0:
		print " + Dependencies OK."
	else:
		print "*+ EXITING. Please see error messages."
		sys.exit()
dep_check()
print 'all done'

import meica_figures
import sphinx_files
import rst_files

__version__ = "v2.5 beta8"
startdir = '/Users/gutierrezbe/Documents/NIFTI/WE'
nsmprage = '/Users/gutierrezbe/Documents/NIFTI/WE/20140528_094651MPRAGE1isoG2s013a1001_ns_at.nii.gz'
outprefix = 'prefix'
setname = 'meica.WE_TE1234label'
MNI = True
min_component_number = 20
min_variance_explained = 90
montage_threshold = 1.96 #threshold for accepted components onto anatomcial
corr_threshold = 1.96 #threshold for correlation for ROI



os.chdir(startdir)
os.system('rm -rf %s/sphinx' % startdir)
os.system('mkdir %s/sphinx' % startdir)
os.system('mkdir %s/sphinx/_build' % startdir)
os.system('mkdir %s/sphinx/_static' % startdir)
os.system('mkdir %s/sphinx/_templates' % startdir)



components = meica_figures.file_parse('%s/%s_ctab.txt' % (startdir,outprefix))
maps = meica_figures.collect_data(anatomical = nsmprage, overlay = '%s/%s_mefl.nii.gz' % (startdir,outprefix),
	threshold_map = '%s/%s/TED/feats_OC2.nii' % (startdir,setname))
accept, reject, middle, ignore = meica_figures.split_components('%s/%s_ctab.txt' % (startdir,outprefix), components)

os.system('rm -rf %s/png_dump' % startdir)
os.mkdir('%s/png_dump' % startdir)
os.system('mv warning.png %s/png_dump' % startdir)
os.chdir('%s/png_dump' % startdir)


#make figures
print('++ making figures')
meica_figures.kr_vs_component('%s/%s_ctab.txt' % (startdir, outprefix))
meica_figures.kappa_vs_rho_plot(accept, reject, middle, ignore)
meica_figures.tsnr('%s/%s_tsoc.nii.gz' % (startdir,outprefix),'%s/%s_medn.nii.gz' % (startdir,outprefix))
print('++ this set of figures may take a while')
meica_figures.montage(maps, accept, montage_threshold, 0.8, startdir, setname, outprefix, Axial = 0, Sagital = 0, Coronal = 1)
if nsmprage != '':
	meica_figures.coreg(setname,nsmprage[len(startdir)+1:])
	if MNI == True:
		meica_figures.correlation(startdir, setname, nsmprage, threshold = corr_threshold)

os.chdir('%s/sphinx' % startdir)

#set up sphinx documentation
sphinx_files.conf(__version__)
sphinx_files.make_bat()
sphinx_files.make_file()

#make .rst files for sphinx to use
print('++ occupying sphinx directory with .rst files')
os.system('rm -f index.rst')#remove sphinx version of .rst
rst_files.diagnostics_rst(nsmprage)
rst_files.index_rst()
rst_files.intro_rst()
rst_files.analysis_rst(accept, reject, middle, ignore, nsmprage, montage_threshold, '%s/%s_ctab.txt' % (startdir,outprefix),
	min_component_number, min_variance_explained)
if nsmprage != '' and MNI == True:
	rst_files.correlation_rst()



#run sphinx build
os.system('make html')
os.system('make latex')
