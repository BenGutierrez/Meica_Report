#!/usr/bin/env python
"""
Gutierrez, B.  Generates Meica report form.

"""
__version__ = "v2.5 beta8"

import commands
import subprocess
import argparse
import sys
import ast
import os
import re

def MNI_check(MNI, User_ROI, ROI_def, ROI_att, ROI_ref):
	if not MNI:
		if User_ROI or ROI_att or ROI_def or ROI_ref:
			print '*+ MNI needs to be specified if options User_ROI, ROI_att, ROI_def, or ROI_ref are specified'
			sys.exit()
		corr = False
	else:
		if not (User_ROI or ROI_att or ROI_def or ROI_ref):
			print '++ !!MNI specified but not User_ROI, ROI_att, ROI_def, or ROI_ref.  No seed based correlation will be computed!!'
			corr = False
		else:
			corr = True
	return corr

def seed_split(ROI):
	ROI = re.sub('[(]','[',ROI)
	ROI = re.sub('[)]',']',ROI)
	ROI = ast.literal_eval(ROI)
	List = []
	for i in range(len(ROI)):
		List.append(ROI[i])
	if len(ROI) == 1:
		List = [List]
	return List

def file_check(anat, startdir, prefix, setname, MNI):
	fails = 0
	if not os.path.isfile(anat) and anat !='': 
		print '*+ Can\'t find the specified anantomical'
		fails += 1
	if not os.path.isfile('%s/%s_ctab.txt' % (startdir,prefix)):
		print '*+ Can\'t find "%s/%s_ctab.txt" check starting directory and file\'s existance.\n' % (startdir,prefix)
		fails += 1
	if not os.path.isfile('%s/%s_mefl.nii.gz' % (startdir,prefix)):
		print '*+ Can\'t find "%s/%s_mefl.nii.gz" check starting directory and file\'s existance.\n' % (startdir,prefix)
		fails += 1
	if not os.path.isfile('%s/%s_tsoc.nii.gz' % (startdir,prefix)):
		print '*+ Can\'t find "%s/%s_tsoc.nii.gz" check directory and file\'s existance.\n' % (startdir,prefix)
		fails += 1
	if not os.path.isfile('%s/%s_medn.nii.gz' % (startdir,prefix)):
		print '*+ Can\'t find "%s/%s_medn.nii.gz" check directory and file\'s existance.\n' % (startdir,prefix)
		fails += 1
	if not os.path.isfile('%s/%s/ocv_uni_vr.nii.gz' % (startdir,setname)):
		print '*+ Can\'t find "%s/%s/ocv_uni_vr.nii.gz check directory and file\'s existance.\n' % (startdir,setname)
		fails += 1
	if not os.path.isfile('%s/%s/TED/meica_mix.1D' % (startdir,setname)):
		print '*+ Can\'t find "%s/%s/TED/meica_mix.1D" check directory and file\'s existance.\n' % (startdir,setname)
		fails += 1
	if MNI == True:
		if not os.path.isfile('%s/%s/eBvrmask.nii.gz' % (startdir,setname)):
			print '*+ Can\'t find "%s/%s/eBvrmask.nii.gz check directory and file\'s existance.\n' % (startdir,setname)
			fails += 1
		if not os.path.isfile('%s/%s/TED/feats_OC2.nii' % (startdir,setname)):
			print '*+ Can\'t find "%s/%s/TED/feats_OC2.nii" check directory and file\'s existance.\n' % (startdir,setname)
			fails += 1
	if fails != 0:
		print "*+ EXITING. Please see error messages."
		sys.exit()

#Run dependency check
def dep_check():
	print '++ Checking system for dependencies...'
	fails = 0
	numpy_installed = 0
	sphinx_installed = 0
	try:
		import numpy
		numpy_installed = 1
	except:
		print "*+ Can't import Numpy! Please check Numpy installation for this Python version."
		fails += 1
	try:
		import matplotlib
	except:
		print "*+ Can't import Matplotlib! Please check Matplotlib installation."
		fails += 1
	try:
		import nibabel
	except:
		print "*+ Can't import Nibabel! Please check Nibabel installation."
		fails += 1
	try:
		import parse
		sphinx_installed = 1
	except:
		print "*+ Can't import Parse! Please check Parse installation."
		fails += 1
	try:
		import sphinx
	except:	
		print "*+ Can't import Sphinx! Please check Sphinx installation."
		fails += 1
	if sphinx_installed:
		print " + Sphinx version: %s" % (sphinx.__version__)
		if float('.'.join(sphinx.__version__.split('.')[0:2]))<1.2:
			fails += 1
			print "*+ Sphinx version is too old! Please upgrade to Sphinx >=1.2.x!"
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

parser = argparse.ArgumentParser()
parser.add_argument('-mt','--montage_threshold', dest = 'montage_threshold', help = 'normalized z-score value for thresholding' + 
	' the accepted component w/ anatomical activation map', type = float, default = 1.96)
parser.add_argument('-ct', '--correlation_threshold', dest = 'corr_threshold', help = 'normalized z-score value for thresholding' + 
	' the seed voxel map', type = float, default = 1.96)
parser.add_argument('-min_c', '--min_comp', dest = 'min_component_number', help = 'minimum total component number before warning raised in report', type = int, default = 20)
parser.add_argument('-min_v', '--min_var', dest = 'min_variance_explained', help = 'minimum variance explained before warning raised in report', type = int, default = 85)
parser.add_argument('-ax', '--axial' , dest = 'Axial', help = 'option to add axial image from montages', action = 'store_true')
parser.add_argument('-sag', '--sag' , dest = 'Sagital', help = 'option to add coronal image from montages', action = 'store_true')
parser.add_argument('-cor', '--coronal', dest = 'Coronal', help = 'option to add coronal image from montages', action = 'store_true')
parser.add_argument('--ROI' , dest = 'User_ROI', help = 'ex: "--ROI \'(0,0,0),(0,-53,26)\'"   MNI coordinates for seed voxel correlation computation', default = '[]')
parser.add_argument('--alpha' , dest = 'alpha', help = 'transparency value for montage overlay', type = float, default = 0.8)
parser.add_argument('--ROI_def', dest = 'ROI_default', help = 'specified if default mode seed voxel analysis to be shown in report. voxels already specified.', action = 'store_true')
parser.add_argument('--ROI_att', dest = 'ROI_attention', help = 'specified if attention network seed voxel analysis to be shown in report. voxels already specified.', action = 'store_true')
parser.add_argument('--ROI_ref', dest = 'ROI_refference', help = 'specified if refference network seed voxel analysis to be shown in report. voxels already specified.', action = 'store_true')
parser.add_argument('--anat', dest = 'anat', help = 'anatomical specified in meica.py (optional)', default = '')
parser.add_argument('--prefix', dest = 'prefix', help = 'prefix specified in meica.py', default = '')
parser.add_argument('--setname', dest = 'setname', help = 'directory meica.py creates', default = '')
parser.add_argument('--MNI', dest = 'MNI', help = 'specified if seed voxel coorelation wanted and MNI option was specified in meica.py', action = 'store_true')
args = parser.parse_args()

dep_check()

import meica_figures
import sphinx_files
import rst_files

reportdir = os.path.dirname(sys.argv[0])
p = subprocess.Popen('pwd', stdout = subprocess.PIPE, stderr = subprocess.PIPE)
startdir, err = p.communicate()
startdir = startdir[:-1]
setname = args.setname
User_ROI = seed_split(args.User_ROI)
prefix = args.prefix
file_check(args.anat, startdir, prefix, setname, args.MNI)

if args.ROI_default == True:#default mode MNI coordinates
	ROI_default = 	[[0,-53,26,'pC'],
					[0,52,6,'mPFC'],
					[-48,-62,36,'LatPar_L',],
					[46,-62,32,'Latpar_R'],
					[-24,-22,-20,'HF_L'],
					[24,-20,-22,'HF_R']]
else:
	ROI_default = []

if args.ROI_attention == True:#attention network MNI coordinates
	ROI_attention = [[-38,-4,48,'FEF_L'],
					[40,-4,48,'FEF_R'],
					[-24,-58,52,'IPS_L'],
					[22,-58,54,'IPS_R'],
					[-56,-60,-2,'MT+_L'],
					[54,-58,-4,'MT+_R']]
else:
	ROI_attention = []
				
if args.ROI_refference == True:#refference network MNI coordinates
	ROI_refference =[[-36,-25,57,'Mot_L'],
					[36,-25,57,'Mot_R'],
					[-43,-26,12,'Aud_L'],
					[43,-26,12,'Aud_R'],
					[-30,-88,0,'Vis_L'],
					[30,-88,0,'Vis_R']]
else:
	ROI_refference = []

corr = MNI_check(args.MNI, User_ROI, ROI_default, ROI_attention, ROI_refference)#checks to make sure MNI is True if ROI's are specified
if args.MNI == True:
	if ROI_default != []:
		meica_figures.check_ROI(ROI_default,startdir,setname, 'Default mode')#check default MNI within bounds
	if ROI_attention != []:
		meica_figures.check_ROI(ROI_attention,startdir,setname, 'Attention network')#check attention MNI within bounds
	if ROI_refference != []:
		meica_figures.check_ROI(ROI_refference,startdir,setname, 'Refference netowrk')#check refference MNI within bounds
	if User_ROI != []:
		meica_figures.check_ROI(User_ROI,startdir,setname, 'User specified')#check User_ROI MNI within bounds

os.chdir(startdir)

subprocess.call('mkdir %s/sphinx' % startdir, shell = True)#make directories
subprocess.call('mkdir %s/sphinx/_build' % startdir, shell = True)
subprocess.call('mkdir %s/sphinx/_static' % startdir, shell = True)
subprocess.call('mkdir %s/sphinx/_templates' % startdir, shell = True)

components = meica_figures.file_parse('%s/%s_ctab.txt' % (startdir,prefix))#collect components from ctab
maps = meica_figures.collect_data(args.anat,'%s/%s_mefl.nii.gz' % (startdir,prefix),'%s/%s/TED/feats_OC2.nii' % (startdir,setname))#collect nifti data
accept, reject, middle, ignore = meica_figures.split_components('%s/%s_ctab.txt' % (startdir,prefix), components)#seperate components into their respective bins

subprocess.call('mkdir %s/png_dump' % startdir, shell = True)
subprocess.call('cp %s/warning.png %s/png_dump' % (reportdir,startdir), shell = True)
os.chdir('%s/png_dump' % startdir)

#make figures
print('++ making figures')
# meica_figures.kr_vs_component('%s/%s_ctab.txt' % (startdir, prefix))#make kappa and rho vs component figure
# meica_figures.kappa_vs_rho_plot(accept, reject, middle, ignore)#make kappa vs rho figure
# meica_figures.tsnr('%s/%s_tsoc.nii.gz' % (startdir,prefix),'%s/%s_medn.nii.gz' % (startdir,prefix))#create tsnr figures
print('++ this set of figures may take a while')
# meica_figures.montage(maps, accept, args.montage_threshold, args.alpha, startdir, setname, prefix, args.Axial, args.Sagital, args.Coronal)#create activation montage
if args.anat != '':
	# meica_figures.coreg(startdir,setname,args.anat)#create corregistration figure
	if args.MNI == True:
		if len(accept) > 3:
			# meica_figures.correlation(startdir, setname, args.anat, ROI_default, ROI_attention, ROI_refference, User_ROI, args.corr_threshold)#create correlation for ROIs
			corr = True
		else:
			print '++ not enough degrees of freedom to compute standard error'
	else:
		print'++ cannot compute correlation, MNI coordinates not specified'
else:
	print '++ no anatomical specified, cannot create coregistration or correlation maps'

os.chdir('%s/sphinx' % startdir)

#set up sphinx documentation
sphinx_files.conf(__version__)
sphinx_files.make_bat()
sphinx_files.make_file()

#make .rst files for sphinx to use
print('++ occupying sphinx directory with .rst files')
rst_files.diagnostics_rst(args.anat)
rst_files.index_rst(corr)
rst_files.intro_rst()
rst_files.analysis_rst(accept, reject, middle, ignore, args.anat, args.montage_threshold, '%s/%s_ctab.txt' % (startdir,prefix),
	args.min_component_number, args.min_variance_explained)
if args.anat != '' and args.MNI == True and (len(ROI_default)>0 or len(ROI_attention)>0 or len(ROI_refference)>0 or len(User_ROI)>0):
	rst_files.correlation_rst(ROI_default,ROI_attention,ROI_refference,User_ROI)

#run sphinx build
subprocess.call('make html', shell = True)
subprocess.call('make latex', shell = True)
