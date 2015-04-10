#!/usr/bin/env python
"""
Gutierrez, B.  Generates Meica report form.

"""
__version__ = "v2.5 beta10"

import commands
import subprocess
import argparse
import sys
import ast
import os
import re

ROI_default = 	[[0,-53,26,'pC'],#default mode MNI coordinates
				[0,52,6,'mPFC'],
				[-48,-62,36,'LatPar_L',],
				[46,-62,32,'Latpar_R'],
				[-24,-22,-20,'HF_L'],
				[24,-20,-22,'HF_R']]

ROI_attention = [[-38,-4,48,'FEF_L'],#attention network MNI coordinates
				[40,-4,48,'FEF_R'],
				[-24,-58,52,'IPS_L'],
				[22,-58,54,'IPS_R'],
				[-56,-60,-2,'MT+_L'],
				[54,-58,-4,'MT+_R']]

ROI_reference = [[-36,-25,57,'Mot_L'],#reference network MNI coordinates
				[36,-25,57,'Mot_R'],
				[-43,-26,12,'Aud_L'],
				[43,-26,12,'Aud_R'],
				[-30,-88,0,'Vis_L'],
				[30,-88,0,'Vis_R']]

def dep_check():
	print '++ INFO: Checking system for dependencies...'
	fails = 0
	numpy_installed = 0
	sphinx_installed = 0

	try:
		import numpy
		numpy_installed = 1
	except:
		print "++ Error: Can't import Numpy! Please check Numpy installation for this Python version."
		fails += 1
	try:
		import matplotlib
	except:
		print "++ Error: Can't import Matplotlib! Please check Matplotlib installation."
		fails += 1
	try:
		import nibabel
	except:
		print "++ Error: Can't import Nibabel! Please check Nibabel installation."
		fails += 1
	try:
		import sphinx
		sphinx_installed = 1
	except:	
		print "++ Error: Can't import Sphinx! Please check Sphinx installation."
		fails += 1
	if sphinx_installed:
		print "++ INFO: Sphinx version: %s" % (sphinx.__version__)
		if float((sphinx.__version__)[0:3])<1.2:
			fails += 1
			print "++ Error: Sphinx version is too old! Please upgrade to Sphinx >=1.2.x!"
	if numpy_installed:
		print "++ INFO: Numpy version: %s" % (numpy.__version__)
		if float((numpy.__version__)[0:3])<1.5:
			fails += 1
			print "++ Error: Numpy version is too old! Please upgrade to Numpy >=1.5.x!"
		import numpy.__config__ as nc
		if nc.blas_opt_info == {}:
			fails += 1
			print "++ Error: Numpy is not linked to BLAS! Please check Numpy installation."

	afnicheck = commands.getstatusoutput("3dinfo")
	if afnicheck[0]!=0:
		print "++ Error: Can't run AFNI binaries. Make sure AFNI is on the path!"
		fails += 1
	elif not afnicheck[1].__contains__('Alternate Alternative Usage'):
		print "++ Error: This seems like an old version of AFNI. Please upgrade to latest version of AFNI."
		fails += 1
	if fails == 0:
		print "++ INFO: Dependencies OK."
	else:
		print "++ INFO: EXITING. Please see above error messages."
		sys.exit()

def MNI_check(MNI, User_ROI, ROI_def, ROI_att, ROI_ref):
	if not MNI:
		if User_ROI or ROI_att or ROI_def or ROI_ref:
			print '++ Error: MNI needs to be specified if options User_ROI, ROI_att, ROI_def, or ROI_ref are specified'
			sys.exit()
		corr = False

	else:
		if not (User_ROI or ROI_att or ROI_def or ROI_ref):
			print '++ Error: !!MNI specified but not User_ROI, ROI_att, ROI_def, or ROI_ref.  No seed based correlation will be computed!!'
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

def file_check(anat, startdir, TED, setname, MNI, reportdir, coreg_anat, coreg):
	fails = 0

	if not os.path.isfile('%s/meica_figures.py' % reportdir):
		print '++ Error: Can\'t find meica_figures.py in %s.  The files from the meica_report directory cannot be seperated into different directories' % reportdir
		fails += 1

	if not os.path.isfile('%s/rst_files.py' % reportdir):
		print '++ Error: Can\'t find rst_files.py in %s.  The files from the meica_report directory cannot be seperated into different directories' % reportdir
		fails += 1

	if not os.path.isfile('%s/sphinx_files.py' % reportdir):
		print '++ Error: Can\'t find sphinx_files.py in %s.  The files from the meica_report directory cannot be seperated into different directories' % reportdir
		fails += 1

	if not os.path.isfile('%s/warning.png' % reportdir):
		print '++ Error: Can\'t find warning.png in %s.  This means that picture warning flags for too few components or low variance will not appear.' % reportdir

	if not os.path.isfile(anat) and anat != '': 
		print '++ Error: Can\'t find the specified anantomical'
		fails += 1

	if not os.path.isfile('%s/comp_table.txt' % TED):
		print '++ Error: Can\'t find "%s/comp_table.txt" check directory and file\'s existance.' % TED
		fails += 1

	if not os.path.isfile('%s/betas_OC.nii' % TED):
		print '++ Error: Can\'t find "%s/betas_OC.nii" check directory and file\'s existance.' % TED
		fails += 1

	if not os.path.isfile('%s/ts_OC.nii' % TED):
		print '++ Error: Can\'t find "%s/ts_OC.nii" check directory and file\'s existance.' % TED
		fails += 1

	if not os.path.isfile('%s/dn_ts_OC.nii' % TED):
		print '++ Error: Can\'t find "%s/dn_ts_OC.nii" check directory and file\'s existance.' % TED
		fails += 1

	if not os.path.isfile('ocv_uni_vr.nii.gz') and (coreg and not os.path.isfile('%s/%s' % (setname,coreg_anat))):
		print '++ Error: Can\'t find "%s/ocv_uni_vr.nii.gz check directory and file\'s existance.' % setname
		fails += 1

	if not os.path.isfile('%s/meica_mix.1D' % TED):
		print '++ Error: Can\'t find "%s/meica_mix.1D" check directory and file\'s existance.' % TED
		fails += 1

	if coreg_anat != '':	
		if not os.path.isfile('%s/%s' % (setname,coreg_anat)):
			print '++ Error: Can\'t find %s in %s. Please check that specified file exists' % (coreg_anat,setname)
			fails += 1

	if MNI == True:
		if not os.path.isfile('eBvrmask.nii.gz'):
			print '++ Error: Can\'t find "%s/eBvrmask.nii.gz check directory and file\'s existance.' % setname
			fails += 1

		if not os.path.isfile('%s/feats_OC2.nii' % TED):
			print '++ Error: Can\'t find "%s/feats_OC2.nii" check directory and file\'s existance.' % TED
			fails += 1

	if fails != 0:
		print "++ INFO: EXITING. Please see error messages."
		sys.exit()
#Run dependency check

def path_names(setname, startdir, TED, anat):
	if setname == '' or setname == None:
		print '++ Error: Need to specify the option -setname to run meica_report.py'
		sys.exit()

	if not os.path.isdir(setname):
		print '++ Error: -setname argument "%s" not found' % setname
		sys.exit()

	if startdir == '':#make sure paths are in correct form and remove wildcards and shortcuts
		startdir = os.getcwd()

	startdir = os.path.abspath(os.path.expanduser(startdir))
	setname  = os.path.abspath(os.path.expanduser(setname))
	TED      = os.path.abspath(os.path.expanduser(TED))

	if anat != '':
		anat = os.path.abspath(os.path.expanduser(anat))

	return(setname, startdir, TED, anat)


if __name__=='__main__':
	parser = argparse.ArgumentParser('Options')
	Required = parser.add_argument_group('Required arguments')
	Required.add_argument('-setname', dest = 'setname', help = 'Directory meica.py creates.  Will be of the form "meica.foo"')
	lab = parser.add_argument_group('File and directory labels')
	lab.add_argument('-label', dest = 'label', help = 'Label to tag directory for all output files, default is "Report"  ', default = 'Report')
	lab.add_argument('-anat', dest = 'anat', help = 'Anatomical specified in meica.py (optional)', default = '')
	lab.add_argument('-dir', dest = 'startdir', help = 'Directory to place report directory.  Default is current directory', default = '')
	lab.add_argument('-TED', dest = 'TED', help = 'Directory containing all files from tedana.py processing steps.  Input files are taken automatically from this directory. Default is "-label/TED"', default = 'TED')
	lab.add_argument('-motion', dest = 'motion_file',help='file containing motion (do not include path). file must be in -setname.  default is motion.1D in -setname.',default='motion.1D')
	lab.add_argument('-overwrite', dest = 'overwrite', help = 'If -overwrite specified and -label directory already exists, will overwrite', action = 'store_false')
	options = parser.add_argument_group('Report options')
	options.add_argument('-ax' , dest = 'Axial', help = 'Add axial images to activation montages', action = 'store_true')
	options.add_argument('-sag', dest = 'Sagittal', help = 'Add sagittal images to activation montages', action = 'store_true')
	options.add_argument('-cor', dest = 'Coronal', help = 'Add coronal images to activation montages', action = 'store_true')
	options.add_argument('-coreg', dest = 'coreg', help = argparse.SUPPRESS, action = 'store_true')#'If specified, redering corregistration.  Need anatomical'---may be buggy
	options.add_argument('-coreg_anat', dest = 'coreg_anat', help = argparse.SUPPRESS, default = '')#'Optional name of anatomical to corregister with.  MUST be in the -\'setname directory\'.  If none specified, will attempt to guess name'---may be buggy
	options.add_argument('-flood', dest = 'flood', help = 'Tells clustering algoithim (flood fill) how many voxels above threshold need to be clustered together (in 26 possible directions) to be kept in activation map. Specify as "0" if you want no clustering', type = int, default = 10)
	options.add_argument('-contrast', dest = 'contrast', help = 'Give contrast to greyscale images in montage.  Ex: "5" will give values in the 5-95 percentile of values.  Default = "2"',type = int, default = 2)
	options.add_argument('-show_ROI', dest = 'show', help = 'Shows prespecified MNI coordinates for Default mode, attention network, and reference network for seed based correlation.  Will NOT make report if specified.', action = 'store_true')
	options.add_argument('-ROI' , dest = 'User_ROI', help = 'ex: "--ROI \'(0,0,0),(0,-53,26)\'"   MNI coordinates for seed voxel correlation computation', default = '[]')
	options.add_argument('-ROI_def', dest = 'ROI_default', help = 'If specified default mode network seed voxel analysis to be shown in report. voxels already specified in code.', action = 'store_true')
	options.add_argument('-ROI_att', dest = 'ROI_attention', help = 'If specified attention network seed voxel analysis to be shown in report. voxels already specified in code.', action = 'store_true')
	options.add_argument('-ROI_ref', dest = 'ROI_reference', help = 'If specified reference network seed voxel analysis to be shown in report. voxels already specified in code.', action = 'store_true')
	options.add_argument('-MNI', dest = 'MNI', help = 'If specified will compute seed voxel correlation needs MNI option specified in meica.py', action = 'store_true')
	options.add_argument('-mt', dest = 'montage_threshold', help = 'Normalized z-score value for thresholding the accepted component w/ anatomical activation map, default = 1.96', type = float, default = 1.96)
	options.add_argument('-ct', dest = 'corr_threshold', help = 'Normalized z-score value for thresholding the seed voxel map, default = 1.96', type = float, default = 1.96)
	options.add_argument('-latex', dest = 'latex', help = 'If specified, will use pdfLatex to build tex version of meica report.', action = 'store_true')
	options.add_argument('-min_v', dest = 'min_variance_explained', help = 'Minimum variance explained before warning raised in report, default = 85', type = int, default = 85)
	options.add_argument('-min_c', dest = 'min_component_number', help = 'Minimum total component number before warning raised in report, default = 20', type = int, default = 20)
	options.add_argument('-alpha' , dest = 'alpha', help = 'Transparency value for montage overlay', type = float, default = 0.8)
	options.add_argument('-title', dest = 'title', help = 'Title of ME-ICA report.  Will be shown on browser tab for easier identification between multiple reports.', default = 'Your ME-ICA Report!')
	args = parser.parse_args()

	if args.show:
		print '++ INFO: Default mode netowrk seed MNI coordinates:\n' + str(ROI_default).replace('],','],\n')
		print '\n++ INFO: Attention network seed MNI cooridnates:\n' + str(ROI_attention).replace('],','],\n')
		print '\n++ INFO: Reference network seed MNI coordinates:\n' + str(ROI_reference).replace('],','],\n')
		print '++ INFO: To create report, do not specify -show_ROI'
		sys.exit()

	dep_check()
	import meica_figures
	import sphinx_files
	import rst_files

	meica_txt=[]
	reportdir = os.path.abspath(os.path.dirname(sys.argv[0]))
	label = 'meica.' + args.label
	figures = 'Report_Figures'

	setname, startdir, TED, anat = path_names(args.setname, args.startdir, args.TED, args.anat)
	User_ROI = seed_split(args.User_ROI)
	os.chdir(setname)

	file_check(anat, startdir, TED, setname, args.MNI, reportdir, args.coreg_anat, args.coreg)


	if os.path.isdir('%s/%s' % (startdir,label)) and args.overwrite:
		print '++ Error: %s directory already exits and -overwrite not specified' % label
		sys.exit()
	if os.path.isdir('%s/%s' % (startdir,label)):
		subprocess.call('rm -rf %s/%s' % (startdir,label), shell = True)

	ctab = '%s/comp_table.txt' % TED
	tsoc = '%s/%s/axialized_nifti/tsoc.nii.gz' % (startdir,label)
	medn =  '%s/%s/axialized_nifti/medn.nii.gz' % (startdir,label)
	mefl = '%s/betas_OC.nii' % TED

	if os.path.isfile('%s/feats_OC2.nii' % TED):
		feats = '%s/feats_OC2.nii' % TED
	else:
		feats = ''

	if not args.ROI_default:
		ROI_default = []
	if not args.ROI_attention:
		ROI_attention = []
	if not args.ROI_reference:
		ROI_reference = []

	corr = MNI_check(args.MNI, User_ROI, ROI_default, ROI_attention, ROI_reference)#checks to make sure MNI is True if ROI's are specified
	if args.MNI:
		if ROI_default != []:
			meica_figures.check_ROI(ROI_default,setname, TED,'Default mode')#check default MNI within bounds
		if ROI_attention != []:
			meica_figures.check_ROI(ROI_attention,setname, TED,'Attention network')#check attention MNI within bounds
		if ROI_reference != []:
			meica_figures.check_ROI(ROI_reference,setname, TED,'Reference netowrk')#check reference MNI within bounds
		if User_ROI != []:
			meica_figures.check_ROI(User_ROI,setname, TED,'User specified')#check User_ROI MNI within bounds

	subprocess.call('mkdir %s/%s' % (startdir,label), shell = True)#make directories
	subprocess.call('mkdir %s/%s/axialized_nifti' % (startdir,label), shell = True)
	if corr:
		subprocess.call('3daxialize -overwrite -prefix %s/%s/axialized_nifti/betas_hik_OC.nii %s/betas_hik_OC.nii' % (startdir,label,TED), shell = True)

	maps = meica_figures.collect_data(startdir,label,TED,anat,mefl,feats)#collect nifti data
	accept, reject, middle, ignore = meica_figures.components(TED)


	subprocess.call('mkdir %s/%s/%s' % (startdir,label,figures), shell = True)
	subprocess.call('cp %s/warning.png %s/%s/%s' % (reportdir,startdir,label,figures), shell = True)
	os.chdir('%s/%s/%s' % (startdir,label,figures))

	#make figures
	print '++ INFO: Making figures'
	print '\t Making Kappa vs component number plot'
	meica_figures.kr_vs_component(ctab)#make kappa and rho vs component figure
	print '\t Making Kappa vs Rho plot'
	meica_txt.append(meica_figures.kappa_vs_rho_plot(accept, reject, middle, ignore,ctab,startdir,label,figures))#make kappa vs rho figure
	print '\t Making TSNR plots'
	meica_txt.append(meica_figures.tsnr(tsoc,medn,startdir,label))#create tsnr figures
	print '\t Making motion plots'
	if os.path.isfile('%s/%s' % (setname,args.motion_file)):
		meica_txt.append(meica_figures.motion(startdir,label,figures,setname,args.motion_file))
	else:
		meica_txt.append("Max head displacement in any one dirrection:   %s\nTR of Max Head displacement:   %s\nMax rate of head motion:   %s\nTR of max head motion rate:   %s" % (' ',' ',' ',' '))

	print('\t Making component images.  This set of figures may take awhile...')
	meica_figures.montage(maps, accept, args.montage_threshold, args.alpha, TED, args.Axial, args.Sagittal, args.Coronal, args.flood, args.contrast)#create activation montage
	if anat != '':
		if args.coreg:
			meica_figures.coreg(startdir, setname, label, figures, anat, args.coreg_anat)#create corregistration figure
		if args.MNI and corr:
			if len(accept) > 3:
				meica_figures.correlation(startdir, label, figures, anat, ROI_default, ROI_attention, ROI_reference, User_ROI, args.corr_threshold)#create correlation for ROIs
				if args.ROI_reference + args.ROI_attention + args.ROI_default + len(args.User_ROI) != 0:
					corr = True
			else:
				print '++ Error: Not enough degrees of freedom to compute standard error for correlation'
		else:
			print'++ Error: Cannot compute correlation, MNI coordinates not specified'
	else:
		print '++ Error: No anatomical specified, cannot create coregistration or correlation maps'

	os.chdir('%s/%s' % (startdir,label))

	# set up sphinx documentation
	sphinx_files.conf(__version__)
	sphinx_files.make_bat()
	sphinx_files.make_file()

	#make .rst files for sphinx to use to generate the report
	print('++ INFO: Occupying sphinx directory with .rst files')
	rst_files.diagnostics_rst(anat,args.coreg,figures)
	rst_files.index_rst(corr, args.title)
	rst_files.intro_rst()
	rst_files.analysis_rst(accept, reject, middle, ignore, anat, args.montage_threshold, ctab,
		args.min_component_number, args.min_variance_explained, figures, setname, args.motion_file, args.Axial, args.Sagittal, args.Coronal)
	if anat != '' and args.MNI and (len(ROI_default)> 0 or len(ROI_attention)> 0 or len(ROI_reference)> 0 or len(User_ROI)>0):
		rst_files.correlation_rst(ROI_default,ROI_attention,ROI_reference,User_ROI,figures)
	ofh = open("meica_report.txt","w")
	ofh.write("\n".join(meica_txt) + "\n")
	ofh.close()

	# run sphinx build
	subprocess.call('make html', shell = True)
	subprocess.call('make latex', shell = True)

	if args.latex:
		subprocess.call('make latexpdf', shell = True)

	subprocess.call('mv %s/%s/_build/* %s/%s' % (startdir,label,startdir,label), shell = True)
	subprocess.call('mv %s/%s/_static/* %s/%s' % (startdir,label,startdir,label), shell = True)
	subprocess.call('rm -rf _*', shell = True)
	subprocess.call('mkdir %s/%s/sphinx_files' % (startdir,label), shell = True)
	subprocess.call('mv %s/%s/*.rst %s/%s/sphinx_files/' % (startdir,label,startdir,label), shell = True)
	subprocess.call('mv %s/%s/Makefile %s/%s/sphinx_files' % (startdir,label,startdir,label), shell = True)
	subprocess.call('mv %s/%s/make.bat %s/%s/sphinx_files' % (startdir,label,startdir,label), shell = True)
	subprocess.call('mv %s/%s/conf.py %s/%s/sphinx_files' % (startdir,label,startdir,label), shell = True)
