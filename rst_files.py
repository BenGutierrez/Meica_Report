#!/usr/bin/env python
"""
Gutierrez, B.  Generates the rst files that sphinx reads to make the eport form
"""
import numpy as np

"""
make analysis.rst file
"""
def analysis_rst(accept, reject, middle, ignore, nsmprage, threshold, ctab, min_component_number, min_variance_explained):
	sl = []
	size = len(str(accept.shape[0] + reject.shape[0] + middle.shape[0] + ignore.shape[0]))
	sl.append('Component Visualization')
	sl.append('=======================')
	line = header(ctab)
	for i in range(len(line)):
		if ':' in line[i][:]:
			index = line[i][:].index(':')
		if line[i][1] == 'T' and line[i][7] == 'c' and int(line[i][index+1:]) < min_component_number:
			sl.append('%s !!!! |warning| !!!! number of components is below %s\n'  % (line[i][1:-1], min_component_number))
		if line[i][1] == 'D' and float(line[i][index+1:]) < min_variance_explained:
			sl.append('%s !!!! |warning| !!!! variance explained below %s\n'  % (line[i][1:-1], min_variance_explained))
		else:
			sl.append('%s' % line[i][1:])
	sl.append('\n+----------------+------------------+-------------------------+\n' +
				'|                | %  Total Vairance| %  Total Variance(norm) |\n' + 
				'+================+==================+=========================+')
	sl.append('| **Accepted**   |        %s    |          %s         |'   % (digit_length(sum(accept[:,3]),6), digit_length(sum(accept[:,4]),6)))
	sl.append('+----------------+------------------+-------------------------+\n' +
				'| **Rejected**   |        %s    |          %s         |' % (digit_length(sum(reject[:,3]),6), digit_length(sum(reject[:,4]),6)))
	sl.append('+----------------+------------------+-------------------------+\n' +
				'|**Middle kappa**|        %s    |          %s         |' % (digit_length(sum(middle[:,3]),6), digit_length(sum(middle[:,4]),6)))
	sl.append('+----------------+------------------+-------------------------+\n' +
				'| **Ignored**    |        %s    |          %s         |' % (digit_length(sum(ignore[:,3]),6), digit_length(sum(ignore[:,4]),6)))
	sl.append('+----------------+------------------+-------------------------+\n')

	sl.append('Graphs')
	sl.append('++++++')
	sl.append('.. image:: .. /png_dump/kappa_vs_rho.png')
	sl.append('	:width: 49%')
	sl.append('.. image:: .. /png_dump/kappa_rho_vs_components.png')
	sl.append('	:width: 49%\n')
	sl.append('The size of the scatter plot points is linearly related to the percent variance of that component.\n')
	sl.append('=============  =============  =================  =============')
	sl.append('# of Accepted  # of Rejected  # of Middle kappa  # of Ignored ')
	sl.append('=============  =============  =================  =============')
	sl.append('     %s             %s              %s                 %s    ' % 
		(len(accept), len(reject), len(middle), len(ignore)))
	sl.append('=============  =============  =================  =============\n\n')
	if nsmprage != '':
		sl.append('Accepted Components with anatomical')
		sl.append('+++++++++++++++++++++++++++++++++++')	
		sl.append('The following images are the thresholded components from the accepted bin of meica.py.  The threshold was set to %s.' % threshold)
		N = 0	

		for i in accept[:,0]:
			sl.append('\nComponent %s' % int(i))
			sl.append('----------'+ '-'*len(str(int(i))) + '\n')
			sl.append('.. image:: ../png_dump/Accepted_Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
			sl.append('	:scale: 75%')
			sl.append('	:align: left\n\n')
			sl.append('.. image:: .. /png_dump/FFT_Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
			sl.append('	:width: 49%\n')
			sl.append('=============  =============  =============  =================')
			sl.append('     kappa         rho         %s Variance     %s Variance(norm)' % ('%','%'))
			sl.append('=============  =============  =============  =================')
			sl.append('%s       %s         %s           %s       ' % 
				(digit_length(accept[N,1],8), digit_length(accept[N,2],7), digit_length(accept[N,3],3), digit_length(accept[N,4],3)))
			sl.append('=============  =============  =============  =================\n')
			N += 1

	sl.append('\nAccepted Components\n' + '+++++++++++++++++++')
	sl.append('The following images are the thresholded components from the accepted bin of meica.py.  The threshold was set to %s.' % threshold)
	N = 0	
	for i in accept[:,0]:
		sl.append('\nComponent %s' % int(i))
		sl.append('----------'+ '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n\n')
		sl.append('.. image:: .. /png_dump/FFT_Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:width: 49%\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(accept[N,1],8), digit_length(accept[N,2],7), digit_length(accept[N,3],3), digit_length(accept[N,4],3)))
		sl.append('=============  =============  =============  =================\n')
		N += 1

	sl.append('\nRejected Components\n' + '+++++++++++++++++++')
	sl.append('The following images are the grey scale components from the rejected bin of meica.py.')
	N = 0
	for i in reject[:,0]:
		sl.append('\nComponent %s' % int(i))
		sl.append('----------' + '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n\n')
		sl.append('.. image:: .. /png_dump/FFT_Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:width: 49%\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(reject[N,1],8), digit_length(reject[N,2],7), digit_length(reject[N,3],3), digit_length(reject[N,4],3)))
		sl.append('=============  =============  =============  =================\n')
		N += 1

	sl.append('\nMiddle Components\n' + '+++++++++++++++++++')
	sl.append('The following images are the grey scale components from the middle kappa bin of meica.py.')
	N = 0
	for i in middle[:,0]:
		sl.append('\nComponent %s' % int(i))
		sl.append('----------' + '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n\n')
		sl.append('.. image:: .. /png_dump/FFT_Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:width: 49%\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(middle[N,1],8), digit_length(middle[N,2],7), digit_length(middle[N,3],4), digit_length(middle[N,4],4)))
		sl.append('=============  =============  =============  =================\n')
		N += 1

	sl.append('\nIgnore Components\n' + '+++++++++++++++++++')
	sl.append('The following images are the grey scale components from the ignore bin of meica.py.')
	N = 0
	for i in ignore[:,0]: 
		sl.append('\nComponent %s' % int(i))
		sl.append('----------' + '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_'  + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n\n')
		sl.append('.. image:: .. /png_dump/FFT_Component_' + (size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:width: 49%\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(ignore[N,1],8), digit_length(ignore[N,2],7), digit_length(ignore[N,3],4), digit_length(ignore[N,4],4)))
		sl.append('=============  =============  =============  =================\n\n')
		N += 1
	sl.append('.. |warning| image:: ../png_dump/warning.png')
	sl.append('    		 :align: middle')
	sl.append('    		 :alt: warning')

	ofh = open("analysis.rst","w")
	ofh.write("\n".join(sl) + "\n")
	ofh.close()

"""
make diagnostics.rst file
"""
def diagnostics_rst(nsmprage):
	sl = []
	sl.append('Preliminary Diagnostics\n' + '==========================\n')
	sl.append('The signal to noise ratio (TSNR) for a dataset is defined as the mean over the standard deviation of the dataset.  ' +
		'Meica.py denoises the BOLD time signal which should increase TSNR. This will be seen as the medn dataset having a greater ' + 
		'TSNR than the tsoc dataset.\n')
	sl.append('TSNR\n' + '+++++++')
	sl.append('Absolute TSNR of the medn NIFTI1 dataset:\n')
	sl.append('.. image:: .. /png_dump/medn_tsnr.png')
	sl.append('	:scale: 99%')
	sl.append('	:align: left\n')
	sl.append('Absolute TSNR of the tsoc NIFTI1 dataset:\n')
	sl.append('.. image:: .. /png_dump/tsoc_tsnr.png')
	sl.append('	:scale: 99%')
	sl.append('	:align: left\n')
	sl.append('Ratio of the TSNR of the medn NIFTI1 dataset to the TSNR of the tsoc NIFTI1 dataset:\n')
	sl.append('.. image:: .. /png_dump/tsnr_ratio.png')
	sl.append('	:scale: 99%')
	sl.append('	:align: left\n')
	sl.append('Histograms of the TSNR of the medn NIFTI1 dataset and the TSNR of the tsoc NIFTI1 dataset:\n')
	sl.append('.. image:: .. /png_dump/medn_tsnr_hist.png')
	sl.append('	:width: 49%')
	sl.append('.. image:: .. /png_dump/tsoc_tsnr_hist.png')
	sl.append('	:width: 49%\n')
	sl.append('Histogram of the ratio of the TSNR of the medn dataset to the TSNR of the tsoc NIFTI1 dataset:\n')
	sl.append('.. image:: .. /png_dump/tsnr_ratio_hist.png')
	sl.append('	:width: 49%\n')
	
	if nsmprage != '':
		sl.append('Coregistration\n' + '+++++++++++++++')
		sl.append('Coregistration of the subject to the anatomical provided:\n')
		sl.append('.. image:: .. /png_dump/coregistration.png')
		sl.append('	:scale: 75%')
	ofh = open("diagnostics.rst","w")
	ofh.write("\n".join(sl) + "\n")
	ofh.close()

"""
make index.rst file
"""
def index_rst(cor):
	sl = []
	sl.append('Welcome to your ME-ICA Report!')
	sl.append('==============================')
	sl.append('The program meica.py was created to form an algorithmic method for performing '
		+ 'independent component analysis on multi-echo data and then algorithmically deciding '
		+ 'which components represent BOLD-like phenomena.\n')

	sl.append('The following content is a report that has taken information provided by '
		+ 'meica.py and summarizes a few of the results.\n')
	sl.append('This report form was created by the Section on Functional Imaging Methods in the NIMH.  '
		+ 'The creators of this report form are Benjamin Gutierrez, Prantik Kundu, Daniel Handwerker, '
		+ 'Javier Gonzalez-Castillo, Souheil Inati, and Peter Bandettini.\n')
	sl.append('Contents:\n')
	sl.append('.. toctree::')
	sl.append('   :maxdepth: 2\n')
	sl.append('   intro\n' + '   diagnostics') 
	if cor:
		sl.append('   correlation')
	sl.append('   analysis')
	sl.append('\n\n\nSearch')
	sl.append('======\n')
	sl.append('* :ref:`search`')

	ofh = open("index.rst","w")
	ofh.write("\n".join(sl) + "\n")
	ofh.close()

"""
make intro.rst file
"""
def intro_rst():
	sl = []
	sl.append('Intro\n' + '=====')
	sl.append('This report collects information from your ME-ICA analysis and displays '
		 + 'several ways to interpret the output of meica.py.  The Preliminary Diagnostics '
		 + 'section contains a view of the TSNR of the "denoised" BOLD time series after: basic '
		 + 'preprocessing, T2* weighted averaging of echoes (i.e. "optimal combination"), and ICA denoising, the TSNR of the '
		 + '"raw" BOLD time series dataset after: basic preprocessing, and a ratio of these two TSNR maps.  '
		 + 'also in this section is a view of the corregistration if an anatomical was specified when using meica.py\n\n')
	sl.append('The Correlation section exists if the MNI option was used when running meica.py and contains the correlation maps '
		 + 'of the whole brain with several different seed voxels.  These seed voxels were chosen based on networks that may be of '
		 + 'particular interest to the viewer, i.e. seeds for the Default network.\n\n')
	sl.append('In the Component Visualization section there are two different ways to visualize the information from the accepted components. '
		 + 'The first is a thresholding of the bold response overlayed onto the anatomical image '
		 + 'that was used for correlation (if one was used).  If no anatomical specified, this last option does not exist.  The other way to view '
		 + 'the accepted component, which will always be present, are the non-thresholded grey scale images without an anatomical underlay.  '
		 + 'For the rejected, middle kappa, and ignored components, the grey scale image is always displayed as well.\n')
	sl.append('Requirements for generating this report form:\n\n' + '* matplotlib\n\n' + '* numpy\n\n' + '* nibabel\n\n'  + '* parse\n\n' + '* sphinx\n')
	ofh = open("intro.rst","w")
	ofh.write("\n".join(sl) + "\n")
	ofh.close()


def correlation_rst(ROI_default=[], ROI_attention=[], ROI_refference=[], User_ROI=[]):
	ROI = [ROI_default,ROI_attention,ROI_refference]
	sl = []
	sl.append('Correlation\n' + '==========================\n')
	sl.append('In this section, several seeds were chosen to represent a few networks that are common in the literature.  ' +
	  'These seeds were chosen using MNI coordinates derived by Koene R. A. Van Dijk [1]_.  Z scores were capped at 5 and -5 for the below figures.\n')
	sl.append('.. [1] **Koene R. A. Van Dijk, Trey Hedden, Archana Venkataraman, Karleyton C. Evans, Sara W. Lazar and Randy L. Buckner.**  ' +
		'Intrinsic Functional Connectivity As a Tool For Human Connectomics: Theory, Properties, and Optimization. *J Neurophysiol 103:297-321, 2010*\n')
	for j in range(3):
		if j == 0 and ROI_default != []:
			sl.append('Default network\n' + '++++++++++++++++++++')
			sl.append('Below is our seed based correlation for several voxels thought to be part of the default network.\n')
		if j == 1 and ROI_attention != []:
			sl.append('Attention network\n' + '+++++++++++++++++++++')
			sl.append('Below is our seed based correlation for several voxels thought to be part of the attention network.\n')
		if j == 2 and ROI_refference != []:
			sl.append('Refference network\n' + '+++++++++++++++++++++')
			sl.append('Below is our seed based correlation for several voxels thought to be part of the Refeence network.\n')
		for i in range(len(ROI[j])):
			title = ROI[j][i][3]
			if '_R' in title:
				title = 'Right ' + title[:-2]
			if '_L' in title:
				title = 'Left ' + title[:-2]
			sl.append('%s' % title)
			sl.append('-'*(len(title)+1))
			sl.append('MNI coordinates(%s,%s,%s)\n' % (ROI[j][i][0], ROI[j][i][1], ROI[j][i][2]))
			sl.append('.. image:: .. /png_dump/%s_seed.png\n' % ROI[j][i][3])
			sl.append('.. image:: .. /png_dump/%s_correlation.png' % ROI[j][i][3])
			sl.append('	:scale: 75%\n')
	if User_ROI != []:
		sl.append('User ROI\n' + '+++++++++++++++')
		for i in range(len(User_ROI)):
			sl.append('MNI coordinates(%s,%s,%s)' % (User_ROI[i][0], User_ROI[i][1], User_ROI[i][2]))
			sl.append('------------------------------')
			sl.append('.. image:: .. /png_dump/(%s,%s,%s)_seed.png\n' % (User_ROI[i][0], User_ROI[i][1], User_ROI[i][2]))
			sl.append('.. image:: .. /png_dump/(%s,%s,%s)_correlation.png' % (User_ROI[i][0], User_ROI[i][1], User_ROI[i][2]))
			sl.append('	:scale: 75%\n')

	ofh = open("correlation.rst","w")
	ofh.write("\n".join(sl) + "\n")
	ofh.close()
"""
Controls length of floats.  Needs to be done because spacing important for sphinx.
This function is not designed to work if length argument is less than the number of digits before the decimal.
"""
def digit_length(n,length):
	n_string = str(n)
	if len(n_string) < length:
		n_string = n_string + ' '*(length - len(n_string))
	elif len(n_string) > length:
		n_string = str(round(n,length))
	return(n_string)

"""
Pulls the header from ctab
"""
def header(ctab):
	txt_file = open(str(ctab))
	txt_file.seek(0)
	line = np.zeros(11, dtype = object)
	for i in range(2):
		txt_file.readline()
	for i in range(11):
		line[i] = txt_file.readline()
	return(line)
