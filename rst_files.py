#!/usr/bin/env python
"""
Gutierrez, B.  Generates the analysis.rst file
"""
import numpy as np

"""
make analysis.rst file
"""
def analysis_rst(accept, reject, middle, ignore, threshold, ctab, min_component_number, min_variance_explained):
	sl = []
	size = len(str(accept.shape[0] + reject.shape[0] + middle.shape[0] + ignore.shape[0]))
	sl.append('Component Visualization')
	sl.append('=======================')
	line = header(ctab)
	for i in range(len(line)):
		if ':' in line[i][:]:
			index = line[i][:].index(':')
		if line[i][1] == 'T' and line[i][7] == 'c' and int(line[i][index+1:]) < min_component_number:
			sl.append('%s !!!! |warning| !!!! number of components is below %s\n'  % (line[i][1:-1],min_component_number))
		if line[i][1] == 'D' and float(line[i][index+1:]) < min_variance_explained:
			sl.append('%s !!!! |warning| !!!! variance explained below %s\n'  % (line[i][1:-1],min_variance_explained))
		else:
			sl.append('%s' % line[i][1:])
	sl.append('\n+----------------+------------------+-------------------------+\n' +
				'|                | %  Total Vairance| %  Total Variance(norm) |\n' + 
				'+================+==================+=========================+')
	sl.append('| **Accepted**   |        %s    |          %s         |' % (digit_length(sum(accept[:,3]),6),digit_length(sum(accept[:,4]),6)))
	sl.append('+----------------+------------------+-------------------------+\n' +
				'| **Rejected**   |        %s    |          %s         |' % (digit_length(sum(reject[:,3]),6),digit_length(sum(reject[:,4]),6)))
	sl.append('+----------------+------------------+-------------------------+\n' +
				'|**Middle kappa**|        %s    |          %s         |' % (digit_length(sum(middle[:,3]),6),digit_length(sum(middle[:,4]),6)))
	sl.append('+----------------+------------------+-------------------------+\n' +
				'| **Ignored**    |        %s    |          %s         |' % (digit_length(sum(ignore[:,3]),6),digit_length(sum(ignore[:,4]),6)))
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
	sl.append('=============  =============  =================  =============\n')
	sl.append('Accepted Components with anatomical')
	sl.append('+++++++++++++++++++++++++++++++++++')	
	sl.append('The following images are the thresholded components from the accepted bin of meica.py.  The threshold was set to %s.' % threshold)
	N = 0	
	for i in accept[:,0]:
		sl.append('\nComponent %s' % int(i))
		sl.append('----------'+ '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Accepted_Component_'+(size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance     %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(accept[N,1],8),digit_length(accept[N,2],7),digit_length(accept[N,3],3),digit_length(accept[N,4],3)))
		sl.append('=============  =============  =============  =================\n')
		N += 1

	sl.append('Accepted Components')
	sl.append('+++++++++++++++++++')
	sl.append('The following images are the thresholded components from the accepted bin of meica.py.  The threshold was set to %s.' % threshold)
	N = 0	
	for i in accept[:,0]:
		sl.append('\nComponent %s' % int(i))
		sl.append('----------'+ '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_'+(size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(accept[N,1],8),digit_length(accept[N,2],7),digit_length(accept[N,3],3),digit_length(accept[N,4],3)))
		sl.append('=============  =============  =============  =================\n')
		N += 1

	sl.append('Rejected Components\n' + '+++++++++++++++++++')
	sl.append('The following images are the grey scale components from the rejected bin of meica.py.')
	N = 0
	for i in reject[:,0]:
		sl.append('\nComponent %s' % int(i))
		sl.append('----------' + '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_'+(size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(reject[N,1],8),digit_length(reject[N,2],7),digit_length(reject[N,3],3),digit_length(reject[N,4],3)))
		sl.append('=============  =============  =============  =================')
		N += 1

	sl.append('\nMiddle Components\n' + '+++++++++++++++++++')
	sl.append('The following images are the grey scale components from the middle kappa bin of meica.py.')
	N = 0
	for i in middle[:,0]:
		sl.append('\nComponent %s' % int(i))
		sl.append('----------' + '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_' +(size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(middle[N,1],9),digit_length(middle[N,2],7),digit_length(middle[N,3],4),digit_length(middle[N,4],4)))
		sl.append('=============  =============  =============  =================')
		N += 1

	sl.append('\nIgnore Components\n' + '+++++++++++++++++++')
	sl.append('The following images are the grey scale components from the ignore bin of meica.py.')
	N = 0
	for i in ignore[:,0]: 
		sl.append('\nComponent %s' % int(i))
		sl.append('----------' + '-'*len(str(int(i))) + '\n')
		sl.append('.. image:: ../png_dump/Component_'  +(size - len(str(int(i))))*'0' + '%s.png' % int(i))
		sl.append('	:scale: 75%')
		sl.append('	:align: left\n')
		sl.append('=============  =============  =============  =================')
		sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
		sl.append('=============  =============  =============  =================')
		sl.append('%s       %s         %s           %s       ' % 
			(digit_length(ignore[N,1],8),digit_length(ignore[N,2],7),digit_length(ignore[N,3],4),digit_length(ignore[N,4],4)))
		sl.append('=============  =============  =============  =================\n\n')
		N += 1
	sl.append('.. |warning| image:: ../png_dump/warning.png')
	sl.append('    		 :align: middle')
	sl.append('    		 :alt: warning')

	ofh = open("analysis.rst","w")
	ofh.write("\n".join(sl)+"\n")
	ofh.close()

"""
make diagnostics.rst file
"""
def diagnostics_rst(nsmprage):
	sl = []
	sl.append('Preliminary Diagnostics\n' + '==========================\n')
	sl.append('The signal to noise ratio (TSNR) for a dataset is defined as the mean over the standard deviationof the dataset.  ' +
		'Meica.py denoises the BOLD time signal which should increase TSNR. This will be seen as the medn dataset has a greater ' + 
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
	ofh.write("\n".join(sl)+"\n")
	ofh.close()

"""
make index.rst file
"""
def index_rst():
	sl = []
	sl.append('Welcome to your ME-ICA Report!')
	sl.append('==============================')
	sl.append('The program meica.py was created to form an algorithmic method for performing '
		+ 'independent component analysis on multi-echo data and then algorithmically deciding '
		+ 'which components represent BOLD-like phenomena.\n')

	sl.append('The following content is a report that has taken information provided by '
		+ 'meica.py and summerizes a few of the results.\n')
	sl.append('This report form was created by the Section on Functional Imaging Methods in the NIMH.  '
		+ 'The creators of this report form are Benjamin Gutierrez, Prantik Kundu, Daniel Handwerker, '
		+ 'Javier Gonzalez-Castillo, Souheil Inati, and Peter Bandettini.\n')
	sl.append('Contents:\n')
	sl.append('.. toctree::')
	sl.append('   :maxdepth: 2\n')
	sl.append('   intro\n' + '   diagnostics\n' + '   correlation\n' + '   analysis')
	sl.append('\n\n\nSearch')
	sl.append('======\n')
	sl.append('* :ref:`search`')

	ofh = open("index.rst","w")
	ofh.write("\n".join(sl)+"\n")
	ofh.close()

"""
make intro.rst file
"""
def intro_rst():
	sl = []
	sl.append('Intro\n' + '=====')
	sl.append('This report collects information from your ME-ICA analysis and displays '
		 + 'several ways to interpret the output of meica.py.  The preliminary diagnostics '
		 + 'section contains a view of the TSNR of the medn dataset, the "denoised" BOLD time series after: basic '
		 + 'preprocessing, T2* weighted averaging of echoes (i.e. "optimal combination"), and ICA denoising, produced by meica.py '
		 + 'and also the TSNR ratio of the medn dataset over tsoc dataset, "raw" BOLD time series dataset after: basic preprocessing '
		 + ' and T2* weighted averaging of echoes (i.e. "optimal combination")\n\n')
	sl.append('This report also includes two different ways to visualize the information from the accepted components. '
		 + 'The first is a thresholding of the bold responce overlayed onto the anatomical image '
		 + 'that was used for correlation.  If this is not present, then it is because this option '
		 + 'was not included when this report was generated.  The other option, which should always '
		 + 'be present are the non-thresholded accepted components without the anatomical visible.  '
		 + 'For the rejected, middle kappa, and ignored components, the latter grey scale image is '
		 + 'displayed as well.\n')
	sl.append('Requirements for generating this report form:\n\n' + '* matplotlib\n\n' + '* numpy\n\n' + '* nibabel\n\n' + '* scipy\n\n' +'* sphinx\n\n' + '* parse\n')
	ofh = open("intro.rst","w")
	ofh.write("\n".join(sl)+"\n")
	ofh.close()


def correlation_rst():
	sl = []
	sl.append('Correlation\n' + '==========================\n')
	sl.append('In this section, several seeds were choosen to represent a few networks that are common in the literature.  ' +
	  'These seeds for the default were taken at MNI coordinates (0,-53,26) and (0,52,6) [1]_ which coorespond to Post cingulate ' +
	   'cortex and the medial prefrontal cortex respectively.  Z scores were capped at 5 and -5 for the below figures.\n')
	sl.append('.. [1] **Koene R. A. Van Dijk, Trey Hedden, Archana Venkataraman, Karleyton C. Evans, Sara W. Lazar and Randy L. Buckner.** ' +
		'Intrinsic Functional Connectivity As a Tool For Human Connectomics: Theory, Properties, and Optimization. *J Neurophysiol 103:297-321, 2010*\n')
	sl.append('Default mode\n' + '+++++++++++++++')
	sl.append('Below is our seed based correlation for several voxels thought to be part of the default mode.\n')
	sl.append('MNI (0,-53,26)\n' + '--------------')
	sl.append('.. image:: .. /png_dump/Pc_correlation.png')
	sl.append('	:scale: 75%\n')
	sl.append('MNI (0,52,6)\n' + '--------------')
	sl.append('.. image:: .. /png_dump/mPFC_correlation.png')
	sl.append('	:scale: 75%')
	ofh = open("correlation.rst","w")
	ofh.write("\n".join(sl)+"\n")
	ofh.close()
"""
Controls length of floats.  Needs to be done because spaceing important for sphinx.
This function is not designed to work if length is less than the number of digits before the decimal.
"""
def digit_length(n,length):
	n_string = str(n)
	if len(n_string) < length:
		n_string = n_string + ' '*(length-len(n_string))
	elif len(n_string)> length:
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

