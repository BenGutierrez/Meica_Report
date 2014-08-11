#!/usr/bin/env python
"""
Gutierrez, B.  Generates the analysis.rst file
"""
import numpy as np
from parse import parse

def analysis_rst(accept, reject, middle, ignore, threshold, ctab):
	size = len(str(accept.shape[0] + reject.shape[0] + middle.shape[0] + ignore.shape[0]))
	file = open("analysis.rst","w")
	file.write('Analysis\n' + '========\n')
	line = header(ctab)
	for i in range(len(line)):
		file.write('%s \n' % line[i][1:])
	file.write('\n')
	file.write('Graphs\n' + '++++++\n')
	file.write('.. image:: .. /png_dump/kappa_vs_rho.png\n')
	file.write('	:width: 49%\n')
	file.write('.. image:: .. /png_dump/kappa_rho_vs_components.png\n')
	file.write('	:width: 49%\n\n')
	file.write('=============  =============  =============  =============\n')
	file.write('  #Accepted      #Rejected    #Middle Kappa    #Ignored   \n')
	file.write('=============  =============  =============  =============\n')
	file.write('     %s             %s              %s             %s     \n' % 
		(len(accept), len(reject), len(middle), len(ignore)))
	file.write('=============  =============  =============  =============\n\n')
	file.write('Accepted Components\n' + '+++++++++++++++++++\n')
	file.write('The following images are the thresholded components from the accepted bin of meica.py.  ')
	file.write('The threshold was set to %s.\n' % threshold)
	N = 0	
	for i in accept[:,0]:
		file.write('\nComponent %s\n' % int(i))
		file.write('----------'+ '-'*len(str(int(i))) + '\n\n')
		file.write('.. image:: ../png_dump/Component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n\n')
		file.write('=============  =============  =============  =============\n')
		file.write('     kappa         Rho           %%Var        %%Var(norm)\n')
		file.write('=============  =============  =============  =============\n')
		file.write('%s       %s         %s           %s       \n' % 
			(digit_length(accept[N,1],8),digit_length(accept[N,2],7),digit_length(accept[N,3],3),digit_length(accept[N,4],3)))
		file.write('=============  =============  =============  =============\n')
		N = N + 1

	file.write('\nRejected Components\n' + '+++++++++++++++++++\n')
	file.write('The following images are the grey scale components from the rejected bin of meica.py.\n')
	N = 0
	for i in reject[:,0]:
		file.write('\nComponent %s\n' % int(i))
		file.write('----------' + '-'*len(str(int(i))) + '\n\n')
		file.write('.. image:: ../png_dump/Component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n\n')
		file.write('=============  =============  =============  =============\n')
		file.write('     kappa         Rho           %%Var        %%Var(norm)\n')
		file.write('=============  =============  =============  =============\n')
		file.write('%s       %s         %s           %s       \n' % 
			(digit_length(reject[N,1],8),digit_length(reject[N,2],7),digit_length(reject[N,3],3),digit_length(reject[N,4],3)))
		file.write('=============  =============  =============  =============\n')
		N = N + 1

	file.write('\nMiddle Components\n' + '+++++++++++++++++++\n')
	file.write('The following images are the grey scale components from the middle kappa bin of meica.py.\n')
	N = 0
	for i in middle[:,0]:
		file.write('\nComponent %s\n' % int(i))
		file.write('----------' + '-'*len(str(int(i))) + '\n\n')
		file.write('.. image:: ../png_dump/Component_' +(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n\n')
		file.write('=============  =============  =============  =============\n')
		file.write('     kappa         Rho           %%Var        %%Var(norm)\n')
		file.write('=============  =============  =============  =============\n')
		file.write('%s       %s         %s           %s       \n' % 
			(digit_length(middle[N,1],9),digit_length(middle[N,2],7),digit_length(middle[N,3],4),digit_length(middle[N,4],4)))
		file.write('=============  =============  =============  =============\n')
		N = N + 1

	file.write('\nIgnore Components\n' + '+++++++++++++++++++\n')
	file.write('The following images are the grey scale components from the ignore bin of meica.py.\n')
	N = 0
	for i in ignore[:,0]: 
		file.write('\nComponent %s\n' % int(i))
		file.write('----------' + '-'*len(str(int(i))) + '\n\n')
		file.write('.. image:: ../png_dump/Component_'  +(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n\n')
		file.write('=============  =============  =============  =============\n')
		file.write('     kappa         Rho           %%Var        %%Var(norm)\n')
		file.write('=============  =============  =============  =============\n')
		file.write('%s       %s         %s           %s       \n' % 
			(digit_length(ignore[N,1],8),digit_length(ignore[N,2],7),digit_length(ignore[N,3],4),digit_length(ignore[N,4],4)))
		file.write('=============  =============  =============  =============\n')
		N = N + 1

def diagnostics_rst():
	file = open("diagnostics.rst","w")
	file.write('Preliminary Diagnostics\n' + '==========================\n\n')
	file.write('TSNR\n' + '+++++++\n')
	file.write('Absolute TSNR of the medn NIFTI1 dataset:\n\n')
	file.write('.. figure:: .. /png_dump/medn_tsnr.png\n')
	file.write('	:width: 100%\n')
	file.write('	:align: left\n\n')
	file.write('Ratio of the TSNR of the medn NIFTI1 dataset to the TSNR of the tsoc NIFTI1 dataset:\n\n')
	file.write('.. figure:: .. /png_dump/tsnr_ratio.png\n')
	file.write('	:width: 100%\n')
	file.write('	:align: left\n\n')
	file.write('Histograms of the TSNR of the medn NIFTI1 dataset and the ratio of the TSNr of the medn dataset to the TSNR of the tsoc NIFTI1 dataset:\n\n')
	file.write('.. image:: .. /png_dump/medn_tsnr_hist.png\n')
	file.write('	:width: 49%\n')
	file.write('.. image:: .. /png_dump/tsnr_ratio_hist.png\n')
	file.write('	:width: 49%\n')




def index_rst():
	file = open("index.rst")
	file.write('Welcome to your ME-ICA Report!\n')
	file.write('==================================\n')
	file.write('This report collects information from your ME-ICA analysis and displays'
		+ 'several ways to interpret the output of meica.py.\n\n')
	file.write('Requirements: matplotlib, numpy, nibabel, scipy, sphinx, and parse\n\n')
	file.write('Contents:\n\n')
	file.write('.. toctree::\n')
	file.write('   :maxdepth: 2\n\n')
	file.write('   intro\n' + '   diagnostics\n' + '   analysis\n')
	file.write('\n\n\nSearch\n')
	file.write('======\n\n')
	file.write('* :ref:`search`')


def digit_length(n,length):
	n_string = str(n)
	if len(n_string) < length:
		n_string = n_string + ' '*(length-len(n_string))
	elif len(n_string)> length:
		n_string = round(n,length)
	return(n_string)

def header(ctab):
	txt_file = open(str(ctab))
	txt_file.seek(0)
	line = np.zeros(11, dtype = object)
	for i in range(2):
		txt_file.readline()
	for i in range(11):
		line[i] = txt_file.readline()
	return(line)

