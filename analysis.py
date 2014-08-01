#!/usr/bin/env python
"""
Gutierrez, B.  Generates the analysis.rst file
"""
import numpy as np

def analysis_rst(accept, reject, middle, ignore, threshold):
	size = len(str(accept.shape[0] + reject.shape[0] + middle.shape[0] + ignore.shape[0]))
	file = open("analysis.rst","w")
	file.write('Analysis\n' + '========\n')
	file.write('Graphs\n' + '++++++\n')
	file.write('.. image:: .. /png_dump/kappa_vs_rho.png\n')
	file.write('	:width: 50%\n\n')

	file.write('Accepted Components\n' + '+++++++++++++++++++\n')
	file.write('The following images are the thresholded components from the accepted bin of meica.py.  ')
	file.write('The threshold was set to %s.\n' % threshold)
	N = 0	
	for i in accept[:,0]:
		file.write('\nComponent %s\n' % int(i))
		file.write('----------'+ '-'*len(str(int(i))) + '\n\n')
		file.write('.. image:: ../png_dump/axial_component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n')
		file.write('.. image:: ../png_dump/sagital_component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n')
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
		file.write('.. image:: ../png_dump/axial_component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n')
		file.write('.. image:: ../png_dump/sagital_component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n')
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
		file.write('.. image:: ../png_dump/axial_component_' +(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n')
		file.write('.. image:: ../png_dump/sagital_component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n')
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
	for l in ignore[:,0]: 
		file.write('\nComponent %s\n' % int(i))
		file.write('----------' + '-'*len(str(int(i))) + '\n\n')
		file.write('.. image:: ../png_dump/axial_component_'  +(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n')
		file.write('.. image:: ../png_dump/sagital_component_'+(size - len(str(int(i))))*'0' + '%s.png\n' % int(i))
		file.write('	:scale: 75%\n')
		file.write('	:align: left\n\n')
		file.write('=============  =============  =============  =============\n')
		file.write('     kappa         Rho           %%Var        %%Var(norm)\n')
		file.write('=============  =============  =============  =============\n')
		file.write('%s       %s         %s           %s       \n' % 
			(digit_length(ignore[N,1],8),digit_length(ignore[N,2],7),digit_length(ignore[N,3],4),digit_length(ignore[N,4],4)))
		file.write('=============  =============  =============  =============\n')
		N = N + 1

def digit_length(n,length):
	n_string = str(n)
	if len(n_string) < length:
		n_string = n_string + ' '*(length-len(n_string))
	elif len(n_string)> length:
		n_string = round(n,length)
	return(n_string)



