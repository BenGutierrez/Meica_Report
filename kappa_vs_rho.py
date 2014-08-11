#!/usr/bin/env python
"""
Gutierrez, B.  Parses ctab txt file created by the program meica.py and 
plots kappa vs rho while distinguishing the bins the components lie in.
The percent variance is 
"""

"""
Parse the component file from meica.py and return a list of lists containing the
component numbers from the accepted, rejected, middle, and ignore bins.
"""
from parse import parse
import numpy as np
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

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

"""
plot kappa vs rho and represent %varaince by the size of the markers
"""
def kappa_vs_rho_plot(accept,reject,middle,ignore):
	plt.figure(2)
	trial_1 = plt.scatter(1,1,c = 'r',marker = 'o')
	trial_2 = plt.scatter(1,1,c = 'b',marker = '^')
	trial_3 = plt.scatter(1,1,c = 'g',marker = 'v')
	trial_4 = plt.scatter(1,1,c = 'c',marker = '*')
	plt.close(2)
	plt.title('ME-ICA Analysis, ' + r'$\kappa$' + ' vs ' + r'$\rho$',fontsize = 14)
	ACC = plt.scatter(accept[:,1], accept[:,2],c = 'r',marker = 'o', s = 50*accept[:,4]) 
	REJ = plt.scatter(reject[:,1], reject[:,2],c = 'b',marker = '^', s = 50*reject[:,4])
	MID = plt.scatter(middle[:,1], middle[:,2],c = 'g',marker = 'v', s = 50*middle[:,4])
	IGN = plt.scatter(ignore[:,1], ignore[:,2],c = 'c',marker = '*', s = 50*ignore[:,4])
	plt.legend((trial_1,trial_2,trial_3,trial_4),('Accepted','Rejected','Middle',
		'Ignore'),scatterpoints = 1, loc = 'upper right',markerscale = 2)
	plt.gca().xaxis.set_major_locator( MaxNLocator(nbins = 5))
	plt.tick_params(axis = 'x',which = 'both',top = 'off')
	plt.tick_params(axis = 'y',which = 'both',right = 'off')
	plt.xlabel(r'$\kappa$',fontsize = 15)
	plt.ylabel(r'$\rho$',fontsize = 15)

def kr_vs_component(comp_table_title):
	components = np.loadtxt(str(comp_table_title))
	plt.figure()
	plt.title('ME-ICA Analysis, ' + r'$\kappa$' + ' and ' + r'$\rho$' + ' vs Component Rank',fontsize = 14)
	plt.xlabel(r'$\kappa$' ', ' + r'$\rho$' ,fontsize = 15)
	plt.xlabel('Component Rank' ,fontsize = 15)
	kappa = plt.plot(components[:,0],components[:,1])
	rho = plt.plot(components[:,0],components[:,2])
	plt.legend((r'$\kappa$',r'$\rho$'))





