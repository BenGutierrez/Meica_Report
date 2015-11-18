#!/usr/bin/env python
__version__="1.0"

"""
Gutierrez, B.  Generates overview Meica report form for many me-ica reports.

Introduction:
    By putting in many of the txt files generated inside of Meica.Report directories its is possible to quickly
    grab many possibly interesting statistics and place them into a single csv file for analysis.  This is the main function
    of this script. However for "fun" I have also used sphinx to display the Kappa vs Rho plots of each me-ica run on an html
    page along with the TSNR figures and the ability to create simple x vs y plots along with trendlines of all of the statistics
    in the csv file.
Usage:
    This script must be run after running meica_report.py.  A typical call of this script may look something like the following:

        python ~/total_report.py -pattern_1 SBJ01_S01/Meica_Analysis/Task0?/meica.Report/meica_report.txt

    By giving patterns this program can grab many datasets of similar organizational structure quickly.  In the above example the 
    wild card ? will allow any directories of the form Task0_ to be searched for the file at the end of the path 
    meica.Report/meica_report.txt

    In the end, the csv file meica_report.csv will be in the directory meica.Meta_Report unless another name for this directory is specified.
    The html version of the report can be opened by using firefox/meica.Meta_Report/html/index.html

Options:
    -h,          --help        Provide a help menu with all of the options.
    -pattern_1               First path to meica_report.txt files(the options -patter_2,...,-pattern_13 can also be used but are hidden).
    -dest                    Path to place the directory meica.Meta_Report into.
    -label                     Alternate name for the directory meica.Meta_Report directory for all outputs. Default is meica.Meta_Report
    -plots                     This argument is rather strange in order to be more versatile.  It takes in the names of statistics and 
                            makes x vs y plots for example, "Variance Explained by ICA vs Rejected Variance, # Accepted Components vs
                            Variance Explained by ICA" will give the two x vs y figures described in the example argument above. Default
                            is 'Variance Explained by ICA vs Median Meica Denoised TSNR, Variance Explained by ICA vs Median Optimally 
                            Combined TSNR, Variance Explained by ICA vs Median Denoised Over Optimally Combined TSNR Ratio, Accepted 
                            Variance vs Variance Explained by ICA'.
    -view_component            Will display given component number across all patterns.
    -var_component             Will display the highest variance accepted component if specified.


"""
from matplotlib.colors                 import LinearSegmentedColormap
from mpl_toolkits.axes_grid1         import ImageGrid

import matplotlib.pyplot     as plt
import scipy.stats             as stats
import numpy                 as np

import commands
import subprocess
import argparse
import glob
import sys
import ast
import csv
import os

"""
Plotting function designed to be able to easily accept arguments to make plots with a variety of source data, also adds a trend line.
Argument names should be straight forward.
"""
def plot_scaffold(x_data,y_data,title,xlab,ylab,file_name):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_data,y_data)
    ymin=slope*min(x_data)+intercept
    ymax=slope*max(x_data)+intercept
    plt.figure()
    line=plt.plot([min(x_data),max(x_data)],[ymin,ymax],ls='--',c='k',label=r'$y=%sx+%s,R^2=%s$' % (round(slope,2),round(intercept,2),round(r_value**2,3)))
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.legend(loc = 'upper right')
    scatter=plt.scatter(x_data,y_data)
    plt.savefig(file_name,bbox_inches='tight')
    plt.close()

"""
make index.rst file that serves as the table of contents for the sphinx html page generation.
"""
def index_rst(num_comp,var_comp):
    sl = []
    title = 'Total Report'
    sl.append(title)
    sl.append('='* (len(title)+1))
    
    sl.append('Contents:\n')
    sl.append('.. toctree::')
    sl.append('   :glob:   ')
    sl.append('   :maxdepth: 1\n')
    sl.append('   Kappa*')
    sl.append('   General*')
    sl.append('   TSNR*')
    if num_comp != False:
        sl.append('   Number*') 
    if var_comp != False:
        sl.append('   Variance*')
    sl.append('\n\n\nSearch')
    sl.append('======\n')
    sl.append('* :ref:`search`')

    ofh = open("index.rst","w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()

"""
This function creates the csv file from all of the input txt files.  There is quite a bit of semi-intelligent labeling down from the paths provided.
I make no claims that this self labeling works in all input methods.  Hopefully it is stable however and it was done for convienience not necessity.
"""
def table(pattern):
    variables = len(pattern)*[0]
    start1=''
    start0=''

    # Fill list with all of the data taken from input files
    # -----------------------------------------------------
    for k in range(len(pattern)):
        variables[k] = len(pattern[k])*[0]
        for i in range(len(pattern[k])):
            variables[k][i]= [0]*21
        for i in range(len(pattern[k])):
            meica_txt=open(pattern[k][i],'r')
            for j in np.arange(len(variables[k][i])-2)+2:
                temp = meica_txt.readline()
                if temp[temp.index(':')+1:-1]!=len(temp[temp.index(':')+1:-1])*' ':
                    variables[k][i][j]=temp[temp.index(':')+1:-1]
                else:
                    variables[k][i][j]='   NA'

            # Make the first label based on the differences inside an input pattern
            # ---------------------------------------------------------------------
            if len(pattern[k])>1:
                if start1=='':
                    if i != 0:
                        differences=[j for j in range(len(pattern[k][0])) if pattern[k][i][j] != pattern[k][0][j]]
                    else:
                        differences=[j for j in range(len(pattern[k][1])) if pattern[k][i][j] != pattern[k][1][j]]
                    start1=differences[0]
                    end1=differences[0]
                    while pattern[k][i][start1]!='/' and start1!=0:
                        start1 = start1-1
                    while pattern[k][i][end1]!='/' and end1!= len(pattern[k][i])-1:
                        end1 = end1+1
                variables[k][i][1]=pattern[k][i][start1+1:end1].replace('_',' ')
            else:
                variables[k][i][1]=''

        # Make the second label based on the differences across input patterns
        # --------------------------------------------------------------------
        if len(pattern)>1:
            if start0 == '':
                if k != 0:
                    differences=[j for j in range(len(pattern[k][0])) if pattern[k][0][j] != pattern[0][0][j]]
                else:
                    differences=[j for j in range(len(pattern[k][0])) if pattern[k][0][j] != pattern[1][0][j]]
                start0=differences[0]
                end0=differences[0]
                while pattern[k][0][start0]!='/' and start0!=0:
                    start0 = start0-1
                while pattern[k][0][end0]!='/' and end0!= len(pattern[k][0])-1:
                    end0 = end0+1
            for i in range(len(pattern[k])):
                variables[k][i][0]=pattern[k][i][start0+1:end0].replace('_',' ')
        else:
            for i in range(len(pattern[k])):
                variables[k][i][0]=''

    # write csv file and column labels
    # --------------------------------
    with open('meica_report.csv', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerows([['Label','Experimental label','# Accepted Components', '# Rejected Compoents','# Middle Components','# Ignore Components',
            'Accepted Variance','Rejected Variance','Middle Kappa Variance', 'Ignore Variance','Variance Explained by ICA','Largest Variance Accepted Component',
            'Highest Kappa Rejected Component','Highest Kappa Middle Kappa Component','Median Meica Denoised TSNR','Median Optimally Combined TSNR',
            'Median Denoised Over Optimally Combined TSNR Ratio','Max Head Displacement In Any One Dirrection','Time of Max Head Displacement (TR)',
            'Max Rate of Head Motion','Time of Max Head Motion Rate (TR)']])
        for i in range(len(variables)):
            writer.writerows(variables[i])

    return variables

"""
Function to create the file Kappa_plots.rst which sphinx uses to create the html page for containing all
Kappa vs Rho plots.
"""
def kappa_rho_plot(variables,pattern):
    sl = []
    sl.append('Kappa vs Rho plots')
    sl.append('==================')
    N=0
    for k in range(len(variables)):
        title = '%s' % (variables[k][0][0])
        sl.append(title)
        sl.append('+'*len(title))
        for i in range(len(variables[k])):
            N_ = str(N)
            while len(N_) < len(str(num_report)):
                N_ = '0' + N_

            title = '%s, %s' % (variables[k][i][0],variables[k][i][1])
            sl.append(title)
            sl.append('-'*len(title))
            sl.append('\n.. image:: ./Report_Figures/kappa_vs_rho%s.png' % N_)
            sl.append('\t:width: 49%\n')
            sl.append('.. image:: ./Report_Figures/kappa_rho_vs_components%s.png' % N_)
            sl.append('\t:width: 49%\n')
            N += 1

    ofh = open("Kappa_plots.rst","w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()

"""
Function to create the file General_plots.rst which sphinx uses to create an html page containing all of 
the plots created by the function plot_scaffold
"""
def generic_plots(plot_names):
    sl = []
    sl.append('Figures spanning multiple meica runs')
    sl.append('====================================')
    N=0
    for k in range(len(plot_names)):
        title = plot_names[k].replace('_',' ')[:-4]
        sl.append(title)
        sl.append('-'*len(title))
        sl.append('\n.. image:: ./Report_Figures/%s.png' % plot_names[k])
        sl.append('\t:width: 49%\n')

    ofh = open("General_plots.rst","w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()

"""
Function to create the file Number_comp.rst which sphinx uses to create an html page containing all of the 
user specified components for all input patterns
"""
def num_comp_rst(num_comp,variables,num_report):
    sl = []
    title = 'Component number %s' % num_comp
    sl.append(title)
    sl.append('='*len(title))
    N=0
    for k in range(len(variables)):
        for i in range(len(variables[k])):
            N_ = str(N)
            while len(N_) < len(str(num_report)):
                N_ = '0' + N_

            title = '%s , %s' % (variables[k][i][0],variables[k][i][1])
            sl.append(title)
            sl.append('-'*len(title))
            sl.append('\n.. image:: ./Report_Figures/Num_Component_%s_%s.png' % (num_comp,N_))
            sl.append('\t:width: 95%')
            sl.append('\t:align: left\n')
            N+= 1

    ofh = open("Number_comp.rst","w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()

"""
Function to create the file Variance_comp.rst which sphinx uses to create an html page containing the 
highest variance accepted component for all input patterns
"""
def var_comp_rst(variables,num_report):
    sl = []
    title = 'Highest Variance Accepted Components'
    sl.append(title)
    sl.append('='*len(title))
    N=0
    for k in range(len(variables)):
        for i in range(len(variables[k])):
            N_ = str(N)
            while len(N_) < len(str(num_report)):
                N_ = '0' + N_
            var_comp = str(int(variables[k][i][11]))
            while len(var_comp) < len(str(len(glob.glob('%s/Report_Figures/Axial_GS_Component_*' % os.path.dirname(os.path.expanduser(pattern[k][i])))))):
                var_comp = '0' + var_comp
            title = '%s, %s, Component %s' % (variables[k][i][0],variables[k][i][1],var_comp)
            sl.append(title)
            sl.append('-'*len(title))
            sl.append('\n.. image:: ./Report_Figures/Var_Component_%s_%s.png' % (var_comp,N_))
            sl.append('\t:width: 95%')
            sl.append('\t:align: left\n')
            N+=1
    ofh = open("Variance_comp.rst","w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()

"""
Function to create the TSNR figures from all patterns.  These images cannot be so simply pulled from
the original runs since the contrast need to be normalized across all of the runs.  Also creates the file
TSNR.rst for sphink to use for an HTML page showing the TSNR images.
"""
def TSNR(pattern,variables):
    # Create custom color map
    # -----------------------
    cdict = {'red':     ((0.0,  0.0, 0.0),
                 (0.2, 0.6, 0.6),
                 (0.35, 0.9, 0.9),
                 (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),

         'green':   ((0.0, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (0.65, 0.9, 0.9),
                     (0.9, 0.5, 0.5),
                     (1.0, 0.0, 0.0)),

          'blue':    ((0.0, 0.0, 0.0),
                      (1.0, 0.0, 0.0))}

    GYR = LinearSegmentedColormap('GYR', cdict)
    sl=[]
    sl.append('TSNR Figures')
    sl.append('============')
    N=0
    fig_num = 0
    for k in range(len(pattern)):
        fig_num = fig_num + len(pattern[k])

    print "++ INFO: Creating TSNR figures ..."
    for k in range(len(pattern)):
        maxs=[]
        mins=[]

        # Take the max and min TSNR values from each me-ica run
        # -----------------------------------------------------
        for i in range(len(pattern[k])):
            mins.append(np.loadtxt(os.path.dirname(os.path.expanduser(pattern[k][i]))+ '/axialized_nifti/tsnr_thresholds.txt')[0])
            maxs.append(np.loadtxt(os.path.dirname(os.path.expanduser(pattern[k][i]))+ '/axialized_nifti/tsnr_thresholds.txt')[1])
        mean_max=np.mean(maxs)
        mean_min=np.mean(mins)
        for i in range(len(pattern[k])):
            medn_tsnr = np.load(os.path.dirname(os.path.expanduser(pattern[k][i])) + '/axialized_nifti/medn_tsnr.npy')
            tsoc_tsnr = np.load(os.path.dirname(os.path.expanduser(pattern[k][i])) + '/axialized_nifti/tsoc_tsnr.npy')
            background = np.zeros(medn_tsnr[:,:,0].shape)
            fig = plt.figure(figsize = (12,2),dpi=150)

            # create the TSNR images.  ImageGrid very nicely arranges images, but is very slow
            # --------------------------------------------------------------------------------
            grid = ImageGrid(fig, 111, nrows_ncols=(2,10),cbar_location='right',add_all=True,cbar_pad=0.05,axes_pad=0.0,share_all=True,cbar_mode="single")
            for j in range(10):
                grid[j].imshow(background, cmap = 'Greys_r')
                plot = grid[j].imshow(medn_tsnr[:,:,j],vmin=mean_min,vmax=mean_max,cmap=GYR)
                grid[j].axes.get_xaxis().set_ticks([])
                grid[j].axes.get_yaxis().set_ticks([])
                grid[j+10].imshow(background, cmap = 'Greys_r')
                grid[j+10].imshow(tsoc_tsnr[:,:,j],vmin=mean_min,vmax=mean_max,cmap=GYR)
                grid[j+10].axes.get_xaxis().set_ticks([])
                grid[j+10].axes.get_yaxis().set_ticks([])
            grid[0].set_ylabel('medn TSNR',fontsize=9)
            grid[10].set_ylabel('tsoc TSNR',fontsize=9)
            grid.cbar_axes[0].colorbar(plot)
            grid.cbar_axes[0].tick_params(labelsize=9) 
            plt.savefig('./Report_Figures/tsnr_figure_%s.png' % N)
            plt.close()

            title = variables[k][i][0]+ ', ' + variables[k][i][1]
            sl.append(title)
            sl.append('-'*len(title))
            sl.append('\n.. image:: ./Report_Figures/tsnr_figure_%s.png' % N)
            sl.append('\t:width: 95%')
            sl.append('\t:align: left\n')
            N+=1
            print '\t TSNR figure %s out of %s created'  % (N,fig_num)

    ofh = open("TSNR.rst","w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()

if __name__=='__main__':

    # Parse input arguments
    # ---------------------
    parser = argparse.ArgumentParser('Options')
    Required = parser.add_argument_group('Required arguments')
    Required.add_argument('-pattern_1', dest = 'patt1', help = "Pathway pattern from current directory (or gloabl) to find meica report forms. Runs "
        +    "from a single pattern will be grouped together")
    lab = parser.add_argument_group('File and directory labels')
    lab.add_argument('-dest', dest = 'dest', help = 'Directory to place final report form directory, default is "."', default = '.')
    lab.add_argument('-label', dest = 'label', help = 'Label to call final report form directory, default is "meica.Meta_Report"', default = 'meica.Meta_Report')
    opt = parser.add_argument_group('Optinal arguments')
    opt.add_argument('-pattern_2', dest = 'patt2', help = 'Pathway pattern from current directory to find meica report forms', default = None)
    opt.add_argument('-plots', dest='plots', help = "This argument is rather strange in order to be more versatile.  It takes in the names of statistics and "
        +   "makes x vs y plots for example, 'Variance Explained by ICA vs Rejected Variance, # Accepted Components vs Variance Explained by ICA' will give "
        +    "the two x vs y figures described in the example argument above. Default is 'Variance Explained by ICA vs Median Meica Denoised TSNR, Variance "
        +    "Explained by ICA vs Median Optimally Combined TSNR, Variance Explained by ICA vs Median Denoised Over Optimally Combined TSNR Ratio, Accepted "
        +     "Variance vs Variance Explained by ICA.         options: # Accepted Components, # Rejected Compoents, # Middle Components, # Ignore Components,Accepted "
        +    "Variance, Rejected Variance, Middle Kappa Variance, Ignore Variance, Variance Explained by ICA, Largest Variance Accepted Component, Highest "
        +    "Kappa Rejected Component, Highest Kappa Middle Kappa Component, Median Meica Denoised TSNR, Median Optimally Combined TSNR, Median Denoised "
        +    "Over Optimally Combined TSNR Ratio", default = "Variance Explained by ICA vs Median Meica Denoised TSNR, Variance Explained by ICA vs Median "
        +    "Optimally Combined TSNR, Variance Explained by ICA vs Median Denoised Over Optimally Combined TSNR Ratio, Accepted Variance vs Variance Explained by ICA")
    opt.add_argument('-view_component', dest = 'num_comp', help = 'If given a number, will show that component from each meica report', default=False)
    opt.add_argument('-var_component', dest = 'var_comp', help = 'If specified, will show the highest variance accepted component from each report', action = 'store_true')
    opt.add_argument('-high_var', dest = 'high_var', help = 'Displays the highest variance component from each report', action = 'store_true')
    opt.add_argument('-pattern_3', dest = 'patt3', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_4', dest = 'patt4', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_5', dest = 'patt5', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_6', dest = 'patt6', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_7', dest = 'patt7', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_8', dest = 'patt8', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_9', dest = 'patt9', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_10', dest = 'patt10', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_11', dest = 'patt11', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_12', dest = 'patt12', default = None, help=argparse.SUPPRESS)
    opt.add_argument('-pattern_13', dest = 'patt13', default = None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    print "-- ME Meta Report Utility version %s --" % __version__

    # Make sure required inputs are provided
    # --------------------------------------
    if args.patt1 == None:
        print "*+ no pattern_1 given.  exiting"
        sys.exit()

    # Use module glob to use wildcards to expand given patterns
    # ---------------------------------------------------------
    patt=[args.patt1,args.patt2,args.patt3,args.patt4,args.patt5,args.patt6,args.patt7,args.patt8,args.patt9,args.patt10,args.patt11,args.patt12,args.patt13]
    pattern=[]
    pattern_labels = []
    for i in range(len(patt)):
        if patt[i] != None:
            pattern.append(patt[i])

    for i in range(len(pattern)):#make groups of patterns and expand paths to absolute paths for ease
        pattern[i] = glob.glob(os.path.expanduser(pattern[i]))
        for j in range(len(pattern[i])):
            pattern[i][j] = os.path.abspath(pattern[i][j])
        pattern[i][:]=sorted(pattern[i][:])

    # Make output directory
    # ---------------------
    os.chdir(args.dest)
    if os.path.isdir('%s/%s' % (args.dest,args.label)):
        subprocess.call('rm -rf %s/%s' % (args.dest,args.label),shell=True)
    subprocess.call('mkdir %s/%s' % (args.dest,args.label),shell=True)
    os.chdir(args.dest + '/' + args.label)
    subprocess.call('mkdir Report_Figures',shell=True)
    subprocess.call('mkdir Report_data',shell=True)
    
    print "++ INFO: Creating meica_report.csv"
    variables = table(pattern)

    # copy images over from me_ica runs to meica.Meta_Report directory
    # ----------------------------------------------------------------
    print "++ INFO: Copying figures into %s/Report_Figures" % args.dest
    N=0
    num_report=0
    for k in range(len(pattern)):
            num_report = num_report + len(pattern[k])
    num_comp=args.num_comp
    if args.num_comp != False:
        while len(str(num_comp)) < len(str(num_report)):
                num_comp = '0' + num_comp

    for k in range(len(variables)):
        for i in range(len(variables[k])):
            N_ = str(N)
            while len(N_) < len(str(num_report)):
                N_ = '0' + N_

            title = '%s, %s' % (variables[k][i][0],variables[k][i][1])
            subprocess.call('cp %s/Report_Figures/kappa_vs_rho.png ./Report_Figures/kappa_vs_rho%s.png' % (os.path.dirname(os.path.expanduser(pattern[k][i])),N_),shell=True)
            subprocess.call('cp %s/Report_Figures/kappa_rho_vs_components.png ./Report_Figures/kappa_rho_vs_components%s.png' % 
                (os.path.dirname(os.path.expanduser(pattern[k][i])),N_),shell=True)

            if args.num_comp:
                subprocess.call('cp  %s/Report_Figures/Axial_GS_Component_%s.png ./Report_Figures/Num_Component_%s_%s.png' % 
                    (os.path.dirname(os.path.expanduser(pattern[k][i])),num_comp,num_comp,N_),shell=True)

            if args.var_comp:
                var_comp = str(int(variables[k][i][11]))
                while len(var_comp)< len(str(len(glob.glob('%s/Report_Figures/Axial_GS_Component_*.png' % os.path.dirname(os.path.expanduser(pattern[k][i])))))):
                    var_comp = '0' + var_comp

                subprocess.call('cp  %s/Report_Figures/Axial_GS_Component_%s.png ./Report_Figures/Var_Component_%s_%s.png' % (os.path.dirname(os.path.expanduser(pattern[k][i])),
                    var_comp,var_comp,N_),shell=True)

            N += 1
    # Define where in the csv table the variables are
    # -----------------------------------------------
    print "++ INFO: Creating figures from -plot argument"
    plot_names=[]
    plot_labels = [['# Accepted Components',2],['# Rejected Compoents',3],['# Middle Components',4],['# Ignore Components',5],['Accepted Variance',6],
        ['Rejected Variance',7],['Middle Kappa Variance',8], ['Ignore Variance',9],['Variance Explained by ICA',10],['Largest Variance Accepted Component',11],
        ['Highest Kappa Rejected Component',12],['Highest Kappa Middle Kappa Component',13],['Median Meica Denoised TSNR',14],['Median Optimally Combined TSNR',15],
        ['Median Denoised Over Optimally Combined TSNR Ratio',16]]
    all_plots=args.plots

    # Goes through the argument plots and creates x vs y plots from csv file data
    # ---------------------------------------------------------------------------
    for j in range(len(plot_labels)):
        while ('%s vs ' % plot_labels[j][0]) in all_plots:
            x_data=[]
            y_data=[]
            for l in range(len(plot_labels)):
                #print '%s,%s' % (j,l)
                if ('%s vs %s' % (plot_labels[j][0],plot_labels[l][0][:15])) in all_plots:
                    for k in range(len(variables)):
                        for i in range(len(variables[k])):
                            x_data.append(float(variables[k][i][plot_labels[j][1]]))
                            y_data.append(float(variables[k][i][plot_labels[l][1]]))

                    plot_names.append(('%s vs %s' % (plot_labels[j][0],plot_labels[l][0])).replace(' ','_'))
                    all_plots=all_plots.replace('%s vs %s' % (plot_labels[j][0],plot_labels[l][0]),'')
                    plot_scaffold(x_data,y_data,'%s vs %s' % (plot_labels[j][0],plot_labels[l][0]),plot_labels[j][0],plot_labels[l][0],('%s vs %s' % 
                        (plot_labels[j][0],plot_labels[l][0])).replace(' ','_')+'.png')
                    break

    # Create .rst files
    # -----------------
    print "++ INFO: Creating .rst files"
    generic_plots(plot_names)
    if args.num_comp != False:
        num_comp_rst(num_comp,variables,num_report)
    if args.var_comp != False:
        var_comp_rst(variables,num_report)
    kappa_rho_plot(variables,pattern)
    index_rst(args.num_comp,args.var_comp)
    TSNR(pattern,variables)

    # Steals files necessary for sphinx from meica_report.py generated directory
    # --------------------------------------------------------------------------
    subprocess.call('mv *.png Report_Figures', shell = True)
    subprocess.call('cp %s/sphinx_files/make.bat ./' % os.path.dirname(pattern[0][0]),shell=True)
    subprocess.call('cp %s/sphinx_files/Makefile ./' % os.path.dirname(pattern[0][0]),shell=True)
    subprocess.call('cp %s/sphinx_files/conf.py ./' % os.path.dirname(pattern[0][0]),shell=True)
    subprocess.call('cp -r %s/html/_static ./' % os.path.dirname(pattern[0][0]),shell=True)
    subprocess.call('make html', shell = True)

    # Sphinx
    # ------
    subprocess.call('mv _build/* ./', shell = True)
    subprocess.call('rm -rf _*', shell = True)
    subprocess.call('mkdir sphinx_files', shell = True)
    subprocess.call('mv *.rst sphinx_files/', shell = True)
    subprocess.call('mv Makefile sphinx_files', shell = True)
    subprocess.call('mv conf.py sphinx_files', shell = True)
