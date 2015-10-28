#!/usr/bin/env python
"""
Gutierrez, B.  Generates the rst files that sphinx reads to make the report form
"""
import numpy as np
import os

"""
make analysis.rst file
"""
def dynamic_analysis_rst(accept, reject, middle, ignore, ctab, outputDir, motion_file,TED_dir):
    figures = 'Report_Figures'
    sl = []
    size = len(str(accept.shape[0] + reject.shape[0] + middle.shape[0] + ignore.shape[0]))
    sl.append('Dynamic Analysis')
    sl.append('=======================')
    with open("%s/comp_table.txt" % (TED_dir), 'r') as original: ctab_txt = original.read()
    ctab_txt = ctab_txt.split('\n')
    for i in range(len(ctab_txt)-1):
        if '#' in ctab_txt[i]:
            sl.append(ctab_txt[i] + '\n')
    sl.append('\nPlease remember that the denoised time series is accepted AND ignored bin and components in neither rejected nor middle kappa are considered BOLD-like')
    sl.append('\n+----------------+------------------+-------------------------+\n' +
          '|                | %  Total Variance| %  Total Variance(norm) |\n' + 
          '+================+==================+=========================+')
    sl.append('| **Accepted**   |        %s    |          %s         |' % (digit_length(sum(ctab[accept,3]),6), digit_length(sum(ctab[accept,4]),6)))
    sl.append('+----------------+------------------+-------------------------+')
    sl.append('| **Rejected**   |        %s    |          %s         |' % (digit_length(sum(ctab[reject,3]),6), digit_length(sum(ctab[reject,4]),6)))
    sl.append('+----------------+------------------+-------------------------+')
    sl.append('|**Middle kappa**|        %s    |          %s         |' % (digit_length(sum(ctab[middle,3]),6), digit_length(sum(ctab[middle,4]),6)))
    sl.append('+----------------+------------------+-------------------------+')
    sl.append('| **Ignored**    |        %s    |          %s         |' % (digit_length(sum(ctab[ignore,3]),6), digit_length(sum(ctab[ignore,4]),6)))
    sl.append('+----------------+------------------+-------------------------+\n')

    sl.append('In the scatter plots below, you can click on any point to see the component time series.  ' + 
            'The size of the scatter points in the plot furthest to the right coorespond to the variance of the component.  ' +
            'The green components are accepted, the red are everything else.\n')
    sl.append('.. bokeh-plot:: %s/bokeh_plot.py' % figures)
    sl.append('\t:source-position: none \n')

    ofh = open("%s/Dynamic.rst" % outputDir,"w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()
def analysis_rst(accept, reject, middle, ignore, ctab, outputDir, motion_file, TED_dir):
    figures =  'Report_Figures'
    sl = []
    size = len(str(ctab.shape[0]))
    sl.append('Component Visualization')
    sl.append('=======================')
    with open("%s/comp_table.txt" % (TED_dir), 'r') as original: ctab_txt = original.read()
    ctab_txt = ctab_txt.split('\n')
    for i in range(len(ctab_txt)-1):
        if '#' in ctab_txt[i]:
            sl.append(ctab_txt[i] + '\n')
    sl.append('\nPlease remember that the denoised time series is accepted AND ignored bin and components in neither rejected nor middle kappa are considered BOLD-like')
    sl.append('\n+----------------+------------------+-------------------------+\n' +
          '|                | %  Total Variance| %  Total Variance(norm) |\n' + 
          '+================+==================+=========================+')
    sl.append('| **Accepted**   |        %s    |          %s         |' % (digit_length(sum(ctab[accept,3]),6), digit_length(sum(ctab[accept,4]),6)))
    sl.append('+----------------+------------------+-------------------------+')
    sl.append('| **Rejected**   |        %s    |          %s         |' % (digit_length(sum(ctab[reject,3]),6), digit_length(sum(ctab[reject,4]),6)))
    sl.append('+----------------+------------------+-------------------------+')
    sl.append('|**Middle kappa**|        %s    |          %s         |' % (digit_length(sum(ctab[middle,3]),6), digit_length(sum(ctab[middle,4]),6)))
    sl.append('+----------------+------------------+-------------------------+')
    sl.append('| **Ignored**    |        %s    |          %s         |' % (digit_length(sum(ctab[ignore,3]),6), digit_length(sum(ctab[ignore,4]),6)))
    sl.append('+----------------+------------------+-------------------------+\n')

    sl.append('Graphs')
    sl.append('++++++\n')
    #sl.append('.. bokeh-plot:: %s/bokeh_plot.py' % figures)
    #sl.append('\t:source-position: none \n')
    sl.append('.. image:: %s/kappa_vs_rho.png' % figures)
    sl.append('\t:width: 49%')
    sl.append('.. image:: %s/kappa_rho_vs_components.png' % figures)
    sl.append('\t:width: 49%\n')
    sl.append('The size of the scatter plot points is linearly related to the percent variance of that particular component.\n')
    if os.path.isfile('%s' % (motion_file)):
        sl.append('\nThe below figures show subject motion.  The plot on the right is the derivative of the 6 plots on the left with the euclidean norm applied in order to combine the data:\n')
        sl.append('.. image:: %s/motion_plot.png' % figures)
        sl.append('\t:width: 49%')
        sl.append('.. image:: %s/motion_rate.png' % figures)
        sl.append('\t:width: 49%\n')
    sl.append('=============  =============  =================  =============')
    sl.append('# of Accepted  # of Rejected  # of Middle kappa  # of Ignored ')
    sl.append('=============  =============  =================  =============')
    sl.append('     %s             %s              %s                 %s    ' % 
        (len(accept), len(reject), len(middle), len(ignore)))
    sl.append('=============  =============  =================  =============\n\n')
    sl.append('\nAccepted Components\n' + '+++++++++++++++++++')
    sl.append('The following images are the accepted components from meica ouput')
    
    N = 0
    for i in accept:
        sl.append('\nComponent %s: Var %s' % (int(i),digit_length(ctab[accept[N],3],2)))
        sl.append('-----------------------------\n')
        sl.append('=============  =============  =============  =================')
        sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
        sl.append('=============  =============  =============  =================')
        sl.append('%s       %s         %s           %s       ' % 
            (digit_length(ctab[accept[N],1],8), digit_length(ctab[accept[N],2],7), digit_length(ctab[accept[N],3],3), digit_length(ctab[accept[N],4],3)))
        sl.append('=============  =============  =============  =================\n')
        sl.append('.. image:: %s/Axial_GS_Component_'    % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Sagittal_GS_Component_' % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Coronal_GS_Component_'  % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('\n.. image:: %s/TimeSeries_'          % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        sl.append('.. image:: %s/FFT_Component_'         % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        N += 1

    sl.append('\nRejected Components\n' + '+++++++++++++++++++')
    sl.append('The following images are the rejected components from meica ouput')
    N = 0
    for i in reject:
        sl.append('\nComponent %s: Var %s' % (int(i),digit_length(ctab[reject[N],3],2)))
        sl.append('-----------------------------\n')
        sl.append('=============  =============  =============  =================')
        sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
        sl.append('=============  =============  =============  =================')
        sl.append('%s       %s         %s           %s       ' % 
            (digit_length(ctab[reject[N],1],8), digit_length(ctab[reject[N],2],7), digit_length(ctab[reject[N],3],3), digit_length(ctab[reject[N],4],3)))
        sl.append('=============  =============  =============  =================\n')
        sl.append('.. image:: %s/Axial_GS_Component_'    % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Sagittal_GS_Component_' % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Coronal_GS_Component_'  % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('\n.. image:: %s/TimeSeries_'          % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        sl.append('.. image:: %s/FFT_Component_'         % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        N += 1

    sl.append('\nMiddle Components\n' + '+++++++++++++++++++')
    sl.append('The following images are the middle kappa components from meica ouput')
    N = 0
    for i in middle:
        sl.append('\nComponent %s: Var %s' % (int(i),digit_length(ctab[middle[N],3],2)))
        sl.append('-----------------------------\n')
        sl.append('=============  =============  =============  =================')
        sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
        sl.append('=============  =============  =============  =================')
        sl.append('%s       %s         %s           %s       ' % 
            (digit_length(ctab[middle[N],1],8), digit_length(ctab[middle[N],2],7), digit_length(ctab[middle[N],3],4), digit_length(ctab[middle[N],4],4)))
        sl.append('=============  =============  =============  =================\n')
        sl.append('.. image:: %s/Axial_GS_Component_'    % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Sagittal_GS_Component_' % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Coronal_GS_Component_'  % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('\n.. image:: %s/TimeSeries_'          % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        sl.append('.. image:: %s/FFT_Component_'         % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        N += 1

    sl.append('\nIgnore Components\n' + '+++++++++++++++++++')
    sl.append('The following images are the ignored components from meica ouput.  These are kept with in the denoised timeseries for increased variance')
    N = 0
    for i in ignore: 
        sl.append('\nComponent %s: Var %s' % (int(i),digit_length(ctab[ignore[N],3],2)))
        sl.append('-----------------------------\n')
        sl.append('=============  =============  =============  =================')
        sl.append('     kappa         rho         %s Variance    %s Variance(norm)' % ('%','%'))
        sl.append('=============  =============  =============  =================')
        sl.append('%s       %s         %s           %s       ' % 
            (digit_length(ctab[ignore[N],1],8), digit_length(ctab[ignore[N],2],7), digit_length(ctab[ignore[N],3],4), digit_length(ctab[ignore[N],4],4)))
        sl.append('=============  =============  =============  =================\n\n')
        sl.append('.. image:: %s/Axial_GS_Component_'    % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Sagittal_GS_Component_' % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('.. image:: %s/Coronal_GS_Component_'  % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 95%\n')
        sl.append('\n.. image:: %s/TimeSeries_'          % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        sl.append('.. image:: %s/FFT_Component_'         % figures + str(i).zfill(size) + '.png'); sl.append('\t:width: 49%\n')
        N += 1


    ofh = open("%s/Analysis.rst" % outputDir,"w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()
"""
make diagnostics.rst file
"""

def diagnostics_rst(outputDir):
    figures = "Report_Figures"
    sl = []
    sl.append('Preliminary Diagnostics\n' + '==========================\n')
    sl.append('The signal to noise ratio (TSNR) for a dataset is defined as the mean over the standard deviation of the dataset.  ' +
        'Meica.py denoises the BOLD time signal which should increase TSNR. This will be seen as the medn dataset having a greater ' + 
        'TSNR than the tsoc dataset.\n')
    sl.append('TSNR\n' + '+++++++')
    sl.append('Absolute TSNR of the medn NIFTI1 dataset, Accepted and Ignored components:\n')
    sl.append('.. image::  %s/medn_tsnr.png' % figures)
    sl.append('\t:width: 95%')
    sl.append('\t:align: left\n')
    sl.append('Absolute TSNR of the tsoc NIFTI1 dataset, all components:\n')
    sl.append('.. image:: %s/tsoc_tsnr.png' % figures)
    sl.append('\t:width: 95%')
    sl.append('\t:align: left\n\n')
    sl.append('Ratio of the TSNR of the medn NIFTI1 dataset to the TSNR of the tsoc NIFTI1 dataset:\n')
    sl.append('.. image:: %s/tsnr_ratio.png' % figures)
    sl.append('\t:width: 95%')
    sl.append('\t:align: left\n')
    sl.append('Histograms of the TSNR of the medn NIFTI1 dataset and the TSNR of the tsoc NIFTI1 dataset:\n')
    sl.append('.. image:: %s/medn_tsnr_hist.png' % figures)
    sl.append('\t:width: 49%')
    sl.append('.. image:: %s/tsoc_tsnr_hist.png' % figures)
    sl.append('\t:width: 49%\n')
    sl.append('Histogram of the ratio of the TSNR of the medn dataset to the TSNR of the tsoc NIFTI1 dataset:\n')
    sl.append('.. image:: %s/tsnr_ratio_hist.png' % figures)
    sl.append('\t:width: 49%\n')
    ofh = open("%s/Diagnostics.rst" % outputDir,"w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()
"""
make index.rst file
"""
def index_rst(outputDir):
    title = outputDir.split('/')[-2]
    sl = []
    sl.append(title)
    sl.append('='* (len(title)+1))
    sl.append('The program meica.py was created to form an algorithmic method for performing '
        + 'independent component analysis on multi-echo data and then algorithmically deciding '
        + 'which components represent BOLD-like phenomena.\n')
    sl.append('The following content is a report that has taken information provided by '
        + 'meica.py and summarizes a few of the results.\n')
    sl.append('This report form was created by the Section on Functional Imaging Methods in the NIMH.  '
        + 'The creators of this report form are Benjamin Gutierrez, Daniel Handwerker, '
        + 'Javier Gonzalez-Castillo, Prantik Kundu, Souheil Inati, and Peter Bandettini.\n')
    sl.append('Contents:\n')
    sl.append('.. toctree::')
    sl.append('   :maxdepth: 2\n')
    sl.append('   Intro')
    sl.append('   Diagnostics') 
    sl.append('   Dynamic')
    sl.append('   Analysis')
    sl.append('\n\n\nSearch')
    sl.append('======\n')
    sl.append('* :ref:`search`')

    ofh = open("%s/index.rst" % outputDir,"w")
    ofh.write("\n".join(sl) + "\n")
    ofh.close()
"""
make intro.rst file
"""
def intro_rst(outputDir):
    sl = []
    sl.append('Intro\n' + '=====')
    sl.append('This report collects information from your ME-ICA analysis and displays '
         + 'several ways to interpret the output of meica.py.  The Preliminary Diagnostics '
         + 'section contains a view of the TSNR of the "denoised" BOLD time series after: basic '
         + 'preprocessing, T2* weighted averaging of echoes (i.e. "optimal combination"), and ICA denoising, the TSNR of the '
         + '"raw" BOLD time series dataset after: basic preprocessing, and a ratio of these two TSNR maps.\n\n')

    sl.append('''There are two ways to visualize your data in ths Report, the first is the dynamic report.
  The dynamic report uses bokeh to create multiple plots which can be interacted with in order to go through the your data.
  The static report contains the same information, but cannot be interacted with.\n''')
    sl.append('Requirements for generating this report form:\n\n' + '* matplotlib\n\n' + '* numpy\n\n' + '* nibabel\n\n' + '* sphinx\n' + '* bokeh\n')
    ofh = open("%s/Intro.rst" % outputDir,"w")
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

