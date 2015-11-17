#!/usr/bin/env python
"""
Gutierrez, B.  Generates Meica report form.

"""
__version__ = "3.0"

from multiprocessing import cpu_count
import subprocess
import argparse
import sys
import os

path = os.path.abspath(os.path.dirname(__file__))
print("++ INFO [Main]: Using Report library located in: %s" % path)
if not path in sys.path:
    sys.path.insert(1, path)
del path

def dep_check():
    print("++ INFO [Main]: Checking for dependencies....")
    fails   = 0
    modules = set(["numpy", "matplotlib","argparse","scipy","bokeh","nibabel","sphinx"])

    for m in modules:
        try:
            __import__(m)
        except ImportError:
            fails += 1
            print("++ ERROR [Main]: Can't import Module %s. Please install." % m)
                
    if fails == 0:
        print(" +              All Dependencies are OK.")
    else:
        print(" +              All dependencies not available. Please install according to above error messages.")
        print(" +              Program exited.")
        sys.exit()

if __name__=='__main__':
    print("------------------------------------")
    print("-- SFIM ME-ICA Report version %s --" % __version__)
    print("------------------------------------")
    dep_check()
    
    
    from argparse import ArgumentParser
    import meica_figures
    import rst_files
    import sphinx_files
    import nibabel as ni
    import numpy   as np
    import pandas  as pd
    # Parse input parameters
    # ----------------------
    parser = ArgumentParser()
    parser.add_argument("-t","--TED_dir",          dest='TED_dir',        help="Path to meica output TED directory",                        type=str, default=None)
    parser.add_argument("-o","--out_dir",          dest='out_dir',        help="Output directory to output report, default='./meica.Report",type=str, default='./meica.Report')
    parser.add_argument(     "--motion",           dest='motion',         help="Path to motion.1D file",                                    type=str, default=None)
    parser.add_argument(     "--ncpus",            dest='Ncpus',          help='Number of cpus available. Default will be #available/2',    type=int, default=None)
    parser.add_argument(     "--overwrite",        dest='overwrite',      help="overwrite files previous created", action='store_true')
    options = parser.parse_args()
    
    if (options.Ncpus is None) or (options.Ncpus > cpu_count()):
        Ncpu = int(cpu_count()/2)
    else:
        Ncpu = int(options.Ncpus) 

    print("++ INFO [Main]: Output dir = %s" % os.path.abspath(os.path.expanduser(options.out_dir)))
    print("++ INFO [Main]: Overwrite old files? %s" % options.overwrite)
    print("++ INFO [Main]: Number of CPUs to use: %d" % Ncpu)
    # Control all necessary inputs are available
    # ------------------------------------------
    if options.TED_dir is None:
        print("++ Error: No information about TED directory given.")
        sys.exit()
    # Control for existence of files and directories
    # ----------------------------------------------
    fails = 0
    if not os.path.exists(options.TED_dir) and (not os.path.isdir(options.TED_dir)):
        print("++ Error: TED directory [%s] does not exists." % options.TED_dir)
        sys.exit()
    if (os.path.exists(options.out_dir)) and (not options.overwrite):
        print("++ Error: Output directory exits [%s], will not overwrite existing files without overwrite option." % options.out_dir)
        sys.exit()
    if not os.path.exists(options.TED_dir + '/' + 'accepted.txt'):
        print("++ Warning: File with accepted components [%s] not found. Attempting to source component bins by comp_table." % options.TED_dir + '/' + 'accepted.txt')
        source_comp_table = True
    else:
        source_comp_table = False
    if not os.path.exists(options.TED_dir + '/' + 'comp_table.txt'):
        print("++ Error: Component table [%s] does not exists." % (options.TED_dir + '/' + 'comp_table.txt'))
        fails += 1
    if not os.path.exists(options.TED_dir + '/' + 'meica_mix.1D'):
        print("++ Error: ME-ICA component timeseries [%s] does not exists." % options.TED_dir + '/' + 'meica_mix.1D')
        fails += 1
    if not os.path.exists(options.TED_dir + '/' + 'ts_OC.nii'):
        print("++ Error: Optimally Combined timeseries [%s] does not exists." % options.TED_dir + '/' + 'ts_OC.nii')
        fails += 1
    if not os.path.exists(options.TED_dir + '/' + 'dn_ts_OC.nii'):
        print("++ Error: ME-ICA Denoised timeseries [%s] does not exists." % options.TED_dir + '/' + 'dn_ts_OC.nii')
        fails += 1
    if not os.path.exists(options.TED_dir + '/' + 'betas_OC.nii') and not os.path.exists(glob.glob(options.TED_dir + '/' + '*.ICA.Zmaps.nii')[0]):
        print("++ Error: ME-ICA Denoised timeseries [%s] does not exists." % options.TED_dir + '/' + 'betas_OC.nii')
        print("++ Error: ME-ICA Denoised timeseries [%s] does not exists.  Only worry about this if you are Javier" % glob.glob(options.TED_dir + '/' + '*.ICA.Zmaps.nii')[0])
        fails += 1
    if fails != 0:
        sys.exit()
    
    # Set all output paths
    # --------------------
    print("++ INFO [Main]: Setting up directories.")
    outputDir = os.path.abspath(options.out_dir)
    figureDir = os.path.abspath(options.out_dir)
    if options.overwrite:
        subprocess.call('rm -rf %s' % outputDir, shell=True)
        subprocess.call('mkdir %s' % outputDir,  shell=True)
    if not (os.path.exists(outputDir)):
        subprocess.call('mkdir %s' % outputDir,  shell=True)
    subprocess.call('mkdir %s/Report_Figures' % outputDir,    shell=True)
    subprocess.call('mkdir %s/sphinx_files' % outputDir,      shell=True)
    subprocess.call('mkdir %s/axialized_nifti' % outputDir, shell=True)
    subprocess.call('cp %s/bokeh_plot.py %s/Report_Figures/bokeh_plot.py' % (os.path.dirname(__file__), outputDir),shell=True)
    # =========  report starts   =========
    # ====================================
    
    # load data
    # ---------
    print("++ INFO [Main]: Loading data...")
    
    ctab_unordered          = np.loadtxt(options.TED_dir + '/' + 'comp_table.txt')
    Denoised_components_ts  = np.loadtxt(options.TED_dir + '/' + 'meica_mix.1D')
    Denoised_ts             = ni.load(options.TED_dir + '/' + 'dn_ts_OC.nii').get_data()
    OptCombine_ts           = ni.load(options.TED_dir + '/' + 'ts_OC.nii').get_data()
    try:
        Denoised_components = ni.load(options.TED_dir + '/' + 'betas_OC.nii').get_data()
    except:
        Denoised_components = ni.load(glob.glob(options.TED_dir + '/' + '*.ICA.Zmaps.nii')[0]).get_data()
    
    
    # parse and order component table columns
    # ---------------------------------------
    ctab = np.zeros((ctab_unordered.shape[0],5))
    print("++ INFO [Main]: Parsing comp_table.txt")
    with open("%s/comp_table.txt" % (options.TED_dir), 'r') as original: ctab_txt = original.read()
    ctab_txt = ctab_txt.split('\n')
    N = 0
    while '#' not in ctab_txt[-2-N][0]:
        N += 1
    ctab_columns = ctab_txt[-2 -N].split()

    if ctab_columns[0] == '#':
        ctab_columns = ctab_columns[1:]

    for i in range(len(ctab_columns)):
        if ctab_columns[i] == 'Kappa':
            ctab[:,1] = ctab_unordered[:,i]
        if ctab_columns[i] == 'Rho':
            ctab[:,2] = ctab_unordered[:,i]
        if ctab_columns[i] == '%%Var' or ctab_columns[i] == 'Var':
            ctab[:,3] = ctab_unordered[:,i]
        if ctab_columns[i] == '%%Var(norm)':
            ctab[:,4] = ctab_unordered[:,i]

    if np.max(ctab[:,4] == 0):
        ctab[:,4] = ctab[:,3]

    ctab = ctab[ctab[:,1].argsort()[::-1]]
    ctab[:,0] = np.arange(ctab.shape[0])
    
    # set variables
    # -------------
    meica_txt = []
    if source_comp_table == True:
        with open(options.TED_dir + '/' + 'comp_table.txt') as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "ACC" in lines[i]:
                accepted = ((lines[i]).split(' ')[1]).split(',')
            if "REJ" in lines[i]:
                rejected = ((lines[i]).split(' ')[1]).split(',')
            if "MID" in lines[i]:
                middle_kappa = ((lines[i]).split(' ')[1]).split(',')
            if "IGN" in lines[i]:
                ignored = ((lines[i]).split(' ')[1]).split(',')
                
    else:
        accepted = np.loadtxt(options.TED_dir + '/' + 'accepted.txt', delimiter = ',', dtype = 'int', ndmin = 1)
        rejected = np.loadtxt(options.TED_dir + '/' + 'rejected.txt', delimiter = ',', dtype = 'int', ndmin = 1)
        if os.path.exists(options.TED_dir + '/' + 'midk_rejected.txt'):
            middle_kappa = np.loadtxt(options.TED_dir + '/' + 'midk_rejected.txt', delimiter = ',',dtype = 'int', ndmin = 1)
        else:
            middle_kappa = []
        Nc = Denoised_components.shape[3]
        ignored = np.array([],dtype = 'int')
        for i in range(Nc):
            if i not in np.concatenate((accepted,rejected,middle_kappa)):
                ignored = np.append(ignored,[i])

    print("++ INFO [Main]: Creating report figures.")
    meica_figures.kr_vs_component(ctab,outputDir)
    meica_txt.append(meica_figures.kappa_vs_rho_plot(accepted, rejected, middle_kappa, ignored,ctab,outputDir,options.TED_dir))
    meica_txt.append(meica_figures.tsnr(OptCombine_ts, Denoised_ts, outputDir))
    meica_figures.montage_control(options.TED_dir, outputDir, Denoised_components, Denoised_components_ts, Ncpu)

    if os.path.isfile(options.motion):
        meica_txt.append(meica_figures.motion(outputDir,options.motion))
        motion_file = options.motion
    else:
        meica_txt.append("Max head displacement in any one dirrection:   %s\nTR of Max Head displacement:   %s\nMax rate of head motion:   %s\nTR of max head motion rate:   %s" % (' ',' ',' ',' '))
        motion_file = ''
    
    if options.TED_dir[-1] == '/':
        TED_dir = options.TED_dir[:-1]
    else:
        TED_dir = options.TED_dir
    newline = "#!/usr/bin/env python \n" + ("TED_dir = '%s' \n" % os.path.abspath(TED_dir)) + ("outputDir,figures = '%s','%s'\n" % (outputDir,'Report_Figures'))
    with open("%s/Report_Figures/bokeh_plot.py" % outputDir, 'r') as original: bokeh_plot = original.read()
    with open("%s/Report_Figures/bokeh_plot.py" % outputDir, 'w') as modified: modified.write(newline + bokeh_plot)
    
    # set up sphinx documentation
    # ---------------------------
    
    sphinx_files.conf(__version__,outputDir)
    sphinx_files.make_file(outputDir)
    
    # make .rst files for sphinx to use to generate the report
    # --------------------------------------------------------
    rst_files.diagnostics_rst(outputDir)
    rst_files.index_rst(outputDir)
    rst_files.intro_rst(outputDir)
    rst_files.analysis_rst(accepted, rejected, middle_kappa, ignored, ctab, outputDir, motion_file, options.TED_dir)
    rst_files.dynamic_analysis_rst(accepted, rejected, middle_kappa, ignored, ctab, outputDir, motion_file, options.TED_dir)

    ofh = open("meica_report.txt","w")
    ofh.write("\n".join(meica_txt) + "\n")
    ofh.close()

    # run sphinx build
    # ----------------
    os.chdir(outputDir)
    subprocess.call('mkdir _static', shell = True)
    subprocess.call('make html'    , shell = True)
    subprocess.call('make latex'   , shell = True)
    
    #if options.latex:
    #    subprocess.call('make latexpdf', shell = True)

    subprocess.call('mv %s/_build/* %s' % (outputDir, outputDir), shell = True)
    #subprocess.call('rm -rf _*', shell = True)
    subprocess.call('mv %s/*.rst %s/sphinx_files/' % (outputDir, outputDir), shell = True)
    subprocess.call('mv %s/Makefile %s/sphinx_files' % (outputDir, outputDir), shell = True)
    subprocess.call('mv %s/conf.py %s/sphinx_files' % (outputDir, outputDir), shell = True)
    
