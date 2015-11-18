#!/usr/bin/env python 

__version__ ="0.5"
import subprocess
import sys
import os
import nibabel as ni
import seaborn as sns
import numpy   as np

from glob                import glob
from scipy.stats         import rankdata
from bokeh.models        import ColumnDataSource, HoverTool, CustomJS, TapTool, Plot, Range1d
from bokeh.models.widgets import DataTable, TableColumn, NumberFormatter, BooleanFormatter, CheckboxEditor
from bokeh.models.widgets import Tabs, Panel, Select
from bokeh.document      import Document
from bokeh.plotting      import figure, gridplot, output_file, show, hplot, vplot
from scipy.stats         import scoreatpercentile
from bokeh.models.glyphs import ImageURL

# Load Inputs into memory
# =======================
accepted_components = np.loadtxt(TED_dir + '/accepted.txt',delimiter=',').astype('int')
ctab_unordered      = np.loadtxt(TED_dir + '/comp_table.txt')
comp_timeseries     = np.loadtxt(TED_dir + '/meica_mix.1D')
Nt, Nc              = comp_timeseries.shape
path                = TED_dir + '/' +'betas_OC.nii'
report_figures      = outputDir + '/' + 'Report_Figures'

for i in range(Nc):
    N = str(i).zfill(len(str(Nc)))
    subprocess.call('convert %s/Axial_GS_Component_%s.png %s/Sagittal_GS_Component_%s.png %s/Coronal_GS_Component_%s.png -pointsize 20 ' % (report_figures,N,report_figures,N,report_figures,N) +
                    'label:\'Component %s\' -append %s/AxSagCor_Component_%s.png' % (N,report_figures,N),shell=True)

subprocess.call('convert %s/Axial_GS_Component_XX.png %s/Sagittal_GS_Component_XX.png %s/Coronal_GS_Component_XX.png -pointsize 20 ' % (report_figures,report_figures,report_figures) +
                    'label:\'Component XX\' -append %s/AxSagCor_Component_XX.png' % report_figures ,shell=True)


p              = subprocess.Popen(['file','%s/AxSagCor_Component_XX.png' % report_figures],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
string,error   = p.communicate()

width , height = (str(string).split(',')[1]).split('x')

p         = subprocess.Popen(['3dinfo','-tr', path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)# retrieve TR
TR, err   = p.communicate()
comp_ffts = np.zeros(comp_timeseries.shape)

ICAmap_paths        = glob(outputDir +'/Report_Figures/'+ 'AxSagCor*.png')
ICAmap_default_path = outputDir + '/Report_Figures/' + 'AxSagCor_Component_XX.png'

# Compute FFTS
# =======================
for i in range(comp_timeseries.shape[1]):
    fft = abs(np.fft.fft(comp_timeseries[:,i]))
    comp_ffts[:,i] = np.fft.fftshift(fft)
    
freq      = np.fft.fftfreq(comp_timeseries.shape[0],float(TR[:-2]))
freq      = np.fft.fftshift(freq)
freq_axis = freq[np.where(freq == 0)[0][0]:][::-1]

# Setting up output file
# ======================
output_file("%s/Bokeh_plot.html" % (report_figures))

# Correcting comp_table columns
# =============================
ctab = np.zeros((ctab_unordered.shape[0],5))
with open("%s/comp_table.txt" % (TED_dir), 'r') as original: ctab_txt = original.read()
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

# Split components into correct bins
# ==================================
meica_txt = []
accepted_components_mask     = np.zeros((Nc,), dtype = bool)
rejected_components_mask     = np.zeros((Nc,), dtype = bool)
middle_kappa_components_mask = np.zeros((Nc,), dtype = bool)
ignored_components_mask      = np.zeros((Nc,), dtype = bool)

with open(TED_dir + '/' + 'comp_table.txt') as f:
    lines = f.readlines()
for i in range(len(lines)):
    if "ACC" in lines[i]:
        accepted_components = ((lines[i]).split(' ')[1]).split(',')
        accepted_components_mask[np.array(((lines[i]).split(' ')[1]).split(',')    ,dtype=int)] = True
    if "REJ" in lines[i]:
        rejected_components = ((lines[i]).split(' ')[1]).split(',')
        rejected_components_mask[np.array(((lines[i]).split(' ')[1]).split(',')    ,dtype=int)] = True
    if "MID" in lines[i]:
        middle_kappa_components = ((lines[i]).split(' ')[1]).split(',')
        middle_kappa_components_mask[np.array(((lines[i]).split(' ')[1]).split(','),dtype=int)] = True
    if "IGN" in lines[i]:
        ignored_components = ((lines[i]).split(' ')[1]).split(',')
        ignored_components_mask[np.array(((lines[i]).split(' ')[1]).split(',')     ,dtype=int)] = True
# Reading the different features of interest
# ==========================================
Nc,Nf = ctab.shape
kappa = ctab[:,1]
rho   = ctab[:,2]
var   = ctab[:,3]
cID   = ctab[:,0]
ratio = ctab[:,1]/ctab[:,2]

loc_by_kappa = Nc - rankdata(kappa)
loc_by_rho   = Nc - rankdata(rho)
loc_by_var   = Nc - rankdata(var)
loc_by_ratio = Nc - rankdata(ratio)

fica_psel                      = np.zeros((Nc,))
fica_psel[accepted_components_mask]     = 3
fica_psel[ignored_components_mask]      = 2
fica_psel[middle_kappa_components_mask] = 1
component_colormap             = { "3.0" : "#0000ff",  "2.0" : "#00ffff",  "1.0" : "#00ff00",      "0.0" : "#ff0000",}
component_color                = [component_colormap[str(x)] for x in fica_psel]
component_status_labels        = { "3.0" : "Accepted", "2.0" : "Ignored",  "1.0" : "Middle Kappa", "0.0" : "Rejected"}
component_status               = [component_status_labels[str(x)] for x in fica_psel]

Source = ColumnDataSource(data = dict(cID = cID, 
				      kappa = kappa, loc_by_kappa = loc_by_kappa, 
			 	      rho = rho, loc_by_rho = loc_by_rho, 
				      var = var, loc_by_var = loc_by_var, 
				      comp_color  = component_color,
				      comp_status = component_status,
				      ratio = ratio, loc_by_ratio = loc_by_ratio))

# ==============================================================================
#                                 FEATURE TABLE

comp_table_columns = [
    TableColumn(field="cID"        ,title="ID"      ,formatter=NumberFormatter(format='0o')),
    TableColumn(field="comp_status",title="Status"  ,editor=CheckboxEditor()),
    TableColumn(field="var"        ,title="Variance",formatter=NumberFormatter(format='0.00')),
    TableColumn(field="kappa"      ,title="Kappa"   ,formatter=NumberFormatter(format='0.00')),
    TableColumn(field="rho"        ,title="Rho"     ,formatter=NumberFormatter(format='0.00')),
    TableColumn(field="ratio"      ,title="Ka/Rh"   ,formatter=NumberFormatter(format='0.00'))
]

comp_table_DTABLE = DataTable(source=Source,columns=comp_table_columns,width=1350, height=250, editable=False, selectable=True, sortable=True)

# ==============================================================================
#                                 FEATURE PLOTS
# Feaute Plots Tools
# ==================
TOOLS      = "tap,box_zoom,reset"
HoverKappa = HoverTool(
tooltips   =[("Component", "@cID"),("Kappa", "@kappa"),("Rho", "@rho"),("Variance", "@var"),("Ratio", "@ratio")])

HoverRho   = HoverTool(
tooltips   =[("Component", "@cID"),("Kappa", "@kappa"),("Rho", "@rho"),("Variance", "@var"),("Ratio", "@ratio")])

HoverVar   = HoverTool(
tooltips   =[("Component", "@cID"),("Kappa", "@kappa"),("Rho", "@rho"),("Variance", "@var"),("Ratio", "@ratio")])

HoverRatio = HoverTool(
tooltips   =[("Component", "@cID"),("Kappa", "@kappa"),("Rho", "@rho"),("Variance", "@var"),("Ratio", "@ratio")])

HoverKvsR  = HoverTool(
tooltips   =[("Component", "@cID"),("Kappa", "@kappa"),("Rho", "@rho"),("Variance","@var") ,("Ratio", "@ratio")])

# Feature Plots
# =============
sp_kappa = figure(tools=[TOOLS, HoverKappa],width=325, height=250, y_axis_label='Kappa', toolbar_location='left')
sp_kappa.circle('loc_by_kappa','kappa',size=5,color='comp_color',source=Source)
sp_legend = figure( title="", x_axis_type=None, y_axis_type=None, x_range=[-100, 100], y_range=[-100, 100],
            min_border=0, outline_line_color="#ffffff", width=325, height=250, toolbar_location="left", tools=[TOOLS])
sp_legend.circle([-30,-30,-30,-30], [25,8,-9,-26], color =  ["#0000ff" ,"#ff0000" ,"#00ff00"     ,"#00ffff"], radius=5)
sp_legend.text(  [-10,-10,-10,-10], [24,7,-10,-27], text =  ["Accepted","Rejected","Middle Kappa","Ignored"],    text_font_size="12pt", text_align="left", text_baseline="middle",text_font_style = "normal")

sp_kappa.yaxis.axis_label_text_font_size = "12pt"
sp_tab_kappa  = Panel(child=sp_kappa,  title='Sorted by Kappa')
sp_tab_legend = Panel(child=sp_legend, title='Legend')
sp_tabs_kappa = Tabs(tabs=[sp_tab_kappa,sp_tab_legend])
 
sp_rho   = figure(tools=[TOOLS, HoverRho],width=325, height=250, y_axis_label='Rho', toolbar_location=None,x_range=sp_kappa.x_range)
sp_rho.circle('loc_by_rho','rho',size=5,color='comp_color',source=Source)
sp_rho.yaxis.axis_label_text_font_size = "12pt"
sp_tab_rho  = Panel(child=sp_rho, title='Sorted by Rho')
sp_tabs_rho = Tabs(tabs=[sp_tab_rho])

sp_var   = figure(tools=[TOOLS,HoverVar],width=325, height=250, y_axis_label='Variance', toolbar_location=None,x_range=sp_kappa.x_range)
sp_var.circle('loc_by_var','var',size=5,color='comp_color',source=Source)
sp_var.yaxis.axis_label_text_font_size = "12pt"
sp_tab_var  = Panel(child=sp_var, title='Sorted by Variance')
sp_tabs_var = Tabs(tabs=[sp_tab_var])

sp_ratio = figure(tools=[TOOLS,HoverRatio],width=325, height=250, y_axis_label='K/R Ratio', toolbar_location=None,x_range=sp_kappa.x_range)
sp_ratio.circle('loc_by_ratio','ratio',size=5,color='comp_color',source=Source)
sp_ratio.yaxis.axis_label_text_font_size = "12pt"

sp_kvr = figure(tools=[TOOLS,HoverKvsR],width=325, height=250, y_axis_label='rho', x_axis_label='kappa', toolbar_location=None)
sp_kvr.circle('kappa','rho',size=((var)/np.max(var))*20 + 2,color='comp_color',source=Source)
sp_kvr.xaxis.axis_label_text_font_size = "12pt"
sp_kvr.yaxis.axis_label_text_font_size = "12pt"

sp_tab_ratio = Panel(child=sp_ratio, title='Kappa/Rho Ratio')
sp_tab_kvr   = Panel(child=sp_kvr, title='Kappa vs. Rho')
sp_tabs_left = Tabs(tabs=[sp_tab_ratio,sp_tab_kvr])

# ==============================================================================

# ==============================================================================
#                          TIME SERIES PLOTS

# Load Default Data on Plots
# ==========================
default_ts_x  = range(Nt)
default_ts_y  = np.zeros((Nt,))
default_fft_x = freq_axis
default_fft_y = np.zeros((Nt,))

# Generate Plots
# ==============
sp_ts = figure(tools=[],toolbar_location=None, width=680,height=200, x_axis_label='Time [TR]', 
               title='Component Timeseries',x_range=(0,Nt), title_text_font_size='12pt')
sp_ts.xaxis.axis_label_text_font_size = "12pt"

sp_fft = figure(tools=[],toolbar_location=None, width=680,height=200, x_axis_label='Frequency', 
                title='Component Spectrum',x_range=(min(default_fft_x), max(default_fft_x)), title_text_font_size='12pt')
sp_fft.xaxis.axis_label_text_font_size = "12pt"

# Generate Data Sources for interactivity
# =======================================
timeseries_to_display = ColumnDataSource(data=dict(x=default_ts_x, y=default_ts_y))
available_timeseries  = ColumnDataSource(data=dict(x=default_ts_x, y=comp_timeseries.T,comp_color  = component_color, cID=cID))
sp_ts.line(x = 'x',y = 'y', source=timeseries_to_display, line_width=3)

ffts_to_display = ColumnDataSource(data=dict(x=default_fft_x, y=default_fft_y))
available_ffts  = ColumnDataSource(data=dict(x=default_fft_x, y=comp_ffts.T,comp_color  = component_color,cID=cID))
sp_fft.line(x = 'x',y = 'y',source=ffts_to_display, line_width=3)
# ==============================================================================

# ==============================================================================
#                    BRAIN MAP PLOTS
# Convert Input maps into usable mosaics
# ======================================
available_ICAmaps    = ColumnDataSource(data=dict(urls=ICAmap_paths,cID=cID))
ICAmap_to_display    = ColumnDataSource(data=dict(x=[0], y=[0], w=[int(width)], h=[int(height)], url=[ICAmap_default_path]))
xdr                  = Range1d(start=-(int(width)/2), end=(int(width)/2))
ydr                  = Range1d(start=-(int(height)/2), end=(int(height)/2))
ICAmapFigure         = figure(tools=[],title="ICA maps", x_range=xdr, y_range=ydr,width=1200, height=int((int(height)*9)/10), x_axis_type=None, y_axis_type=None, toolbar_location=None, title_text_font_size ='12pt')
ICAmapImg            = ImageURL(url="url", x="x", y="y", w="w", h="h", anchor="center")
ICAmapFigure.add_glyph(ICAmap_to_display,ICAmapImg)
ICAmapFigure.outline_line_color='#ffffff'

# ==============================================================================

# ==============================================================================
#                   JAVA SCRIPT INTERACTIVITY
 
update_ts = CustomJS(args=dict(timeseries_to_display=timeseries_to_display, 
                               comp_ts=available_timeseries, 
                               ffts_to_display=ffts_to_display, 
                               comp_fft=available_ffts,
                               ICApaths=available_ICAmaps,
                               ICAmap_to_display=ICAmap_to_display), 
       code="""
         var c            = cb_obj.get('selected')['1d'].indices
         
         var data2disp_ts = timeseries_to_display.get('data')
         x2disp_ts        = data2disp_ts['x'];
         y2disp_ts        = data2disp_ts['y'];
         var comp_ts      = comp_ts.get('data');
         ts_x             = comp_ts['x'];
         ts_y             = comp_ts['y'];
         for (i = 0; i < x2disp_ts.length; i++) {
            y2disp_ts[i]  = ts_y[c][i];
         }
         
         var data2disp_fft = ffts_to_display.get('data');
         x2disp_fft        = data2disp_fft['x'];
         y2disp_fft        = data2disp_fft['y'];
         var comp_fft      = comp_fft.get('data');
         fft_x             = comp_fft['x'];
         fft_y             = comp_fft['y'];
         for (i=0; i < x2disp_fft.length; i++) {
             y2disp_fft[i] = fft_y[c][i]
         }
         
         var ICA2display  = ICAmap_to_display.get('data');
         url2display      = ICA2display['url'];
         var availICAurls = ICApaths.get('data');
         allICAurls       = availICAurls['urls'];
         url2display[0]   = allICAurls[c];
         
         ICAmap_to_display.trigger('change');
         timeseries_to_display.trigger('change');
         ffts_to_display.trigger('change');
    """)

# Additional Code to link with time series
kappa_taptool          = sp_kappa.select(type=TapTool)
kappa_taptool.callback = update_ts 
rho_taptool            = sp_rho.select(type=TapTool)
rho_taptool.callback   = update_ts 
var_taptool            = sp_var.select(type=TapTool)
var_taptool.callback   = update_ts 
ratio_taptool          = sp_ratio.select(type=TapTool)
ratio_taptool.callback = update_ts 
kvr_taptool            = sp_kvr.select(type=TapTool)
kvr_taptool.callback   = update_ts


# ==============================================================================

# ==============================================================================
#                       GRAPH   LAYOUT
top_left  = hplot(sp_tabs_kappa, sp_tabs_rho)
top_right = hplot(sp_tabs_var, sp_tabs_left)
top       = hplot(top_left, top_right)
middle    = hplot(sp_ts, sp_fft)
pl        = vplot(comp_table_DTABLE,top, middle,ICAmapFigure)
p         = hplot(pl)
show(p)
# ==============================================================================
