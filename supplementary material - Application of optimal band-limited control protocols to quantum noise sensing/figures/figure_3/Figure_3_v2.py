# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:47:31 2017

@author: virginia
"""
from __future__ import division, print_function
import sys
sys.path.append('/home/virginia/Dropbox/Simple fidelity experiments/Core modules')
sys.path.append('C:\\Users\\Sandeep\\Dropbox\\Simple fidelity experiments\\Core modules')
from Slepian import DPSS
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from scipy.optimize import curve_fit

from PyPDF2 import  PdfFileReader, PdfFileWriter


def pdfMerger(pdfs_in, pdf_out, path=''):
    ''' Merge multiple PDF files together '''
    output = PdfFileWriter()
    
    main_pdf = PdfFileReader(open(pdfs_in[0],"rb")).getPage(0)
    
    for i in range(1, len(pdfs_in)):
        pdf = PdfFileReader(open(pdfs_in[1],"rb")).getPage(0)
        main_pdf.mergePage(pdf)
    
    output.addPage(main_pdf)
    output.write(open(path+pdf_out,"wb"))
    
def gradient_fill(x, y, fill_color=None, ax=None, **kwargs):
    """
    Plot a line with a linear alpha gradient filled beneath it.

    Parameters
    ----------
    x, y : array-like
        The data values of the line.
    fill_color : a matplotlib color specifier (string, tuple) or None
        The color for the fill. If None, the color of the line will be used.
    ax : a matplotlib Axes instance
        The axes to plot on. If None, the current pyplot axes will be used.
    Additional arguments are passed on to matplotlib's ``plot`` function.

    Returns
    -------
    line : a Line2D instance
        The line plotted.
    im : an AxesImage instance
        The transparent gradient clipped to just the area beneath the curve.
    """
    if ax is None:
        ax = plt.gca()
    x = np.array(x); y = np.array(y)
    line, = ax.plot(x, y, **kwargs)
    if fill_color is None:
        fill_color = line.get_color()

    zorder = line.get_zorder()
    alpha = line.set_alpha(0.0)
    alpha = 0.9#1.0 if alpha is None else alpha

    z = np.empty((100, 1, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)
    z[:,:,:3] = rgb
    z[:,:,-1] = np.linspace(0, alpha, 100)[::-1,None]

    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min()*0.01, y.max()*0.95
    im = ax.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder, interpolation='nearest')

    xy = np.column_stack([x, y])
    xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    clip_path = Polygon(xy, facecolor='none', edgecolor='none', closed=True)
    ax.add_patch(clip_path)
    im.set_clip_path(clip_path)

    #ax.autoscale(True)
    return line, im

def mm2inch(mm):
    return mm * 0.0393701
    

# General rc parameter    
matplotlib.rcParams['figure.figsize'] = (mm2inch(80), mm2inch(60))
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.size'] = 1.5
matplotlib.rcParams['ytick.major.size'] = 1.5
font = {'family' : 'Arial',
        'size'   : 8}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 0.5
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)



fig = plt.figure()
grid = gridspec.GridSpec(2, 1, height_ratios=[1.0, 1.0])

ax1 = plt.subplot(grid[0, 0])
ax2 = plt.subplot(grid[1, 0])




###############################################################################################
#
#       TIGHT LAYOUT
#
###############################################################################################


plt.tight_layout()



###############################################################################################
#
#       Absolute positioning of subplots
#
###############################################################################################

for ax in [ax1, ax2]:
    pos1 = ax.get_position()
    pos2 = [pos1.x0 + 0.025, pos1.y0 + 0.025, pos1.width, pos1.height]
    ax.set_position(pos2)
    

for ax in [ax1]:
    pos1 = ax.get_position()
    h = 0.90; y1  = pos1.height*(1 - h)
    pos2 = [pos1.x0, pos1.y0 + y1, pos1.width, pos1.height*h*1.05]
    ax.set_position(pos2)

for ax in [ax2]:
    pos1 = ax.get_position()
    h = 1.1
    pos2 = [pos1.x0, pos1.y0, pos1.width, pos1.height*h]
    ax.set_position(pos2)



###############################################################################################
#
#       ax1 PLOTS: CALIBRATION
#
###############################################################################################


expect_c = 'Gray'; data_c = 'Green'

expect = np.loadtxt('assets/data/expected_fidelities_016.txt', skiprows=1)[12:]
subbin_data = np.loadtxt('assets/data/20161031016_subbinned_data.txt')[12:-1, :]
alphas = np.loadtxt('assets/data/alpha_values_016.txt', skiprows=1)[2:]

n = len(alphas)

data_noisy = [subbin_data[i+1+i*10:i+1+i*10+10, 0] for i in range(n)]
simulation = [expect[i+1+i*10:i+1+i*10+10] for i in range(n)]
data_clean = [subbin_data[i*11+i, 0] for i in range(n)]




max_peaks = 5


x = range(subbin_data.shape[0])

ax1.plot(expect, color=expect_c, lw=0.75, label='Simulation')
ax1.errorbar(x, subbin_data[:, 0], subbin_data[:, 1:].T, fmt ='o', markersize=1.25, 
             color=data_c, lw=0.75, capsize=1.5, label='Measurements')
             
ax1.errorbar(x[::11], subbin_data[::11, 0], subbin_data[::11, 1:].T, fmt ='s', markersize=1.0, 
             color='k', lw=0.75, capsize=1.5)


baseline = np.median(data_clean)             
#gradient_fill([x[0], x[-1]], [baseline, baseline], fill_color='k', ax=ax1, color='k')

ax1.set_xlim([0, (max_peaks)*11])

y_arr = [0.93-i*0.25 for i in range(max_peaks+1)]
mean_line = np.zeros((max_peaks+1, 2))




dBs = np.log10(1+alphas[:max_peaks])*10


for i in range(max_peaks+1):
    ax1.plot([i*11, i*11], [0, 1.1], 'k--', dashes=(1.5, 0.75), alpha=0.25, lw=0.5)

    mean = np.mean(data_noisy[i])
    ax1.plot([i*11, (i+1)*11], [mean, mean], color='k', lw=0.75, alpha=0.25)
    
    mean_line[i, :] = ((i*11+(i+1)*11)/2., mean)
    
    if i != max_peaks:
        if i == 0:
            label = r'$\alpha=%.4f $' % (dBs[i]) 
            ax1.text(2.5 + i*11, 0.95, label, fontsize=6)
        elif i != max_peaks-1:
            label = r'$%.4f $' % (dBs[i])
            ax1.text(5.5 + i*11, 0.95, label, fontsize=6) 
        else:
            label = r'$%.4f \, \mathrm{dB} $' % (dBs[i])
            ax1.text(3 + i*11, 0.95, label, fontsize=6) 
            
            
            
    
    #ax1.fill_between([i*11, (i+1)*11], [mean, mean], 0, color='DarkGreen', alpha=0.10)
    #ax1.arrow(11*i + 0.17, y_arr[i], 10.5, 0, fc='k', ec='k', lw=0.75, head_length=0.75, head_width=0.05, 
    #      overhang=0.2, length_includes_head=True, alpha=0.25 if i != 1 else 0.75)
         
#ax1.text(12, 0.75, r'$\varphi = 0 \rightarrow 2\pi$', fontsize=6)         


ax1.plot(mean_line[:, 0], mean_line[:, 1], 'kx-', lw=0.75, markersize=3.0, mew=1.0, label='Mean signal')
ax1.plot([0, 55], [mean_line[-2][1], mean_line[-2][1]], lw=0.5, color='k', ls=':', dashes=(2,1))
#ax1.plot([i*11, (i+1)*11], [mean, mean], color='DarkGreen', lw=0.75, label='Mean signal')
        

# Set ticklabel positions to match peaks in simulation
#x_ticks = [np.argmax(expect[i*11:11*(i+1)])+11*i for i in range(max_peaks)]
x_ticks = [1, 5, 10, 12, 16, 21, 23, 27, 32, 34, 39, 43, 45, 50, 54]

ax1.set_xticks(x_ticks)
#xticklabels = ['{:.1f}'.format(dBs[i]*1e3) for i in range(max_peaks)]
xticklabels = [r'$0$', '', r'$2\pi$', r'$0$', '', r'$2\pi$', r'$0$', '', r'$2\pi$', r'$0$', '', r'$2\pi$', r'$0$', '', r'$2\pi$']

ax1.set_xticklabels(xticklabels)

all_ticks = ax1.xaxis.get_major_ticks()

for tick in all_ticks[2::3]:
    tick.label1.set_horizontalalignment('right')
for tick in all_ticks[0::3]:
    tick.label1.set_horizontalalignment('left')

ax1.set_xlabel(r'Modulation phase (rad)', labelpad=1)


ax1.tick_params(top='off', right='off', pad=2.0)
ax1.set_ylabel(r'$P\,(\uparrow_z)$', labelpad=3.5)

ax1.set_ylim([0, 1.05])

ax1.text(48, 0.12, 'Baseline', fontsize=6, color='k')

legend  = ax1.legend(loc='center right', frameon=True, fontsize=7, numpoints=1, handlelength=1.0)
legend.get_frame().set_facecolor('w')
legend.get_frame().set_edgecolor('w')


###############################################################################################
#
#       ax2 PLOTS: SIGNAL
#
###############################################################################################





c1 = 'Gray'; c2 = 'RoyalBlue'
mk = 'o'

# -------------- Read in data -----------------------------------------------------------------
d1 = np.loadtxt('assets/data/20161031014_subbinned_data.txt')[2:, :]
d2 = np.loadtxt('assets/data/20161031015_subbinned_data.txt')[:-2, :][::-1, :]

f = np.arange(4, 101, 2)



# -------------- Plots ------------------------------------------------------------------------

ax2.fill_between(f, d1[:, 0], d2[:, 0], color=c1, lw=0, alpha=0.15)

ax2.plot(f, d1[:, 0],  color=c1, lw=0.5, alpha=0.5)
ax2.plot(f, d2[:, 0],  color=c1, lw=0.5, alpha=0.5)

ax2.errorbar(f, d1[:, 0], d1[:, 1:].T, color=c1, fmt='o', capsize=1.0, lw=0.75, alpha=0.75,
             markersize=1.0, mec=c1)
ax2.errorbar(f, d2[:, 0], d2[:, 1:].T, color=c1, fmt='o', capsize=1.0, lw=0.75, alpha=0.75, 
             label='Individual runs', markersize=1.0, mec=c1)

mean = (d1[:, 0] + d2[:, 0]) / 2
mean_err = np.sqrt((d1[:, 1:]**2).T + (d2[:, 1:]**2).T)

ax2.plot(f, mean,  color=c2, lw=0.5, alpha=0.5)

ax2.errorbar(f, mean, yerr=mean_err, color=c2, fmt='o', capsize=1.25, lw=0.75, 
             label='Average', markersize=1.0, mec=c2)
             
             

ax2.tick_params(top='off', right='off', pad=2.5)

ax2.set_ylim([-0.1, 0.85])
ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])

ax2.set_xlim([0, 100])
ax2.set_xticks([0, 20, 40, 60, 80, 100])

ax2.set_ylabel(r'$P\,(\uparrow_z)$', labelpad=3.5)
ax2.set_xlabel(r'Band shift $\omega_s$ (kHz)', labelpad=2.5)

ax2.legend(loc='upper left', frameon=False, fontsize=7, ncol=2, columnspacing=1.0,
           numpoints=1, handlelength=0.5, handleheight=-0.5, bbox_to_anchor=(0.20, 1.05))





###############################################################################################
#
#       ax2 inset PLOTS: Offset scans
#
###############################################################################################

offset_values = np.loadtxt('assets/data/offset_values_004.txt', skiprows=1)/1e3/21.6 # divided by pi area


data10 = np.loadtxt('assets/data/20161031004_subbinned_data.txt')
data50 = np.loadtxt('assets/data/20161031005_subbinned_data.txt')
data100 = np.loadtxt('assets/data/20161031006_subbinned_data.txt')


def fit_phase(x,  b, c=np.pi/2 * 0.9):
        return np.cos(c*x+b)**2 

popt10, pcov10 = curve_fit(fit_phase, offset_values, data10[:, 2], p0=[np.pi/2])
popt50, pcov50 = curve_fit(fit_phase, offset_values, data50[:, 2], p0=[np.pi/2])
popt100, pcov100 = curve_fit(fit_phase, offset_values, data100[:, 2], p0=[np.pi/2])


data = [data10, data50, data100]
popt = [popt10, popt50, popt100]

popt = [[np.pi/2*0.92], [np.pi*0.91], [np.pi*1.05]]



pos1 = ax2.get_position()
ax2_inset = fig.add_axes([pos1.x0 + pos1.width*0.50, pos1.y0+pos1.height*0.05, pos1.width*0.38, 
                          pos1.height*0.3])
                        
                  
colours = ['DarkGreen', 'DarkBlue', 'DarkBlue']
symbols = ['o', 'x']
delta = [10, 50]

for i in range(2):
    ax2_inset.plot(offset_values, fit_phase(offset_values, *popt[i]), color=colours[i], 
                                            lw=0.5)
                                            
    ax2_inset.plot(offset_values, data[i][:, 0], symbols[i], markersize=1.5, 
             color=colours[i], lw=0.75,  mec=colours[i], mfc='None' if i == 0 else colours[i],
            label=r'$\omega_s %i\,\mathrm{kHz}$' % delta[i])


ax2_inset.tick_params(top='on', right='on', labeltop='on', labelright='on', labelleft='off', 
                labelbottom='off', left='off', bottom='off', pad=0.9, labelsize=6)
                

ax2_inset.tick_params(axis='x', pad=-0.5)

ax2_inset.set_ylabel(r'$P\,(\uparrow_z)$', fontsize=7)   
ax2_inset.yaxis.set_label_coords(1.28, 0.5)             
                
ax2_inset.set_yticks([0.0,  0.5,  1.0])

ax2_inset.set_xlim([-1, 1])
ax2_inset.set_xticks([-1, -0.5,  0, 0.5, 1])
ax2_inset.set_xticklabels([r'$-\pi$', r'$-\pi/2$', 0, r'$\pi/2$', r'$\pi$']) 

ax2_inset.set_xlabel(r'Offset area', fontsize=7)
ax2_inset.xaxis.set_label_coords(0.780, 1.67)

ax2_inset.legend(loc='center left', frameon=False, fontsize=6, handlelength=1.0, 
                 bbox_to_anchor=(-0.50, 0.5), handletextpad=0.15, labelspacing=0.25)



fig.text(0.13, 0.940, 'a', fontsize=8, fontweight='bold')
fig.text(0.13, 0.48, 'b', fontsize=8, fontweight='bold')


fig_version = '2'




fig.savefig('Figure_3_v{}.pdf'.format(fig_version))
plt.close()






