# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 13:43:15 2017

@author: virginia
"""

from __future__ import division, print_function
import sys
sys.path.append('/home/virginia/Dropbox/Simple fidelity experiments/Core modules')
from Slepian import DPSS
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


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
    

def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if not exponent:
        exponent = int(np.floor(np.log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if not precision:
        precision = decimal_digits

    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)   


def mm2inch(mm):
    return mm * 0.0393701
    

# General rc parameter    
matplotlib.rcParams['figure.figsize'] = (mm2inch(80), mm2inch(75))
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


grid = gridspec.GridSpec(2, 1, height_ratios=[0.4, 0.6])

ax0 = plt.subplot(grid[0, 0]); ax0.axis('off')


bottom_grid = gridspec.GridSpecFromSubplotSpec(2,2, grid[1, 0], width_ratios=[1, 1])

ax1 = plt.subplot(bottom_grid[0, 0])
ax2 = plt.subplot(bottom_grid[0, 1])
ax3 = plt.subplot(bottom_grid[1, 0])
ax4 = plt.subplot(bottom_grid[1, 1])




plt.tight_layout()

for ax in [ax1, ax2, ax3, ax4]:
    pos1 = ax.get_position()
    x1 = pos1.width*(1.06 - 1.005) if (ax in [ax2, ax4]) else 0
    pos2 = [pos1.x0 - x1 - 0.01, pos1.y0 + 0.02, pos1.width*1.06, pos1.height*0.98]
    ax.set_position(pos2)
    
for ax in [ax2, ax4]:
    pos1 = ax.get_position()
    w = 1.21
    x1 = pos1.width*(w - 1)
    pos2 = [pos1.x0 - x1, pos1.y0, pos1.width*w, pos1.height]
    ax.set_position(pos2)
    
for ax in [ax2, ax4]:
    pos1 = ax.get_position()
    x1 = 0.047
    pos2 = [pos1.x0 + x1, pos1.y0, pos1.width, pos1.height*1.11]
    ax.set_position(pos2)
    
    
for ax in [ax1, ax2]:
    pos1 = ax.get_position()
    pos2 = [pos1.x0, pos1.y0-0.02, pos1.width, pos1.height]
    if True: ax.set_position(pos2)
    


for ax in [ax1, ax3]:
    pos1 = ax.get_position()
    pos2 = [pos1.x0, pos1.y0, pos1.width*0.60, pos1.height]
    if True: ax.set_position(pos2)
        
for ax in [ax1]:
    pos1 = ax.get_position()
    pos2 = [pos1.x0, pos1.y0+0.02, pos1.width, pos1.height]
    ax.set_position(pos2)


ax2_pos = ax2.get_position()
ax4_pos = ax4.get_position()
ax4.set_position([ax4_pos.x0, ax4_pos.y0, ax4_pos.width, ax2_pos.y0 - ax4_pos.y1 + ax4_pos.height*2])

ax2.axis('off')


colours = ['Gray', 'RoyalBlue', 'DarkOrange', 'DarkRed', 'RoyalBlue']
band_colour = (0, 1, 1, 20./255) 
a_start = 0.4; a_step = 0.30; 


# Time-domain pulses
n_points = 120; n_zeros = 5; zeros = np.zeros(n_zeros)
orders = [1, 2, 3, 4]
n_orders = len(orders)
dpss = DPSS(n_points, max(orders), max(orders))

x = np.arange(n_points)
ones = np.ones(n_points)

p_offset = 3.25
f_offset = 1.25


lw = 0.75
n_fft = 350
ones_fft = np.ones(n_fft-1)

# Units for frequency axis:
# total pulse duration: tau;  delt_t = tau/n_points;  f_max = 1/ 2*delt_t = n_points/ 2*tau
# Normalise f-axis (multiply by tau): f_max = n_points / 2

x_fft = np.linspace(0, n_points/2, n_fft)[1:]
band_end = 2*np.pi*max(orders) #  = 2*pi * NW / tau
x_band_end = np.argmin(np.abs(x_fft-band_end))


bumps = [3.2, 4.8, 6.4, 8.0]

text_x0 = 8.75; text_y0 = 0.5


pos_ax4 = ax4.get_position()

subplt_spacing = 0.125 # relative to subplot height

subplt_height = pos_ax4.height / n_orders 


ylim = [[1e-14, 5e4], [0.5e-11, 5], [1e-11, 5], [1e-11, 5]]; 
fill = ylim[3][1]*1e-4

# Loop through all orders of flat-top and slepian pulses
ax1.text(-20, 3.5, r'$k$', clip_on=False)
square_txt_pos = [1.3, -3.3, -9.0, -10.7]
slep_txt_pos = [0.3, -3.3, -6.7, -10.7]
for i in range(n_orders):

    o = orders[i]; a_curve = a_start+a_step*i if a_start+a_step*i < 1 else 1
    a_fill = a_curve / 2.5

    square_pulse = np.array(np.array_split(np.ones(n_points-2*n_zeros), o))
    square_pulse[1::2] *= -1
    square_pulse = np.concatenate((zeros, np.concatenate(square_pulse), zeros))
    
    invert = 1#-1 if i == 2 else 1
    ax1.plot(np.ones(n_points)*i*(-p_offset), ls='-', alpha=0.25, color='k', lw=0.5)    
    ax1.plot(square_pulse*invert - i*p_offset, color=colours[0], alpha=a_curve, lw=lw)
    ax1.fill_between(x, square_pulse*invert - i*p_offset, ones*(-i*p_offset), color=colours[0], alpha=a_fill, lw=0.1)
    
    #ax1.text(106, square_txt_pos[i], r'${}$'.format(o-1), fontsize=6)
    
    slepian = dpss[:, o-1]
    slepian /= np.max(np.abs(slepian))

    ax3.plot(np.ones(n_points)*i*(-p_offset), ls='-', alpha=0.25, color='k', lw=0.5)
    ax3.plot(slepian - i*p_offset, color=colours[1], alpha=a_curve, lw=lw)
    ax3.fill_between(x, slepian - i*p_offset, ones*(-i*p_offset), color=colours[1], alpha=a_fill, lw=0.1)
    #ax3.text(106, slep_txt_pos[i], r'${}$'.format(o-1), fontsize=6)
    
    fft_square = (np.abs(np.fft.fft(square_pulse, n=4096))**2)[1:n_fft]
    fft_slepian = (np.abs(np.fft.fft(slepian, n=4096))**2)[1:n_fft]

    fft_square /= np.max(fft_square)
    fft_slepian /= np.max(fft_slepian)
    
    ax4_inset = fig.add_axes([pos_ax4.x0, pos_ax4.y0 + subplt_height*(n_orders-1-i + subplt_spacing) - 0.003*i, 
                              pos_ax4.width, subplt_height*(1-subplt_spacing)])
    
    
    ax4_inset.semilogy(x_fft, fft_square, color=colours[0], alpha=0.9, lw=lw)
    ax4_inset.fill_between(x_fft, fft_square, fft_slepian, 
                           where=np.logical_and(fft_square>fft_slepian,  x_fft > band_end), 
                            color='Grey', alpha=0.25, lw=0.0, clip_on=True)    
    
    ax4_inset.semilogy(x_fft, fft_slepian, color=colours[1], alpha=0.9, lw=lw, clip_on=True if i != n_orders-1 else True)
    
    # Spectral concentration within band
    lambda_slep = 1 - np.sum(fft_slepian[:x_band_end]) / np.sum(fft_slepian)
    lambda_square = 1 - np.sum(fft_square[:x_band_end]) / np.sum(fft_square)
    
    
    ax4_inset.text(40, 0.5e-1, r'$L_{%d} =$' % (o-1,) + ' %s ' % sci_notation(lambda_square/lambda_slep ), fontsize=6)
    
    #ax4_inset.text(0.2, 0.0001*fill, r'$k = {}$'.format(o-1), fontsize=7)
    
    ax4_inset.set_ylim(ylim[i])

    for spine in ['top', 'bottom', 'right']:
        ax4_inset.spines[spine].set_visible(False)
    ax4_inset.tick_params(which='both', right='off', top='off', bottom='off',
                          labelbottom='off', pad=1.0)
    ax4_inset.tick_params(which='minor', length=1.0, width=0.25)
    if i != 0:
        ax4_inset.set_yticks([1e-8, 1e-2])
    else:
        ax4_inset.set_yticks([1e-10, 1e-2])
    ax4_inset.patch.set_alpha(0.0)
    
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax4_inset.transAxes, color='k', clip_on=False, lw=0.5)
    if i != n_orders - 1:        
        ax4_inset.plot((-d, +d), (0, 0), **kwargs)        # top-left diagonal
    if i != 0:
        ax4_inset.plot((-d, +d), (1, 1), **kwargs)        # top-left diagonal
        
        
    
ax4.fill_between([0, band_end], [-n_orders*f_offset-0.5, -n_orders*f_offset-0.5], [f_offset, f_offset], color=band_colour) 
ax4.plot([band_end, band_end], [-n_orders*f_offset-0.5, +n_orders*f_offset-0.5],  color='Grey', 
         dashes=(1.2,0.75), lw=0.70, alpha=0.5) 

  

for ax in [ax1, ax3]:
    ax.tick_params(top='off', right='off')
    ax.set_ylim([-i*p_offset-1.5, 3.3])
    ax.set_yticks([-i*p_offset for i in range(n_orders)][::-1])
    ax.set_yticklabels([i for i in range(n_orders)][::-1])
    ax.set_xticks([0,  60,  120])
    ax.set_xticklabels([0,  0.50,  1.00])
    
    
#ax1.tick_params(labelbottom=False)

for ax in [ax4]:
    ax.spines['left'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    ax.tick_params(which='both', top='off', right='off', left='off', labelleft='off', labelright='off')
    ax.set_xlim([0, x_fft[-1]])
    
    ax.set_ylim([-i*f_offset-0.5, f_offset])
    ax.set_yticks([-i*f_offset for i in range(n_orders)][::-1])
    ax.set_yticklabels([i+1 for i in range(n_orders)])
    
  


        
ax2.tick_params(labelbottom=False)


ax3.set_ylabel(r'Control sequence $\Omega_k(t)$ (a.u.)')
ax3.yaxis.set_label_coords(-0.21, 1.05)
ax3.set_xlabel(r'Time, $t$ ($\tau$)', labelpad=1.9)

ax4.set_ylabel(r'Filter function $F_{\Omega}(\omega)$ (a.u.)')
#ax4.yaxis.set_label_coords(1.21, 1.10)
ax4.yaxis.set_label_coords(-0.17, 0.5)
ax4.set_xlabel(r'Frequency, $\omega$ $(1/\tau)$', labelpad=1.9)

fig.text(0.01, 0.943, 'a', fontsize=8, fontweight='bold')
fig.text(0.5, 0.943, 'b', fontsize=8, fontweight='bold')
fig.text(0.1, 0.535, 'c', fontsize=8, fontweight='bold')
fig.text(0.1, 0.27, 'd', fontsize=8, fontweight='bold')
fig.text(0.5, 0.12, 'e', fontsize=8, fontweight='bold')



svg_version = '7'
plt_version = '5'
fig_version = '5'



if True:    # Export pdf from Inkscape if there were any changes
    import os
    cmd = 'inkscape assets/svg_schematic_v{v}.svg --export-pdf=assets/schematic_v{v}.pdf'.format(v=svg_version)
    os.system(cmd)

plot = 'assets/plot_v{}.pdf'.format(plt_version)
schematic = 'assets/schematic_v{}.pdf'.format(svg_version)
fig.savefig(plot)

# Combine plot with schematic
pdfMerger([plot, schematic], 'Figure_0_v{}.pdf'.format(fig_version), path='')

plt.close()