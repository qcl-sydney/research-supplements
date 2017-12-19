# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 08:31:14 2017

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
from scipy.signal import hilbert

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
    



def mm2inch(mm):
    return mm * 0.0393701
    

# General rc parameter    
matplotlib.rcParams['figure.figsize'] = (mm2inch(80), mm2inch(110))
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
grid = gridspec.GridSpec(4, 1, height_ratios=[0.4, 0.4, 0.4, 0.4])

ax1 = plt.subplot(grid[0, 0])
ax2 = plt.subplot(grid[1, 0])
ax3 = plt.subplot(grid[2, 0])
ax4 = plt.subplot(grid[3, 0])

plt.tight_layout()

###############################################################################################
#
#       absolute positioning of subplots
#
###############################################################################################

for ax in [ax1, ax2, ax3, ax4]:
    pos1 = ax.get_position()
    x1 = pos1.width*(0.95 - 1)
    pos2 = [pos1.x0 - x1 + 0.005, pos1.y0 + 0.005, pos1.width*0.98, pos1.height]
    ax.set_position(pos2)
    
for ax in [ax1, ax2]:
    pos1 = ax.get_position()
    y1 = 0.090 if ax == ax2 else 0.025
    pos2 = [pos1.x0, pos1.y0 + y1, pos1.width, pos1.height*0.86]
    ax.set_position(pos2)

for ax in [ax3]:
    pos1 = ax.get_position()
    pos2 = [pos1.x0, pos1.y0 + 0.035, pos1.width, pos1.height*1.3]
    ax.set_position(pos2)
    
for ax in [ax4]:
    pos1 = ax.get_position()
    pos2 = [pos1.x0, pos1.y0 + 0.005, pos1.width, pos1.height*1.08]
    ax.set_position(pos2)




###############################################################################################
#
#       ax1 & ax2 plots (band-shifted Slepians)
#
###############################################################################################
colours = ['k', 'Green', 'RoyalBlue']; symbols = ['o', 'x', 's']

mod_freqs = [5, 15, 25]



# ---------------- Read in data ---------------------------------------------------------------
am_data = []; ssb_data = []
for mf in mod_freqs:
    am_data.append(np.loadtxt('assets/data/AM Slepians {} kHz.txt'.format(mf)))
    ssb_data.append(np.loadtxt('assets/data/SSB Slepians {} kHz.txt'.format(mf)))

# --------------- Create the Slepian waveforms for the insets ---------------------------------
dpss = DPSS(500, 2, 6)[:,1]
t = np.linspace(0, 1.08e-3, num=500)


# -------------- Plots ------------------------------------------------------------------------
mk = 2; mew = 0.75; lw = 0.75
#h_pos = [[0.065, 0.342, 0.618], [0.065, 0.342, 0.618]]  # horizontal positions for Slepian pulses
h_pos = [[0.065, 0.342, 0.618], [0.015, 0.282, 0.558]]  # horizontal positions for Slepian pulses
text_pos = [[0.10, 9.2, 19.2], [0.1, 9.25, 19.5]]
axes = [ax1, ax2]

for i in range(3):
    data = [am_data[i], ssb_data[i]]; c = colours[i]; 
    m = [np.argmax(data[0][1:,0] == 0), np.argmax(data[1][1:,0] == 0)]
    ssb = ssb_data[i]; c = colours[i]; m_ssb = np.argmax(ssb[1:,0] == 0)
    
    for j in range(2):
        
        ax = axes[j]; s = data[j]; mi = m[j]
        #ax.grid(True, ls='-', alpha=0.25)
        ax.plot([mod_freqs[i], mod_freqs[i]], [0, 0.25], '--', color='Gray', dashes=(3,1), lw=0.5)
        ax.plot(s[:,3]/1e3, 1-s[:,4], ls='-', color=c, lw=lw)    
        ax.plot(s[:mi,0]/1e3, s[:mi,1], symbols[i], color=c, markersize=mk, 
             mfc='none' if symbols[i] != 'x' else c, mec=c, mew=mew, label=r'$\omega_s = %i \ \mathrm{kHz}$ ' % mod_freqs[i])
        

        #ax2.text(text_pos[1][i], 0.11, r'$f_{\mathrm{mod}} = %i \ \mathrm{kHz}$' % mod_freqs[i] , fontsize=6)
        ax1.text(text_pos[0][i], 0.125, r'$\omega_{s} = %i \, \mathrm{kHz}$' % mod_freqs[i] , fontsize=6)   
        
        
        # Add inset with Slepian pulse
        pos1 = ax.get_position(); y1 = 0.015 if j == 1 else 0.025; h1 = 0.35 if j == 1 else 0.3
        pos2 = [pos1.x0  + h_pos[j][i], pos1.y0 + pos1.height/2 + 0.005 + y1,  pos1.width*0.17 , pos1.height*h1] 
        ax_inset = fig.add_axes(pos2)
        if j ==0:
            invert = 1 if i != 0 else 1
            pulse = invert*dpss*np.cos(2*np.pi*mod_freqs[i]*1000*t)
        else:
            dpss_h = np.imag(hilbert(dpss))
            pulse = dpss*np.cos(2*np.pi*mod_freqs[i]*1000*t) - dpss_h*np.sin(2*np.pi*mod_freqs[i]*1000*t)
        ax_inset.plot(pulse, color=c, lw=0.5)
        ax_inset.axis('off')


#ax1.legend(loc='upper left', frameon=False, fontsize=7, ncol=3, bbox_to_anchor=(-0.00, 1.27),
#           handlelength=1.6, handletextpad=0.5, columnspacing=2, numpoints=3)


for ax in [ax1, ax2]:
    ax.tick_params(top='off', right='off', pad=2.5)
    ax.set_ylim([0, 0.20])
    ax.set_yticks([0.0, 0.1, 0.2])
ax1.tick_params(labelbottom='off')


props = dict(facecolor='w', lw=0.5, boxstyle='square, pad=0.15')

eq = False

if eq:
    y = ax1.get_ylim()[1]+0.0112#0.2112
    # 18.47, 0.2112
    ax1.text(15.00, y, r'Cosine:  $\Omega_{\mathrm{mod}}(t) = \Omega(t) \cos(2\pi f_{\mathrm{mod}}t)$', 
             clip_on=False,  fontsize=6)

    # 9.165, 0.2112
    ax2.text(6.80, y, r'SSB:  $\Omega_{\mathrm{mod}}(t) = \Omega(t) \cos(2\pi f_{\mathrm{mod}}t) \pm \mathcal{H}[\Omega(t)] \sin(2\pi f_{\mathrm{mod}}t)$', 
             clip_on=False, fontsize=6)

ax1.text(0.0, 0.207, 'Cosinusoidal modulation', alpha=0.75, fontsize=6, clip_on=False)
ax2.text(0.0, 0.207, 'Single-sideband modulation', alpha=0.75, fontsize=6, clip_on=False)

a_c = np.array([0.8, 0.8, 0.8])

ax2.arrow(0.5, 0.075, 4.0, 0, fc=a_c, ec=a_c, head_length=0.30, head_width=0.01, lw=0.75)
ax2.arrow(0.5, 0.055, 14.0, 0, fc=a_c, ec=a_c, head_length=0.30, head_width=0.01, lw=0.75)
ax2.arrow(0.5, 0.035, 24.0, 0, fc=a_c, ec=a_c, head_length=0.30, head_width=0.01, lw=0.75)

ax2.text(2.0, 0.08, r'$\omega_s$', fontsize=7)

#ax1.text(29, 0.13, 'Cosine', rotation=90, fontsize=7, color='Gray')
#ax2.text(29, 0.11, 'SSB', rotation=90, fontsize=7, color='Gray')


###############################################################################################
#
#       ax3  plots (spectrum reconstruction)
#
###############################################################################################

s_c = 'firebrick'

axes = ['X', 'Y', 'Z']
data = []
# --------------- Read in data ----------------------------------------------------------------
for i in range(3):
    data.append(np.loadtxt('assets/data/ThreeAxes_{}_Bayes_dephasing.txt'.format(axes[i])))
expect1 = np.array([2.220446049250313081e-16,
4.966755649881771095e-02,
1.796403833988369225e-01,
2.885864183078730161e-01,
3.118931560209533549e-01,
3.118086965652168141e-01,
3.046769380042700925e-01,
2.337537135528695575e-01,
1.307519932270677643e-01,
5.852239577515694346e-02,
2.437104106931253789e-02,
7.785252100278361453e-03,
1.398191439469953501e-03,
9.716551234240089485e-05,
2.503733753655978944e-08])

expect2 = np.loadtxt('assets/data/Nov30_SLEP_NW3_dephasing_Bayes_Expected.txt', skiprows=1)                           
    

# --------------- Plot ----------------------------------------------------------------

mk = 2.0; mew = 0.75
colours = ['k', 'k', 'k', s_c]; mfc = 'rosybrown'
line_colours = ['k', 'k', 'k']
symbols = ['s', 'o', 'D']; mfc_s = ['none', 'DimGray', 'none']



break_axis = True


if break_axis:
    
    pos1 = ax3.get_position()
    ax3.axis('off')
    
    y1 = 0.47
    ax3_1 = fig.add_axes([pos1.x0, pos1.y0, pos1.width, pos1.height*y1])
    ax3_2 = fig.add_axes([pos1.x0, pos1.y0+pos1.height*(1-y1), pos1.width, pos1.height*y1])
    
    ax3_1.patch.set_alpha(0.0)
    ax3_2.patch.set_alpha(0.0)
    
    ax3_2.tick_params(top='off', right='off', bottom='off', labelbottom='off', pad=2.5)
    ax3_1.tick_params(top='off', right='off', pad=2.5)
    
    ax3_2.spines['bottom'].set_visible(False)
    ax3_1.spines['top'].set_visible(False)
    
    ax3_1.set_ylim([0, 0.29])
    ax3_1.set_yticks([0.0, 0.1, 0.2])
    ax3_1.set_xlim([0, 3])
    ax3_1.set_xticks([0, 1, 2, 3])
    
    ax3_2.set_ylim([0.71, 1.0])
    ax3_2.set_yticks([0.8, 0.9, 1.0])
    
    d = .010  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax3_2.transAxes, color='k', clip_on=False, lw=0.5)
    ax3_2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax3_2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    
    kwargs.update(transform=ax3_1.transAxes)  # switch to the bottom axes
    ax3_1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax3_1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    
    for i in range(3):  
        ax3_2.plot(data[0][:, 0]/1e3, data[i][:, 1], symbols[i], markersize=mk, mec=colours[i], mfc=mfc_s[i], 
                 mew=mew, color=colours[i], label=r'$P\, (\uparrow_{})$'.format(axes[i].lower()))
        ax3_2.plot(data[0][:, 0]/1e3, data[i][:, 1], '-', color=line_colours[i], lw=0.75, alpha=0.20)
    
    ax3_1.plot(data[0][:, 0]/1e3, data[i][:, 1], symbols[i], markersize=mk, mec=colours[i], mfc='none', 
                 mew=mew, color=colours[i])
    ax3_1.plot(data[0][:, 0]/1e3, data[i][:, 1], '-', color=colours[i], lw=0.75, alpha=0.25)
    ax3_1.plot(data[0][:, 0]/1e3, expect2, color=colours[3], lw=0.75, alpha=0.75, label=r'$\mathcal{S}_{\mathrm{target}}$')

    
    reconstructed = ( 1 + data[0][:, 1] - data[1][:, 1] -(1 - data[2][:, 1]) ) / 2 
    
    ax3_1.plot(data[0][:, 0]/1e3, reconstructed, 'o', color=colours[3], markersize=2.5,
                mew=1.0, mec=colours[3], mfc=mfc, label=r'$\mathcal{S}_{\mathrm{reconstructed}}$')
    #ax3_2.plot(data[0][:, 0]/1e3, reconstructed, 'o', color=colours[3], markersize=2.5,
    #            mew=1.0, mec=colours[3], mfc=mfc, label=r'$P(\uparrow_z | \beta_{\Omega})$')
    #ax3_2.plot(data[0][:, 0]/1e3, expect2, color=colours[3], lw=0.75, alpha=0.75, label=r'$\hat{P}(\uparrow_z | \beta_{\Omega})$')

    

    ax3_1.set_ylabel(r'$P \, (\uparrow_{x, y ,z})$')
    ax3_1.yaxis.set_label_coords(-0.094, 1.1)

    
    ax3_1.set_xlabel(r'$\omega_{\mathrm{sid}}/2\pi$ (kHz)', labelpad=1.0)

    #ax3_1.text(1.28, 0.27, r'$\beta_z^{\mathrm{rms}} = 0.5 \Omega_{x}$', clip_on=False)

    ax3_1.legend(loc='upper right', frameon=False, fontsize=8, handlelength=1.5,
                  handletextpad=0.1, bbox_to_anchor=(0.99, 1.54)
                  )


    ax3_2.legend(loc='upper left', ncol=5, fontsize=7, handlelength=1.20, frameon=False,
                     columnspacing=0.75, handletextpad=0.05, bbox_to_anchor=(0.16, 1.12))





###############################################################################################
#
#       ax4  plots (spectrum reconstruction)
#
###############################################################################################

# ------------- Read in the data --------------------------------------------------------------

multitaper = np.loadtxt('assets/data/multitaper_estimate.txt')
bayes = np.loadtxt('assets/data/bayes_estimate')

b_scale = 10.75
s_scale = 3.5

amp = 3.3329923e-4
true_spectrum = np.zeros(13*20); true_spectrum[:7*20] = amp
f = np.arange(13*20) / 20

# alternatively: tones
t_f = np.arange(0.1, 7.1, 0.1)

bandwidth_mean = 2.8223299736137117




# ------------- Plots -------------------------------------------------------------------------

#ax4.plot(f, true_spectrum, color='k', lw=0.75, alpha=0.7)
#ax4.fill_between(f, true_spectrum, 0, color='k', lw=1.5, alpha=0.05)

for tf in t_f:
    ax4.plot([tf, tf], [0, amp], color='k', lw=0.5, alpha=0.25)

ax4.plot(multitaper[:,0], multitaper[:, 1], '^-', color='k', lw=0.75,
             label=r'Multitaper', markersize=3, mec='k', mfc='none', mew=1)
             
ax4.plot(bayes[:, 0]/2/np.pi/1e3, bayes[:, 1], 'o-', color='SeaGreen', lw=0.75,
                  label=r'Bayesian', markersize=2, mec='SeaGreen', mew=0.5)


x = 2; y = 4.75; y1 = 0.1
ax4.plot([x, x+bandwidth_mean], [y*1e-4, y*1e-4], color='k', lw=0.5)
ax4.plot([x, x], [(y-y1)*1e-4, (y+y1)*1e-4], color='k', lw=0.5)
ax4.plot([x + bandwidth_mean, x + bandwidth_mean], [(y-y1)*1e-4, (y+y1)*1e-4], color='k', lw=0.5)
ax4.text(x+0.55, (y+0.1)*1e-4, r'$2\pi NW/\tau$', fontsize=6)

        
ax4.set_xlim([0,12])
ax4.set_ylim([0,0.0006])
#ax.set_yticks([0.0, 0.0002, 0.0004, 0.0006])
ax4.set_yticks([0.0, 0.0002, 0.0004, 0.0006])
ax4.set_yticklabels([0.0, 2.0, 4.0, 6.0])

ax4.tick_params(top='off', right='off', pad=2.5)
#ax.set_yticklabels([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
ax4.legend(loc='upper left', fontsize=7, frameon=False, ncol=2, bbox_to_anchor=(-0.02, 1.22),
           columnspacing=0.5, handletextpad=0.25, handlelength=1.5)
ax4.text(0.2, 0.00005, 'Engineered spectrum', fontsize=7, color='DimGray')

#ax4.grid(alpha=0.1, ls='-')



###############################################################################################
#
#       ax4_inset  plots (schematic filter functions)
#
###############################################################################################

pos1 = ax4.get_position()
pos2 = [pos1.x0 + pos1.width*0.68, pos1.y0 + pos1.height*0.4, pos1.width*0.3, pos1.height*0.55]
ax4_inset = fig.add_axes(pos2)

filters = np.loadtxt('assets/data/FlatTop_1stOrderFilters.txt')
w = filters[:, 0]

c = 'Olive'

for i in range(9):
    ax4_inset.fill_between(w/2/np.pi/1e3, filters[:, i+1], lw=0, color=c, alpha=0.05)
    ax4_inset.plot(w/2/np.pi/1e3, filters[:, i+1], lw=0.5, color=c, alpha=0.5)


ax4_inset.arrow(1.5, 1.85, 7.4, 0, lw=0.75, color='DimGray', head_length=0.5, head_width=0.1)
ax4_inset.text(4.8, 1.95, r'$\omega_{s}$', fontsize=6)


ax4_inset.tick_params(top='off', right='off', labelsize=6, pad=2.5)

ax4_inset.set_xlim(0, 12)
ax4_inset.set_xticks([0, 4, 8, 12])

ax4_inset.set_ylim([0, 2.5])
ax4_inset.set_yticks([0, 1, 2])

ax4_inset.set_xlabel(r'$\omega/2\pi$ (kHz)', labelpad=1.05)
ax4_inset.set_ylabel(r'$F_{k=1}(\omega)$', labelpad=1.05)




###############################################################################################
#
#       axes label for all plots
#
###############################################################################################


ax2.set_xlabel(r'Tone frequency, $\omega_{\mathrm{sid}}/2\pi$ (kHz)', labelpad=1)
ax2.set_ylabel(r'Probability bright, $P\,(\uparrow_z)$')
ax2.yaxis.set_label_coords(-0.1, 1.05)


ax4.set_xlabel(r'Frequency $\omega/2\pi$ (kHz)', labelpad=1.75)
ax4.set_ylabel(r'$\hat{S}(\omega)$ ($10^{-4}$rad$^2$/Hz)')






fig.text(0.155, 0.95, 'a', fontsize=8, fontweight='bold')
fig.text(0.155, 0.767, 'b', fontsize=8, fontweight='bold')
fig.text(0.155, 0.55, 'c', fontsize=8, fontweight='bold')
fig.text(0.155, 0.235, 'd', fontsize=8, fontweight='bold')






svg_version = '3'
plt_version = '4'
fig_version = '4'



if False:    # Export pdf from Inkscape if there were any changes
    import os
    cmd = 'inkscape assets/svg_schematic_v{v}.svg --export-pdf=assets/schematic_v{v}.pdf'.format(v=svg_version)
    os.system(cmd)

plot = 'assets/plot_v{}.pdf'.format(plt_version)
#schematic = 'assets/schematic_v{}.pdf'.format(svg_version)
fig.savefig('Figure_2_v{}.pdf'.format(fig_version))

# Combine plot with schematic
#pdfMerger([plot, schematic], 'Figure_2_v{}.pdf'.format(fig_version), path='')

plt.close()

