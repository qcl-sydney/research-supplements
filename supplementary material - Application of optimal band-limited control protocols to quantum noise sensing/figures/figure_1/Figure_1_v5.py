# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 10:41:39 2017

@author: virginia
"""

from __future__ import division, print_function
import sys
sys.path.append('/home/virginia/Dropbox/Multitaper paper/Figures')


import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter
import matplotlib
import matplotlib.pyplot as plt

from Slepian import DPSS
import numpy as np


# ------------------------ Axes formatting for log-scales --------------------------------------------
def xaxis_major_formatter(x, pos):
    format_str = '%.0f' % x
    return format_str
major_formatter = FuncFormatter(xaxis_major_formatter)

def xaxis_minor_formatter(x, pos):
    if x < 1:
        format_str = '' if x != 0.5 else '%.1f' % x
    elif x > 1 and x < 10:
        format_str = '' if x != 2 and x != 5 else '%.0f' % x
    else:
        format_str = '' if x != 20 and x != 50 else '%.0f' % x
    
    return format_str
minor_formatter = FuncFormatter(xaxis_minor_formatter)


def mm2inch(mm):
    return mm * 0.0393701
    

# General rc parameter    
matplotlib.rcParams['figure.figsize'] = (mm2inch(90), mm2inch(105))
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.size'] = 1.5
matplotlib.rcParams['ytick.major.size'] = 1.5
matplotlib.rcParams['xtick.minor.size'] = 1.0
matplotlib.rcParams['ytick.minor.size'] = 1.0
font = {'family' : 'Arial',
        'size'   : 8}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 0.5
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)



fig = plt.figure()
grid = gridspec.GridSpec(4, 1, height_ratios=[1.0, 0.5, 0.5, 0.5])

top_grid = gridspec.GridSpecFromSubplotSpec(1, 2, grid[0, 0])

ax1 = plt.subplot(top_grid[0, 0])
ax1_2 = plt.subplot(top_grid[0, 1])
ax2 = plt.subplot(grid[1, 0])
ax3 = plt.subplot(grid[2, 0])
ax4 = plt.subplot(grid[3, 0])







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

for ax in [ax1, ax1_2,  ax2, ax3, ax4]:
    pos1 = ax.get_position()
    x1 = 0.035 if ax != ax1_2 else - 0.011
    w1 = 0.90 if ax != ax1_2 else 0.90
    pos2 = [pos1.x0 + x1, pos1.y0 + 0.02,  pos1.width*w1, pos1.height*1] 
    ax.set_position(pos2) 

for ax in [ax1, ax1_2]:
    pos1 = ax.get_position()
    w1 = 1.1 if ax == ax1 else 1
    pos2 = [pos1.x0, pos1.y0 + 0.02,  pos1.width*w1, pos1.height*0.95] 
    ax.set_position(pos2)

y1 = iter([0.10, 0.055, 0])
for ax in [ax2, ax3, ax4]:
    pos1 = ax.get_position()
    pos2 = [pos1.x0, pos1.y0 - next(y1) * 0.65 - 0.0,  pos1.width, pos1.height*1.30] 
    ax.set_position(pos2) 
    

###############################################################################################
#
#       ax1 plots
#
###############################################################################################

dpss = DPSS(100, 4, 4)

c = 'DarkGreen'#(60/255, 100/255, 40/255)
n = 'DarkRed'
    
slepian = dpss[:,3]
ft = np.abs(np.fft.fft(slepian, n=4096) )**2

f_axis = np.linspace(0, 100./2, 180)

ax1.fill_between(f_axis, ft[:180], 0, color=c, alpha=0.15, lw=0)
ax1.plot(f_axis, ft[:180], color=c, lw=.75)
    
ax1.plot([16.5, 16.5], [0, 30], color=n)
ax1.arrow(16, 25, -8, 0.0, head_width=2, head_length=3, fc='Gray', ec='Gray', lw=0.75)
ax1.arrow(17, 25, 8, 0.0, head_width=2, head_length=3, fc='Gray', ec='Gray', lw=0.75)
ax1.text(13, 33, r'$\omega_{\mathrm{sid}}$', fontsize=10) 


ax1.tick_params(top='off', right='off', pad=2.5)

ax1.set_ylim([0, 98])
ax1.set_yticks([0, 40, 80])
ax1.set_yticklabels([0.0, 1.0, 2.0])

ax1.set_xlim([0, 50])
ax1.set_xticks([0, 10, 20, 30, 40, 50])
#ax1.set_xticklabels([0, 4, 8, 12, 16])

# ----------- Axes label ----------------------------------------------------------------------
ax1.set_ylabel(r'$F_{\Omega}(\omega)$ (a.u.)')
ax1.set_xlabel(r'$\omega$ ($1/\tau$)', labelpad=2)

###############################################################################################
#
#       ax1_inset plots
#
###############################################################################################

# shrink ax1 first, then create the inset axis, then position it (please don't judge me on this)

pos1 = ax1.get_position()
pos2 = [pos1.x0, pos1.y0,  pos1.width, pos1.height*1.0] 
ax1.set_position(pos2)

pos1 = ax1.get_position()
w1 = 0.60
x1 = pos1.width*(1-w1)/2
pos2 = [pos1.x0+x1, pos1.y0 + pos1.height/2 + 0.005,  pos1.width*w1, pos1.height*0.45] 
ax1_inset = fig.add_axes(pos2) 

ax1_inset.patch.set_alpha(0.0)


for spine in ['top', 'left', 'right', 'bottom']:
    ax1_inset.spines[spine].set_visible(False)
ax1_inset.tick_params(top='off', left='off', right='off', labelbottom='off', labelleft='off', bottom='off')



# ----------- Data ----------------------------------------------------------------------------
t = np.linspace(0, 150, num=100)

# delt_t = 150./100 -> f_max = 1/2 delt_t = 100 / 2*150

w_sid = np.pi/16.5
noise = np.cos(w_sid*t + 1.3*np.pi)*0.25
slepian = slepian/np.max(np.abs(slepian))*0.9


# ----------- Plots ---------------------------------------------------------------------------
ax1_inset.plot(t, noise - 1.6, color=n, alpha=1, lw=0.75)  
#ax1_inset.fill_between(t, noise-1.5, np.ones(100)*-1.5, color=n, alpha=0.15, lw=0)
   
ax1_inset.fill_between(t, slepian, 0, color=c, alpha=0.15, lw=0)
ax1_inset.plot(t, slepian, color=c, lw=0.75)


# ---------- Tick params and ticks ------------------------------------------------------------
ax1_inset.tick_params(top='off', right='off', labelsize=6)
#ax1_inset.spines['right'].set_visible(False)
#ax1_inset.spines['top'].set_visible(False)

ax1_inset.set_xlim([0, 150])
ax1_inset.set_xticks([0, 37.5, 75, 112.5, 150])
#ax1_inset.set_xticklabels([0,  0.5, 1.0])


ax1_inset.set_ylim([-3, 1])
ax1_inset.set_yticks([-3, -1, 1])
ax1_inset.set_yticklabels([0, 1, 2])
 

# ----------- Axes label ----------------------------------------------------------------------
#ax1_inset.set_xlabel(r'Time, $t$', labelpad=2.5, fontsize=7)

ax1_inset.text(65, -1.3, r'+', fontsize=10)
ax1_inset.text(-15, 0.3, r'$\Omega(t)$', fontsize=6)
ax1_inset.text(-15, -2.7, r'$\beta_{\Omega}(t) = \alpha \cos(\omega_{\mathrm{sid}}t + \varphi)$', fontsize=6)
ax1_inset.text(160, -0.75, r'$t$', fontsize=7)

ax1_inset.arrow(-5, 0.00, 155, 0, head_length=7.5, head_width=0.20, clip_on=False,
                fc='k', ec='k', lw=0.5, overhang=0.2)



###############################################################################################
#
#       ax1_2 plots
#
###############################################################################################


# ----------- Read in data --------------------------------------------------------------------
s_3_1 = np.loadtxt('assets/data/Slepian k=1 NW 3.txt')
s_9_1 = np.loadtxt('assets/data/Slepian k=1 NW 9.txt')
s_15_1 = np.loadtxt('assets/data/Slepian k=1 NW 15.txt')

piTime = 1.09e-3

data = [s_3_1, s_9_1, s_15_1]
nws = [3, 9, 15]

colours = ['k', '#2ca02c', '#1f77b4' ]
symbols = ['o', 'x', 's']
mk = 2.5

pos1 = ax1_2.get_position()
inset_pos = [0.19, 0.13, 0.06]
text_pos = [2.4, 3.9, 3.9]; y0 = 0.18; y1 = 0.045
text = ['NW=3', '9', '15']

x_pulse = np.arange(100)

for i in range(3):
    s = data[i]; c = colours[i]; m = np.argmax(s[1:,1] == 0)
    ax1_2.plot(s[:,3]/1e3, 1-s[:,4], ls='-', color=c, lw=0.75)    
    #ax1_2.plot(s[:m,0]*piTime, 1-s[:m,2], ls='-', color=c, lw=0.75)
    ax1_2.plot(s[:m,0]/1e3, 1-s[:m,1], symbols[i], color=c, markersize=mk,
               mec=c, mew=0.75, mfc='none' if i == 2 else c, label=r'$NW={}$'.format(nws[i]))

    dpss = DPSS(100, 2, nws[i])[:, 1]
    dpss /= np.sum(np.abs(dpss))
    inset = fig.add_axes([pos1.x0+0.65*pos1.width, pos1.y0+inset_pos[i],
                           pos1.width*0.3, pos1.height*0.3])
    
    inset.axis('off')
    inset.fill_between(x_pulse, dpss, color=c, lw=0.0, alpha=0.10)
    inset.plot(dpss, color=c, lw=0.75)
    inset.plot(68, 0.015, symbols[i], color=c, markersize=mk,
               mec=c, mew=0.75, mfc='none' if i == 2 else c)
    inset.plot(88, 0.015, symbols[i], color=c, markersize=mk,
               mec=c, mew=0.75, mfc='none' if i == 2 else c)
    inset.set_ylim([-0.05, 0.05])
    
    ax1_2.text(text_pos[i], y0-i*y1, r'${}$'.format(text[i]), fontsize=7)
    
ax1_2.set_ylim([0, 0.20])
ax1_2.set_yticks([0.0, 0.1,  0.2])

ax1_2.set_xlim([0, 6])
ax1_2.set_xticks([0, 2, 4, 6])

ax1_2.tick_params(top='off', left='off', labelleft='off', labelright='on', pad=2.5)

#ax1_2.legend(loc='upper left', frameon=False, fontsize=7, ncol=2, handlelength=1.2,
#             handletextpad=0.05, columnspacing=0.5, bbox_to_anchor=(-0.05, 1.07), numpoints=1)

ax1_2.set_xlabel(r'$\omega_{\mathrm{sid}}/2\pi$ (kHz)', labelpad=2)
ax1_2.set_ylabel(r'$P\,(\uparrow_z)$')
ax1_2.yaxis.set_label_coords(1.25, 0.5)





###############################################################################################
#
#       ax2, ax3 & ,ax4 plots
#
###############################################################################################


p_c = 'Grey'
s_c = 'RoyalBlue'

# ---------------- Read in data -----------------------------------------------
s_6_1 = np.loadtxt('assets/data/Slepian NW=small order 1 wide.txt')
s_6_3 = np.loadtxt('assets/data/Slepian NW=small order 3 wide.txt')
s_6_5 = np.loadtxt('assets/data/Slepian NW=6 order 5 wide.txt')

prse_1 = np.loadtxt('assets/data/PRSE repetitions 1.txt')
prse_2 = np.loadtxt('assets/data/PRSE repetitions 2.txt')
prse_3 = np.loadtxt('assets/data/PRSE repetitions 3.txt')    


slepian_data = [s_6_1, s_6_3, s_6_5]; prse_data = [prse_1, prse_2, prse_3]


# ----------- Time-domain pulses  -------------------------------------------------------------
dpss_5 = DPSS(100, 6, 6)[:,5];  dpss_3 = DPSS(100, 4, 4)[:,3]; dpss_1 = DPSS(100, 2, 2)[:,1];

orders = [1, 3, 5]
t = np.arange(100)

square_pulse_1 = np.ones(100); square_pulse_1[50:] *= -1
square_pulse_1[:5] = 0; square_pulse_1[-5:] = 0

square_pulse_2 = np.ones(100); square_pulse_2[25:] *= -1; square_pulse_2[50:75] *= -1
square_pulse_2[:2] = 0; square_pulse_2[-2:] = 0

square_pulse_3 = np.ones(100); square_pulse_3[16:] *= -1; square_pulse_3[32:49] *= -1; 
square_pulse_3[65:81] *= -1; 
square_pulse_3[:2] = 0; square_pulse_3[-2:] = 0

square_pulses = [square_pulse_1, square_pulse_2, square_pulse_3]
slepians = [dpss_1, dpss_3, dpss_5]




colours = ['RoyalBlue', 'DarkOrange', 'DarkGreen']; symbols = ['x', '^', '^']
band_colour = (0, 1, 1, 20./255)

piTime = 1.08e-3

mk = 2.5

bumps = [2.62, 5.381, 8.1]
band_ends = [2.0, 3.52, 5.1]
hws = [0.1, 0.1, 0.1]

label_positions = [band_ends[i] + 0.05 for i in range(3)]
orders = [1, 3, 5]

axes = iter([ax2, ax3, ax4])

arrow_c = 'Crimson'

for i in range(0, 3):
    s = slepian_data[i]; c = colours[i]; m = np.argmax(s[1:,1] == 0); ax = next(axes)
    
    delete_points = np.where(abs(s[:m,1] - s[:m,2]) > 0.03 )
    s = np.delete(s, delete_points, axis=0)
    
        
    l1 = ('DPSS pulse', 'Flat-top pulse') if i == 0 else ('', '')
    
    ax.fill_between([0.0, band_ends[i]], [2, 2], [0.0001, 0.0001], color=band_colour)
    #ax.semilogy(s[:m,0]/1e3, s[:m,2], ls='-', color=s_c, lw=0.75)
    ax.semilogy(s[:,3]/1e3, 1 - s[:,4], ls='-', color=s_c, lw=0.75)
    ax.semilogy(s[:m,0]/1e3, s[:m,1], symbols[0], color=s_c, markersize=mk,
            mew=0.75, label=l1[0])
    
    p = prse_data[i]
    m = np.where(p[1:,0] == 0)[0][0]
    #ax.semilogy(p[:m,0]/1e3, 1-p[:m,2], ls='-', color=p_c, lw=0.75)
    ax.semilogy(p[:,3]/1e3, 1 - p[:,4], ls='-', color=p_c, lw=0.75)
    ax.semilogy(p[:m,0]/1e3, 1-p[:m,1], symbols[1], color=p_c, markersize=mk, mfc='none',
            mew=0.75, mec=p_c, label=l1[1])
            
    ax.plot([band_ends[i], band_ends[i]], [0.0001, 2.0], '--', color='Grey', dashes=(1.2, 0.75), lw=0.7, alpha=0.5)
            
    if i != 2:
        ax.arrow(bumps[i], 0.15, 0.0, -0.06, lw=1.0, fc=arrow_c, ec=arrow_c, head_width=hws[i], head_length=0.025,
                 zorder=100)
    else:
        ax.arrow(bumps[i], 0.011, 0.0, 0.007, lw=1.0, fc=arrow_c, ec=arrow_c, head_width=hws[i], head_length=0.01,
                 zorder=100, clip_on=False)
                 
    ax.text(label_positions[i], 0.20, 
            (r'$k={order}$'+'\n'+'$NW={NW}$').format(order=orders[i], NW=orders[i]+1),
            fontsize=6, alpha=0.8)
    
    ax.tick_params(which="both", top='off', right='off', pad=2.5)
    
    pos1 = ax.get_position()
    pos2 = [pos1.x0 + pos1.width/2 + 0.157, pos1.y0 + pos1.height/2 - 0.013,  pos1.width*0.282 , pos1.height*0.55] 
    ax_inset = fig.add_axes(pos2)
    
    slepian = slepians[i]; square_pulse = square_pulses[i]
    slepian /= np.max(np.abs(slepian))
    
    ax_inset.fill_between(t, slepian + 2.2, np.ones(100)*2.2, color=s_c, alpha=0.25, lw=0)
    ax_inset.plot(slepian + 2.2, color=s_c, lw=0.75)
    
    ax_inset.fill_between(t, square_pulse, 0, color=p_c, alpha=0.25, lw=0)
    ax_inset.plot(square_pulse, color=p_c, lw=0.75)
    
    ax_inset.tick_params(top='off', right='off', labelbottom='off', labelleft='off', labelsize=6)

    
    ax_inset.set_ylim([-2, 4])
    
    ax_inset.set_xlabel(r'$t$', labelpad=1.5, fontsize=7)
    ax_inset.xaxis.set_label_coords(0.9, -0.08)
    ax_inset.set_ylabel(r'$\Omega(t)$', labelpad=1, fontsize=7)
    
    #ax_inset.axis('off')
    
    ax_inset.set_xticks([0, 50, 100])
    ax_inset.set_xticklabels([0.0, 0.5, 1.0])
    ax_inset.set_yticks([-2, 1, 4])
    ax_inset.set_yticklabels([0, 3, 6])
    ax_inset.patch.set_alpha(0.0)        
    
        

ax2.tick_params(which='both', labelbottom='off', pad=2.5)
ax3.tick_params(which='both', labelbottom='off', pad=2.5)


for ax in [ax2, ax3, ax4]:
    ax.set_xlim([0.0, 10])
    ax.set_ylim([0.003, 2.0])
    ax.set_yticks([1e-2, 1e-1, 1])
    ax.set_yticklabels([0.01, 0.1, 1.0])
    ax.tick_params(which='minor', labelsize=6)
    ax.tick_params('y', pad=2)

#ax4.set_ylim([0.003, 3])

ax4.xaxis.set_minor_formatter(minor_formatter)
ax4.xaxis.set_major_formatter(major_formatter)



x0 = 0.2; x1 = 1.4; text_w = 2.0
y0 = 1.0; text_y = 0.7; y1 = 0.2

xend = x0+2*x1+text_w

ax4.text(x0+x1+0.1, 0.7, 'Target band', alpha=0.75, fontsize=6)


ax4.plot([x0, x0+x1], [1.0, 1.0], 'k', lw=0.5)
ax4.plot([x0+x1+text_w, xend], [1.0, 1.0], 'k', lw=0.5)

ax4.plot([x0, x0], [y0-y1, y0+y1], 'k', lw=0.5)
ax4.plot([xend, xend], [y0-y1, y0+y1], 'k', lw=0.5)


ax4.set_ylabel(r'Probability bright, $P\,(\uparrow_z)$')
ax4.yaxis.set_label_coords(-0.105, 1.5)

ax4.set_xlabel(r'Tone frequency, $\omega_{\mathrm{sid}}/2\pi$ (kHz)', labelpad=1.1)


                  
ax2.plot(-1, -1, 'k', lw=0.75, label=r'Expected fidelity, $\mathcal{F}_{\mathrm{av}}$')
#ax2.plot(-1, -1, 'k-', lw=0.75, label='Simulation')

ax2.legend(loc='upper center', frameon=False, fontsize=7, handlelength=1.3, ncol=4,
                  handletextpad=0.25, bbox_to_anchor=(0.5, 1.28), labelspacing=0.25,
                  columnspacing=0.90)





fig.text(0.135, 0.965, 'a', fontsize=8, fontweight='bold')
fig.text(0.555, 0.965, 'b', fontsize=8, fontweight='bold')

fig.text(0.135, 0.585, 'c', fontsize=8, fontweight='bold')
fig.text(0.135, 0.405, 'd', fontsize=8, fontweight='bold')
fig.text(0.135, 0.22, 'e', fontsize=8, fontweight='bold')






fig_version = '5'



plt.savefig('Figure_1_v{}.pdf'.format(fig_version))
plt.close()


