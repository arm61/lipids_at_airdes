
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns; sns.set('paper', palette='colorblind')
import numpy as np
from refnx import analysis
import sys
from scipy.stats.mstats import mquantiles


# In[4]:


n = len(sys.argv)
import os
cwd = os.getcwd()
figures_dir = cwd + '/../../reports/figures/'
output_dir = '../../output/'
lab = sys.argv[-1].lower()
sp1 = int(sys.argv[1][-10:-8])
sp2 = int(sys.argv[3][-10:-8])
sp3 = int(sys.argv[5][-10:-8])
sp4 = int(sys.argv[7][-10:-8])

sps = [sp1, sp2, sp3, sp4]

def get_value(file, sp):
    f = open(output_dir + lab + '/' + file + sp + '.txt', 'r')
    for line in f:
        k = line
    l = k.split('$')[1].split('^')[0]
    m = k.split('$')[1].split('^')[1].split('+')[1].split('}')[0]
    n = k.split('$')[1].split('^')[1].split('+')[1].split('-')[1][:-1]
    return (float(l), float(m), float(n))

tails = np.zeros((4))
sols = np.zeros((4))
tails_n = np.zeros_like(sols)
sols_n = np.zeros_like(sols)
tails_p = np.zeros_like(sols)
sols_p = np.zeros_like(sols)

for i in range(0, len(sps)):
    tails[i] = get_value('tail', str(sps[i]))[0]
    sols[i] = get_value('solh', str(sps[i]))[0]
    tails_n[i] = get_value('tail', str(sps[i]))[1]
    sols_n[i] = get_value('solh', str(sps[i]))[1]
    tails_p[i] = get_value('tail', str(sps[i]))[2]
    sols_p[i] = get_value('solh', str(sps[i]))[2]


# In[9]:


#plotting reflectometry
def plotref(data, gs, offset):
    ax = plt.subplot(gs)
    ax.errorbar(data[0], data[1] * offset, yerr=data[2] * offset, linestyle='', marker='s', markersize=7, 
                markeredgecolor='k', markerfacecolor='k', ecolor='k')
    ax.plot(data[0], data[3] * offset, linewidth=4)
    for i in range(4, data.shape[0]):
        ax.plot(data[0], data[i] * offset, color='k', linewidth=2, alpha=0.005)
    ax.set_ylabel(r'$Rq^4$/Å$^{-4}$')
    ax.set_yscale('log')
    ax.set_xlabel(r'$q$/Å$^{-1}$')
    
def plotsld(data, gs, offset, label):
    ax = plt.subplot(gs)
    z = data[0] - 10
    true_sld = data[1]
    ax.plot(z, true_sld + offset, linewidth=4)
    for i in range(2, data.shape[0]):
        sld = data[i]
        ax.plot(z, sld + offset, color='k', linewidth=2, alpha=0.005)
    ax.set_xlabel(r'$z$/Å')
    ax.set_ylabel(r'SLD/$10^{-6}$Å$^{-2}$')
    ax.text(0.80, 0.05, '(' + label + ')', fontsize=44, transform=ax.transAxes)
    return plt

def plothist(tohist, gs, name):
    ax = plt.subplot(gs)
    a = mquantiles(tohist, prob=[0.025, 0.5, 0.975])
    ax.hist(tohist, bins=50, histtype='stepfilled')
    ax.set_ylabel('PDF({}-$V_h$)'.format(name))
    ax.set_xlabel('{}-$V_h$/Å$^3$'.format(name))
    ax.set_xticks([a[0], a[1], a[2]])
    ax.set_ylim([0, 15000])
    ax.set_yticks([])
    ax.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])

def plotgraph(gs, sp, tails, sols, n1, n2, p1, p2, x1, x2, label):
    def make_patch_spines_invisible(ax):
        ax.set_frame_on(True)
        ax.patch.set_visible(False)
        for sp in ax.spines.values():
            sp.set_visible(False)


    host = plt.subplot(gs)


    par1 = host.twinx()

    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    #par2.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    #make_patch_spines_invisible(par2)
    # Second, show the right spine.
    #par2.spines["right"].set_visible(True)
    p11, = host.plot(sp, tails, c='#0173B2', marker='s', ls='', ms=15)
    p21, = par1.plot(sp, sols, c='#029E73', marker='o', ls='', ms=15)
    
    host.set_xlim(x1, x2)
    host.set_ylim(np.min(tails)-0.5, np.max(tails)+0.5)
    a = 1
    if len(np.arange(np.floor(np.min(tails)-0.5), np.max(tails)+1.5, 1)) > 5:
        a = 2
    host.set_yticks(np.arange(np.floor(np.min(tails)-0.5), np.max(tails)+a+0.5, a))
    par1.set_ylim(np.min(sols)-5, np.max(sols)+5)
    host.set_xticks(np.arange(np.min(sp), np.max(sp)+5, 5))
    #par2.set_ylim(42, 48)

    host.set_xlabel(r'Surface Pressure/mNm$^{-1}$')
    host.set_ylabel(r'$d_t$/Å')
    par1.set_ylabel(r'$\phi_h$/$\times 10^2$')

    host.yaxis.label.set_color(p11.get_color())
    par1.yaxis.label.set_color(p21.get_color())
    host.text(0.88, 0.07, '(' + label + ')', fontsize=44, transform=host.transAxes)

    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p11.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p21.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)
    
mpl.rcParams['axes.labelsize']=44
mpl.rcParams['xtick.labelsize']=32
mpl.rcParams['ytick.labelsize']=32
mpl.rcParams['grid.linestyle'] = ''
mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.facecolor'] = 'w'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.edgecolor'] = 'k'
mpl.rcParams['xtick.bottom'] = True
mpl.rcParams['ytick.left'] = True

l = [1, 10, 100, 1000, 10000]
fig = plt.figure(figsize=(20, 7.5))
gs = mpl.gridspec.GridSpec(1, 3) 
k = 0
for i in range(0, int(n)-4, 2):
    data = np.loadtxt(sys.argv[i+1])
    plotref(data, gs[0, 0:2], l[k])
    k += 1
f = 0
for i in range(1, int(n)-3, 2):
    data = np.loadtxt(sys.argv[i+1])
    plotsld(data, gs[0, 2], f, str(sys.argv[n-2]))
    f += 5
plt.tight_layout()
plt.savefig('{}{}_all_data.pdf'.format(figures_dir, sys.argv[n-1]))
plt.close()

fig = plt.figure(figsize=(20, 5.5))
gs = mpl.gridspec.GridSpec(1, 2) 
j = analysis.curvefitter.load_chain(sys.argv[n-3])
plothist(j[:, :, 4].flatten(), gs[0, 0], sys.argv[n-1])
plotgraph(gs[0, 1], sps, tails, sols, tails_n, sols_n, tails_p, sols_p, np.min(sps)-2, np.max(sps)+2, str(sys.argv[n-2]))
plt.tight_layout()
plt.savefig('{}{}_other_data.pdf'.format(figures_dir, sys.argv[n-1]))
plt.close()
#plt.show()

