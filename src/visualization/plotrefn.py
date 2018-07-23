
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns; sns.set('paper', palette='colorblind')
import numpy as np
from refnx import analysis
import sys
from scipy.stats.mstats import mquantiles


# In[ ]:


n = len(sys.argv)
import os
cwd = os.getcwd()
figures_dir = cwd + '/../../reports/figures/'


# In[6]:


#plotting reflectometry
def plotref(data, gs, offset):
    ax = plt.subplot(gs)
    ax.errorbar(data[0], data[1] * offset, yerr=data[2] * offset, linestyle='', marker='s', markersize=5, 
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
        ax.plot(z, sld + offset, linewidth=2, alpha=0.05)
    ax.text(0.80, 0.9, '(' + label + ')', fontsize=30, transform=ax.transAxes)
    ax.set_xlabel(r'$z$/Å')
    ax.set_ylabel(r'SLD/$10^{-6}$Å$^{-2}$')
    return plt

def plothist(tohist, gs, label, name):
    ax = plt.subplot(gs)
    a = mquantiles(tohist, prob=[0.025, 0.5, 0.975])
    ax.hist(tohist, density=True, bins=20)
    ax.set_ylabel('PDF({}-$V_h$)'.format(name))
    ax.set_xlabel('{}-$V_h$/Å$^3$'.format(name))
    ax.set_xticks([a[0], a[1], a[2]])
    ax.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])

mpl.rcParams['axes.labelsize']=28
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18

l = [1, 10, 100, 1000, 10000]
fig = plt.figure(figsize=(15, 6))
gs = mpl.gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
k = 0
for i in range(0, int(n)-4, 2):
    data = np.loadtxt(sys.argv[i+1])
    plotref(data, gs[0, 0], l[k])
    k += 1
f = 0
for i in range(1, int(n)-3, 2):
    data = np.loadtxt(sys.argv[i+1])
    plotsld(data, gs[0, 1], f, sys.argv[n-2])
    f += 5
plt.tight_layout()
plt.savefig('{}{}_all_data.png'.format(figures_dir, sys.argv[n-1]), dpi=600)
plt.close()
#plt.show()

