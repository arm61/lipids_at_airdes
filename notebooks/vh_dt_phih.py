#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Standard libraries to import
from __future__ import division
import numpy as np 
import scipy
from scipy.stats.mstats import mquantiles
import matplotlib as mpl
import matplotlib.pyplot as plt

# The refnx library
import refnx
from refnx.reflect import structure, ReflectModel, SLD
from refnx.dataset import ReflectDataset
from refnx.analysis import Transform, CurveFitter, Objective, GlobalObjective, Parameter

# The custom class to constain the monolayer model. 
import sys
sys.path.insert(0, '../src/models')
import mol_vol as mv
sys.path.insert(0, '../src/tools')
import helper


# In[2]:


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
mpl.rcParams['legend.fontsize'] = 32


# In[3]:


figures_dir = '../reports/figures/'


# In[11]:


fig = plt.figure(figsize=(20, 22))
gs = mpl.gridspec.GridSpec(4, 2) 
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[1, 0])
ax3 = plt.subplot(gs[2, 0])
ax4 = plt.subplot(gs[3, 0])
ax5 = plt.subplot(gs[:2, 1])
ax6 = plt.subplot(gs[2:, 1])
dlpc_vh = np.loadtxt('{}{}_vh.txt'.format(figures_dir, 'dlpc'))
weights = np.ones_like(dlpc_vh)/float(
    len(dlpc_vh))
ax1.hist(dlpc_vh, bins=50, histtype='stepfilled', 
         color='k', weights=weights)
ax1.text(0.02, 0.96, '(a)', fontsize=44, transform=ax1.transAxes, ha='left', va='top')
ax1.set_ylabel('PDF(DLPC-$V_h$)')
ax1.set_xlabel('DLPC-$V_h$/Å$^3$')
a = mquantiles(dlpc_vh, prob=[0.025, 0.5, 0.975])
ax1.set_xticks([a[0], a[1], a[2]])
ax1.set_xlim([np.min(dlpc_vh)-0.01, 
              np.max(dlpc_vh)+0.01])
ax1.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])
dmpc_vh = np.loadtxt('{}{}_vh.txt'.format(figures_dir, 'dmpc'))
weights = np.ones_like(dmpc_vh)/float(
    len(dmpc_vh))
ax2.hist(dmpc_vh, bins=50, histtype='stepfilled', 
         color='k', weights=weights)
ax2.text(0.02, 0.96, '(b)', fontsize=44, transform=ax2.transAxes, ha='left', va='top')
ax2.set_ylabel('PDF(DMPC-$V_h$)')
ax2.set_xlabel('DMPC-$V_h$/Å$^3$')
a = mquantiles(dmpc_vh, prob=[0.025, 0.5, 0.975])
ax2.set_xticks([a[0], a[1], a[2]])
ax2.set_xlim([np.min(dmpc_vh)-0.01, 
              np.max(dmpc_vh)+0.01])
ax2.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])
dppc_vh = np.loadtxt('{}{}_vh.txt'.format(figures_dir, 'dppc'))
weights = np.ones_like(dppc_vh)/float(
    len(dppc_vh))
ax3.hist(dppc_vh, bins=50, histtype='stepfilled', 
         color='k', weights=weights)
ax3.text(0.02, 0.96, '(c)', fontsize=44, transform=ax3.transAxes, ha='left', va='top')
ax3.set_ylabel('PDF(DPPC-$V_h$)')
ax3.set_xlabel('DPPC-$V_h$/Å$^3$')
a = mquantiles(dppc_vh, prob=[0.025, 0.5, 0.975])
ax3.set_xticks([a[0], a[1], a[2]])
ax3.set_xlim([np.min(dppc_vh)-0.01, 
              np.max(dppc_vh)+0.01])
ax3.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])
dmpg_vh = np.loadtxt('{}{}_vh.txt'.format(figures_dir, 'dmpg'))
weights = np.ones_like(dmpg_vh)/float(
    len(dmpg_vh))
ax4.hist(dmpg_vh, bins=50, histtype='stepfilled', 
         color='k', weights=weights)
ax4.text(0.02, 0.96, '(d)', fontsize=44, transform=ax4.transAxes, ha='left', va='top')
ax4.set_ylabel('PDF(DMPG-$V_h$)')
ax4.set_xlabel('DMPG-$V_h$/Å$^3$')
a = mquantiles(dmpg_vh, prob=[0.025, 0.5, 0.975])
ax4.set_xticks([a[0], a[1], a[2]])
ax4.set_xlim([np.min(dmpg_vh)-0.01, 
              np.max(dmpg_vh)+0.01])
ax4.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])
dlpc_t = np.loadtxt('{}{}_tailplot.txt'.format(figures_dir, 'dlpc'))
dmpc_t = np.loadtxt('{}{}_tailplot.txt'.format(figures_dir, 'dmpc'))
dppc_t = np.loadtxt('{}{}_tailplot.txt'.format(figures_dir, 'dppc'))
dmpg_t = np.loadtxt('{}{}_tailplot.txt'.format(figures_dir, 'dmpg'))
ax5.plot(dlpc_t[0], dlpc_t[1], c='k', ls='-', marker='o', ms=20, mfc='none', mew=4)
ax5.plot(dmpc_t[0], dmpc_t[1], c='k', ls='--', marker='s', ms=20, mfc='none', mew=4)
ax5.plot(dppc_t[0], dppc_t[1], c='k', ls='-.', marker='v', ms=20, mfc='none', mew=4)
ax5.plot(dmpg_t[0], dmpg_t[1], c='k', ls=':', marker='x', ms=20, mew=4)
ax5.text(0.98, 0.96, '(e)', fontsize=44, transform=ax5.transAxes, ha='right', va='top')
ax5.set_xlabel(r'Surface Pressure/mNm$^{-1}$')
ax5.set_ylabel(r'$d_t$/Å')
ax5.set_ylim([5, 20])
ax5.set_xticks(range(15, 45, 5))
dlpc_s = np.loadtxt('{}{}_solhplot.txt'.format(figures_dir, 'dlpc'))
dmpc_s = np.loadtxt('{}{}_solhplot.txt'.format(figures_dir, 'dmpc'))
dppc_s = np.loadtxt('{}{}_solhplot.txt'.format(figures_dir, 'dppc'))
dmpg_s = np.loadtxt('{}{}_solhplot.txt'.format(figures_dir, 'dmpg'))
ax6.plot(dlpc_s[0], dlpc_s[1], c='k', ls='-', marker='o', ms=20, mfc='none', mew=4)
ax6.plot(dmpc_s[0], dmpc_s[1], c='k', ls='--', marker='s', ms=20, mfc='none', mew=4)
ax6.plot(dppc_s[0], dppc_s[1], c='k', ls='-.', marker='v', ms=20, mfc='none', mew=4)
ax6.plot(dmpg_s[0], dmpg_s[1], c='k', ls=':', marker='x', ms=20, mew=4)
ax6.text(0.98, 0.96, '(f)', fontsize=44, transform=ax6.transAxes, ha='right', va='top')
ax6.set_xlabel(r'Surface Pressure/mNm$^{-1}$')
ax6.set_ylabel(r'$\phi_h$/$\times 10^{-2}$')
ax6.set_ylim([35, 85])
ax6.set_xticks(range(15, 45, 5))
plt.tight_layout()
plt.savefig('{}vh_dt_phih.pdf'.format(figures_dir))
plt.close()


# In[ ]:




