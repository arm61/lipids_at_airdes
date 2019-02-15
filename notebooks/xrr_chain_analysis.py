#!/usr/bin/env python
# coding: utf-8

# # X-ray reflectometry analysis
# 
# This is a custom Python analysis notebook for the analysis of the chain output by the `lipid_xrr` notebook. 
# 
# It is first necessary to import the necessary modules for the analysis.

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


# These are parameters to make the plots pretty.

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


# When running the `Makefile` in the top directory of this ESI, a this notebook is converted to a Python script and running for four different lipids, each at four surface pressures. The necessary variables are assigned here. 

# In[4]:


# The type of lipid being investigated
lipid = sys.argv[1]
length = int(sys.argv[2])
sp1 = sys.argv[3]
sp2 = sys.argv[4]
sp3 = sys.argv[5]
sp4 = sys.argv[6]
label = sys.argv[7]
sps = [sp1, sp2, sp3, sp4]


# Here we assign the directories that contain the data, as well as where the figures and analysis outputs should be stored. If you directory structure does not match that in the GitHub repository these should be adapted. 

# In[5]:


# Relative directory locations
data_dir = '../data/processed/{}_'.format(lipid)
figures_dir = '../reports/figures/'
analysis_dir = '../output/'


# In order for the analysis to be exactly reproducible the same package versions must be used. The conda packaging manager, and pip, can be used to ensure this is the case. The versions of refnx and scipy used original are:
# 
# ```
# refnx.version.full_version = 0.0.17
# scipy.version.version = 1.1.0
# ```

# In[6]:


refnx.version.full_version, scipy.version.version


# ## Setup of the processing of the MCMC chain
# 
# For details see the `lipid_xrr` notebook.

# In[19]:


# Reading datasets into refnx format
dataset1 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp1)), 0.6)
dataset2 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp2)), 0.6)
dataset3 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp3)), 0.6)
dataset4 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp4)), 0.6)

datasets = [dataset1, dataset2, dataset3, dataset4]

if lipid == 'dmpg':
    head = {'C': 8, 'H': 12, 'O': 10, 'Na': 1, 'P': 1}
else:
    head = {'C': 10, 'H': 18, 'O': 8, 'N': 1, 'P': 1}
tail = {'C': length * 2, 'H': length * 4 + 2}

head_sl = mv.get_scattering_length(head)
tail_sl = mv.get_scattering_length(tail)

thick_heads = 11.057
if lipid == 'dlpc':
    vols = [330., 667.]
if lipid == 'dmpg':
    vols = [330., 779.]
if lipid == 'dmpc':
    vols = [330., 779.]
if lipid == 'dppc':
    vols = [330., 891.]

tail_length = 1.54 + 1.265 * length

lipid1 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}1'.format(lipid))
lipid2 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}2'.format(lipid))
lipid3 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}3'.format(lipid))
lipid4 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}4'.format(lipid))

air = SLD(0, 'air')
des = SLD(10.8, 'des')

structure_lipid1 = air(0, 0) | lipid1 | des(0, 3.3)
structure_lipid2 = air(0, 0) | lipid2 | des(0, 3.3)
structure_lipid3 = air(0, 0) | lipid3 | des(0, 3.3)
structure_lipid4 = air(0, 0) | lipid4 | des(0, 3.3)

lipid1.head_mol_vol.setp(vary=True, bounds=(vols[0]*0.8, vols[0]*1.2))
lipid1.tail_mol_vol.setp(vary=True, bounds=(vols[1]*0.8, vols[1]*1.2))
lipid1.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid1.rough_head_tail.constraint = structure_lipid1[-1].rough
lipid1.rough_preceding_mono.constraint = structure_lipid1[-1].rough
lipid1.phih.constraint = 1 - (lipid1.head_mol_vol /  lipid1.tail_mol_vol) * (
    lipid1.thick_tails / lipid1.thick_heads)
lipid1.thick_heads.setp(vary=True, bounds=(6, 20))
structure_lipid1[-1].rough.setp(vary=True, bounds=(2.5, 6))

lipid2.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid2.rough_head_tail.constraint = structure_lipid2[-1].rough
lipid2.rough_preceding_mono.constraint = structure_lipid2[-1].rough
lipid2.thick_heads.setp(vary=True, bounds=(6, 20))
lipid2.phih.constraint = 1 - (lipid2.head_mol_vol / lipid2.tail_mol_vol) * (
    lipid2.thick_tails / lipid2.thick_heads)
structure_lipid2[-1].rough.setp(vary=True, bounds=(2.5, 6))

lipid3.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid3.rough_head_tail.constraint = structure_lipid3[-1].rough
lipid3.rough_preceding_mono.constraint = structure_lipid3[-1].rough
lipid3.thick_heads.setp(vary=True, bounds=(6, 20))
lipid3.phih.constraint = 1 - (lipid3.head_mol_vol / lipid3.tail_mol_vol) * (
    lipid3.thick_tails / lipid3.thick_heads)
structure_lipid3[-1].rough.setp(vary=True, bounds=(2.5, 6))

lipid4.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid4.rough_head_tail.constraint = structure_lipid4[-1].rough
lipid4.rough_preceding_mono.constraint = structure_lipid4[-1].rough
lipid4.thick_heads.setp(vary=True, bounds=(6, 20))
lipid4.phih.constraint = 1 - (lipid4.head_mol_vol / lipid4.tail_mol_vol) * (
    lipid4.thick_tails / lipid4.thick_heads)
structure_lipid4[-1].rough.setp(vary=True, bounds=(2.5, 6))

lipids = [lipid1, lipid2, lipid3, lipid4]
structures = [structure_lipid1, structure_lipid2, structure_lipid3, structure_lipid4]
lipids = mv.set_constraints(lipids, vary_tails=True)

model_lipid1 = ReflectModel(structure_lipid1)
model_lipid1.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid1.bkg.setp(dataset1.y[-1], vary=False)

model_lipid2 = ReflectModel(structure_lipid2)
model_lipid2.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid2.bkg.setp(dataset2.y[-1], vary=False)

model_lipid3 = ReflectModel(structure_lipid3)
model_lipid3.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid3.bkg.setp(dataset3.y[-1], vary=False)

model_lipid4 = ReflectModel(structure_lipid4)
model_lipid4.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid4.bkg.setp(dataset4.y[-1], vary=False)

models = [model_lipid1, model_lipid2, model_lipid3, model_lipid4]

objective1 = Objective(model_lipid1, dataset1, transform=Transform('YX4'))
objective2 = Objective(model_lipid2, dataset2, transform=Transform('YX4'))
objective3 = Objective(model_lipid3, dataset3, transform=Transform('YX4'))
objective4 = Objective(model_lipid4, dataset4, transform=Transform('YX4'))

global_objective = GlobalObjective([objective1, objective2, objective3, objective4])


# The chain is read in by refnx, and processed to assigned it to the global objective. 

# In[25]:


chain = refnx.analysis.load_chain('{}/{}/chain.txt'.format(analysis_dir, lipid))

processed_chain = refnx.analysis.process_chain(global_objective, chain)


# The global objective is printed to check it is accurate.

# In[26]:


print(global_objective)


# Using the probability distribution functions from the processed chain, the PDFs for the tail layer thickness and the solvent content of the headgroup are defined. 

# In[33]:


tail1 = processed_chain[2].chain
tail2 = processed_chain[7].chain
tail3 = processed_chain[10].chain
tail4 = processed_chain[13].chain

solh1 = 1 - (processed_chain[4].chain / processed_chain[3].chain) * (
    tail1 / processed_chain[1].chain)
solh2 = 1 - (processed_chain[4].chain / processed_chain[3].chain) * (
    tail2 / processed_chain[1].chain)
solh3 = 1 - (processed_chain[4].chain / processed_chain[3].chain) * (
    tail3 / processed_chain[1].chain)
solh4 = 1 - (processed_chain[4].chain / processed_chain[3].chain) * (
    tail4 / processed_chain[1].chain)


# The reflectometry and SLD profile are then plotted. 

# In[35]:


fig = plt.figure(figsize=(20, 7.5))
gs = mpl.gridspec.GridSpec(1, 3)
lines = [':', '-.', '--', '-']

for i, dataset in enumerate(datasets):
    choose = global_objective.pgen(ngen=100)
    ax1 = plt.subplot(gs[0, 0:2])
    ax2 = plt.subplot(gs[0, 2])
    ax1.errorbar(dataset.x, dataset.y*(dataset.x)**4 * 10**(i-1), 
                 yerr=dataset.y_err*(dataset.x)**4 * 10**(i-1), 
                 linestyle='', marker='s', markersize=7, markeredgecolor='k', 
                 markerfacecolor='k', ecolor='k')
    for pvec in choose:
        global_objective.setp(pvec)
        ax1.plot(dataset.x, models[i](dataset.x, x_err=dataset.x_err)*(dataset.x)**4 * 10**(i-1), 
                 linewidth=4, color='k', ls=lines[i], alpha=0.1)
        zs, sld = structures[i].sld_profile()
        ax2.plot(zs, sld + i*5, color='k', ls=lines[i], linewidth=2, 
                 alpha=0.1)
    ax1.plot(dataset.x, models[i](dataset.x, x_err=dataset.x_err)*(dataset.x)**4 * 10**(i-1), 
             linewidth=4, color='k', ls=lines[i], label = sps[i] + ' mNm$^{-1}$')
    ax1.set_ylabel(r'$Rq^4$/Å$^{-4}$')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$q$/Å$^{-1}$')
    ax2.set_xlabel(r'$z$/Å')
    ax2.set_ylabel(r'SLD/$10^{-6}$Å$^{-2}$')
ax1.legend(bbox_to_anchor=(0., 1.02, 1.57, .102), loc=3,
                ncol=4, mode="expand", borderaxespad=0., frameon=False)
ax2.text(0.80, 0.05, '(' + label + ')', fontsize=44, transform=ax2.transAxes)
plt.tight_layout()
plt.savefig('{}{}_ref_sld.pdf'.format(figures_dir, lipid))
plt.close()


# The PDFs for the head volumes of the lipid and the variation of the tail thickness and solvent content with surface pressure are plotted. 

# In[38]:


fig = plt.figure(figsize=(20, 5.5))
gs = mpl.gridspec.GridSpec(1, 2) 
ax1 = plt.subplot(gs[0, 0])
a = mquantiles(processed_chain[4].chain.flatten(), prob=[0.025, 0.5, 0.975])
weights = np.ones_like(processed_chain[4].chain.flatten())/float(
    len(processed_chain[4].chain.flatten()))
ax1.hist(processed_chain[4].chain.flatten(), bins=50, histtype='stepfilled', 
         color='k', weights=weights)
ax1.set_ylabel('PDF({}-$V_h$)'.format(lipid.upper()))
ax1.set_xlabel('{}-$V_h$/Å$^3$'.format(lipid.upper()))
ax1.set_xticks([a[0], a[1], a[2]])
ax1.set_xlim([np.min(processed_chain[4].chain.flatten())-0.01, 
              np.max(processed_chain[4].chain.flatten())+0.01])
ax1.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])
ax2 = plt.subplot(gs[0, 1])
ax3 = ax2.twinx()
ax2.plot([float(sp1), float(sp2), float(sp3), float(sp4)], 
         [np.average(tail1), np.average(tail2), np.average(tail3), np.average(tail4)], 
          c='k', marker='x', ls='', ms=20)
ax3.plot([float(sp1), float(sp2), float(sp3), float(sp4)],
         np.array([np.average(solh1), np.average(solh2), np.average(solh3), np.average(solh4)]) * 100, 
         c='k', marker='o', ls='', ms=20)

k = np.array([np.average(tail1), np.average(tail2), np.average(tail3), np.average(tail4)])
ax2.set_ylim([np.min(k)-0.2, np.max(k)+0.2])
l = np.array([np.average(solh1), np.average(solh2), np.average(solh3), np.average(solh4)]) * 100
ax3.set_ylim([np.min(l)-2, np.max(l)+2])
ax2.set_xlabel(r'Surface Pressure/mNm$^{-1}$')
ax2.set_ylabel(r'$d_t$/Å')
ax3.set_ylabel(r'$\phi_h$/$\times 10^{-2}$')

ax2.yaxis.label.set_color('k')
ax3.yaxis.label.set_color('k')
ax2.text(0.88, 0.07, '(' + label + ')', fontsize=44, transform=ax2.transAxes)
ax2.tick_params(axis='y', colors='k')
ax3.tick_params(axis='y', colors='k')
plt.tight_layout()
plt.savefig('{}{}_vh_dt_phi.pdf'.format(figures_dir, lipid))
plt.close()


# Each of the variables is output to a text file, so that they may be easily imported into the final document if necessary. 

# In[40]:


lab = ['scale{}'.format(sp1), 'head', 'tail{}'.format(sp1), 'vt', 'vh', 'rough{}'.format(sp1), 
       'scale{}'.format(sp2), 'tail{}'.format(sp2), 'rough{}'.format(sp2), 
       'scale{}'.format(sp3), 'tail{}'.format(sp3), 'rough{}'.format(sp3), 
       'scale{}'.format(sp4), 'tail{}'.format(sp4), 'rough{}'.format(sp4)]

alpha = 0.05

for i in range(0, len(processed_chain)):
    len_to = int(processed_chain[i].chain.size/5000)
    f_out = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab[i]), 'w')
    stat, p = scipy.stats.shapiro(processed_chain[i].chain[::len_to])
    if p > alpha:
        print('{} - normal'.format(lab[i]))
        a = mquantiles(processed_chain[i].chain, prob=[0.025, 0.5])
        k = [a[1], a[1] - a[0]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        f_out.write(helper.latex_sym(q, e))
        f_out.close()
    else:
        print('{} - not normal'.format(lab[i]))
        a = mquantiles(processed_chain[i].chain, prob=[0.025, 0.5, 0.975])
        k = [a[1], a[1] - a[0], a[2] - a[1]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        w = '{:.2f}'.format(k[2])
        f_out.write(helper.latex_asym(q, e, w))
        f_out.close()
    
lab = ['solh{}'.format(sp1), 'solh{}'.format(sp2), 'solh{}'.format(sp3), 'solh{}'.format(sp4)]
kl = [solh1, solh2, solh3, solh4]
for i in range(0, len(lab)):
    len_to = int(kl[i].size/5000)
    f_out = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab[i]), 'w')
    stat, p = scipy.stats.shapiro(kl[i][::len_to])
    if p > alpha:
        print('{} - normal'.format(lab[i]))
        a = mquantiles(kl[i], prob=[0.025, 0.5])
        k = [a[1]*100, (a[1] - a[0])*100]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        f_out.write(helper.latex_sym(q, e))
        f_out.close()
    else:
        print('{} - not normal'.format(lab[i]))
        a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
        k = [a[1]*100, (a[1] - a[0])*100, (a[2] - a[1])*100]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        w = '{:.2f}'.format(k[2])
        f_out.write(helper.latex_asym(q, e, w))
        f_out.close()


# The corner plots for each of the surface pressures is produced, these are presented in the ESI.

# In[43]:


# plotting pdfs
import corner

mpl.rcParams['axes.labelsize']=22
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.edgecolor'] = 'k'


label=['$V_t$/Å$^3$', '$V_h$/Å$^3$', '$d_h$/Å', '$d_t$/Å', r'ϕ$_h/\times10^{-2}$', 'σ$_{t,h,s}$/Å']

new_flat = np.zeros((processed_chain[0].chain.size, 6))

new_flat[:, 0] = list(processed_chain[3].chain.flatten())
new_flat[:, 1] = list(processed_chain[4].chain.flatten())
new_flat[:, 3] = list(processed_chain[2].chain.flatten())
new_flat[:, 5] = list(processed_chain[5].chain.flatten())
new_flat[:, 2] = list(processed_chain[1].chain.flatten())
new_flat[:, 4] = list(solh1.flatten() * 100)

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}1_all_corner.pdf'.format(figures_dir, lipid))
plt.close()


new_flat = np.zeros((processed_chain[0].chain.size, 6))

new_flat[:, 0] = list(processed_chain[3].chain.flatten())
new_flat[:, 1] = list(processed_chain[4].chain.flatten())
new_flat[:, 3] = list(processed_chain[7].chain.flatten())
new_flat[:, 5] = list(processed_chain[8].chain.flatten())
new_flat[:, 2] = list(processed_chain[1].chain.flatten())
new_flat[:, 4] = list(solh2.flatten() * 100)

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}2_all_corner.pdf'.format(figures_dir, lipid))
plt.close()


new_flat = np.zeros((processed_chain[0].chain.size, 6))

new_flat[:, 0] = list(processed_chain[3].chain.flatten())
new_flat[:, 1] = list(processed_chain[4].chain.flatten())
new_flat[:, 3] = list(processed_chain[10].chain.flatten())
new_flat[:, 5] = list(processed_chain[11].chain.flatten())
new_flat[:, 2] = list(processed_chain[1].chain.flatten())
new_flat[:, 4] = list(solh3.flatten() * 100)

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}3_all_corner.pdf'.format(figures_dir, lipid))
plt.close()


new_flat = np.zeros((processed_chain[0].chain.size, 6))

new_flat[:, 0] = list(processed_chain[3].chain.flatten())
new_flat[:, 1] = list(processed_chain[4].chain.flatten())
new_flat[:, 3] = list(processed_chain[13].chain.flatten())
new_flat[:, 5] = list(processed_chain[14].chain.flatten())
new_flat[:, 2] = list(processed_chain[1].chain.flatten())
new_flat[:, 4] = list(solh4.flatten() * 100)                      

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}4_all_corner.pdf'.format(figures_dir, lipid))
plt.close()


# ## Bibliography
# 
# 1. Andrew Nelson, Stuart Prescott, Isaac Gresham, & Andrew R. McCluskey. (2018, August 3). refnx/refnx: v0.0.17 (Version v0.0.17). Zenodo. http://doi.org/10.5281/zenodo.1345464

# In[ ]:




