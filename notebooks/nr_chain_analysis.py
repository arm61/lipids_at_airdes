#!/usr/bin/env python
# coding: utf-8

# # NR reflectometry analysis
# 
# This is a custom Python analysis notebook for the analysis of the chain output by the `lipid_nr` notebook. 
# 
# It is first necessary to import the necessary modules for the analysis.

# In[ ]:


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

# In[ ]:


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

# In[ ]:


# The type of lipid being investigated
lipid = sys.argv[1]
length = int(sys.argv[2])
sp = sys.argv[3]
label = sys.argv[4]


# Here we assign the directories that contain the data, as well as where the figures and analysis outputs should be stored. If you directory structure does not match that in the GitHub repository these should be adapted. 

# In[ ]:


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

# In[ ]:


refnx.version.full_version, scipy.version.version


# ## Setup of the processing of the MCMC chain
# 
# For details see the `lipid_xrr` notebook.

# In[ ]:


# Reading datasets into refnx format
dataset1 = ReflectDataset('{}nr_h_sp_{}.dat'.format(data_dir, sp))
dataset2 = ReflectDataset('{}nr_hd_sp_{}.dat'.format(data_dir, sp))

datasets = [dataset1, dataset2]

head = {'C': 10, 'H': 18, 'O': 8, 'N': 1, 'P': 1}
tail = {'C': length * 2, 'D': length * 4 + 2}

head_sl = mv.get_scattering_length(head, neutron=True)
tail_sl = mv.get_scattering_length(tail, neutron=True)

solvent_sld = [0.43, 3.15]
super_sld = [0, 0]
thick_heads = 13.1117
tail_length = 1.54 + 1.265 * length
vols = [200.497, 891.]

lipid_1 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                      reverse_monolayer=True, name='{}_1'.format(lipid))
lipid_2 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                      reverse_monolayer=True, name='{}_2'.format(lipid))

# build the structures
air = SLD(0, '')
des_1 = SLD(solvent_sld[0], '')
des_2 = SLD(solvent_sld[1], '')

structure_lipid_1 = air(0, 0) | lipid_1 | des_1(0, 0)
structure_lipid_2 = air(0, 0) | lipid_2 | des_2(0, 0)

def get_value(file):
    f = open(analysis_dir + lipid + '/' + file + '.txt', 'r')
    for line in f:
        k = line
    if '^' in k:
        l = k.split('$')[1].split('^')[0]
    else:
        l = k.split('$')[1].split('\\pm')[0]
    return float(l)

lipid_1.head_mol_vol.setp(get_value('vh'), vary=False)
lipid_1.tail_mol_vol.setp(get_value('vt'), vary=False)
lipid_1.thick_tails.setp(get_value('tail{}'.format(sp)), vary=True, bounds=(5, tail_length))
lipid_1.rough_head_tail.constraint = structure_lipid_1[-1].rough
lipid_1.rough_preceding_mono.constraint = structure_lipid_1[-1].rough
lipid_1.phih.constraint = 1 - (lipid_1.head_mol_vol * lipid_1.thick_tails / 
                               (lipid_1.tail_mol_vol * lipid_1.thick_heads))
lipid_1.thick_heads.setp(get_value('head'), vary=False)
structure_lipid_1[-1].rough.setp(get_value('rough{}'.format(sp)), vary=True, bounds=(2.5, 6))

lipid_2.head_mol_vol.constraint = lipid_1.head_mol_vol
lipid_2.tail_mol_vol.constraint = lipid_1.tail_mol_vol
lipid_2.thick_tails.constraint = lipid_1.thick_tails
lipid_2.rough_head_tail.constraint = structure_lipid_1[-1].rough
lipid_2.rough_preceding_mono.constraint = structure_lipid_1[-1].rough
lipid_2.phih.constraint = lipid_1.phih
lipid_2.thick_heads.constraint = lipid_1.thick_heads
structure_lipid_2[-1].rough.constraint = structure_lipid_1[-1].rough

model_lipid1 = ReflectModel(structure_lipid_1)
model_lipid1.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid1.bkg.setp(dataset1.y[-2], vary=False)

model_lipid2 = ReflectModel(structure_lipid_2)
model_lipid2.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid2.bkg.setp(dataset2.y[-2], vary=False)

models = [model_lipid1, model_lipid2]
structures = [structure_lipid_1, structure_lipid_2]

# building the global objective
objective_n1 = Objective(model_lipid1, dataset1, transform=Transform('YX4'))
objective_n2 = Objective(model_lipid2, dataset2, transform=Transform('YX4'))

global_objective = GlobalObjective([objective_n1, objective_n2])


# The chain is read in by refnx, and processed to assigned it to the global objective. 

# In[ ]:


chain = refnx.analysis.load_chain('{}/{}/{}_chain_neutron.txt'.format(analysis_dir, lipid, sp))

processed_chain = refnx.analysis.process_chain(global_objective, chain)


# The global objective is printed to check it is accurate.

# In[ ]:


print(global_objective)


# The reflectometry and SLD profile are then plotted. 

# In[ ]:


fig = plt.figure(figsize=(20, 7.5))
gs = mpl.gridspec.GridSpec(1, 3)
lines = ['--', '-']
contrast = ['hDES', 'hdDES']

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
        ax2.plot(zs, sld + i*5, color='k', ls=lines[i], linewidth=2, alpha=0.1)
    ax1.plot(dataset.x, models[i](dataset.x, x_err=dataset.x_err)*(dataset.x)**4 * 10**(i-1), 
                 linewidth=4, color='k', ls=lines[i], label = contrast[i])
    ax1.set_ylabel(r'$Rq^4$/Å$^{-4}$')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$q$/Å$^{-1}$')
    ax2.set_xlabel(r'$z$/Å')
    ax2.set_ylabel(r'SLD/$10^{-6}$Å$^{-2}$')
ax1.legend(bbox_to_anchor=(0., 1.02, 1.57, .102), loc=3,
                ncol=2, mode="expand", borderaxespad=0., frameon=False)
ax2.text(0.80, 0.87, '(' + label + ')', fontsize=44, transform=ax2.transAxes)
plt.tight_layout()
plt.savefig('{}{}_{}n_ref_sld.pdf'.format(figures_dir, lipid, sp))
plt.close()


# In[ ]:


lab = ['scale{}'.format(sp), 'tail{}'.format(sp), 'rough{}'.format(sp), 'scalea{}'.format(sp)]

alpha = 0.05

for i in range(0, len(processed_chain)):
    len_to = int(processed_chain[i].chain.size/5000)
    total_pearsons = open('{}{}/{}_neutron.txt'.format(analysis_dir, lipid, lab[i]), 'w')
    stat, p = scipy.stats.shapiro(processed_chain[i].chain[::len_to])
    if p > alpha:
        print('{} - normal'.format(lab[i]))
        a = mquantiles(processed_chain[i].chain.flatten(), prob=[0.025, 0.5])
        k = [a[1], a[1] - a[0]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        total_pearsons.write(helper.latex_sym(q, e))
        total_pearsons.close()
    else:
        print('{} - not normal'.format(lab[i]))
        a = mquantiles(processed_chain[i].chain.flatten(), prob=[0.025, 0.5, 0.975])
        k = [a[1], a[1] - a[0], a[2] - a[1]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])        
        w = '{:.2f}'.format(k[2])
        total_pearsons.write(helper.latex_asym(q, e, w))
        total_pearsons.close()
    
lab2 = ['solh{}'.format(sp)]
kl = 1 - ((lipid_1.head_mol_vol.value * processed_chain[1].chain.flatten()) / (lipid_1.tail_mol_vol.value * lipid_1.thick_heads.value))
kl = kl * 100
for i in range(0, len(lab2)):
    len_to = int(kl.size/5000)
    total_pearsons = open('{}{}/{}_neutron.txt'.format(analysis_dir, lipid, lab2[i]), 'w')
    stat, p = scipy.stats.shapiro(kl[::len_to])
    if p > alpha:
        print('{} - normal'.format(lab2))
        a = mquantiles(kl.flatten(), prob=[0.025, 0.5])
        k = [a[1], a[1] - a[0]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        total_pearsons.write(helper.latex_sym(q, e))
        total_pearsons.close()
    else:
        print('{} - not normal'.format(lab2))
        a = mquantiles(kl.flatten(), prob=[0.025, 0.5, 0.975])
        k = [a[1], a[1] - a[0], a[2] - a[1]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])        
        w = '{:.2f}'.format(k[2])
        total_pearsons.write(helper.latex_asym(q, e, w))
        total_pearsons.close()


# The corner plots for each of the surface pressures is produced, these are presented in the ESI.

# In[ ]:


# plotting pdfs
import corner

mpl.rcParams['axes.labelsize']=22
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.edgecolor'] = 'k'


label=['$d_t$/Å', r'ϕ$_h/\times10^{-2}$', 'σ$_{t,h,s}$/Å']

new_flat = np.zeros((processed_chain[0].chain.size, 3))

new_flat[:, 0] = processed_chain[1].chain.flatten()
new_flat[:, 1] = (1 - ((get_value('vh') * processed_chain[1].chain.flatten()) / (
    get_value('head') * get_value('vt')))) * 100
new_flat[:, 2] = processed_chain[2].chain.flatten()

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}_{}n_all_corner.pdf'.format(figures_dir, lipid, sp))
plt.close()


# ## Bibliography
# 
# 1. Andrew Nelson, Stuart Prescott, Isaac Gresham, & Andrew R. McCluskey. (2018, August 3). refnx/refnx: v0.0.17 (Version v0.0.17). Zenodo. http://doi.org/10.5281/zenodo.1345464

# In[ ]:




