
# coding: utf-8

# # X-ray reflectometry analysis
# 
# This is a custom Python analysis notebook for analysing XRR data using the class `VolMono`, as defined in `src/models/mol_vol.py`, and the refnx [1] package. 
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


# When running the `Makefile` in the top directory of this ESI, a this notebook is converted to a Python script and running for four different lipids, each at four surface pressures. The necessary variables are assigned here. 

# In[ ]:


# The type of lipid being investigated
lipid = sys.argv[1]
length = int(sys.argv[2])
sp1 = sys.argv[3]
sp2 = sys.argv[4]
sp3 = sys.argv[5]
sp4 = sys.argv[6]
label = sys.argv[7]


# Here we assign the directories that contain the data, as well as where the figures and analysis outputs should be stored. If you directory structure does not match that in the GitHub repository these should be adapted. 

# In[ ]:


# Relative directory locations
data_dir = '../data/processed/{}/'.format(lipid)
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


# ## Setup of the analysis
# 
# The experimental datafiles are then read in and redefined such that all data after $q = 0.6 Å^{-1}$ is ignored as this is considered background. 

# In[ ]:


# Reading datasets into refnx format
dataset1 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp1)), 0.6)
dataset2 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp2)), 0.6)
dataset3 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp3)), 0.6)
dataset4 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp4)), 0.6)

datasets = [dataset1, dataset2, dataset3, dataset4]


# The scattering lengths for the head and tail components are defined.

# In[ ]:


if lipid == 'dmpg':
    head = {'C': 8, 'H': 12, 'O': 10, 'Na': 1, 'P': 1}
else:
    head = {'C': 10, 'H': 18, 'O': 8, 'N': 1, 'P': 1}
tail = {'C': length * 2, 'H': length * 4 + 2}

head_sl = mv.get_scattering_length(head)
tail_sl = mv.get_scattering_length(tail)


# Initial 'guesses' for a series of parameters are defined.

# In[ ]:


thick_heads = 11.057
chain_tilt = 0.792674
if lipid == 'dlpc':
    vols = [330., 667.]
if lipid == 'dmpg':
    vols = [330., 779.]
if lipid == 'dmpc':
    vols = [330., 779.]
if lipid == 'dppc':
    vols = [330., 891.]


# The length of the carbon tail ($t_l$) is defined based on the Tanford equation.
# 
# $$ t_l = 1.54 + 1.265(n-1) $$
# 
# where $n$ is the length of the carbon chain (e.g. 16 for DPPC).

# In[ ]:


tail_length = 1.54 + 1.265 * length


# The `VolMono` class objects are defined. 

# In[ ]:


lipid1 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                    reverse_monolayer=True, name='{}1'.format(lipid))
lipid2 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                    reverse_monolayer=True, name='{}2'.format(lipid))
lipid3 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                    reverse_monolayer=True, name='{}3'.format(lipid))
lipid4 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                    reverse_monolayer=True, name='{}4'.format(lipid))


# A series of structures for each surface pressure is defined. 

# In[ ]:


air = SLD(0, 'air')
des = SLD(10.8, 'des')

structure_lipid1 = air(0, 0) | lipid1 | des(0, 3.3)
structure_lipid2 = air(0, 0) | lipid2 | des(0, 3.3)
structure_lipid3 = air(0, 0) | lipid3 | des(0, 3.3)
structure_lipid4 = air(0, 0) | lipid4 | des(0, 3.3)


# The variables, and bounds for each of the four surface pressures are setup

# In[ ]:


lipid1.head_mol_vol.setp(vary=True, bounds=(vols[0]*0.8, vols[0]*1.2))
lipid1.tail_mol_vol.setp(vary=True, bounds=(vols[1]*0.8, vols[1]*1.2))
lipid1.tail_length.setp(vary=False)
lipid1.cos_rad_chain_tilt.setp(vary=True, bounds=(0.01, 0.99))
lipid1.rough_head_tail.constraint = structure_lipid1[-1].rough
lipid1.rough_preceding_mono.constraint = structure_lipid1[-1].rough
lipid1.phih.constraint = 1 - (lipid1.head_mol_vol /  lipid1.tail_mol_vol) * (
    lipid1.cos_rad_chain_tilt * lipid1.tail_length / lipid1.thick_heads)
lipid1.thick_heads.setp(vary=True, bounds=(6, 20))
structure_lipid1[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid2.cos_rad_chain_tilt.setp(vary=True, bounds=(0.01, 0.99))
lipid2.rough_head_tail.constraint = structure_lipid2[-1].rough
lipid2.rough_preceding_mono.constraint = structure_lipid2[-1].rough
lipid2.thick_heads.setp(vary=True, bounds=(6, 20))
lipid2.phih.constraint = 1 - (lipid2.head_mol_vol / lipid2.tail_mol_vol) * (
    lipid2.cos_rad_chain_tilt * lipid2.tail_length / lipid2.thick_heads)
structure_lipid2[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid3.cos_rad_chain_tilt.setp(0.57, vary=True, bounds=(0.01, 0.99))
lipid3.rough_head_tail.constraint = structure_lipid3[-1].rough
lipid3.rough_preceding_mono.constraint = structure_lipid3[-1].rough
lipid3.thick_heads.setp(vary=True, bounds=(6, 20))
lipid3.phih.constraint = 1 - (lipid3.head_mol_vol / lipid3.tail_mol_vol) * (
    lipid3.cos_rad_chain_tilt * lipid3.tail_length / lipid3.thick_heads)
structure_lipid3[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid4.cos_rad_chain_tilt.setp(0.57, vary=True, bounds=(0.01, 0.99))
lipid4.rough_head_tail.constraint = structure_lipid4[-1].rough
lipid4.rough_preceding_mono.constraint = structure_lipid4[-1].rough
lipid4.thick_heads.setp(vary=True, bounds=(6, 20))
lipid4.phih.constraint = 1 - (lipid4.head_mol_vol / lipid4.tail_mol_vol) * (
    lipid4.cos_rad_chain_tilt * lipid4.tail_length / lipid4.thick_heads)
structure_lipid4[-1].rough.setp(vary=True, bounds=(2.5, 6))


# The three lipids are constrained such that the tail and head volumes and head thickness are kept constant.

# In[ ]:


lipids = [lipid1, lipid2, lipid3, lipid4]
structures = [structure_lipid1, structure_lipid2, structure_lipid3, structure_lipid4]
lipids = mv.set_constraints(lipids)


# Each model is then associated with a dataset

# In[ ]:


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


# The global objective fitting object is defined and the fitting and MCMC performed. 

# In[ ]:


objective1 = Objective(model_lipid1, dataset1, transform=Transform('YX4'))
objective2 = Objective(model_lipid2, dataset2, transform=Transform('YX4'))
objective3 = Objective(model_lipid3, dataset3, transform=Transform('YX4'))
objective4 = Objective(model_lipid4, dataset4, transform=Transform('YX4'))

global_objective = GlobalObjective([objective1, objective2, objective3, objective4])


# ## Analysis 
# 
# The differential evolution algorithm is used to find optimal parameters, before the MCMC algorithm probes the parameter space for 1000 steps.

# In[ ]:


fitter = CurveFitter(global_objective)
res = fitter.fit('differential_evolution', seed=1)

fitter.sample(200, random_state=1)
fitter.sampler.reset()
res = fitter.sample(1000, nthin=1, random_state=1, f='{}/{}/chain.txt'.format(analysis_dir, lipid))
flatchain = fitter.sampler.flatchain


# The `global_objective` is printed containing information about the models.

# In[ ]:


print(global_objective)


# Using the probability distribution functions, the PDFs for the tail layer thickness and the solvent content of the headgroup are defined. 

# In[ ]:


tail1 = flatchain[:, 2] * lipid1.tail_length.value
tail2 = flatchain[:, 7] * lipid2.tail_length.value
tail3 = flatchain[:, 10] * lipid3.tail_length.value
tail4 = flatchain[:, 13] * lipid4.tail_length.value

solh1 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail1 / flatchain[:, 1])
solh2 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail2 / flatchain[:, 1])
solh3 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail3 / flatchain[:, 1])
solh4 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail4 / flatchain[:, 1])


# The reflectometry and SLD profile are then plotted. 

# In[ ]:


fig = plt.figure(figsize=(20, 7.5))
gs = mpl.gridspec.GridSpec(1, 3)
colorblind = ["#0173B2", "#DE8F05", "#029E73", "#D55E00"]

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
                 linewidth=4, color=colorblind[i], alpha=0.1)
        zs, sld = structures[i].sld_profile()
        ax2.plot(zs, sld + i*5, color=colorblind[i], linewidth=2, alpha=0.1)
    ax1.set_ylabel(r'$Rq^4$/Å$^{-4}$')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$q$/Å$^{-1}$')
    ax2.set_xlabel(r'$z$/Å')
    ax2.set_ylabel(r'SLD/$10^{-6}$Å$^{-2}$')
ax2.text(0.80, 0.05, '(' + label + ')', fontsize=44, transform=ax2.transAxes)
plt.tight_layout()
plt.savefig('{}{}_ref_sld.pdf'.format(figures_dir, lipid))
plt.close()


# The PDFs for the head volumes of the lipid and the variation of the tail thickness and solvent content with surface pressure are plotted. 

# In[ ]:


fig = plt.figure(figsize=(20, 5.5))
gs = mpl.gridspec.GridSpec(1, 2) 
ax1 = plt.subplot(gs[0, 0])
a = mquantiles(flatchain[:, 4], prob=[0.025, 0.5, 0.975])
weights = np.ones_like(flatchain[:, 4])/float(len(flatchain[:, 4]))
ax1.hist(flatchain[:, 4], bins=50, histtype='stepfilled', weights=weights)
ax1.set_ylabel('PDF({}-$V_h$)'.format(lipid.upper()))
ax1.set_xlabel('{}-$V_h$/Å$^3$'.format(lipid.upper()))
ax1.set_xticks([a[0], a[1], a[2]])
ax1.set_xlim([np.min(flatchain[:, 4])-0.01, np.max(flatchain[:, 4])+0.01])
ax1.set_xticklabels(['{:.1f}'.format(a[0]), '{:.1f}'.format(a[1]), '{:.1f}'.format(a[2])])
ax2 = plt.subplot(gs[0, 1])
ax3 = ax2.twinx()
ax2.plot([sp1, sp2, sp3, sp4], 
         [np.average(tail1), np.average(tail2), np.average(tail3), np.average(tail4)], 
          c="#0173B2", marker='s', ls='', ms=15)
ax3.plot([sp1, sp2, sp3, sp4], 
         np.array([np.average(solh1), np.average(solh2), np.average(solh3), np.average(solh4)]) * 100, 
         c="#DE8F05", marker='o', ls='', ms=15)
ax2.set_xlabel(r'Surface Pressure/mNm$^{-1}$')
ax2.set_ylabel(r'$d_t$/Å')
ax3.set_ylabel(r'$\phi_h$/$\times 10^2$')

ax2.yaxis.label.set_color("#0173B2")
ax3.yaxis.label.set_color("#DE8F05")
ax2.text(0.88, 0.07, '(' + label + ')', fontsize=44, transform=ax2.transAxes)
ax2.tick_params(axis='y', colors="#0173B2")
ax3.tick_params(axis='y', colors="#DE8F05")
plt.tight_layout()
plt.savefig('{}{}_vh_dt_phi.pdf'.format(figures_dir, lipid))
plt.close()


# Each of the variables is output to a text file, so that they may be easily imported into the final document if necessary. 

# In[ ]:


lab = ['scale{}'.format(sp1), 'head{}'.format(sp1), 'angle{}'.format(sp1), 'vt', 'vh', 'rough{}'.format(sp1), 
       'scale{}'.format(sp2), 'angle{}'.format(sp2), 'rough{}'.format(sp2), 
       'scale{}'.format(sp3), 'angle{}'.format(sp3), 'rough{}'.format(sp3), 
       'scale{}'.format(sp4), 'angle{}'.format(sp4), 'rough{}'.format(sp4)]
for i in range(0, flatchain.shape[1]):
    f_out = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab[i]), 'w')
    a = mquantiles(flatchain[:, i], prob=[0.025, 0.5, 0.975])
    if 'angle' in lab[i]:
        c = np.rad2deg(np.arccos(a))
        k = [c[1], c[0] - c[1], c[1] - c[2]]
        q = '{:.2f}'.format(k[0])
        w = '{:.2f}'.format(k[1])
        e = '{:.2f}'.format(k[2])
        f_out.write(helper.latex_asym(q, e, w))
    elif 'sol' in lab[i]:
        k = [a[1]*100, (a[1] - a[0])*100, (a[2] - a[1])*100]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        w = '{:.2f}'.format(k[2])
        f_out.write(helper.latex_asym(q, e, w))
    else:
        k = [a[1], a[1] - a[0], a[2] - a[1]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        w = '{:.2f}'.format(k[2])
        f_out.write(helper.latex_asym(q, e, w))
    f_out.close()
            
lab = ['tail{}'.format(sp1), 'tail{}'.format(sp2), 'tail{}'.format(sp3), 'tail{}'.format(sp4)]
kl = [tail1, tail2, tail3, tail4]
for i in range(0, len(lab)):
    f_out = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    f_out.write(helper.latex_asym(q, e, w))
    f_out.close()
    
lab = ['solh{}'.format(sp1), 'solh{}'.format(sp2), 'solh{}'.format(sp3), 'solh{}'.format(sp4)]
kl = [solh1, solh2, solh3, solh4]
for i in range(0, len(lab)):
    f_out = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    k = [a[1]*100, (a[1] - a[0])*100, (a[2] - a[1])*100]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    f_out.write(helper.latex_asym(q, e, w))
    f_out.close()


# The corner plots for each of the surface pressures is produced, these are presented in the ESI.

# In[ ]:


# plotting pdfs
import corner

mpl.rcParams['axes.labelsize']=22
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.edgecolor'] = 'k'


label=['$V_t$/Å$^3$', '$V_h$/Å$^3$', '$d_h$/Å', 'θ$_t$/°', r'ϕ$_h/\times10^{-2}$', 'σ$_{t,h,s}$/Å']

new_flat = np.zeros((flatchain.shape[0], 6))

new_flat[:, 0] = list(flatchain[:, 3].flatten())
new_flat[:, 1] = list(flatchain[:, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, 2].flatten())))
new_flat[:, 5] = list(flatchain[:, 5].flatten())
new_flat[:, 2] = list(flatchain[:, 1].flatten())
new_flat[:, 4] = list(solh1 * 100)

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}1_all_corner.pdf'.format(figures_dir, lipid))
plt.close()

new_flat = np.zeros((flatchain.shape[0], 6))

new_flat[:, 0] = list(flatchain[:, 3].flatten())
new_flat[:, 1] = list(flatchain[:, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, 7].flatten())))
new_flat[:, 5] = list(flatchain[:, 8].flatten())
new_flat[:, 2] = list(flatchain[:, 1].flatten())
new_flat[:, 4] = list(solh2 * 100)

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}2_all_corner.pdf'.format(figures_dir, lipid))
plt.close()

new_flat = np.zeros((flatchain.shape[0], 6))

new_flat[:, 0] = list(flatchain[:, 3].flatten())
new_flat[:, 1] = list(flatchain[:, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, 10].flatten())))
new_flat[:, 5] = list(flatchain[:, 11].flatten())
new_flat[:, 2] = list(flatchain[:, 1].flatten())
new_flat[:, 4] = list(solh3 * 100)


plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}3_all_corner.pdf'.format(figures_dir, lipid))
plt.close()


new_flat = np.zeros((flatchain.shape[0], 6))

new_flat[:, 0] = list(flatchain[:, 3].flatten())
new_flat[:, 1] = list(flatchain[:, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, 13].flatten())))
new_flat[:, 5] = list(flatchain[:, 14].flatten())
new_flat[:, 2] = list(flatchain[:, 1].flatten())
new_flat[:, 4] = list(solh4 * 100)
                      

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}4_all_corner.pdf'.format(figures_dir, lipid))
plt.close()


# ## Bibliography
# 
# 1. Andrew Nelson, Stuart Prescott, Isaac Gresham, & Andrew R. McCluskey. (2018, August 3). refnx/refnx: v0.0.17 (Version v0.0.17). Zenodo. http://doi.org/10.5281/zenodo.1345464
