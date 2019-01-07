#!/usr/bin/env python
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


thick_heads = 12.057
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


lipid1 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}1'.format(lipid))
lipid2 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}2'.format(lipid))
lipid3 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}3'.format(lipid))
lipid4 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                    reverse_monolayer=True, name='{}4'.format(lipid))


# A series of structures for each surface pressure is defined. 

# In[ ]:


air = SLD(0, 'air')
des = SLD(10.8, 'des')

structure_lipid1 = air(0, 0) | lipid1 | des(0, 3.3)
structure_lipid2 = air(0, 0) | lipid2 | des(0, 3.3)
structure_lipid3 = air(0, 0) | lipid3 | des(0, 3.3)
structure_lipid4 = air(0, 0) | lipid4 | des(0, 3.3)


# The variables, and bounds for each of the four surface pressures are setup.

# In[ ]:


lipid1.head_mol_vol.setp(vary=True, bounds=(vols[0]*0.8, vols[0]*1.2))
lipid1.tail_mol_vol.setp(vary=True, bounds=(vols[1]*0.8, vols[1]*1.2))
lipid1.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid1.rough_head_tail.constraint = structure_lipid1[-1].rough
lipid1.rough_preceding_mono.constraint = structure_lipid1[-1].rough
lipid1.phih.constraint = 1 - (lipid1.head_mol_vol /  lipid1.tail_mol_vol) * (
    lipid1.thick_tails / lipid1.thick_heads)
lipid1.thick_heads.setp(vary=True, bounds=(7, 20))
structure_lipid1[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid2.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid2.rough_head_tail.constraint = structure_lipid2[-1].rough
lipid2.rough_preceding_mono.constraint = structure_lipid2[-1].rough
lipid2.thick_heads.setp(vary=True, bounds=(6, 20))
lipid2.phih.constraint = 1 - (lipid2.head_mol_vol / lipid2.tail_mol_vol) * (
    lipid2.thick_tails / lipid2.thick_heads)
structure_lipid2[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid3.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid3.rough_head_tail.constraint = structure_lipid3[-1].rough
lipid3.rough_preceding_mono.constraint = structure_lipid3[-1].rough
lipid3.thick_heads.setp(vary=True, bounds=(6, 20))
lipid3.phih.constraint = 1 - (lipid3.head_mol_vol / lipid3.tail_mol_vol) * (
    lipid3.thick_tails / lipid3.thick_heads)
structure_lipid3[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid4.thick_tails.setp(vary=True, bounds=(5, tail_length))
lipid4.rough_head_tail.constraint = structure_lipid4[-1].rough
lipid4.rough_preceding_mono.constraint = structure_lipid4[-1].rough
lipid4.thick_heads.setp(vary=True, bounds=(6, 20))
lipid4.phih.constraint = 1 - (lipid4.head_mol_vol / lipid4.tail_mol_vol) * (
    lipid4.thick_tails / lipid4.thick_heads)
structure_lipid4[-1].rough.setp(vary=True, bounds=(2.5, 6))


# The three lipids are constrained such that the tail and head volumes and head thickness are kept constant.

# In[ ]:


lipids = [lipid1, lipid2, lipid3, lipid4]
structures = [structure_lipid1, structure_lipid2, structure_lipid3, structure_lipid4]
lipids = mv.set_constraints(lipids, vary_tails=True)


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


# ## Fitting 
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


# ## Bibliography
# 
# 1. Andrew Nelson, Stuart Prescott, Isaac Gresham, & Andrew R. McCluskey. (2018, August 3). refnx/refnx: v0.0.17 (Version v0.0.17). Zenodo. http://doi.org/10.5281/zenodo.1345464

# In[ ]:




