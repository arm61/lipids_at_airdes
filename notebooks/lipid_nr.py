#!/usr/bin/env python
# coding: utf-8

# # Neutron reflectometry analysis
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

# The refnx library, and associated classes
import refnx
from refnx.reflect import structure, ReflectModel, SLD#
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


lipid = sys.argv[1]
length = int(sys.argv[2])
sp = sys.argv[3]


# Here we assign the directories that contain the data, as well as where the figures and analysis outputs should be stored. If you directory structure does not match that in the GitHub repository these should be adapted. 

# In[ ]:


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
# The experimental datafiles are then read in.

# In[ ]:


# Reading datasets into refnx format
dataset_1 = ReflectDataset('{}nr_h_sp_{}.dat'.format(data_dir, sp))
dataset_2 = ReflectDataset('{}nr_hd_sp_{}.dat'.format(data_dir, sp))


# The scattering lengths for the head and tail components are defined.

# In[ ]:


head = {'C': 10, 'H': 18, 'O': 8, 'N': 1, 'P': 1}
tail = {'C': length * 2, 'D': length * 4 + 2}

head_sl = mv.get_scattering_length(head, neutron=True)
tail_sl = mv.get_scattering_length(tail, neutron=True)


# Initial 'guesses' for a series of parameters are defined.

# In[ ]:


solvent_sld = [0.43, 3.15]
super_sld = [0, 0]
thick_heads = 13.1117
tail_length = 1.54 + 1.265 * length
vols = [200.497, 891.]


# The `VolMono` class objects are defined. 

# In[ ]:


lipid_1 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                      reverse_monolayer=True, name='{}_1'.format(lipid))
lipid_2 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, vols, 
                      reverse_monolayer=True, name='{}_2'.format(lipid))


# A series of structures for each surface pressure is defined. 

# In[ ]:


# build the structures
air = SLD(0, '')
des_1 = SLD(solvent_sld[0], '')
des_2 = SLD(solvent_sld[1], '')

structure_lipid_1 = air(0, 0) | lipid_1 | des_1(0, 0)
structure_lipid_2 = air(0, 0) | lipid_2 | des_2(0, 0)


# A function is defined to allow floats to be obtained from previously defined output text files. 

# In[ ]:


def get_value(file):
    f = open(analysis_dir + lipid + '/' + file + '.txt', 'r')
    for line in f:
        k = line
    if '^' in k:
        l = k.split('$')[1].split('^')[0]
    else:
        l = k.split('$')[1].split('\\pm')[0]
    return float(l)


# The variables, and bounds for both surface pressures are setup.

# In[ ]:


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


# Each model is then associated with a dataset

# In[ ]:


model_lipid_1 = ReflectModel(structure_lipid_1)
model_lipid_1.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid_1.bkg.setp(dataset_1.y[-2], vary=False)

model_lipid_2 = ReflectModel(structure_lipid_2)
model_lipid_2.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid_2.bkg.setp(dataset_2.y[-2], vary=False)


# The global objective fitting object is defined and the fitting and MCMC performed. 

# In[ ]:


# building the global objective
objective_n1 = Objective(model_lipid_1, dataset_1, transform=Transform('YX4'))
objective_n2 = Objective(model_lipid_2, dataset_2, transform=Transform('YX4'))

global_objective = GlobalObjective([objective_n1, objective_n2])


# ## Fitting 
# 
# The differential evolution algorithm is used to find optimal parameters, before the MCMC algorithm probes the parameter space for 1000 steps.

# In[ ]:


# A differential evolution algorithm is used to obtain an best fit
fitter = CurveFitter(global_objective)
# A seed is used to ensure reproduciblity
res = fitter.fit('differential_evolution', seed=1)
# The first 200*200 samples are binned
fitter.sample(200, random_state=1)
fitter.sampler.reset()
# The collection is across 5000*200 samples
# The random_state seed is to allow for reproducibility
res = fitter.sample(1000, nthin=1, random_state=1, f='{}{}/{}_chain_neutron.txt'.format(analysis_dir, lipid, sp))
flatchain = fitter.sampler.flatchain


# The `global_objective` is printed containing information about the models.

# In[ ]:


#print total objective
print(global_objective)


# ## Bibliography
# 
# 1. Andrew Nelson, Stuart Prescott, Isaac Gresham, & Andrew R. McCluskey. (2018, August 3). refnx/refnx: v0.0.17 (Version v0.0.17). Zenodo. http://doi.org/10.5281/zenodo.1345464

# In[ ]:




