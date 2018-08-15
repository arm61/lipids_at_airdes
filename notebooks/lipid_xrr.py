
# coding: utf-8

# Import the modules used within the analysis.

# In[2]:


# Standard libraries to import
from __future__ import division
import numpy as np 
import scipy
from scipy.stats.mstats import mquantiles

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

# The type of lipid being investigated
lipid = sys.argv[1]
length = int(sys.argv[2])
sp1 = sys.argv[3]
sp2 = sys.argv[4]
sp3 = sys.argv[5]
sp4 = sys.argv[6]

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


# The experimental datafiles are then read in and redefined such that all data after q = 0.6 Å<sup>-1</sup> is ignored as this is considered background. 

# In[ ]:


# Reading datasets into refnx format
dataset3 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp1)), 0.6)
dataset4 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp2)), 0.6)
dataset5 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp3)), 0.6)
dataset6 = helper.data_cutoff(ReflectDataset('{}xrr_sp_{}.dat'.format(data_dir, sp4)), 0.6)


# Defining the scattering lengths for the lipid head and tails

# In[ ]:


if lipid == 'dmpg':
    head = {'C': 8, 'H': 12, 'O': 10, 'Na': 1, 'P': 1}
else:
    head = {'C': 10, 'H': 18, 'O': 8, 'N': 1, 'P': 1}
tail = {'C': length * 2, 'H': length * 4 + 2}

head_sl = mv.get_scattering_length(head)
tail_sl = mv.get_scattering_length(tail)


# Defining initial values for the system variables

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


# Defining the length of the carbon tail ($t_l$) based on the Tanford equation.
# 
# $$ t_l = 1.54 + 1.265(n-1) $$
# 
# where $n$ is the length of the carbon chain (e.g. 16 for DPPC).

# In[ ]:


tail_length = 1.54 + 1.265 * length


# Initialise the lipid objects

# In[ ]:


lipid3 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                   reverse_monolayer=True, name='{}1'.format(lipid))
lipid4 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                   reverse_monolayer=True, name='{}2'.format(lipid))
lipid5 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                   reverse_monolayer=True, name='{}3'.format(lipid))
lipid6 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                   reverse_monolayer=True, name='{}4'.format(lipid))


# Building a structure of each of the lipid systems

# In[ ]:


air = SLD(0, 'air')
des = SLD(10.8, 'des')

structure_lipid3 = air(0, 0) | lipid3 | des(0, 3.3)
structure_lipid4 = air(0, 0) | lipid4 | des(0, 3.3)
structure_lipid5 = air(0, 0) | lipid5 | des(0, 3.3)
structure_lipid6 = air(0, 0) | lipid6 | des(0, 3.3)


# The variables, and bounds are then defined.

# In[ ]:


lipid3.head_mol_vol.setp(vary=True, bounds=(vols[0]*0.8, vols[0]*1.2))
lipid3.tail_mol_vol.setp(vary=True, bounds=(vols[1]*0.8, vols[1]*1.2))
lipid3.tail_length.setp(vary=False)
lipid3.cos_rad_chain_tilt.setp(vary=True, bounds=(0.01, 0.99))
lipid3.rough_head_tail.constraint = structure_lipid3[-1].rough
lipid3.rough_preceding_mono.constraint = structure_lipid3[-1].rough
lipid3.phih.constraint = 1 - (lipid3.head_mol_vol / 
                              lipid3.tail_mol_vol) * (lipid3.cos_rad_chain_tilt * 
                                                      lipid3.tail_length / lipid3.thick_heads)
lipid3.thick_heads.setp(vary=True, bounds=(6, 15))
structure_lipid3[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid4.cos_rad_chain_tilt.setp(vary=True, bounds=(0.01, 0.99))
lipid4.rough_head_tail.constraint = structure_lipid4[-1].rough
lipid4.rough_preceding_mono.constraint = structure_lipid4[-1].rough
lipid4.phih.constraint = 1 - (lipid4.head_mol_vol / lipid4.tail_mol_vol) * (lipid4.cos_rad_chain_tilt * lipid4.tail_length / lipid4.thick_heads)
structure_lipid4[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid5.cos_rad_chain_tilt.setp(0.57, vary=True, bounds=(0.01, 0.99))
lipid5.rough_head_tail.constraint = structure_lipid5[-1].rough
lipid5.rough_preceding_mono.constraint = structure_lipid5[-1].rough
lipid5.phih.constraint = 1 - (lipid5.head_mol_vol / lipid5.tail_mol_vol) * (lipid5.cos_rad_chain_tilt * lipid5.tail_length / lipid5.thick_heads)
structure_lipid5[-1].rough.setp(vary=True, bounds=(2.5, 6))


# In[ ]:


lipid6.cos_rad_chain_tilt.setp(0.57, vary=True, bounds=(0.01, 0.99))
lipid6.rough_head_tail.constraint = structure_lipid6[-1].rough
lipid6.rough_preceding_mono.constraint = structure_lipid6[-1].rough
lipid6.phih.constraint = 1 - (lipid6.head_mol_vol / lipid6.tail_mol_vol) * (lipid6.cos_rad_chain_tilt * lipid6.tail_length / lipid6.thick_heads)
structure_lipid6[-1].rough.setp(vary=True, bounds=(2.5, 6))


# The three lipids are constrained such that the tail and head volumes and head thickness is kept constant.

# In[ ]:


lipids = [lipid3, lipid4, lipid5, lipid6]

lipids = mv.set_constraints(lipids)


# Each model is then associated with a dataset

# In[ ]:


model_lipid3 = ReflectModel(structure_lipid3)
model_lipid3.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid3.bkg.setp(dataset3.y[-1], vary=False)

model_lipid4 = ReflectModel(structure_lipid4)
model_lipid4.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid4.bkg.setp(dataset4.y[-1], vary=False)

model_lipid5 = ReflectModel(structure_lipid5)
model_lipid5.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid5.bkg.setp(dataset5.y[-1], vary=False)

model_lipid6 = ReflectModel(structure_lipid6)
model_lipid6.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid6.bkg.setp(dataset6.y[-1], vary=False)


# The global objective fitting object is defined and the fitting and MCMC performed. 

# In[ ]:


objective3 = Objective(model_lipid3, dataset3, transform=Transform('YX4'))
objective4 = Objective(model_lipid4, dataset4, transform=Transform('YX4'))
objective5 = Objective(model_lipid5, dataset5, transform=Transform('YX4'))
objective6 = Objective(model_lipid6, dataset6, transform=Transform('YX4'))

global_objective = GlobalObjective([objective3, objective4, objective5, objective6])


# In[ ]:


fitter = CurveFitter(global_objective)
res = fitter.fit('differential_evolution', seed=1)

fitter.sample(200, random_state=1)
fitter.sampler.reset()
res = fitter.sample(1000, nthin=1, random_state=1, f='{}{}_chain.txt'.format(analysis_dir, lipid))
flatchain = fitter.sampler.flatchain


# In[ ]:


print(global_objective)


# The reflectometry and scattering length density profiles are printed to a file.

# In[ ]:


def printref(n, dataset, model, objective, analysis_dir):
    file_open = open('{}{}{}_ref.txt'.format(analysis_dir, lipid, n), 'w')
    saved_params = np.array(objective.parameters)
    for i in range(0, len(dataset.x)):
        file_open.write('{} '.format(dataset.x[i]))
    file_open.write('\n')
    for i in range(0, len(dataset.x)):
        file_open.write('{} '.format(dataset.y[i]*(dataset.x[i])**4))
    file_open.write('\n')
    for i in range(0, len(dataset.x)):
        file_open.write('{} '.format(dataset.y_err[i]*(dataset.x[i])**4))
    file_open.write('\n')
    for i in range(0, len(dataset.x)):
        file_open.write('{} '.format((model(dataset.x, x_err=dataset.x_err)[i])*(dataset.x[i])**4))
    file_open.write('\n')
    choose = objective.pgen(ngen=100)
    for pvec in choose:
        objective.setp(pvec)
        calc = model(dataset.x, x_err=dataset.x_err) * np.power(dataset.x, 4)
        for i in range(0, len(dataset.x)):
            file_open.write('{} '.format(calc[i]))
        file_open.write('\n')
    file_open.close()
    
def printsld(n, structure, objective):
    file_open = open('{}{}{}_sld.txt'.format(analysis_dir, lipid, n), 'w')
    z, true_sld = structure.sld_profile()
    for i in range(0, len(z)):
        file_open.write('{} '.format(z[i]))
    file_open.write('\n')
    for i in range(0, len(z)):
        file_open.write('{} '.format(true_sld[i]))
    file_open.write('\n')
    choose = objective.pgen(ngen=100)
    for pvec in choose:
        objective.setp(pvec)
        zs, sld = structure.sld_profile()
        for i in range(0, len(z)):
            file_open.write('{} '.format(sld[i]))   
        file_open.write('\n')
    file_open.close()
    
    
printref(sp1, dataset3, model_lipid3, global_objective, analysis_dir)
printref(sp2, dataset4, model_lipid4, global_objective, analysis_dir)
printref(sp3, dataset5, model_lipid5, global_objective, analysis_dir)
printref(sp4, dataset6, model_lipid6, global_objective, analysis_dir)
printsld(sp1, structure_lipid3, global_objective)
printsld(sp2, structure_lipid4, global_objective)
printsld(sp3, structure_lipid5, global_objective)
printsld(sp4, structure_lipid6, global_objective)


# The tail length and head solvation are determined and the variables are printed to files to allow of inclusion in the paper. 

# In[ ]:


tail3 = flatchain[:, 2] * lipid3.tail_length.value
tail4 = flatchain[:, 7] * lipid4.tail_length.value
tail5 = flatchain[:, 10] * lipid5.tail_length.value
tail6 = flatchain[:, 13] * lipid6.tail_length.value

solh3 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail3 / flatchain[:, 1])
solh4 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail4 / flatchain[:, 1])
solh5 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail5 / flatchain[:, 1])
solh6 = 1 - (flatchain[:, 4] / flatchain[:, 3]) * (tail6 / flatchain[:, 1])


# In[ ]:


lab = ['scale{}'.format(sp1), 'head{}', 'angle{}'.format(sp1), 'vt', 'vh', 'rough{}'.format(sp1), 
       'scale{}'.format(sp2), 'angle{}'.format(sp2), 'rough{}'.format(sp2), 
       'scale{}'.format(sp3), 'angle{}'.format(sp3), 'rough{}'.format(sp3), 
       'scale{}'.format(sp4), 'angle{}'.format(sp4), 'rough{}'.format(sp4)]
for i in range(0, flatchain.shape[1]):
    total_pearsons = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab[i]), 'w')
    a = mquantiles(flatchain[:, i], prob=[0.025, 0.5, 0.975])
    if 'angle' in lab[i]:
        c = np.rad2deg(np.arccos(a))
        k = [c[1], c[0] - c[1], c[1] - c[2]]
        q = '{:.2f}'.format(k[0])
        w = '{:.2f}'.format(k[1])
        e = '{:.2f}'.format(k[2])
        total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    elif 'sol' in lab[i]:
        k = [a[1]*100, (a[1] - a[0])*100, (a[2] - a[1])*100]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        w = '{:.2f}'.format(k[2])
        total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    else:
        k = [a[1], a[1] - a[0], a[2] - a[1]]
        q = '{:.2f}'.format(k[0])
        e = '{:.2f}'.format(k[1])
        w = '{:.2f}'.format(k[2])
        total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()
            
lab2 = ['tail{}'.format(sp1), 'tail{}'.format(sp2), 'tail{}'.format(sp3), 'tail{}'.format(sp4)]
kl = [tail3, tail4, tail5, tail6]
for i in range(0, len(lab2)):
    total_pearsons = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab2[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()
    
lab2 = ['solh{}'.format(sp1), 'solh{}'.format(sp2), 'solh{}'.format(sp3), 'solh{}'.format(sp4)]
kl = [solh3, solh4, solh5, solh6]
for i in range(0, len(lab2)):
    total_pearsons = open('{}{}/{}.txt'.format(analysis_dir, lipid, lab2[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    k = [a[1]*100, (a[1] - a[0])*100, (a[2] - a[1])*100]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()

