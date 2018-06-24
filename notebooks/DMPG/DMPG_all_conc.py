
# coding: utf-8

# In[1]:


# Standard libraries to import
from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
#from matplotlib import rcParams, rc
import seaborn as sns; sns.set('paper', palette='colorblind')
import matplotlib as mpl
from matplotlib import gridspec
from scipy.stats import pearsonr
from scipy.stats.mstats import mquantiles

# The refnx library, and associated classes
import refnx
from refnx.reflect import structure, ReflectModel, SLD
from refnx.dataset import ReflectDataset
from refnx.analysis import Transform, CurveFitter, Objective, GlobalObjective, Parameter

# The custom class to constain the monolayer model. 
import sys
sys.path.insert(0, '/home/arm61/work/writing/articles/lipids_at_airdes/src/models')
import mol_vol as mv

data_dir = '/home/arm61/work/writing/articles/lipids_at_airdes/data/processed/DMPG/'
figures_dir = '/home/arm61/work/writing/articles/lipids_at_airdes/reports/figures/'
analysis_dir = '/home/arm61/work/writing/articles/lipids_at_airdes/output/'


# In[2]:


# Reading datasets into refnx format
dataset4 = ReflectDataset('{}DMPG_Xray_conc4.dat'.format(data_dir))
dataset5 = ReflectDataset('{}DMPG_Xray_conc5.dat'.format(data_dir))


# In[3]:


# Scattering length of the lipid head group 
# (found from summing the electrons in the head group 
# and multiplying by the classical radius of an electron)
head_sl = 4731e-6
# Scattering length of the lipid tail group 
tail_sl = 5985e-6
# Some initial values for the head and tail thicknesses & APM
thick_heads = [13.1117, 11.0571]
tail_length = 1.54 + 1.265 * 13
chain_tilt = [0.792674, 0.79015]
vols = [200.497, 891.]
head_tail_rough = 3.3
tail_air_rough = 5.1


# In[4]:


# set up the chemical context system
dmpg4 = mv.VolMono(head_sl, thick_heads[0], tail_sl, tail_length, chain_tilt[0], vols, 
                  head_tail_rough, tail_air_rough, reverse_monolayer=True, name='dmpg4')
dmpg5 = mv.VolMono(head_sl, thick_heads[1], tail_sl, tail_length, chain_tilt[1], vols, 
                  head_tail_rough, tail_air_rough, reverse_monolayer=True, name='dmpg5')


# In[5]:


# build the structures
air = SLD(0, '')
des = SLD(10.8, '')

structure_dmpg4 = air(0, 0) | dmpg4 | des(0, 0)
structure_dmpg5 = air(0, 0) | dmpg5 | des(0, 0)


# In[6]:


dmpg4.head_mol_vol.setp(vary=True, bounds=(72., 472.))
dmpg4.tail_mol_vol.setp(779., vary=False)
dmpg4.tail_length.setp(vary=False)
dmpg4.cos_rad_chain_tilt.setp(vary=True, bounds=(0.56, 0.96))
dmpg4.rough_head_tail.setp(vary=True, bounds=(1, 12))
dmpg4.rough_preceding_mono.setp(vary=True, bounds=(1,12))
dmpg4.phit.setp(0, vary=True, bounds=(0, 0.2))
dmpg4.phih.setp(0.312487, vary=True, bounds=(0, 0.9))
dmpg4.solventsld.setp(vary=False)
dmpg4.solventsldi.setp(vary=False)
dmpg4.supersld.setp(vary=False)
dmpg4.supersldi.setp(vary=False)
dmpg4.thick_heads.constraint = (dmpg4.head_mol_vol * dmpg4.tail_length * dmpg4.cos_rad_chain_tilt * 
                                (1 - dmpg4.phit)) / (dmpg4.tail_mol_vol * (1 - dmpg4.phih))
structure_dmpg4[-1].rough.setp(vary=False)

dmpg5.tail_mol_vol.setp(779., vary=False)
dmpg5.tail_length.setp(vary=False)
dmpg5.cos_rad_chain_tilt.setp(vary=True, bounds=(0.57, 0.97))
dmpg5.rough_head_tail.setp(vary=True, bounds=(1, 12))
dmpg5.rough_preceding_mono.setp(vary=True, bounds=(1,12))
dmpg5.phit.setp(0.001, vary=True, bounds=(0, 0.1))
dmpg5.phih.setp(0.291706, vary=True, bounds=(0, 0.9))
dmpg5.solventsld.setp(vary=False)
dmpg5.solventsldi.setp(vary=False)
dmpg5.supersld.setp(vary=False)
dmpg5.supersldi.setp(vary=False)
dmpg5.thick_heads.constraint = (dmpg5.head_mol_vol * dmpg5.tail_length * dmpg5.cos_rad_chain_tilt * 
                                (1 - dmpg5.phit)) / (dmpg5.tail_mol_vol * (1 - dmpg5.phih))
structure_dmpg5[-1].rough.setp(vary=False)


# In[7]:


# constraining the head and tail molecular volumes
lipids = [dmpg4, dmpg5]

lipids = mv.set_contraints(lipids)


# In[8]:


# Creating a ReflectModel class object, add setting an initial scale 
model_dmpg4 = ReflectModel(structure_dmpg4)
model_dmpg4.scale.setp(0.9364, vary=True, bounds=(0.005, 10))
# The background for held constant to a value determined from a previous fitting
model_dmpg4.bkg.setp(dataset4.y[-1], vary=False)

# Creating a ReflectModel class object, add setting an initial scale 
model_dmpg5 = ReflectModel(structure_dmpg5)
model_dmpg5.scale.setp(0.9364, vary=True, bounds=(0.005, 10))
# The background for held constant to a value determined from a previous fitting
model_dmpg5.bkg.setp(dataset5.y[-4], vary=False)


# In[9]:


# building the global objective
objective4 = Objective(model_dmpg4, dataset4, transform=Transform('YX4'))
objective5 = Objective(model_dmpg5, dataset5, transform=Transform('YX4'))

global_objective = GlobalObjective([objective4, objective5])


# In[10]:


# A differential evolution algorithm is used to obtain an best fit
fitter = CurveFitter(global_objective)
# A seed is used to ensure reproduciblity
res = fitter.fit('differential_evolution', seed=1)
# The first 200*200 samples are binned
fitter.sample(200, random_state=1)
fitter.sampler.reset()
# The collection is across 5000*200 samples
# The random_state seed is to allow for reproducibility
res = fitter.sample(1000, nthin=1, random_state=2, f='{}dmpg_highconc_chain.txt'.format(analysis_dir))
flatchain = fitter.sampler.flatchain


# In[11]:


#print total objective
print(global_objective)


# In[12]:


head4 = flatchain[:, 2] * dmpg4.tail_length.value * flatchain[:, 1] * (1 - flatchain[:, 5])
head4 = head4 / (dmpg4.tail_mol_vol.value)
a = 1 - flatchain[:, 6]
head4 = np.array(head4) / a

head5 = flatchain[:, 2] * dmpg5.tail_length.value * flatchain[:, 8] * (1 - flatchain[:, 11])
head5 = head5 / (dmpg5.tail_mol_vol.value)
a = 1 - flatchain[:, 12]
head5 = np.array(head5) / a

tail4 = flatchain[:, 1] * dmpg4.tail_length.value
tail5 = flatchain[:, 8] * dmpg5.tail_length.value


# In[13]:


def printref(n, dataset, model, objective, analysis_dir):
    file_open = open('{}dmpg{}_ref.txt'.format(analysis_dir, n), 'w')
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
    file_open = open('{}dmpg{}_sld.txt'.format(analysis_dir, n), 'w')
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
    
printref(4, dataset4, model_dmpg4, global_objective, analysis_dir)
printref(5, dataset5, model_dmpg5, global_objective, analysis_dir)
printsld(4, structure_dmpg4, global_objective)
printsld(5, structure_dmpg5, global_objective)


# In[14]:


# determine the average pearson r coefficient magnitude 
from scipy.stats import pearsonr

sumr = np.array([])

for i in range(0, flatchain.shape[1]-1):
    for j in range(i+1, flatchain.shape[1]):
        sumr = np.append(sumr, np.abs(pearsonr(flatchain[:, i], flatchain[:, j])[0]))
        
total_pearsons = open('{}dmpg_all_pearsons.txt'.format(analysis_dir), 'w')
total_pearsons.write('pearsons {}'.format(np.average(sumr)))
total_pearsons.close()


# In[16]:


lab = ['scale4', 'angle4', 'vh', 'roughh4', 'rought4', 'solt4', 'solh4', 
       'scale5', 'angle5', 'roughh5', 'rought5', 'solt5', 'solh5']
for i in range(0, flatchain.shape[1]):
    total_pearsons = open('{}dmpg/{}.txt'.format(analysis_dir, lab[i]), 'w')
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
    
lab2 = ['head4', 'head5']
kl = [head4, head5]
for i in range(0, len(lab2)):
    total_pearsons = open('{}dmpg/{}.txt'.format(analysis_dir, lab2[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()
        
lab2 = ['tail4', 'tail5']
kl = [tail4, tail5]
for i in range(0, len(lab2)):
    total_pearsons = open('{}dmpg/{}.txt'.format(analysis_dir, lab2[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()
    
total_pearsons = open('{}dmpg/{}.txt'.format(analysis_dir, 'vh'), 'w')
a = mquantiles(flatchain[:, 2], prob=[0.025, 0.5, 0.975])
k = [a[1], a[1] - a[0], a[2] - a[1]]
q = '{:.2f}'.format(k[0])
e = '{:.2f}'.format(k[1])
w = '{:.2f}'.format(k[2])
total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
total_pearsons.close()

