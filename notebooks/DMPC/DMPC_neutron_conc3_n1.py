
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
from refnx.reflect import structure, ReflectModel, SLD#
from refnx.dataset import ReflectDataset
from refnx.analysis import Transform, CurveFitter, Objective, GlobalObjective, Parameter

# The custom class to constain the monolayer model. 
import sys
sys.path.insert(0, '/home/arm61/work/writing/articles/lipids_at_airdes/src/models')
import mol_vol as mv

data_dir = sys.argv[1] + '/data/processed/DMPC/'
figures_dir = sys.argv[1] + '/reports/figures/'
analysis_dir = sys.argv[1] + '/output/'


# In[2]:


# Reading datasets into refnx format
dataset3_n1 = ReflectDataset('{}DMPC_Neutron_conc2_dDMPC_hDES.mft'.format(data_dir))


# In[3]:


# Scattering length of the lipid head group 
# (found from summing the electrons in the head group 
# and multiplying by the classical radius of an electron)
head_sl = [602.7e-6, 602.7e-6]
# Scattering length of the lipid tail group",
tail_sl = [5329.76e-6, 5329.76e-6]
# Solvent SLD from ref [3]
solvent_sld = [0.43, 3.15]
# SLD of air",
super_sld = [0, 0]
# Some initial values for the head and tail thicknesses & APM
thick_heads = [13.1117, 11.0571]
tail_length = 1.54 + 1.265 * 13
chain_tilt = [0.792674, 0.79015]
vols = [200.497, 891.]
head_tail_rough = 3.3
tail_air_rough = 5.1


# In[4]:


# set up the chemical context system
dmpc3_n1 = mv.VolMono(head_sl[0], thick_heads[0], tail_sl[0], tail_length, chain_tilt[0], vols, 
                  head_tail_rough, tail_air_rough, reverse_monolayer=True, name='dmpc3_n1')


# In[5]:


# build the structures
air = SLD(0, '')
des_n1 = SLD(solvent_sld[0], '')

structure_dmpc3_n1 = air(0, 0) | dmpc3_n1 | des_n1(0, 0)


# In[6]:


def get_value(file):
    f = open(analysis_dir + 'dmpc/' + file + '.txt', 'r')
    for line in f:
        k = line
    l = k.split('$')[1].split('^')[0]
    return float(l)


# In[7]:


dmpc3_n1.head_mol_vol.setp(get_value('vh'), vary=False, bounds=(72., 472.))
dmpc3_n1.tail_mol_vol.setp(891., vary=False)
dmpc3_n1.tail_length.setp(vary=False)
dmpc3_n1.rough_head_tail.setp(get_value('roughh5'), vary=False, bounds=(0.1, 20))
dmpc3_n1.rough_preceding_mono.setp(get_value('rought5'), vary=True, bounds=(0.1,12))
dmpc3_n1.phit.setp(get_value('solt5')/100, vary=True, bounds=(0, 0.6))
dmpc3_n1.phih.setp(get_value('solh5')/100, vary=True, bounds=(0, 0.9999))
dmpc3_n1.solventsld.setp(vary=False)
dmpc3_n1.solventsldi.setp(vary=False)
dmpc3_n1.supersld.setp(vary=False)
dmpc3_n1.supersldi.setp(vary=False)
dmpc3_n1.thick_heads.setp(get_value('head5'), vary=False)
dmpc3_n1.cos_rad_chain_tilt.constraint = (dmpc3_n1.tail_mol_vol * (1 - dmpc3_n1.phih) * dmpc3_n1.thick_heads) / (dmpc3_n1.head_mol_vol * dmpc3_n1.tail_length * (1 - dmpc3_n1.phit))
structure_dmpc3_n1[-1].rough.setp(vary=False)
dmpc3_n1.solventsld.setp(solvent_sld[0], vary=False)


# In[8]:


# Creating a ReflectModel class object, add setting an initial scale 
model_dmpc3_n1 = ReflectModel(structure_dmpc3_n1)
model_dmpc3_n1.scale.setp(0.9364, vary=True, bounds=(0.005, 10))
# The background for held constant to a value determined from a previous fitting
model_dmpc3_n1.bkg.setp(dataset3_n1.y[-1], vary=False)


# In[9]:


# building the global objective
objective_n1 = Objective(model_dmpc3_n1, dataset3_n1, transform=Transform('YX4'))


# In[10]:


# A differential evolution algorithm is used to obtain an best fit
fitter = CurveFitter(objective_n1)
# A seed is used to ensure reproduciblity
res = fitter.fit('differential_evolution', seed=1)
# The first 200*200 samples are binned
fitter.sample(200, random_state=1)
fitter.sampler.reset()
# The collection is across 5000*200 samples
# The random_state seed is to allow for reproducibility
res = fitter.sample(1000, nthin=1, random_state=2, f='{}dmpc_highconc_chain_neutron_n1.txt'.format(analysis_dir))
flatchain = fitter.sampler.flatchain


# In[11]:


#print total objective
print(objective_n1)


# In[12]:


a = (dmpc3_n1.tail_mol_vol.value * (1 - flatchain[:, 3]) * dmpc3_n1.thick_heads.value) 
b = (dmpc3_n1.head_mol_vol.value * dmpc3_n1.tail_length.value * (1 - flatchain[:, 2]))
angle3 = a / b


# In[13]:


def printref(n, dataset, model, objective, analysis_dir, choose):
    file_open = open('{}dmpc{}_ref_neutron.txt'.format(analysis_dir, n), 'w')
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
    for pvec in choose:
        objective.setp(pvec)
        calc = model(dataset.x, x_err=dataset.x_err) * np.power(dataset.x, 4)
        for i in range(0, len(dataset.x)):
            file_open.write('{} '.format(calc[i]))
        file_open.write('\n')
    file_open.close()
    
def printsld(n, structure, objective, choose):
    file_open = open('{}dmpc{}_sld_neutron.txt'.format(analysis_dir, n), 'w')
    z, true_sld = structure.sld_profile()
    for i in range(0, len(z)):
        file_open.write('{} '.format(z[i]))
    file_open.write('\n')
    for i in range(0, len(z)):
        file_open.write('{} '.format(true_sld[i]))
    file_open.write('\n')
    for pvec in choose:
        objective.setp(pvec)
        zs, sld = structure.sld_profile()
        for i in range(0, len(z)):
            file_open.write('{} '.format(sld[i]))   
        file_open.write('\n')
    file_open.close()
    
choose = objective_n1.pgen(ngen=100)
printref("3_n1", dataset3_n1, model_dmpc3_n1, objective_n1, analysis_dir, choose)
printsld("3_n1", structure_dmpc3_n1, objective_n1, choose)


# In[15]:


lab = ['scale3', 'rought3', 'solt3', 'solh3']

for i in range(0, flatchain.shape[1]):
    total_pearsons = open('{}dmpc/{}_neutron_n1.txt'.format(analysis_dir, lab[i]), 'w')
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
    
lab2 = ['angle3']
kl = [angle3]
for i in range(0, len(lab2)):
    total_pearsons = open('{}dmpc/{}_neutron_n1.txt'.format(analysis_dir, lab2[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    c = np.rad2deg(np.arccos(a))
    k = [c[1], c[0] - c[1], c[1] - c[2]]
    q = '{:.2f}'.format(k[0])
    w = '{:.2f}'.format(k[1])
    e = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()
    
lab2 = ['tail3']
kl = angle3 * dmpc3_n1.tail_length.value
for i in range(0, len(lab2)):
    total_pearsons = open('{}dmpc/{}_neutron_n1.txt'.format(analysis_dir, lab2[i]), 'w')
    a = mquantiles(kl[i], prob=[0.025, 0.5, 0.975])
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()

