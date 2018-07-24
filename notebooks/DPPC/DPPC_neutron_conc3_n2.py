
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

data_dir = sys.argv[1] + '/data/processed/DPPC/'
figures_dir = sys.argv[1] + '/reports/figures/'
analysis_dir = sys.argv[1] + '/output/'


# In[2]:


# Reading datasets into refnx format
dataset3_n2 = ReflectDataset('{}DPPC_Neutron_conc3_dDPPC_hdDES.mft'.format(data_dir))


# In[3]:


# Scattering length of the lipid head group 
# (found from summing the electrons in the head group 
# and multiplying by the classical radius of an electron)
head_sl = [602.7e-6, 602.7e-6]
# Scattering length of the lipid tail group",
tail_sl = [6129.2e-6, 6129.2e-6]
# Solvent SLD from ref [3]
solvent_sld = [0.43, 3.15]
# SLD of air",
super_sld = [0, 0]
# Some initial values for the head and tail thicknesses & APM
thick_heads = [13.1117, 11.0571]
tail_length = 1.54 + 1.265 * 15
chain_tilt = [0.792674, 0.79015]
vols = [200.497, 891.]
head_tail_rough = 3.3
tail_air_rough = 5.1


# In[4]:


# set up the chemical context system
dppc3_n2 = mv.VolMono(head_sl[1], thick_heads[1], tail_sl[1], tail_length, chain_tilt[1], vols, 
                  head_tail_rough, tail_air_rough, reverse_monolayer=True, name='dppc3_n2')


# In[5]:


# build the structures
air = SLD(0, '')
des_n2 = SLD(solvent_sld[1], '')

structure_dppc3_n2 = air(0, 0) | dppc3_n2 | des_n2(0, 0)


# In[6]:


def get_value(file):
    f = open(analysis_dir + 'dppc/' + file + '.txt', 'r')
    for line in f:
        k = line
    l = k.split('$')[1].split('^')[0]
    return float(l)


# In[7]:


dppc3_n2.head_mol_vol.setp(get_value('vh'), vary=False, bounds=(72., 472.))
dppc3_n2.tail_mol_vol.setp(get_value('vt'), vary=False)
dppc3_n2.tail_length.setp(vary=False)
dppc3_n2.rough_head_tail.constraint = dppc3_n2.solventrough
dppc3_n2.rough_preceding_mono.constraint = dppc3_n2.solventrough
dppc3_n2.solventrough.setp(get_value('rough4'), vary=True, bounds=(2.5, 8.))
dppc3_n2.phih.constraint = 1 - (dppc3_n2.head_mol_vol * dppc3_n2.tail_length * dppc3_n2.cos_rad_chain_tilt / (dppc3_n2.tail_mol_vol * dppc3_n2.thick_heads))
dppc3_n2.solventsld.setp(vary=False)
dppc3_n2.solventsldi.setp(vary=False)
dppc3_n2.supersld.setp(vary=False)
dppc3_n2.supersldi.setp(vary=False)
dppc3_n2.thick_heads.setp(get_value('head4'), vary=False)
dppc3_n2.phit.setp(0, vary=False)
dppc3_n2.cos_rad_chain_tilt.setp(np.cos(np.deg2rad(get_value('angle5'))), vary=True, bounds=(0.01, 0.99))
structure_dppc3_n2[-1].rough.setp(vary=False)
dppc3_n2.solventsld.setp(solvent_sld[1], vary=False)


# In[8]:


# Creating a ReflectModel class object, add setting an initial scale 
model_dppc3_n2 = ReflectModel(structure_dppc3_n2)
model_dppc3_n2.scale.setp(0.9364, vary=True, bounds=(0.005, 10))
# The background for held constant to a value determined from a previous fitting
model_dppc3_n2.bkg.setp(dataset3_n2.y[-2], vary=False)


# In[9]:


# building the global objective
objective_n2 = Objective(model_dppc3_n2, dataset3_n2, transform=Transform('YX4'))


# In[10]:


# A differential evolution algorithm is used to obtain an best fit
fitter = CurveFitter(objective_n2)
# A seed is used to ensure reproduciblity
res = fitter.fit('differential_evolution', seed=1)
# The first 200*200 samples are binned
fitter.sample(200, random_state=1)
fitter.sampler.reset()
# The collection is across 5000*200 samples
# The random_state seed is to allow for reproducibility
res = fitter.sample(1000, nthin=1, random_state=2, f='{}dppc_highconc_chain_neutron_n2.txt'.format(analysis_dir))
flatchain = fitter.sampler.flatchain


# In[11]:


#print total objective
print(objective_n2)


# In[12]:


angle3 = flatchain[:, 1]#(dppc3_n2.head_mol_vol.value * dppc3_n2.tail_length.value * flatchain[:, 1] * 
        #(1 - flatchain[:, 3])) / (dppc3_n2.tail_mol_vol.value * (1 - flatchain[:, 4]))


# In[13]:


lab = ['scale3', 'angle3', 'rought3', 'solh3']

def printref(n, dataset, model, objective, analysis_dir, choose):
    file_open = open('{}dppc{}_ref_neutron.txt'.format(analysis_dir, n), 'w')
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
    file_open = open('{}dppc{}_sld_neutron.txt'.format(analysis_dir, n), 'w')
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
    
choose = objective_n2.pgen(ngen=100)
printref("3_n2", dataset3_n2, model_dppc3_n2, objective_n2, analysis_dir, choose)
printsld("3_n2", structure_dppc3_n2, objective_n2, choose)


# In[14]:


lab = ['scale3', 'angle3', 'rought3']

for i in range(0, flatchain.shape[1]):
    total_pearsons = open('{}dppc/{}_neutron_n2.txt'.format(analysis_dir, lab[i]), 'w')
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
    
lab2 = ['solh3']
kl = 1 - ((dppc3_n2.head_mol_vol.value * flatchain[:, 1] * dppc3_n2.tail_length.value) / (dppc3_n2.tail_mol_vol.value * dppc3_n2.thick_heads.value))
kl = kl * 100
for i in range(0, len(lab2)):
    total_pearsons = open('{}dppc/{}_neutron_n2.txt'.format(analysis_dir, lab2[i]), 'w')
    a = mquantiles(kl, prob=[0.025, 0.5, 0.975])
    c = a
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    w = '{:.2f}'.format(k[1])
    e = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()
    
lab2 = ['tail3']
kl = flatchain[:, 1] * dppc3_n2.tail_length.value
for i in range(0, len(lab2)):
    total_pearsons = open('{}dppc/{}_neutron_n2.txt'.format(analysis_dir, lab2[i]), 'w')
    a = mquantiles(kl, prob=[0.025, 0.5, 0.975])
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()

