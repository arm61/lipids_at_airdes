
# coding: utf-8

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

lipid = sys.argv[1]
length = int(sys.argv[2])
sp = sys.argv[3]

data_dir = '../data/processed/{}/'.format(lipid)
figures_dir = '../reports/figures/'
analysis_dir = '../output/'


# In[ ]:


refnx.version.full_version, scipy.version.version


# In[ ]:


# Reading datasets into refnx format
dataset_1 = ReflectDataset('{}nr_h_sp_{}.dat'.format(data_dir, sp))
dataset_2 = ReflectDataset('{}nr_hd_sp_{}.dat'.format(data_dir, sp))


# In[ ]:


head = {'C': 10, 'H': 18, 'O': 8, 'N': 1, 'P': 1}
tail = {'C': length * 2, 'D': length * 4 + 2}

head_sl = mv.get_scattering_length(head, neutron=True)
tail_sl = mv.get_scattering_length(tail, neutron=True)

solvent_sld = [0.43, 3.15]
super_sld = [0, 0]
thick_heads = 13.1117
tail_length = 1.54 + 1.265 * length
chain_tilt = 0.792674
vols = [200.497, 891.]


# In[ ]:


lipid_1 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                      reverse_monolayer=True, name='{}_1'.format(lipid))
lipid_2 = mv.VolMono(head_sl, thick_heads, tail_sl, tail_length, chain_tilt, vols, 
                      reverse_monolayer=True, name='{}_2'.format(lipid))


# In[ ]:


# build the structures
air = SLD(0, '')
des_1 = SLD(solvent_sld[0], '')
des_2 = SLD(solvent_sld[1], '')

structure_lipid_1 = air(0, 0) | lipid_1 | des_1(0, 0)
structure_lipid_2 = air(0, 0) | lipid_2 | des_2(0, 0)


# In[ ]:


def get_value(file):
    f = open(analysis_dir + lipid + '/' + file + '.txt', 'r')
    for line in f:
        k = line
    l = k.split('$')[1].split('^')[0]
    return float(l)


# In[ ]:


lipid_1.head_mol_vol.setp(get_value('vh'), vary=False)
lipid_1.tail_mol_vol.setp(get_value('vt'), vary=False)
lipid_1.tail_length.setp(vary=False)
lipid_1.rough_head_tail.constraint = structure_lipid_1[-1].rough
lipid_1.rough_preceding_mono.constraint = structure_lipid_1[-1].rough
lipid_1.phih.constraint = 1 - (lipid_1.head_mol_vol * lipid_1.tail_length * lipid_1.cos_rad_chain_tilt / 
                               (lipid_1.tail_mol_vol * lipid_1.thick_heads))
lipid_1.thick_heads.setp(get_value('head'), vary=False)
lipid_1.cos_rad_chain_tilt.setp(np.cos(np.deg2rad(get_value('angle{}'.format(sp)))), vary=True, bounds=(0.001, 0.909))
structure_lipid_1[-1].rough.setp(get_value('rough{}'.format(sp)), vary=True, bounds=(2.5, 6))

lipid_2.head_mol_vol.constraint = lipid_1.head_mol_vol
lipid_2.tail_mol_vol.constraint = lipid_1.tail_mol_vol
lipid_2.tail_length.constraint = lipid_1.tail_length
lipid_2.rough_head_tail.constraint = structure_lipid_1[-1].rough
lipid_2.rough_preceding_mono.constraint = structure_lipid_1[-1].rough
lipid_2.phih.constraint = lipid_1.phih
lipid_2.thick_heads.constraint = lipid_1.thick_heads
lipid_2.cos_rad_chain_tilt.constraint = lipid_1.cos_rad_chain_tilt
structure_lipid_2[-1].rough.constraint = structure_lipid_1[-1].rough


# In[ ]:


model_lipid_1 = ReflectModel(structure_lipid_1)
model_lipid_1.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid_1.bkg.setp(dataset_1.y[-2], vary=False)

model_lipid_2 = ReflectModel(structure_lipid_2)
model_lipid_2.scale.setp(vary=True, bounds=(0.005, 10))
model_lipid_2.bkg.setp(dataset_2.y[-2], vary=False)


# In[ ]:


# building the global objective
objective_n1 = Objective(model_lipid_1, dataset_1, transform=Transform('YX4'))
objective_n2 = Objective(model_lipid_2, dataset_2, transform=Transform('YX4'))

global_objective = GlobalObjective([objective_n1, objective_n2])


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
res = fitter.sample(1000, nthin=1, random_state=1, f='{}/{}/{}_chain_neutron.txt'.format(analysis_dir, lipid, sp))
flatchain = fitter.sampler.flatchain


# In[ ]:


#print total objective
print(global_objective)


# In[ ]:


def printref(n, dataset, model, objective, analysis_dir, choose):
    file_open = open('{}{}{}_ref_neutron.txt'.format(analysis_dir, lipid, n), 'w')
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
    file_open = open('{}{}{}_sld_neutron.txt'.format(analysis_dir, lipid, n), 'w')
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
    
choose = global_objective.pgen(ngen=100)
printref("{}_1".format(sp), dataset_1, model_lipid_1, global_objective, analysis_dir, choose)
printsld("{}_1".format(sp), structure_lipid_1, global_objective, choose)
printref("{}_2".format(sp), dataset_2, model_lipid_2, global_objective, analysis_dir, choose)
printsld("{}_2".format(sp), structure_lipid_2, global_objective, choose)


# In[ ]:


lab = ['scale{}'.format(sp), 'angle{}'.format(sp), 'rough{}'.format(sp), 'scalea{}'.format(sp)]

for i in range(0, flatchain.shape[1]):
    total_pearsons = open('{}{}/{}_neutron.txt'.format(analysis_dir, lipid, lab[i]), 'w')
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
    
lab2 = ['solh{}'.format(sp)]
kl = 1 - ((lipid_1.head_mol_vol.value * flatchain[:, 1] * 
           lipid_1.tail_length.value) / (lipid_1.tail_mol_vol.value * lipid_1.thick_heads.value))
kl = kl * 100
for i in range(0, len(lab2)):
    total_pearsons = open('{}{}/{}_neutron.txt'.format(analysis_dir, lipid, lab2[i]), 'w')
    a = mquantiles(kl, prob=[0.025, 0.5, 0.975])
    c = a
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    w = '{:.2f}'.format(k[1])
    e = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()
    
lab2 = ['tail{}'.format(sp)]
kl = flatchain[:, 1] * lipid_1.tail_length.value
for i in range(0, len(lab2)):
    total_pearsons = open('{}{}/{}_neutron.txt'.format(analysis_dir, lipid, lab2[i]), 'w')
    a = mquantiles(kl, prob=[0.025, 0.5, 0.975])
    k = [a[1], a[1] - a[0], a[2] - a[1]]
    q = '{:.2f}'.format(k[0])
    e = '{:.2f}'.format(k[1])
    w = '{:.2f}'.format(k[2])
    total_pearsons.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
    total_pearsons.close()

