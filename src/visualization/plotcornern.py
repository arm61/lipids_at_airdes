
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns; sns.set('paper', palette='colorblind')
import numpy as np
from refnx import analysis
import sys


# In[3]:


f = sys.argv[1] #'../../output/dlpc_highconc_chain.txt'
tail_length = float(sys.argv[2]) #1.54 + 1.265 * 12
vole = float(sys.argv[3]) #667.)
flatchain = analysis.curvefitter.load_chain(f)
figures_dir = sys.argv[6] + '/reports/figures/'
output_dir = sys.argv[6] + '/output/'
lab = sys.argv[4] #'dlpc'
num = sys.argv[5]
f = open(output_dir + sys.argv[4] + '/vh.txt', 'r')
for line in f:
    k = line
l = k.split('$')[1].split('^')[0]
head_vol = float(l)
def get_value(file):
    f = open(output_dir + lab + '/' + file + '.txt', 'r')
    for line in f:
        k = line
    l = k.split('$')[1].split('^')[0]
    return float(l)


# In[6]:


# plotting pdfs
import corner

mpl.rcParams['axes.labelsize']=22
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14

label=['θ$_t$/°', r'ϕ$_t/\times10^{-2}$', r'ϕ$_h/\times10^{-2}$', 'σ$_t$/Å']

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 4))
a = (vole * (1 - flatchain[:, :, 3]) * get_value('head5')) 
b = (get_value('vh') * tail_length * (1 - flatchain[:, :, 2]))
angle3 = a / b

new_flat[:, 0] = np.rad2deg(np.arccos(angle3.flatten()))
new_flat[:, 1] = flatchain[:, :, 2].flatten()
new_flat[:, 2] = flatchain[:, :, 3].flatten()
new_flat[:, 3] = flatchain[:, :, 1].flatten()


plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}3_neutron_corner{}.png'.format(figures_dir, lab, num), dpi=600)
plt.close()

