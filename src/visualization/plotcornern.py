
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

label=['θ$_t$/°', r'ϕ$_h/\times10^{-2}$', 'σ$_{t,h,s}$/Å']

print(flatchain.shape)

new_flat = np.zeros((flatchain.shape[0]*flatchain.shape[1], 3))

new_flat[:, 0] = np.rad2deg(np.arccos(flatchain[:, :, 1].flatten()))
new_flat[:, 1] = (1 - ((get_value('vh') * flatchain[:, :, 1].flatten() * tail_length) / (get_value('head4') * get_value('vt')))) * 100
new_flat[:, 2] = flatchain[:, :, 2].flatten()
print(new_flat.shape)


plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}3_neutron_corner{}.png'.format(figures_dir, lab, num), dpi=600)
plt.close()

