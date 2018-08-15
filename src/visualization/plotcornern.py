
# coding: utf-8

# In[ ]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns; sns.set('paper', palette='colorblind')
import numpy as np
from refnx import analysis
import sys


# In[ ]:


f = sys.argv[1] #'../../output/dlpc_highconc_chain.txt'
tail_length = 1.54 + 1.265 * int(sys.argv[2]) 
flatchain = analysis.curvefitter.load_chain(f)
figures_dir = '../../reports/figures/'
output_dir = '../../output/'
lab = sys.argv[3] #'dlpc'
num = sys.argv[4]
def get_value(file, sp):
    f = open(output_dir + lab + '/' + file + sp + '.txt', 'r')
    for line in f:
        k = line
    l = k.split('$')[1].split('^')[0]
    return float(l)


# In[ ]:


# plotting pdfs
import corner

mpl.rcParams['axes.labelsize']=22
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.edgecolor'] = 'k'

label=['θ$_t$/°', r'ϕ$_h/\times10^{-2}$', 'σ$_{t,h,s}$/Å']

print(flatchain.shape)

new_flat = np.zeros((flatchain.shape[0]*flatchain.shape[1], 3))

new_flat[:, 0] = np.rad2deg(np.arccos(flatchain[:, :, 1].flatten()))
new_flat[:, 1] = (1 - ((get_value('vh', '') * flatchain[:, :, 1].flatten() * tail_length) / (
    get_value('head', '') * get_value('vt', '')))) * 100
new_flat[:, 2] = flatchain[:, :, 2].flatten()
print(new_flat.shape)


plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}{}_neutron_corner.pdf'.format(figures_dir, lab, num), dpi=600)
plt.close()

