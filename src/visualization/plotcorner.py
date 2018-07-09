
# coding: utf-8

# In[4]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns; sns.set('paper', palette='colorblind')
import numpy as np
from refnx import analysis
import sys


# In[5]:


f = sys.argv[1] #'../../output/dlpc_highconc_chain.txt'
tail_length = float(sys.argv[2]) #1.54 + 1.265 * 12
vole = float(sys.argv[3]) #667.)
flatchain = analysis.curvefitter.load_chain(f)
figures_dir = sys.argv[5] + '/reports/figures/'
lab = sys.argv[4] #'dlpc'


# In[6]:


# plotting pdfs
import corner

mpl.rcParams['axes.labelsize']=22
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14

label=['$V_h$/Å$^3$', '$d_h$/Å', 'θ$_t$/°', r'ϕ$_t/\times10^{-2}$', r'ϕ$_h/\times10^{-2}$', 'σ$_t$/Å', 'σ$_h$/Å']

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 7))

new_flat[:, 0] = list(flatchain[:, :, 2].flatten())
new_flat[:, 2] = list(np.rad2deg(np.arccos(flatchain[:, :, 1].flatten())))
new_flat[:, 3] = list(flatchain[:, :, 5].flatten() * 100)
new_flat[:, 4] = list(flatchain[:, :, 6].flatten() * 100)
new_flat[:, 5] = list(flatchain[:, :, 4].flatten())
new_flat[:, 6] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list((flatchain[:, :, 2].flatten() * tail_length * flatchain[:, :, 1].flatten() * (1 - flatchain[:, :, 5].flatten())))
new_flat[:, 1] = new_flat[:, 1] / (vole)
a = 1 - flatchain[:, :, 6].flatten()
new_flat[:, 1] = new_flat[:, 1] / a


plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}4_all_corner.png'.format(figures_dir, lab), dpi=600)
plt.close()

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 7))

new_flat[:, 0] = list(flatchain[:, :, 2].flatten())
new_flat[:, 2] = list(np.rad2deg(np.arccos(flatchain[:, :, 8].flatten())))
new_flat[:, 3] = list(flatchain[:, :, 11].flatten() * 100)
new_flat[:, 4] = list(flatchain[:, :, 12].flatten() * 100)
new_flat[:, 5] = list(flatchain[:, :, 10].flatten())
new_flat[:, 6] = list(flatchain[:, :, 9].flatten())
new_flat[:, 1] = list((flatchain[:, :, 2].flatten() * tail_length * flatchain[:, :, 8].flatten() * (1 - flatchain[:, :, 11].flatten())))
new_flat[:, 1] = new_flat[:, 1] / (vole)
a = 1 - flatchain[:, :, 12].flatten()
new_flat[:, 1] = new_flat[:, 1] / a


plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}5_all_corner.png'.format(figures_dir, lab), dpi=600)
plt.close()

