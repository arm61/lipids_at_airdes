
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
tail_length = 1.54 + 1.265 * int(sys.argv[2])
flatchain = analysis.curvefitter.load_chain(f)
figures_dir = '../../reports/figures/'
lab = sys.argv[3] #'dlpc'


# In[6]:


# plotting pdfs
import corner

mpl.rcParams['axes.labelsize']=22
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.edgecolor'] = 'k'


label=['$V_t$/Å$^3$', '$V_h$/Å$^3$', '$d_h$/Å', 'θ$_t$/°', r'ϕ$_h/\times10^{-2}$', 'σ$_{t,h,s}$/Å']

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 6))

new_flat[:, 0] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list(flatchain[:, :, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, :, 2].flatten())))
new_flat[:, 5] = list(flatchain[:, :, 5].flatten())
new_flat[:, 2] = list(flatchain[:, :, 1].flatten())
new_flat[:, 4] = list((1 - (flatchain[:, :, 4].flatten() / flatchain[:, :, 3].flatten()) * (
    flatchain[:, :, 2].flatten() * tail_length / flatchain[:, :, 1].flatten())) * 100)

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}2_all_corner.pdf'.format(figures_dir, lab), dpi=600)
plt.close()

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 6))

new_flat[:, 0] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list(flatchain[:, :, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, :, 7].flatten())))
new_flat[:, 5] = list(flatchain[:, :, 8].flatten())
new_flat[:, 2] = list(flatchain[:, :, 1].flatten())
new_flat[:, 4] = list((1 - (flatchain[:, :, 4].flatten() / flatchain[:, :, 3].flatten()) * (
    flatchain[:, :, 7].flatten() * tail_length / flatchain[:, :, 1].flatten())) * 100)

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}3_all_corner.pdf'.format(figures_dir, lab), dpi=600)
plt.close()

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 6))

new_flat[:, 0] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list(flatchain[:, :, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, :, 10].flatten())))
new_flat[:, 5] = list(flatchain[:, :, 11].flatten())
new_flat[:, 2] = list(flatchain[:, :, 1].flatten())
new_flat[:, 4] = list((1 - (flatchain[:, :, 4].flatten() / flatchain[:, :, 3].flatten()) * (
    flatchain[:, :, 10].flatten() * tail_length / flatchain[:, :, 1].flatten())) * 100)


plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}4_all_corner.pdf'.format(figures_dir, lab), dpi=600)
plt.close()


new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 6))

new_flat[:, 0] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list(flatchain[:, :, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, :, 13].flatten())))
new_flat[:, 5] = list(flatchain[:, :, 14].flatten())
new_flat[:, 2] = list(flatchain[:, :, 1].flatten())
new_flat[:, 4] = list((1 - (flatchain[:, :, 4].flatten() / flatchain[:, :, 3].flatten()) * (
    flatchain[:, :, 13].flatten() * tail_length / flatchain[:, :, 1].flatten())) * 100)
                      

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}5_all_corner.pdf'.format(figures_dir, lab), dpi=600)
plt.close()


# In[ ]:


#lab = [0'scale4', 1'angle4', 2'vt', 3'vh', 4'rough4', 5'solh4', 
#       6'scale5', 7'angle5', 8'rough5', 9'solh5']

