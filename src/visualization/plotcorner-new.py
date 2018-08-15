
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


label=['$V_t$/Å$^3$', '$V_h$/Å$^3$', '$d_h$/Å', 'θ$_t$/°', r'ϕ$_h/\times10^{-2}$', 'σ$_{t,h,s}$/Å']

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 6))

new_flat[:, 0] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list(flatchain[:, :, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, :, 2].flatten())))
tail3 = flatchain[:, :, 2].flatten() * tail_length
solh3 = 1 - (flatchain[:, :, 4].flatten() / flatchain[:, :, 3].flatten()) * (tail3 / flatchain[:, :, 1].flatten())
new_flat[:, 4] = list(solh3 * 100)
new_flat[:, 5] = list(flatchain[:, :, 5].flatten())
new_flat[:, 2] = list(flatchain[:, :, 1].flatten())

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}3_all_corner.png'.format(figures_dir, lab), dpi=600)
plt.close()

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 6))

new_flat[:, 0] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list(flatchain[:, :, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, :, 7].flatten())))
tail4 = flatchain[:, :, 7].flatten() * tail_length
solh4 = 1 - (flatchain[:, :, 4].flatten() / flatchain[:, :, 3].flatten()) * (tail4 / flatchain[:, :, 1].flatten())
new_flat[:, 4] = list(solh4 * 100)
new_flat[:, 5] = list(flatchain[:, :, 8].flatten())
new_flat[:, 2] = list(flatchain[:, :, 1].flatten())

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}4_all_corner.png'.format(figures_dir, lab), dpi=600)
plt.close()

new_flat = np.zeros((flatchain.shape[0] * flatchain.shape[1], 6))

new_flat[:, 0] = list(flatchain[:, :, 3].flatten())
new_flat[:, 1] = list(flatchain[:, :, 4].flatten())
new_flat[:, 3] = list(np.rad2deg(np.arccos(flatchain[:, :, 10].flatten())))
tail5 = flatchain[:, :, 10].flatten() * tail_length
solh5 = 1 - (flatchain[:, :, 4].flatten() / flatchain[:, :, 3].flatten()) * (tail5 / flatchain[:, :, 1].flatten())
new_flat[:, 4] = list(solh5 * 100)
new_flat[:, 5] = list(flatchain[:, :, 11].flatten())
new_flat[:, 2] = list(flatchain[:, :, 2].flatten())

plt1 = corner.corner(new_flat, max_n_ticks=3, labels=label)
plt.savefig('{}{}5_all_corner.png'.format(figures_dir, lab), dpi=600)
plt.close()


# In[ ]:


#lab = [0'scale2', 1'angle2', 2'vt2', 3'vh', 4'rough2', 5'solh2', 
#       6'scale3', 7'angle3', 8'rough3', 9'solh3', 
#       10'scale4', 11'angle4', 12'vt4', 13'rough4', 14'solh4', 
#       15'scale5', 16'angle5', 17'rough5', 18'solh5']

