## /notebooks

This contains the Jupyter notebook files, and converted Python scripts, used for the analysis of the X-ray and Neutron reflectometry data.

When the `Makefile` in the top level directory is run `nbconvert` generates the `.py` files from the existing `.ipynb` notebooks. These Python scripts are then run for each lipid for XRR and NR.

#### Contents

- `lipid_nr.ipynb/.py` - notebook and Python script for the fitting of the neutron reflectometry (NR) data, and production of the MCMC chain
- `lipid_xrr.ipynb/.py` - notebook and Python script for the fitting of the X-ray reflectometry (XRR) data, and production of the MCMC chain
- `nr_chain_analysis.ipynb/.py` - notebook and Python script for the analysis of the chain from NR and subsequent plotting of the NR, and SLD, profiles.
- `xrr_chain_analysis.ipynb/.py` - notebook and Python script for the analysis of the chain from XRR and subsequent plotting of the XRR, and SLD, profiles, as well as the figures presenting the head volume PDF and the variation of tail layer thickness and solvent fraction in the head layer. 
