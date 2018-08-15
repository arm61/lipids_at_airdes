Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers
==============================

Electronic supplementary information for "Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers", this ESI aims to be a fully reproducible analysis of the data discussed within the paper. The supplied Makefile while reproduce all analysis, draw figures and compile the paper when run.

Requirements:

- anaconda python
- make

The analysis can be completely reproduced using the following commands:

```
conda create --name paper_env python

source activate paper_env

pip install -r config/requirements.txt

make clean #this will remove all of the output from previous runs

make
```

Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── LICENSE              # MIT License
    ├── README.md            # You are here
    ├── Makefile             # Makefile to outline workflow
    ├── output               # Files and data output by analysis scripts
    ├── config               # requirements.txt file
	├── data        
    │   ├── external
    │   ├── interim
    │   ├── processed        # Processed XRR data in q, r, dr txt
    │   └── raw
    ├── docs
    ├── notebooks            # Notebooks for analysis
    ├── reports              # Paper and ESI
    │   └── figures
    └── src
        ├── data             
        ├── external  
        ├── models           # mol_vol.py custom model for refnx
        ├── tools            # helper.py script
        └── visualization    # Plotting scripts
