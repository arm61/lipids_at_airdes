lipids_at_airdes
==============================

Electronic supplementary information for "Probabilistic determination of the effect of a deep eutectic solvent on the structure of lipid monolayers", this ESI aims to be a fully reproducible analysis of the data discussed within the paper. The supplied Makefile while reproduce all analysis, draw figures and compile the paper when run. 

The analysis can be completely reproduced using the following commands:

```
conda create --name paper_env --file config/spec-file.txt

source activate paper_env

cd src/data/refnx-0.0.15

python setup.py build

python setup.py install

python setup.py test

cd ../../

make
```

Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── LICENSE              # MIT License
    ├── README.md            # You are here
    ├── Makefile             # Makefile to outline workflow
    ├── config               # spec_file.txt for building the conda env
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
        ├── data             # refnx v0.0.15
        ├── external  
        ├── models           # mol_vol.py custom model for refnx
        ├── tools      
        └── visualization    # Plotting scripts
