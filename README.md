lipids_at_airdes
==============================

Electronic supplementary information for "", this ESI aims to be a fully reproducible analysis of the data discussed within the paper. The supplied Makefile while reproduce all analysis, draw figures and compile the paper when run. 

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
    ├── LICENSE
    ├── README.md
    ├── data
    │   ├── external
    │   ├── interim
    │   ├── processed
    │   └── raw
    ├── docs
    ├── notebooks
    ├── reports
    │   └── figures
    └── src
        ├── data
        ├── external
        ├── models
        ├── tools
        └── visualization
