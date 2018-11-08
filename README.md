# ESI for "Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers"

[![arXiv](https://img.shields.io/badge/arXiv-1810.07616-orange.svg)](https://arxiv.org/abs/1810.07616) [![DOI](https://zenodo.org/badge/144010644.svg)](https://zenodo.org/badge/latestdoi/144010644)

[Andrew R. McCluskey](https://orcid.org/0000-0003-3381-5911)&ast;, [Adrian Sanchez-Fernandez](https://orcid.org/0000-0002-0241-1191), [Karen J. Edler](https://orcid.org/0000-0001-5822-0127), [Stephen C. Parker](https://orcid.org/0000-0003-3804-0975), [Andrew J. Jackson](https://orcid.org/0000-0002-6296-0336), [Richard A. Campbell](https://orcid.org/0000-0002-6296-314X), and [Thomas Arnold](https://orcid.org/0000-0001-7196-7831)&ast;.

&ast;[a.r.mccluskey@bath.ac.uk](mailto:a.r.mccluskey@bath.ac.uk) & [tom.arnold@esss.se](mailto:tom.arnold@esss.se)

This is the electronic supplementary information (ESI) associated with the publication "Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers".
This ESI provides full details of the analyses performed in the work and access to an automated analysis workflow, through this we aim to provide better access to analysis reproduciblility.
The [Supplementary Information document](reports/si.pdf) can be found in the [reports](/reports) folder, alongside a preprint copy of the publication.
For more information about reproducible workflows in Python, check out [Tania Allard's talk from Pycon2018](http://bitsandchips.me/Talks/PyCon.html#/title).

## [Data](https://researchdata.bath.ac.uk/id/eprint/548)

The reduced X-ray  and neutron reflectometry data can be obtained from the University of Bath Research Data Archive.

[https://researchdata.bath.ac.uk/id/eprint/548](https://researchdata.bath.ac.uk/id/eprint/548)

The full neutron reflectometry data can be obtained from the ILL Data Archive.

DOI: [10.5291/ILL-DATA.9-13-612](http://doi.org/10.5291/ILL-DATA.9-13-612)

## Analysis

This ESI aims to provide a fully reproducible workflow to the data analysis presented within the paper.

Requirements:

- anaconda or miniconda python
- make
- [REVTeX](https://journals.aps.org/revtex)

The supplied Makefile, will reproduce all of the analysis, plot the figures, and build a preprint version of the paper (`reports/preprint.pdf`) when run. Be aware that the analyses within this work are non-trivial and take many hours to run so **use caution** before re-running.

If you **still** want to re-run all of the analysis, please download the [experimental data](https://researchdata.bath.ac.uk/id/eprint/548), place it in a directory named `data/processed` before running the following commands:

```
conda create --name paper_env python

source activate paper_env

pip install --upgrade pip

pip install -r config/requirements.txt

make clean # this will remove all of the output from previous runs

make
```

## [Figures](/reports/figures)

PDF versions of the figures, can be found in the `reports/figures` directory.

## Acknowledgements

A. R. M. is grateful to the University of Bath and Diamond Light Source for co-funding a studentship (Studentship Number STU0149). The authors thank the European Spallation Source and the University of Bath Alumni Fund for supporting A. S.-F.

## Project Organization

    .
    ├── AUTHORS.md
    ├── LICENSE              # CC-BY-SA-4.0 License
    ├── README.md            # You are here
    ├── Makefile             # Makefile to outline workflow
    ├── output               # Files and data output by analysis scripts
        ├── dlpc             #
        ├── dmpc             #
        ├── dmpg             #
        └── dppc             #
    ├── config               # requirements.txt file
    ├── notebooks            # Notebooks for analysis
    ├── reports              # Paper and ESI
    │   └── figures
    └── src
        ├── models           # mol_vol.py custom model for refnx
        ├── tools            # helper.py script
        └── visualization    # Plotting scripts
