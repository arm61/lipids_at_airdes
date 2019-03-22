# ESI for "Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers"

[![PCCP](https://img.shields.io/badge/publication%20DOI-10.1039/C9CP00203K-yellow.svg)](http://dx.doi.org/10.1039/C9CP00203K) [![DOI](https://zenodo.org/badge/144010644.svg)](https://zenodo.org/badge/latestdoi/144010644) [![arXiv](https://img.shields.io/badge/arXiv-1810.07616-orange.svg)](https://arxiv.org/abs/1810.07616) [![License](https://img.shields.io/github/license/arm61/lipids_at_airdes.svg?color=lightgrey)](https://github.com/arm61/lipids_at_airdes/blob/master/LICENSE)

[Andrew R. McCluskey](https://orcid.org/0000-0003-3381-5911)&ast;, [Adrian Sanchez-Fernandez](https://orcid.org/0000-0002-0241-1191), [Karen J. Edler](https://orcid.org/0000-0001-5822-0127), [Stephen C. Parker](https://orcid.org/0000-0003-3804-0975), [Andrew J. Jackson](https://orcid.org/0000-0002-6296-0336), [Richard A. Campbell](https://orcid.org/0000-0002-6296-314X), and [Thomas Arnold](https://orcid.org/0000-0001-7196-7831)&ast;.

&ast;[a.r.mccluskey@bath.ac.uk](mailto:a.r.mccluskey@bath.ac.uk)/[andrew.mccluskey@diamond.ac.uk](mailto:andrew.mccluskey@diamond.ac.uk) & [tom.arnold@esss.se](mailto:tom.arnold@esss.se)

[![ToCFigure](https://raw.githubusercontent.com/arm61/lipids_at_airdes/master/toc_figure/ToC_figure.png)](http://dx.doi.org/10.1039/C9CP00203K)
*A novel reflectometry analysis method reveals the structure of lipid monolayers at the air-DES interface.*

This is the electronic supplementary information (ESI) associated with the publication "Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers".
This ESI provides full details of the analyses performed in the work and access to an automated analysis workflow, through this we aim to provide better access to analysis reproducibility.
The [Supplementary Information document](reports/si.pdf) can be found in the [reports](/reports) folder, alongside a preprint copy of the publication.
For more information about reproducible workflows in Python, check out [Tania Allard's talk from Pycon2018](http://bitsandchips.me/Talks/PyCon.html#/title).

## [Data](https://doi.org/10.15125/BATH-00548)

The reduced X-ray and neutron reflectometry data can be obtained from the University of Bath Research Data Archive.

DOI: [10.15125/BATH-00548](https://doi.org/10.15125/BATH-00548)

The full neutron reflectometry data can be obtained from the ILL Data Archive.

DOI: [10.5291/ILL-DATA.9-13-612](http://doi.org/10.5291/ILL-DATA.9-13-612)

## Analysis

This ESI aims to provide a fully reproducible workflow to the data analysis presented within the paper.

Requirements:

- anaconda or miniconda python
- make
- [REVTeX](https://journals.aps.org/revtex)

The supplied Makefile, will reproduce all of the analysis, plot the figures, and build a preprint version of the paper (`reports/paper.pdf`) when run. Be aware that the analyses within this work are non-trivial and take many hours to run so **use caution** before re-running.

If you **still** want to re-run all of the analysis, please download the [experimental data zip file](https://doi.org/10.15125/BATH-00548), and unzip it (in the `lipids_at_airdes` directory) using the following command:

```
unzip lipids_at_airdes_data.zip
```

Then run the following commands:

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

## Citations

**Paper**: A. R. McCluskey, A. Sanchez-Fernandez, K. J. Edler, S. C. Parker, A. J. Jackson, R. A. Campbell and T. Arnold, *Phys. Chem. Chem. Phys.*, 2019, **21**(11), 6133–6141.

**ESI**: A. R. McCluskey, A. Sanchez-Fernandez, K. J. Edler, S. C. Parker, A. J. Jackson, R. A. Campbell and T. Arnold, lipids_at_airdes (Version 1.0), http://doi.org/10.5281/zenodo.2577796.

**Data**: A. R. McCluskey, A. Sanchez-Fernandez, K. J. Edler, S. C. Parker, A. J. Jackson, R. A. Campbell and T. Arnold, Dataset for ‘Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers’, https://doi.org/10.15125/BATH-00548.

### BibTeX

```
@article{lipids_at_airdes_paper, 
    title={Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers}, 
    volume={21}, 
    DOI={10.1039/C9CP00203K}, 
    number={11}, 
    journal={Phys. Chem. Chem. Phys.}, 
    author={McCluskey, A. R. and Sanchez-Fernandez, A. and Edler, K. J. and Parker, S. C. and Jackson, A. J. and Campbell, R. A. and Arnold, T.}, 
    year={2019}, 
    pages={6133–6141} 
}

@misc{lipids_at_airdes_esi, 
    title={lipids_at_airdes (Version 1.0)}, 
    url={http://doi.org/10.5281/zenodo.2577796}, 
    author={McCluskey, A. R. and Sanchez-Fernandez, A. and Edler, K. J. and Parker, S. C. and Jackson, A. J. and Campbell, R. A. and Arnold, T.}, 
    year={2019} 
}

@misc{lipids_at_airdes_data, 
    title={Dataset for “Bayesian determination of the effect of a deep eutectic solvent on the structure of lipid monolayers”}, 
    url={https://doi.org/10.15125/BATH-00548}, 
    author={McCluskey, A. R. and Sanchez-Fernandez, A. and Edler, K. J. and Parker, S. C. and Jackson, A. J. and Campbell, R. A. and Arnold, T.}, 
    year={2019} 
}
```

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
