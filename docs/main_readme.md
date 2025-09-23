# rbcodes

A public release of Python codes for astrophysical research by Rongmon Bordoloi. 

## Table of Contents
- [Citation](#-citation)
- [Overview](#overview)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Contents](#contents)
  - [Graphical User Interfaces (GUIs)](#graphical-user-interfaces-guis)
  - [Intergalactic Medium (IGM) Tools](#intergalactic-medium-igm-tools)
  - [Halo Analysis](#halo-analysis)
  - [Statistical Tools](#statistical-tools)
  - [Lensing](#lensing)
  - [Utility Modules](#utility-modules)
  - [Catalog Tools](#catalog-tools)
- [Contributing](#contributing)
- [License](#license)

### 📖 Citation

If you use this package, please cite it using the following DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6079263.svg)](https://doi.org/10.5281/zenodo.6079263)


## Overview
This package is constantly under development and will be periodically updated.

## Installation

```bash
# Recommended: Python 3.10 (also supports 3.9.6)
conda create -n rbcodes-env python=3.10

# Activate the environment
conda activate rbcodes-env

# Install dependencies (use conda to ensure compatibility)
conda install -c conda-forge --file requirements_simple.txt

# Install package in development mode
cd /path/to/rbcodes
pip install -e .
```


### Dependencies
- Core Dependencies (automatically installed):
  - PyQt5 >= 5.15.7
  - numpy >= 1.22.3
  - matplotlib >= 3.5.0
  - astropy >= 5.3.3
  - linetools >= 0.3
  - scipy >= 1.7.3
  - pandas == 1.3.5
  - emcee >= 3.0
  - photutils >= 1.0
  - corner >= 2.0
  - scikit-learn >= 1.5.0
  - tqdm >= 4.65.0

## Contents

### Graphical User Interfaces (GUIs)

#### Active Tools
1. [multispecviewer](/docs/GUIs/multispec/multispec.md): Enhanced spectrum viewer for handling multiple 1D spectra simultaneously with advanced line identification and absorber cataloging capabilities.
2. [rb_spec.py](/docs/GUIs/rb_spec/rb_spec.md): Comprehensive absorption line analysis pipeline
   - Continuum fitting
   - Equivalent width/column density measurements
   - Simple Voigt profile fitting
3. [launch_specgui.py](/docs/GUIs/rb_spec/launch_specgui.md): An interactive GUI wrapper for rb_spec pipeline.
4. [AbsTools](/docs/GUIs/AbsTools/README.md): Complex absorption line analysis GUI
5. [Tutorial for zgui/main.py (PDF)](/docs/GUIs/zgui/Tutorial_for_Emission_Line_Redshift_Estimator_GUI.pdf): Redshift measurement GUI for galaxies
   - Supports 1D and 2D spectra
   - Optimized for JWST NIRCam/Grism spectroscopy
6. [LLSFitter_GUI.py](/docs/GUIs/LLSFitter/LLSFitter.md): A Simple GUI to fit Lyman Limit System column densities.
7. [rb_cont.py](/docs/GUIs/rb_cont.md): Simple interactive continuum fitter
8. [rb_interactive_vpfit_singlet.py](/docs/GUIs/rb_interactive_vpfit_singlet.md): Interactive Voigt profile fitter

#### Deprecated Tools
9. [rb_specgui](/docs/GUIs/rb_specgui/rb_specgui.md): ⚠️ **DEPRECATED** - Advanced 1D spectrum viewer and line identifier (Superseded by multispecviewer)
10. [rb_plot_spec.py](/docs/GUIs/rb_plot_spec.md): ⚠️ **DEPRECATED** - Spectrum plotting and analysis tool

### Intergalactic Medium (IGM) Tools

[Full Documentation](/docs/IGM/IGM_README.md)

1. [compute_EW.py](/docs/IGM/compute_EW_readme.md): Equivalent width and column density calculations
2. [rb_setline.py](/docs/IGM/rb_setline.md): Atomic transition and f-value finder
3. [rb_iter_contfit.py](/docs/IGM/rb_iter_contfit.md): Advanced Iterative polynomial continuum fitting with sigma clipping. Has options to select polynomial order using Bayesian Information Criterion (BIC),
4. [fit_continuum_full_spec.py](/docs/IGM/rb_iter_contfit.md): Use the previous code to fit continuum to a full QSO spectrum. Look at documentation for examples.
5. [rb_specbin.py](/docs/IGM/rb_specbin.md): Spectrum rebinning
6. ransac_contfit.py: Advanced continuum fitting
7. [lens_sep_to_kpc.py](/docs/IGM/lens_sep_to_kpc.md): Sightline separation calculation
8. [LLSFitter.py](/docs/GUIs/LLSFitter/LLSFitter.md): Lyman Limit System column density measurement. 
    
### Halo Analysis
       1) rb_nfw.py      : Compute NFW halo profile
       2) mstar2mhalo.py : Convert stellar mass to Halo mass using Moster et al. 2010.

### Statistical Tools
[Full Documentation](/docs/rbstat/rb_stat_readme.md)

    1. rb_wilsonscore.py: Wilson score confidence intervals
    2. rb_boot.py: Bootstrap function

### Lensing
    1. lens_ang_sep.py: Deflection matrix ray tracing

### Utility Modules
     1) rb_utility.py     :  Several utility functions

      2) rb_x1d_id.py      : This function will read all HST/COS x1d fits files in a folder and print out
                             the header info for raw spectra. It will print filename, exposure type, 
                             and object type entries from the header.

      3) readmultispec.py  : Read IRAF (echelle) spectrum in multispec format from a FITS file. Can read 
                             most multispec formats including linear, log, cubic spline, Chebyshev or 
                             Legendre dispersion spectra. I got this code from https://github.com/kgullikson88/General.

     4) cos_inspec.py      : Routine to quickly inspect HST/COS x1d files.

     5) filter_2d_spec.py  : Custom routine to mask and filter JWST/NIRCam 2D grism spectra from the EIGER survey. [Optimized for a custom data format]

     6) compute_SNR_1d.py  : Compute signal-to-noise ratio of an 1D spectrum and plot the result.
             

### Catalog Tools
    1)rb_search.py.    : Function to do a cone search around any (ra,dec) pointing with respect to a list of  ra,dec entries.
    2) convert_FIRE_coordinates.py:  Custom code only to take Magellan 1D FIRE spectrum and use the header  information
                                          to transform the co-ordinates to J2000 epoch.

## Contributing

Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Submit a pull request with a clear description of your changes.

For major changes, please open an issue first to discuss what you would like to change.
## ![License](https://img.shields.io/badge/license-MIT-green)

This project is licensed under the [MIT License](LICENSE).
