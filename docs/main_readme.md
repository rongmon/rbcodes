# RBCodes

A public release of Python codes for astrophysical research by Rongmon Bordoloi. 

[![DOI](https://zenodo.org/badge/192408573.svg)](https://zenodo.org/badge/latestdoi/192408573)

## Overview
This package is constantly under development and will be periodically updated.

## Installation

### Create Conda Environment
```bash
# Create the environment with Python 3.9.6
conda create -n myenv python=3.9.6

# Activate the environment
conda activate myenv

# Install packages (use conda to ensure compatibility)
conda install -c conda-forge --file requirements_simple.txt
```

### Dependencies
- Core Dependencies: 
  - astropy
  - lmfit
  - scipy
  - numpy
  - matplotlib
  - linetools

## Contents

### Graphical User Interfaces (GUIs)
1. [rb_cont.py](/docs/GUIs/rb_cont.md): Simple interactive continuum fitter
2. [rb_spec.py](/docs/GUIs/rb_spec.md): Absorption line analysis pipeline
   - Continuum fitting
   - Equivalent width/column density measurements
   - Simple Voigt profile fitting
3. [rb_interactive_vpfit_singlet.py](/docs/GUIs/rb_interactive_vpfit_singlet.md): Interactive Voigt profile     fitter
4. [rb_plot_spec.py](/docs/GUIs/rb_plot_spec.md): Spectrum plotting and analysis tool
   - Pan and zoom spectrum
   - Plot absorption lines at various redshifts
   - Equivalent width and Gaussian fitting
5. [AbsTools](/docs/GUIs/AbsTools/README.md): Complex absorption line analysis GUI
6. [rb_specgui](/docs/GUIs/rb_specgui/rb_specgui.md): Advanced 1D spectrum viewer and line identifier
7. [zgui/main.py](/docs/GUIs/zgui.md): Redshift measurement GUI for galaxies
   - Supports 1D and 2D spectra
   - Optimized for JWST NIRCam/Grism spectroscopy

8. [LLSFitter_GUI.py](/docs/GUIs/LLSFitter/LLSFitter.md): A Simple GUI to fit Lyman Limit System column densities.

### Intergalactic Medium (IGM) Tools

[Full Documentation](/docs/IGM/IGM_README.md)

    1. compute_EW.py: Equivalent width and column density calculations
    2. rb_setline.py: Atomic transition and f-value finder
    3. rb_iter_contfit.py: Iterative continuum fitting
    4. rb_specbin.py: Spectrum rebinning
    5. ransac_contfit.py: Advanced continuum fitting
    6. lens_sep_to_kpc.py: Sightline separation calculation
    7 [LLSFitter.py](/docs/GUIs/LLSFitter/LLSFitter.md): Lyman Limit System column density measurement. 
    
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
- Package is under active development
- Contributions and improvements are welcome

## License
[Insert license details]

## Contact
[Insert contact information]