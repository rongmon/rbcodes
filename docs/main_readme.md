# RBCodes

A public release of Python codes for astrophysical research by Rongmon Bordoloi. 

[![DOI](https://zenodo.org/badge/192408573.svg)](https://zenodo.org/badge/latestdoi/192408573)

## Overview
This package is constantly under development and will be periodically updated.

## Installation

### Create Conda Environment
```bash
# Create the environment with Python 3.9.5
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
- Partial Dependencies: PysimpleGUI for some GUIs

## Contents

### Graphical User Interfaces (GUIs)
1. [rb_cont.py](/docs/GUIs/rb_cont.md): Simple interactive continuum fitter
2. [rb_spec.py](/docs/GUIs/rb_spec.md): Absorption line analysis pipeline
   - Continuum fitting
   - Equivalent width/column density measurements
   - Simple Voigt profile fitting
3. [rb_interactive_vpfit_singlet.py](/docs/GUIs/rb_interactive_vpfit_singlet.md): Interactive Voigt profile fitter
4. [rb_plot_spec.py](/docs/GUIs/rb_plot_spec.md): Spectrum plotting and analysis tool
   - Pan and zoom spectrum
   - Plot absorption lines at various redshifts
   - Equivalent width and Gaussian fitting
5. [AbsTools](/docs/GUIs/AbsTools/README.md): Complex absorption line analysis GUI
6. [rb_specgui](/docs/GUIs/rb_specgui/rb_specgui.md): Advanced 1D spectrum viewer and line identifier
7. [zgui/main.py](/docs/GUIs/zgui.md): Redshift measurement GUI for galaxies
   - Supports 1D and 2D spectra
   - Optimized for JWST NIRCam/Grism spectroscopy

### Intergalactic Medium (IGM) Tools

[Full Documentation](/docs/IGM/IGM_README.md)

1. [compute_EW.py](/docs/IGM/compute_EW.md): Equivalent width and column density calculations
2. [rb_setline.py](/docs/IGM/rb_setline.md): Atomic transition and f-value finder
3. [rb_iter_contfit.py](/docs/IGM/rb_iter_contfit.md): Iterative continuum fitting
4. [rb_specbin.py](/docs/IGM/rb_specbin.md): Spectrum rebinning
5. [ransac_contfit.py](/docs/IGM/ransac_contfit.md): Advanced continuum fitting
6. [lens_sep_to_kpc.py](/docs/IGM/lens_sep_to_kpc.md): Sightline separation calculation

### Halo Analysis
1. [rb_nfw.py](/docs/halo/rb_nfw.md): NFW halo profile computation
2. [mstar2mhalo.py](/docs/halo/mstar2mhalo.md): Stellar mass to halo mass conversion

### Statistical Tools
1. [rb_wilsonscore.py](/docs/rbstat/rb_wilsonscore.md): Wilson score confidence intervals
2. [rb_boot.py](/docs/rbstat/rb_boot.md): Bootstrap function

### Lensing
1. [lens_ang_sep.py](/docs/lensing/lens_ang_sep.md): Deflection matrix ray tracing

### Utility Modules
1. [rb_utility.py](/docs/utils/rb_utility.md): Utility functions
2. [rb_x1d_id.py](/docs/utils/rb_x1d_id.md): HST/COS x1d fits file header info
3. [readmultispec.py](/docs/utils/readmultispec.md): IRAF spectrum reader
4. [cos_inspect.py](/docs/utils/cos_inspect.md): HST/COS file inspector
5. [filter_2d_spec.py](/docs/utils/filter_2d_spec.md): JWST/NIRCam 2D grism spectra filter
6. [compute_SNR_1d.py](/docs/utils/compute_SNR_1d.md): Signal-to-noise ratio computation

### Catalog Tools
1. [rb_search.py](/docs/catalog/rb_search.md): Cone search utility
2. [convert_FIRE_coordinates.py](/docs/catalog/convert_FIRE_coordinates.md): Coordinate transformation

## Contributing
- Package is under active development
- Contributions and improvements are welcome

## License
[Insert license details]

## Contact
[Insert contact information]