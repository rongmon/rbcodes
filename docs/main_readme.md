# rbcodes

A public release of Python codes for astrophysical research by Rongmon Bordoloi. 

## Table of Contents
- [Quick Start](#quick-start)
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

## Quick Start

**Want to dive in?** Here are the most commonly used tools:

```bash
# Launch the main spectrum analysis pipeline
launch_specgui

# View and compare multiple spectra
rb_multispec

# Fit Lyman Limit Systems
rb_llsfitter

# Estimate redshifts (esp. for JWST)
rb_zgui
```

For more details, see the individual tool documentation below.

### 📖 Citation

If you use this package, please cite it using the following DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6079263.svg)](https://doi.org/10.5281/zenodo.6079263)


## Overview
This package is constantly under development and will be periodically updated.

## Installation

Quick setup with conda and pip:

```bash
# Create a new conda environment (Python 3.10 recommended)
conda create -n rbcodes-env python=3.10
conda activate rbcodes-env

# Install rbcodes (dependencies installed automatically)
cd /path/to/rbcodes
pip install -e .
```

**Optional**: Pre-install dependencies via conda for better compatibility (not required):
```bash
conda install -c conda-forge --file requirements_simple.txt
```

For more installation options and troubleshooting, see [INSTALLATION.md](../INSTALLATION.md).

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
1. [launch_specgui.py](/docs/GUIs/rb_spec/rb_spec.md): **Main tool** - Comprehensive absorption line analysis pipeline (recommended starting point)
   - Interactive GUI wrapper for rb_spec
   - Continuum fitting
   - Equivalent width/column density measurements
   - Simple Voigt profile fitting
2. [multispecviewer](/docs/GUIs/multispec/multispec.md): Enhanced spectrum viewer for handling multiple 1D spectra simultaneously with advanced line identification and absorber cataloging capabilities.
3. [interactive_continuum_fit.py](/docs/GUIs/interactive_continuum_fit.py): Interactive continuum fitting tool for manual continuum definition
4. [AbsTools](/docs/GUIs/AbsTools/README.md): Complex absorption line analysis GUI with batch processing capabilities
5. [rb_zgui (zgui/main.py)](/docs/GUIs/zgui/Tutorial_for_Emission_Line_Redshift_Estimator_GUI.pdf): Redshift measurement GUI (PDF tutorial)
   - Supports 1D and 2D spectra
   - Optimized for JWST NIRCam/Grism spectroscopy
6. [LLSFitter_GUI.py](/docs/GUIs/LLSFitter/LLSFitter.md): GUI to fit Lyman Limit System column densities
7. [rb_cont.py](/docs/GUIs/rb_cont.md): Interactive continuum fitter

### Intergalactic Medium (IGM) Tools

[Full Documentation](/docs/IGM/IGM_README.md)

1. [compute_EW.py](/docs/IGM/compute_EW_readme.md): Equivalent width and column density calculations with dynamic saturation detection
2. [rb_setline.py](/docs/IGM/rb_setline.md): Atomic transition and f-value finder (supports 12+ line lists)
3. [rb_iter_contfit.py](/docs/IGM/rb_iter_contfit.md): Advanced iterative polynomial continuum fitting with sigma clipping and BIC-based polynomial order selection
4. [fit_continuum_full_spec.py](/docs/IGM/rb_iter_contfit.md): Full spectrum continuum fitting using the iterative polynomial approach
5. [rb_specbin.py](/docs/IGM/rb_specbin.md): Spectrum rebinning utility
6. [ransac_contfit.py](/docs/IGM/): RANSAC-based continuum fitting (advanced alternative to polynomial fitting)
7. [lens_sep_to_kpc.py](/docs/IGM/lens_sep_to_kpc.md): Sightline separation calculations in physical units
8. [LLSFitter.py](/docs/GUIs/LLSFitter/LLSFitter.md): Lyman Limit System column density measurement 
    
### Halo Analysis

1. **rb_nfw.py**: Compute NFW (Navarro-Frenk-White) dark matter halo profiles
2. **mstar2mhalo.py**: Convert stellar mass to halo mass using Moster et al. 2010 relations

### Statistical Tools

[Full Documentation](/docs/rbstat/rb_stat_readme.md)

1. **rb_wilsonscore.py**: Wilson score confidence intervals
2. **rb_boot.py**: Bootstrap resampling function

### Lensing

1. **lens_ang_sep.py**: Deflection matrix ray tracing for gravitational lensing calculations

### Utility Modules

1. **rb_utility.py**: General utility functions (progress reporting, color lists, etc.)
2. **rb_spectrum.py**: Spectrum handling and I/O utilities
3. **rb_x1d_id.py**: Inspect HST/COS x1d FITS file headers (exposure type, object type, etc.)
4. **readmultispec.py**: Read IRAF echelle spectra in multispec format (supports linear, log, cubic spline, Chebyshev, Legendre dispersion)
5. **cos_inspect.py**: Quick inspection utility for HST/COS x1d files
6. **filter_2d_spec.py**: Mask and filter JWST/NIRCam 2D grism spectra (EIGER survey optimized)
7. **compute_SNR_1d.py**: Compute and visualize signal-to-noise ratio of 1D spectra
8. **compute_telescope_offset.py**: Telescope offset calculations
9. **galaxy_qso_pa.py**: Position angle calculations for galaxies and QSOs
10. **rgb_images.py**: RGB image generation from multi-band data
11. **write_ds9_regions.py**: Write DS9 region files

### Catalog Tools

1. **rb_search.py**: Cone search around (RA, Dec) coordinates against a list of objects
2. **convert_FIRE_coordinates.py**: Transform Magellan FIRE spectrum coordinates to J2000 epoch
3. **galaxy_group_finder.py**: Identify galaxy groups and associations

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.rst](../CONTRIBUTING.rst) for detailed guidelines on:
- How to set up a development environment
- The fork-and-pull-request workflow
- Code contribution standards
- Testing requirements

For major changes, please open an issue first to discuss what you would like to change.
## ![License](https://img.shields.io/badge/license-MIT-green)

This project is licensed under the [MIT License](LICENSE).
