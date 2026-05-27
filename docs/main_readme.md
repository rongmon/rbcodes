# rbcodes

A public release of Python codes for astrophysical research by Rongmon Bordoloi. 

## Table of Contents
- [Quick Start](#quick-start)
- [Citation](#-citation)
- [Overview](#overview)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Contents](#contents)
  - [Core Pipeline — rb_spec.py](#core-pipeline--rb_specpy)
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
# Launch the main spectrum analysis pipeline (GUI)
launch_specgui

# View and compare multiple spectra
rb_multispec

# View IFU datacubes and extract spectra
rb_ifuview

# Fit Lyman Limit Systems
rb_llsfitter

# Estimate redshifts (esp. for JWST)
rb_zgui
```

Or use `rb_spec` directly as a Python API:

```python
from rbcodes.GUIs.rb_spec import rb_spec as r

spec = r.from_file('spectrum.fits', filetype='linetools')
spec.shift_spec(zabs=0.5)
spec.slice_spec(1548.2, -1500, 1500, use_vel=True)
spec.fit_continuum(mask=[-200, 200], Legendre=3)
spec.compute_EW(1548.2, vmin=-150, vmax=150)
spec.save_slice('results.json', file_format='json')
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
  - numpy >= 1.22.3, < 1.24  *(upper bound required — numpy 1.24+ breaks compatibility)*
  - matplotlib >= 3.5.0
  - astropy >= 5.3.3
  - linetools >= 0.3
  - scipy >= 1.7.3
  - pandas == 1.3.5  *(pinned — newer versions not yet tested)*
  - emcee >= 3.0
  - photutils >= 1.0
  - corner >= 2.0
  - scikit-learn >= 1.5.0
  - tqdm >= 4.65.0
  - regions >= 0.5
  - pytest >= 6.0
  - pytest-cov >= 2.0

- Optional (install separately):
  - pyds9 — live ds9 bridge for `rb_ifuview` (run `rb_ifuview --install` for instructions)

## Contents

### Core Pipeline — rb_spec.py

[Full Documentation](/docs/GUIs/rb_spec/rb_spec.md)

`rb_spec` is the central Python class in `rbcodes`. It is the programmable backbone behind `launch_specgui` and the recommended entry point for all absorption line analysis — both interactive and scripted.

Key capabilities:
- Load spectra from FITS, ASCII, XSpectrum1D, or numpy arrays
- Shift to absorber rest frame and slice around any transition
- Continuum fitting: interactive spline, interactive masking GUI, iterative Legendre polynomial, RANSAC, BIC-optimized order selection
- Equivalent width and column density (AOD) measurements with full error propagation
- Saturation detection, SNR estimation, velocity centroid and dispersion
- Save/load full analysis state to JSON

```python
from rbcodes.GUIs.rb_spec import rb_spec as r, load_rb_spec_object

# Load and analyse
spec = r.from_file('spectrum.fits', filetype='linetools')
spec.shift_spec(zabs=0.348)
spec.slice_spec(1548.2, -1500, 1500, use_vel=True)
spec.fit_continuum(mask=[-200, 200], Legendre=3)
spec.compute_EW(1548.2, vmin=-150, vmax=150, plot=True)
spec.save_slice('results.json', file_format='json')

# Reload a saved analysis
spec = load_rb_spec_object('results.json')
```

See the [full documentation](/docs/GUIs/rb_spec/rb_spec.md) for the complete API, continuum fitting options, output attributes, and worked examples.

---

### Graphical User Interfaces (GUIs)

[Full Auto-generated Documentation](/docs/GUIs/GUIs_readme.md)

#### Active Tools
1. [launch_specgui](/docs/GUIs/rb_spec/rb_spec.md): **GUI launcher** for `rb_spec` — provides a point-and-click interface around the core pipeline; supports single-spectrum and batch modes
2. [multispecviewer](/docs/GUIs/multispec/multispec.md): Enhanced spectrum viewer for handling multiple 1D spectra simultaneously with advanced line identification and absorber cataloging capabilities.
3. **rb_ifuview** (`rb_ifuview cube.fits`): QFitsView-style IFU datacube viewer — load KCWI, MUSE, or any 3-axis FITS cube; extract 1D spectra interactively; compute moment maps; import ds9 regions; send spectra to `rb_multispec`.
   - Supports Single pixel / Rectangular / Circular / Circular-Annular apertures
   - Variance-weighted and optimal (data/Gaussian) spectral extraction
   - M0 / M1 / M2 moment maps with continuum subtraction and SNR masking
   - Per-dataset state persistence (extractions, image mode, moment params, sky region)
   - Optional ds9 bridge: `rb_ifuview --install` for setup instructions
4. [interactive_continuum_fit.py](/docs/GUIs/interactive_continuum_fit.md): **Recommended continuum fitter** - Interactive tool with polynomial and spline methods, manual masking, and real-time feedback
5. [rb_zgui (zgui/main.py)](/docs/GUIs/zgui/Tutorial_for_Emission_Line_Redshift_Estimator_GUI.pdf): Redshift measurement GUI (PDF tutorial)
   - Supports 1D and 2D spectra
   - Optimized for JWST NIRCam/Grism spectroscopy
6. [LLSFitter_GUI.py](/docs/GUIs/LLSFitter/LLSFitter.md): GUI to fit Lyman Limit System column densities
7. [AbsTools](/docs/GUIs/AbsTools/README.md): ⚠️ **LEGACY** - Complex absorption line analysis GUI with batch processing capabilities (use `launch_specgui` in batch mode instead for new analyses)
8. [rb_cont.py](/docs/GUIs/rb_cont.md): ⚠️ **LEGACY** - Older continuum fitter (use `interactive_continuum_fit.py` instead)

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

[Full Auto-generated Documentation](/docs/halo/halo_readme.md)

Standalone utility modules — available for direct import, not used internally by other package modules.

1. **rb_nfw.py**: Compute NFW (Navarro-Frenk-White) halo profiles — density, circular velocity, escape velocity, and velocity dispersion (Hoeft et al. 2004)
2. **halo_profile.py**: NFW escape velocity, virial radius ↔ mass conversion, and concentration–mass relation (Dutton & Maccio 2014)
3. **mstar2mhalo.py**: Stellar mass → halo mass and virial radius using Moster et al. (2010) SMHM relation
4. **halos.py**: Stellar mass → halo mass using the Behroozi UniverseMachine SMHM relation; also computes R_200 (requires `smhm_med_params.txt` parameter file)

### Statistical Tools

[Full Documentation](/docs/rbstat/rb_stat_readme.md)

1. **rb_wilsonscore.py**: Wilson score confidence intervals
2. **rb_boot.py**: Bootstrap resampling function

### Lensing

[Full Auto-generated Documentation](/docs/lensing/lensing_readme.md)

Standalone utility module — available for direct import.

1. **lens_ang_sep.py**: Gravitational lensing ray-tracing — transports image-plane coordinates to a source plane (or any intermediate redshift plane) using deflection matrices from a lens model. Includes supporting cosmological distance calculations (`cosmic_D`, `ang_D12`). Originally written by Ahmed Shaban (2020), extended by Rongmon Bordoloi (2021).

### Utility Modules

[Full Auto-generated Documentation](/docs/utils/utils_readme.md)

1. **rb_spectrum.py**: Core spectrum I/O and handling — loads FITS, ASCII, and other formats; used internally by `rb_spec` and the GUIs
2. **rb_utility.py**: General utility functions (progress reporting, color lists, etc.)
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

[Full Auto-generated Documentation](/docs/catalog/catalog_readme.md)

Standalone utility modules — available for direct import.

1. **rb_search.py**: Cone search around (RA, Dec) coordinates against a list of objects
2. **convert_FIRE_coordinates.py**: Transform Magellan FIRE spectrum coordinates to J2000 epoch
3. **galaxy_group_finder.py**: Friend-of-Friends galaxy group finder with angular, comoving, and physical distance linking
4. **GalaxyGroupFinder_slow.py**: Reference (slower) FoF implementation — easier to read and verify; prefer `galaxy_group_finder.py` for production use

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.rst](../CONTRIBUTING.rst) for detailed guidelines on:
- How to set up a development environment
- The fork-and-pull-request workflow
- Code contribution standards
- Testing requirements

For major changes, please open an issue first to discuss what you would like to change.
## ![License](https://img.shields.io/badge/license-MIT-green)

This project is licensed under the [MIT License](LICENSE).
