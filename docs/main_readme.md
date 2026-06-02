# rbcodes

Python tools for astrophysical spectroscopy research by Rongmon Bordoloi and collaborators.

## Table of Contents

- [Quick Start](#quick-start)
- [Installation](#installation)
- [Contents](#contents)
  - [Core Pipeline — rb_spec](#core-pipeline--rb_spec)
  - [Graphical User Interfaces](#graphical-user-interfaces)
  - [IGM Tools](#igm-tools)
  - [Halo Analysis](#halo-analysis)
  - [Statistical Tools](#statistical-tools)
  - [Lensing](#lensing)
  - [Utility Modules](#utility-modules)
  - [Catalog Tools](#catalog-tools)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Quick Start

```bash
# Interactive spectrum analysis (GUI)
launch_specgui

# View and compare multiple 1D spectra
rb_multispec

# View IFU datacubes and extract spectra
rb_ifuview

# Fit Lyman Limit Systems
rb_llsfitter

# Estimate redshifts (JWST / grism)
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

---

## Installation

```bash
# 1. Create a conda environment (Python 3.10 recommended)
conda create -n rbcodes-env python=3.10
conda activate rbcodes-env

# 2. Install rbcodes
cd /path/to/rbcodes
pip install -e .
```

**Optional** — pre-install dependencies via conda for better binary compatibility:
```bash
conda install -c conda-forge --file requirements_simple.txt
```

For full installation options, platform notes, and troubleshooting see [INSTALLATION.md](../INSTALLATION.md).

Core dependencies: `numpy`, `astropy`, `matplotlib`, `PyQt5`, `linetools`, `scipy`. See [INSTALLATION.md](../INSTALLATION.md) for the full list and optional packages.

---

## Contents

### Core Pipeline — rb_spec

| | |
|-|-|
| **Command** | `from rbcodes.GUIs.rb_spec import rb_spec as r` |
| **Docs** | [Full Documentation](GUIs/rb_spec/rb_spec.md) |

Central Python class for absorption line spectroscopy — load, shift, slice, fit continuum, measure EW and column density, save results. Backbone behind `launch_specgui`.

---

### Graphical User Interfaces

#### Active Tools

| Tool | Command | Description | Docs |
|------|---------|-------------|------|
| **launch_specgui** | `launch_specgui` | Point-and-click interface for `rb_spec`; single-spectrum and batch modes | [docs](GUIs/rb_spec/rb_spec.md) |
| **rb_multispec** | `rb_multispec` | Multi-spectrum viewer with line identification, redshift overlay, absorber cataloging, vStack, and quick/advanced line fitting | [docs](GUIs/multispec/multispec.md) |
| **rb_ifuview** | `rb_ifuview cube.fits` | IFU datacube viewer — load KCWI/MUSE/generic cubes, extract spectra, compute moment maps, WCS alignment via rb_align, import ds9 regions | [docs](GUIs/ifuview/rb_ifuview.md) |
| **rb_zgui** | `rb_zgui` | Redshift measurement GUI optimized for JWST NIRCam/grism spectroscopy | [PDF tutorial](GUIs/zgui/Tutorial_for_Emission_Line_Redshift_Estimator_GUI.pdf) |
| **rb_zfind** | `rb_zfind` | Semi-automated redshift finder — PCA, template, and picket-fence chi2 search; integrates with rb_zgui | [docs](GUIs/zfind/rb_zfind.md) |
| **rb_llsfitter** | `rb_llsfitter` | GUI to fit Lyman Limit System column densities | [docs](GUIs/LLSFitter/LLSFitter.md) |
| **interactive_continuum_fit** | (import) | Recommended continuum fitter — polynomial and spline methods with interactive masking | [docs](GUIs/interactive_continuum_fit.md) |
| **rb_align** | (import) | Astrometry alignment — align IFU cubes and 2D images to a reference frame with interactive or automated source matching and full WCS fitting | [docs](GUIs/rb_align/rb_align.md) |

**Legacy tools:** `AbsTools` and `rb_cont.py` are superseded by `launch_specgui` (batch mode) and `interactive_continuum_fit` respectively.

---

### IGM Tools

[Full Documentation](IGM/IGM_README.md)

| Module | Description |
|--------|-------------|
| [compute_EW.py](IGM/compute_EW_readme.md) | Equivalent width and column density with saturation detection |
| [rb_setline.py](IGM/rb_setline.md) | Atomic transition and f-value lookup (12+ line lists) |
| [rb_iter_contfit.py](IGM/rb_iter_contfit.md) | Iterative polynomial continuum fitting with sigma clipping and BIC order selection; `fit_continuum_full_spec.py` applies this to a full spectrum |
| [rb_specbin.py](IGM/rb_specbin.md) | Spectrum rebinning utility |
| [ransac_contfit.py](IGM/) | RANSAC-based continuum fitting |
| [lens_sep_to_kpc.py](IGM/lens_sep_to_kpc.md) | Sightline separation in physical units |
| [LLSFitter.py](GUIs/LLSFitter/LLSFitter.md) | Lyman Limit System column density measurement |

---

### Halo Analysis

[Full Documentation](halo/halo_readme.md)

Standalone modules for dark matter halo calculations. Available for direct import.

| Module | Description |
|--------|-------------|
| **rb_nfw.py** | NFW halo profiles — density, circular velocity, escape velocity, dispersion (Hoeft et al. 2004) |
| **halo_profile.py** | NFW escape velocity, virial radius ↔ mass, concentration–mass relation (Dutton & Maccio 2014) |
| **mstar2mhalo.py** | Stellar mass → halo mass and virial radius (Moster et al. 2010 SMHM) |
| **halos.py** | Stellar mass → halo mass via Behroozi UniverseMachine SMHM; computes R_200 |

---

### Statistical Tools

[Full Documentation](rbstat/rb_stat_readme.md)

| Module | Description |
|--------|-------------|
| **rb_wilsonscore.py** | Wilson score confidence intervals |
| **rb_boot.py** | Bootstrap resampling |

---

### Lensing

[Full Documentation](lensing/lensing_readme.md)

| Module | Description |
|--------|-------------|
| **lens_ang_sep.py** | Gravitational lensing ray-tracing — maps image-plane coordinates to source plane using deflection matrices. Includes cosmological distance utilities (`cosmic_D`, `ang_D12`). (A. Shaban 2020, R. Bordoloi 2021) |

---

### Utility Modules

[Full Documentation](utils/utils_readme.md)

| Module | Description |
|--------|-------------|
| **rb_spectrum.py** | Core spectrum I/O — loads FITS, ASCII, and other formats; used internally by `rb_spec` and the GUIs |
| **rb_utility.py** | General utility functions (progress reporting, color lists, etc.) |
| **rb_x1d_id.py** | Inspect HST/COS x1d FITS headers (exposure type, object type) |
| **readmultispec.py** | Read IRAF echelle spectra in multispec format (linear, log, cubic spline, Chebyshev, Legendre) |
| **cos_inspect.py** | Quick inspection utility for HST/COS x1d files |
| **filter_2d_spec.py** | Mask and filter JWST/NIRCam 2D grism spectra (EIGER survey optimized) |
| **compute_SNR_1d.py** | Compute and visualize 1D spectrum SNR |
| **compute_telescope_offset.py** | Telescope offset calculations |
| **galaxy_qso_pa.py** | Position angle calculations for galaxy–QSO pairs |
| **rgb_images.py** | RGB image generation from multi-band data |
| **write_ds9_regions.py** | Write ds9 region files |

---

### Catalog Tools

[Full Documentation](catalog/catalog_readme.md)

| Module | Description |
|--------|-------------|
| **rb_search.py** | Cone search around (RA, Dec) against an object catalog |
| **convert_FIRE_coordinates.py** | Transform Magellan FIRE spectrum coordinates to J2000 |
| **galaxy_group_finder.py** | Friend-of-Friends group finder with angular, comoving, and physical distance linking |
| **GalaxyGroupFinder_slow.py** | Reference FoF implementation (easier to read; use `galaxy_group_finder.py` for production) |

---

## Citation

If you use rbcodes in your research, please cite:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6079263.svg)](https://doi.org/10.5281/zenodo.6079263)

<details>
<summary>BibTeX — version 2.0.0</summary>

```bibtex
@software{bordoloi_2025_15723701,
  author       = {Bordoloi, Rongmon and
                  Liu, Bin and
                  Clark, Sean and
                  Higginson, Jack and
                  Flores, Derick},
  title        = {rongmon/rbcodes: rbcodes v2.0.0},
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v2.0.0},
  doi          = {10.5281/zenodo.6079263},
  url          = {https://doi.org/10.5281/zenodo.6079263}
}
```
</details>

---

## Contributing

Contributions are welcome. See [CONTRIBUTING.rst](../CONTRIBUTING.rst) for guidelines on setting up a development environment, the fork-and-pull-request workflow, and testing requirements. For major changes please open an issue first.

---

## License

This project is licensed under the [MIT License](../LICENSE).
