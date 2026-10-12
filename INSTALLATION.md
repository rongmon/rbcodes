# Installation Guide

Complete instructions for installing rbcodes on your system.

## System Requirements

- **Python**: 3.12
- **OS**: macOS, Linux, or Windows
- **Disk space**: ~500 MB (including dependencies)

## Installation Methods

### Method 1: Using Conda (Recommended)

This is the easiest method as it handles all dependencies correctly.

```bash
# 1. Create a conda environment with Python 3.12
conda create -n rbcodes-env python=3.12
conda activate rbcodes-env

# 2. Clone the repository (if you haven't already)
git clone https://github.com/rongmon/rbcodes.git
cd rbcodes

# 3. Install rbcodes in development mode
pip install -e .
```

### Method 2: Using pip with virtualenv

```bash
# 1. Create a virtual environment
python3.12 -m venv rbcodes-env
source rbcodes-env/bin/activate  # On Windows: rbcodes-env\Scripts\activate

# 2. Clone the repository
git clone https://github.com/rongmon/rbcodes.git
cd rbcodes

# 3. Install rbcodes
pip install -e .
```

### Method 3: Direct pip Installation (for released versions)

Once the package is on PyPI:

```bash
pip install rbcodes
```

## Verify Installation

After installing, run the self-test. It needs no display and takes a few seconds
(up to a minute the first time, while Python and matplotlib build their caches):

```bash
rb_selftest
```

or, from Python:

```python
import rbcodes
rbcodes.selftest()
```

It prints the versions of Python, rbcodes and its dependencies, then checks that every
module imports, that no code uses functions removed from NumPy, SciPy, astropy, pandas
or matplotlib, and that a few calculations (line lists, equivalent widths, NFW halo,
cosmological distances, lens ray tracing, IFU cube reading) give the known answers.
It ends with `ALL CHECKS PASSED` or lists what failed. Use `rb_selftest -v` for full
tracebacks, and include the whole report if you open an issue.

### Full test suite (developers)

From a clone of the repository:

```bash
pip install -e ".[testing]"
pytest
```

## Dependencies

### Core Dependencies (automatically installed)

Lower bounds are the oldest versions that support NumPy 2 and Python 3.12; upper bounds
keep out the next major release until it has been tested. `setup.cfg` is the reference.

```
PyQt5>=5.15.10,<6       # GUI framework
numpy>=2.0,<3           # Numerical computing
scipy>=1.13,<2          # Scientific computing
matplotlib>=3.9,<4      # Plotting
astropy>=6.1,<9         # Astronomy utilities
pandas>=2.2.2,<4        # Data manipulation
linetools>=0.3.2,<0.4   # Spectroscopy tools
emcee>=3.1,<4           # MCMC sampling
photutils>=1.13,<4      # Astronomy photometry
corner>=2.2,<3          # MCMC visualization
scikit-learn>=1.5,<2    # Machine learning
tqdm>=4.65,<5           # Progress bars
regions>=0.9,<1         # ds9 region file I/O (used by rb_ifuview)
ipywidgets>=8.0,<9      # Notebook widgets (interactive continuum fitting)
```

Test dependencies (`pip install -e ".[testing]"`): `pytest`, `pytest-cov`.

The exact set of versions used for testing is in `requirements-py312-lock.txt`; to
reproduce it exactly: `pip install -r requirements-py312-lock.txt && pip install -e . --no-deps`.

### Optional Dependencies

| Package | Purpose | Install |
|---------|---------|---------|
| `pyds9` | Live ds9 bridge for `rb_ifuview` | `rb_ifuview --install` for instructions |

## Troubleshooting

### Installation Issues

#### "ModuleNotFoundError: No module named 'rbcodes'"

**Cause**: Package not installed or wrong environment activated

**Solution**:
```bash
# Make sure you're in the right conda/virtual environment
conda activate rbcodes-env  # or: source rbcodes-env/bin/activate

# Reinstall the package
cd /path/to/rbcodes
pip install -e .
```

#### "error: python setup.py egg_info" failed

**Cause**: Missing setuptools or old version

**Solution**:
```bash
pip install --upgrade setuptools setuptools_scm wheel
pip install -e .
```

#### PyQt5 issues on macOS

**Error**: `NSWindow drag regions should not be empty` or similar

**Solution**: Update PyQt5
```bash
pip install --upgrade PyQt5
```

On M1/M2 Macs, use conda instead:
```bash
conda install -c conda-forge pyqt5
```

#### GUI won't open (blank window or crash)

**Common causes**:
1. Missing Qt backend
2. Display/X11 issues (remote SSH)

**Solutions**:
```bash
# Ensure all GUI dependencies are installed
pip install --upgrade PyQt5 matplotlib

# For remote SSH sessions, use X11 forwarding
ssh -X user@host
```

#### "pandas version 1.3.5" error

**Cause**: Dependency conflict with newer pandas

**Solution**: rbcodes is compatible with pandas 1.3.5. If you need a newer version:
1. Open an issue on GitHub: https://github.com/rongmon/rbcodes/issues
2. Or manually test with newer pandas version:

```bash
pip install pandas==1.4.0  # (at your own risk)
```

#### JWST/NIRCam tools not working

**Cause**: Missing optional data files

**Solution**: Download example data:
```bash
# Download from GitHub
git clone https://github.com/rongmon/rbcodes.git
# Check example-data/ directory
```

### Performance Issues

#### GUIs are slow to start

**Cause**: PyQt5 startup overhead or system resources

**Solution**:
```bash
# Use command line tools instead
python -m rbcodes.IGM.compute_EW  # Example
```

#### Out of memory errors

**Cause**: Processing large spectra or batch files

**Solution**:
1. Process spectra in smaller batches
2. Use spectrum rebinning first: `rb_specbin.py`
3. Increase available RAM or use a more powerful machine

### Platform-Specific Issues

#### Windows

If GUI tools fail to start:

```bash
# Use Python explicitly
python -m rbcodes.GUIs.launch_specgui
```

#### Linux

Ensure X11 is available (for remote machines):
```bash
# Test X11 display
echo $DISPLAY
# If empty, X11 forwarding is not enabled
```

#### macOS

On Apple Silicon (M1/M2/M3) Macs, use conda for better compatibility:

```bash
conda create -n rbcodes-env -c conda-forge python=3.12
conda activate rbcodes-env
pip install -e .
```

## Updating rbcodes

If you installed in development mode (`pip install -e .`):

```bash
cd /path/to/rbcodes
git pull origin master
# No need to reinstall - changes are picked up automatically
```

If you installed normally:

```bash
pip install --upgrade rbcodes
```

## Uninstalling rbcodes

```bash
pip uninstall rbcodes
```

If installed in development mode:

```bash
pip uninstall rbcodes
cd /path/to/rbcodes
python setup.py develop --uninstall
```

## Getting Help

If you encounter issues not covered here:

1. Check [Existing Issues](https://github.com/rongmon/rbcodes/issues)
2. Open a [New Issue](https://github.com/rongmon/rbcodes/issues/new) with:
   - The full output of `rb_selftest`
   - Full error message and traceback
   - Steps to reproduce

3. Check the [Full Documentation](docs/main_readme.md)

## Development Installation

To set up a development environment for contributing:

```bash
# Create environment
conda create -n rbcodes-dev python=3.12
conda activate rbcodes-dev

# Clone and install
git clone https://github.com/rongmon/rbcodes.git
cd rbcodes
pip install -e ".[testing]"

# Run the tests
pytest

# Run tests
pytest
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for more details.
