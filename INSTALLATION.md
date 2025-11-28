# Installation Guide

Complete instructions for installing rbcodes on your system.

## System Requirements

- **Python**: 3.9.6 to 3.10.x
- **OS**: macOS, Linux, or Windows
- **Disk space**: ~500 MB (including dependencies)

## Installation Methods

### Method 1: Using Conda (Recommended)

This is the easiest method as it handles all dependencies correctly.

```bash
# 1. Create a conda environment with Python 3.10
conda create -n rbcodes-env python=3.10
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
python3 -m venv rbcodes-env
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

After installation, verify that everything works:

```bash
# Test 1: Import the package
python -c "import rbcodes; print(rbcodes.__version__)"

# Test 2: Check console entry points
which launch_specgui
which rb_multispec
which rb_llsfitter
which rb_zgui

# Test 3: Launch the main GUI (closes immediately if working)
launch_specgui --help
```

## Dependencies

### Core Dependencies (automatically installed)

```
PyQt5>=5.15.7         # GUI framework
numpy>=1.22.3         # Numerical computing
matplotlib>=3.5.0     # Plotting
astropy>=5.3.3        # Astronomy utilities
linetools>=0.3        # Spectroscopy tools
scipy>=1.7.3          # Scientific computing
pandas==1.3.5         # Data manipulation
emcee>=3.0            # MCMC sampling
photutils>=1.0        # Astronomy photometry
corner>=2.0           # MCMC visualization
scikit-learn>=1.5.0   # Machine learning
tqdm>=4.65.0          # Progress bars
pytest>=6.0           # Testing
```

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
conda create -n rbcodes-env -c conda-forge python=3.10 pyqt5
conda activate rbcodes-env
pip install -e .
```

## Updating rbcodes

If you installed in development mode (`pip install -e .`):

```bash
cd /path/to/rbcodes
git pull origin main
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
   - Python version: `python --version`
   - rbcodes version: `python -c "import rbcodes; print(rbcodes.__version__)"`
   - Full error message and traceback
   - Steps to reproduce

3. Check the [Full Documentation](docs/main_readme.md)

## Development Installation

To set up a development environment for contributing:

```bash
# Create environment
conda create -n rbcodes-dev python=3.10
conda activate rbcodes-dev

# Clone and install
git clone https://github.com/rongmon/rbcodes.git
cd rbcodes
pip install -e .

# Install development tools
pip install pytest pytest-cov tox sphinx

# Run tests
pytest
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for more details.
