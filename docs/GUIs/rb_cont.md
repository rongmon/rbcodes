# Interactive Continuum Fitting Tools

This README provides documentation for the interactive continuum fitting tools in the `rbcodes.GUIs` package:
- `rb_cont.py` - A command-line script for launching the interactive fitter
- `rb_fit_interactive_continuum.py` - The core class for interactive continuum fitting

These tools allow you to interactively fit a continuum to a 1D spectrum by selecting points on a plot and fitting a cubic spline through them.

## Table of Contents
- [Overview](#overview)
- [Usage](#usage)
  - [Command-Line Usage](#command-line-usage)
  - [Python API Usage](#python-api-usage)
  - [Setting Up an Alias](#setting-up-an-alias)
- [Interactive Controls](#interactive-controls)
- [Examples](#examples)
- [Best Practices](#best-practices)
- [Known Issues](#known-issues)

## Overview

The interactive continuum fitter provides a user-friendly interface for fitting the continuum of spectroscopic data. It allows users to:

- Load spectral data from various file formats (FITS, ASCII, pickle)
- Interactively select points along the continuum
- Fit a cubic spline through the selected points
- View the normalized spectrum
- Save the fitted continuum and normalized spectrum

## Usage

### Command-Line Usage

The `rb_cont.py` script can be run directly from the command line:

```bash
python /path/to/rbcodes/GUIs/rb_cont.py filename [filetype]
```

Where:
- `filename` is the path to the spectrum file
- `filetype` (optional) is the type of file ('fits', 'ascii', 'p' for pickle). If not provided, it will be inferred from the file extension.

To view detailed help about controls:

```bash
python /path/to/rbcodes/GUIs/rb_cont.py --help-controls
```

### Python API Usage

You can also use the `rb_fit_interactive_continuum` class directly in your Python code:

```python
import numpy as np
from rbcodes.GUIs.rb_fit_interactive_continuum import rb_fit_interactive_continuum

# Example with wavelength, flux, and error arrays
wave = np.linspace(4000, 5000, 1000)
flux = np.ones_like(wave) + 0.1 * np.sin(wave/100)
error = 0.05 * np.ones_like(wave)

# Initialize the interactive fitter
fitter = rb_fit_interactive_continuum(wave, flux, error)

# After interactive fitting and saving (pressing 'w'), 
# the continuum is available as:
cont = fitter.cont

# Normalize the spectrum
norm_flux = flux / cont
```

Example with loading a spectrum file first:

```python
import numpy as np
from astropy.io import fits
from rbcodes.GUIs.rb_fit_interactive_continuum import rb_fit_interactive_continuum

# Load a FITS file
hdul = fits.open('spectrum.fits')
data = hdul[1].data
wave = data['WAVE']
flux = data['FLUX']
error = data['ERROR']

# Initialize the interactive fitter
fitter = rb_fit_interactive_continuum(wave, flux, error)

# After interactively fitting, get the continuum
cont = fitter.cont
```

### Setting Up an Alias

For convenient access to the continuum fitting tool, you can create an alias in your shell configuration file:

For C shell (`.cshrc` or `.tcshrc`):
```bash
alias rb_cont 'python /path/to/rbcodes/GUIs/rb_cont.py'
```

For Bash (`.bashrc` or `.bash_profile`):
```bash
alias rb_cont='python /path/to/rbcodes/GUIs/rb_cont.py'
```

After setting up the alias, you can simply run:
```bash
rb_cont spectrum.fits
```

## Interactive Controls

The interactive fitter provides the following controls:

### Mouse Controls
- **Left Click**: Select the median flux value within ±2.5 units from the x-coordinate
- **Right Click**: Delete the nearest continuum point

### Keyboard Controls
- **b**: Select a point at the exact (x,y) coordinate of the cursor
- **Enter**: Perform a spline fit using the selected points
- **n**: Show the normalized spectrum
- **w**: Save the continuum (after pressing 'n')
- **h**: Display help screen
- **r**: Reset fit (clear all points)
- **q**: Quit the interactive session

## Examples

### Example Workflow

1. Load a spectrum file:
   ```bash
   rb_cont spectrum.fits
   ```

2. In the interactive plot:
   - Left-click to select points along the continuum
   - Use 'b' for precise point selection where needed
   - Right-click to remove incorrect points
   - Press 'Enter' to fit a cubic spline through the points
   - Press 'n' to view the normalized spectrum
   - Press 'w' to save the continuum
   - Press 'q' to quit

3. The normalized spectrum and continuum will be saved as `spectrum_norm.fits`

### Example Figures

[//]: # (Figure 1: Example of selecting points on a spectrum)
**Figure 1**: Example of selecting points on a spectrum

[//]: # (Figure 2: Fitted continuum on the spectrum)
**Figure 2**: Fitted continuum on the spectrum

[//]: # (Figure 3: Normalized spectrum)
**Figure 3**: Normalized spectrum

## Best Practices

1. **Point Selection**
   - Select points on relatively flat, absorption-free regions of the spectrum
   - Ensure points are distributed across the entire wavelength range
   - Use more points in regions where the continuum has complex curvature
   - Use fewer points in regions with a simple, flat continuum

2. **Spline Fitting**
   - Use at least 4 points for cubic spline fitting
   - Check the fit by pressing 'Enter' frequently as you add points
   - If the fit looks poor, add or remove points and refit

3. **Normalization**
   - Always check the normalized spectrum (press 'n') before saving
   - The normalized spectrum should have a mean value of approximately 1.0 in continuum regions
   - Absorption features should dip below 1.0

4. **File Handling**
   - The original file is not modified
   - A new file is created with '_norm' appended to the original filename
   - For FITS files, the continuum is stored in the FITS file for later use

## Known Issues

1. **Interface Issues**
   - The interactive fitter works best when no toolbar buttons in the plot window are activated
   - On some systems, the plot window may become unresponsive if too many points are selected
   - Occasionally, the keyboard shortcuts may not register if the plot window doesn't have focus

2. **Data Issues**
   - Spectra with large gaps or discontinuities may cause issues with the spline fitting
   - Extremely noisy spectra may benefit from smoothing before continuum fitting
   - Very large files may cause memory issues or slow performance

3. **Compatibility**
   - Requires matplotlib with Qt5Agg backend
   - May not work with older versions of Python (tested with Python 3.6+)
   - Some advanced features require linetools and astropy to be installed

If you encounter any other issues, please report them to the package maintainer.