# MultispecViewer

A tool for visualizing and analyzing multiple spectroscopic datasets simultaneously.

⚠️ **WARNING: UNDER ACTIVE DEVELOPMENT** ⚠️

This tool is still in development. Features and interfaces are subject to change, particularly in the non-classic version. The classic version provides stable functionality while new features are being integrated into the main version.

## Overview

MultispecViewer is a PyQt-based GUI tool that allows users to:

- Load and display multiple FITS spectra simultaneously
- Apply redshifts and identify spectral features
- Perform quick identification of common ionic species
- Adjust view limits interactively
- Catalog multiple absorber systems at different redshifts

## Usage

MultispecViewer can be run from the command line using:

```bash
python /path/to/rbcodes/GUIs/multispecviewer/rb_multispec.py
```

For convenience, you can create an alias:

```bash
alias multispec='python /Users/bordoloi/WORK/python/rbcodes/GUIs/multispecviewer/rb_multispec.py'
```

Then simply run:

```bash
multispec
```

### Command Line Options

```
usage: rb_multispec.py [-h] [-c] [-v]

rb_multispec - A tool for visualizing spectral data

optional arguments:
  -h, --help      show this help message and exit
  -c, --classic   Run the classic version of multispecviewer
  -v, --version   Display version information
```

## Features

### Interface

The MultispecViewer interface consists of several components:

- A toolbar for selecting FITS files
- A main plotting area for displaying spectra
- A redshift input panel for applying and cataloging redshifts
- A message box for displaying status information
- An absorber manager for cataloging multiple systems (new version only)

[Insert screenshot of UI here]

### Loading Data

1. Click the "Select FITS Files" button
2. Navigate to and select one or more FITS format spectral files
3. The spectra will be loaded and displayed in the main plotting area

### Navigation Controls

MultispecViewer supports the following keyboard shortcuts for navigating the display:

| Key | Action |
|-----|--------|
| `x` | Set minimum x-limit to cursor position |
| `X` | Set maximum x-limit to cursor position |
| `t` | Set maximum y-limit to cursor position (per panel) |
| `b` | Set minimum y-limit to cursor position (per panel) |
| `Y` | Manually input y-limits range |
| `[` | Shift view left by one viewport width |
| `]` | Shift view right by one viewport width |
| `o` | Zoom out (x-axis) |
| `r` | Reset view to original state |
| `R` | Clear all line identifications |
| `q` | Quit application |

### Spectral Processing

| Key | Action |
|-----|--------|
| `S` | Increase smoothing (convolution kernel size) |
| `U` | Decrease smoothing (convolution kernel size) |

### Quick Line Identification

Press the following keys with the cursor positioned at a spectral feature to identify it:

| Key | Line |
|-----|------|
| `C` | CIV doublet |
| `M` | MgII doublet |
| `F` | FeII multiplet |
| `6` | OVI doublet |
| `4` | SiIV doublet |
| `8` | NeVIII doublet |
| `2` | Lyman-beta |
| `1` | Lyman-alpha |

Right-clicking anywhere on the spectrum opens a dialog with a complete list of lines from the selected line list at the clicked wavelength.

### Redshift Controls

1. Enter a redshift value in the input field
2. Select a line list from the dropdown menu
3. Choose a color (new version only)
4. Click "Submit" to mark lines at the specified redshift

In the new version, clicking "Catalog" adds the current redshift to the absorber manager for reference.

## Classic vs New Version

### Classic Version

The classic version provides the core functionality in a streamlined interface:
- Multiple spectrum display
- Redshift application
- Line identification
- Basic navigation

[Insert classic version screenshot here]

### New Version

The new version adds several enhancements:
- Dark theme for improved visibility
- Absorber manager for cataloging multiple systems
- Color selection for different redshift systems
- Enhanced metadata display
- Improved performance and error handling

[Insert new version screenshot here]

## Examples

### Loading and Viewing Multiple Spectra

[Insert example screenshot here]

### Identifying Features at Various Redshifts

[Insert example screenshot here]

### Using the Absorber Manager (New Version)

[Insert example screenshot here]

## Development

MultispecViewer is part of the rbcodes package for spectroscopic analysis. It utilizes:
- PyQt5 for the graphical interface
- Matplotlib for plotting
- linetools for spectral data handling
- astropy for astronomical calculations

