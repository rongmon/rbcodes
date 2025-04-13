# rb_specgui

A comprehensive tool for visualizing and analyzing astronomical spectra, with a focus on absorption line identification.

## Overview

The rb_specgui is a Python-based tool designed for astronomers to analyze spectral data, with a particular focus on identifying and characterizing absorption lines in astronomical spectra. The tool provides an interactive interface for visualizing spectra, identifying absorption systems at various redshifts, measuring equivalent widths, and cataloging absorption features.

This tool is particularly useful for:
- Analyzing quasar spectra to identify intervening absorption systems
- Studying Lyman-alpha systems, metal line absorbers, and other absorption features
- Measuring equivalent widths and column densities of absorption lines
- Creating catalogs of absorption systems at various redshifts

## Installation

### Prerequisites

1. Install rbcodes requirements using conda:

```bash
# Create a new conda environment (optional but recommended)
conda create -n rbcodes python=3.9.6
conda activate rbcodes

# Install requirements from requirements_simple.txt
conda install -c conda-forge --file /path/to/rbcodes/requirements_simple.txt
```

Refer to the `requirements_simple.txt` file in the rbcodes repository for specific package details.

### Setup

1. The tool is part of the `rbcodes` distribution, specifically within the `rbcodes.GUIs` package

2. For convenient access from any location, create an alias in your shell configuration file (`.cshrc`):

```csh
# Add to your .cshrc file
alias rb_specgui 'python /path/to/rbcodes/GUIs/rb_specgui.py'
```

For example:
```csh
alias rb_specgui 'python /Users/bordoloi/WORK/python/rbcodes/GUIs/rb_specgui.py'
```

3. After adding the alias, reload your shell configuration:

```bash
source ~/.cshrc
```

4. Now you can run the tool from any directory using:

```bash
rb_specgui -h  # Show help
```

## Usage

The main entry point is `rb_specgui.py`, which can be executed from the command line with various options.

### Command Line Options

```
usage: rb_specgui.py [-h] [-t {ascii,linetools,fits,lt,lt_cont_norm,p}] [-e EFIL] [-n] [-v] [-z REDSHIFT] filename
```

- `filename`: Spectrum file to load
- `-t, --filetype`: File format type (ascii, linetools, fits, lt, lt_cont_norm, p)
- `-e, --error-file`: Optional separate error file
- `-n, --normalize`: Apply continuum normalization if available
- `-v, --verbose`: Enable verbose output
- `-z, --redshift`: Initial redshift (zabs) value (default: 0.0)

### Example Commands

Using the alias:

Load a FITS file with initial redshift of 0.5:
```bash
rb_specgui spectrum.fits -z 0.5
```

Load an ASCII file with separate error file:
```bash
rb_specgui spectrum.txt -t ascii -e error.txt
```

Load a normalized LineTools spectrum:
```bash
rb_specgui spectrum.fits -t lt_cont_norm
```

Get help on available options:
```bash
rb_specgui -h
```

## Workflow Diagram

```
                                 ┌─────────────────┐
                                 │  Load Spectrum  │
                                 │  rb_specgui.py  │
                                 └────────┬────────┘
                                          │
                                          ▼
                          ┌───────────────────────────────┐
                          │      Main GUI Interface       │
                          │   PlotSpec_Integrated.py      │
                          └───────────────┬───────────────┘
                                          │
                        ┌─────────────────┴──────────────────┐
                        │                                     │
                ┌───────┴────────┐              ┌─────────────┴──────────┐
                │                │              │                        │
                ▼                ▼              ▼                        ▼
        ┌──────────────┐ ┌──────────────┐ ┌──────────┐          ┌───────────────┐
        │   Absorber   │ │  Equivalent  │ │  Manual  │          │ Velocity Stack│
        │ Identification│ │Width Analysis│ │Transition│          │   Analysis    │
        │              │ │              │ │Identification│       │ (VStack GUI)  │
        └──────┬───────┘ └──────┬───────┘ └─────┬────┘          └───────┬───────┘
               │                │                │                       │
               ▼                ▼                ▼                       ▼
        ┌──────────────┐ ┌──────────────┐ ┌──────────────┐      ┌───────────────┐
        │ Add Absorber │ │ Measure EW & │ │  Identify    │      │  Categorize   │
        │  to Catalog  │ │Column Density│ │  Absorption  │      │  Transitions  │
        │              │ │              │ │   Systems    │      │               │
        └──────┬───────┘ └──────────────┘ └──────────────┘      └───────┬───────┘
               │                                                         │
               └─────────────────────┬─────────────────────────────────┘
                                     │
                                     ▼
                              ┌──────────────┐
                              │  Save/Load   │
                              │   Catalog    │
                              └──────────────┘