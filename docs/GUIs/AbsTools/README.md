# AbsTools: Absorption Line Analysis Toolbox

AbsTools is an interactive toolkit for analyzing absorption lines in astronomical spectra. It provides a user-friendly graphical interface for continuum fitting, equivalent width measurements, and column density determinations of absorption features in spectra of the Circumgalactic Medium (CGM), Intergalactic Medium (IGM), and Interstellar Medium (ISM).

## Overview

AbsTools is part of the `rbcodes` package under `rbcodes.GUIs.abstools`. This README provides comprehensive instructions on using the toolbox, including the new launcher GUI.

### Features

- **Interactive Continuum Fitting**: Fit polynomials to the spectrum continuum with adjustable orders and masked regions
- **Equivalent Width Measurements**: Calculate equivalent widths and errors with user-defined velocity integration limits
- **Column Density Calculations**: Determine column densities directly from absorption features
- **Intervening Absorption Line Analysis**: Identify and analyze intervening absorption systems
- **Multi-page Support**: Analyze up to 30 transitions simultaneously across 5 tabbed pages
- **Enhanced Save/Load Functionality**: Save your analysis in various formats (JSON, Pickle, PDF, data tables)
- **Improved Error Handling**: Better error recovery and prevention of segmentation faults
- **Signal-Slot Architecture**: Robust communication between application components
- **User-Friendly Launcher**: New GUI for setting up and launching analysis sessions

## Dependencies

AbsTools requires the following libraries:
- Python 3.6 or higher
- NumPy
- Matplotlib
- PyQt5
- Pandas
- Astropy
- linetools (for spectrum handling)

## Quick Start

There are now three ways to start using AbsTools:

### 1. Using the AbsTools Launcher GUI (Recommended)

The new launcher GUI provides the simplest way to get started:

```bash
python abstools_launcher.py
```

This opens a user-friendly interface where you can:
- Load a spectrum file
- Set the redshift
- Select absorption lines to analyze
- Configure analysis parameters
- Load previously saved analysis sessions

### 2. Using the AbsTools API Directly

For more control or scripting, you can use the AbsTools API:

```python
from linetools.spectra.xspectrum1d import XSpectrum1D
from rbcodes.GUIs.abstools import Absorber as A
from rbcodes.GUIs.abstools import Metal_Plot as M

# Read in the 1D spectrum to be analyzed
from pkg_resources import resource_filename
filename = resource_filename('rbcodes', 'example-data/test.fits')
sp = XSpectrum1D.from_file(filename)
wave = sp.wavelength.value
flux = sp.flux.value
error = sp.sig.value

# Specify redshift at which to perform analysis
z = 0.348

# Specify rest-frame wavelengths of absorption lines to analyze
# Common lines: OVI doublet (1031.93, 1037.62), CIV doublet (1548.20, 1550.78)
lines = [1031.93, 1037.62, 1215.67, 1548.20, 1550.78]

# Create an absorber class with the data
absys = A.Absorber(z, wave, flux, error, lines=lines, window_lim=[-2000, 2000])
Abs = absys.ions

# Launch the main GUI
M.Transitions(Abs)

# Optional: if you have pre-identified intervening absorption lines
# M.Transitions(Abs, intervening='intervening_lines.txt')
```

### 3. Loading a Saved Analysis Session

To resume a previously saved analysis:

```python
from rbcodes.GUIs.abstools import Metal_Plot as M

# For JSON files (recommended format)
from rbcodes.GUIs.abstools.json_utils import load_from_json
ions = load_from_json('analysis.json')

# For pickle files
import pickle
with open('analysis.p', 'rb') as f:
    ions = pickle.load(f)

# Launch GUI with loaded data
M.Transitions(ions)
```

## AbsTools Launcher GUI

The new AbsTools Launcher provides a user-friendly interface for setting up analysis sessions.

### Launcher Features

- **Multiple Launch Modes**: Start a new analysis or load a saved session
- **Spectrum Selection**: Browse for and load FITS, ASCII, or HDF5 spectrum files
- **Absorption Line Presets**: Quickly select common absorption lines
- **Parameter Configuration**: Set redshift, velocity window limits, and other parameters
- **Settings Persistence**: Save and recall your settings between sessions
- **Detailed Status Updates**: Track the progress of your analysis setup

### Using the Launcher

1. **Start the Launcher**:
   ```bash
   python abstools_launcher.py
   ```

2. **For New Analysis**:
   - Select the "Load New Spectrum" tab
   - Browse for your spectrum file
   - Enter the redshift
   - Select absorption lines (choose from presets or enter manually)
   - Set velocity window limits
   - Optionally specify an intervening absorbers file
   - Click "Launch AbsTools"

3. **For Saved Analysis**:
   - Select the "Load Saved Analysis" tab
   - Browse for your saved analysis file (.json or .p)
   - Optionally specify an intervening absorbers file
   - Click "Launch AbsTools"

4. **Save Settings**:
   - Click "Save Settings" to store your configuration for future use

## Main Interface Guide

### Main Interface Layout

The AbsTools interface is divided into two columns:
- **Left Side**: Raw spectrum with continuum fitting tools
- **Right Side**: Normalized spectrum with velocity integration tools

Each row represents a different absorption line transition.

### Mouse Controls

| Action | Description |
|--------|-------------|
| Left panel + Left mouse button (double click) | Add wavelength region to continuum fit |
| Left panel + Right mouse button (double click) | Remove wavelength region from continuum fit |
| Right panel + Left mouse button | Set lower velocity integration limit |
| Right panel + Right mouse button | Set upper velocity integration limit |

### Keyboard Controls

| Key | Description |
|-----|-------------|
| Up arrow | Increase polynomial order for continuum fit |
| Down arrow | Decrease polynomial order for continuum fit |
| v | Manually enter regions to mask (left panel) or integration limits (right panel) |
| V | Apply current velocity limits to all subplots |
| m | Measure equivalent width for active subplot |
| M | Measure equivalent width for all subplots |
| 0 | Flag absorber as positive detection |
| 1 | Flag absorber as upper limit |
| 2 | Flag absorber as lower limit |
| t | Toggle text display (equivalent width or column density) |
| q | Exit application |

## Workflow Example

Here's a typical workflow for analyzing absorption lines:

1. **Load Spectrum**: Start AbsTools Launcher and select your spectrum file
2. **Continuum Fitting**:
   - Click on the left panel to select a transition
   - Use left/right mouse buttons to add/remove regions for continuum fitting
   - Adjust polynomial order with up/down arrow keys
   - Areas excluded from the fit appear grayed out

3. **Measure Equivalent Widths**:
   - Click on the right panel to select a transition
   - Set velocity integration limits using left/right mouse buttons
   - Press 'm' to measure the equivalent width
   - For uncertain detections, flag them with keys 0/1/2

4. **Save Results**:
   - Click the "Save" button to open the save dialog
   - Choose to save as PDF (plots), data table, pickle file, or JSON (recommended)

## Working with Saved Analysis Files

AbsTools can save analysis in multiple formats:
- **JSON**: Human-readable, portable format (recommended)
- **Pickle**: Binary Python format
- **Data Table**: ASCII tables of measurements
- **PDF**: Plots of continuum fits and normalized spectra

## Utility Scripts

AbsTools includes several utility scripts for working with analysis data:

### Extracting Measurements

```python
from rbcodes.GUIs.abstools.extract_abstools_measurements import extract_measurements_table

# Extract measurements to a DataFrame
file_path = "Spectrum_Analysis_z_0.348.json"
df = extract_measurements_table(file_path)

# View the measurements
print(f"Measurements for spectrum at z = {df.attrs['redshift']}")
print(df)

# Export to CSV
df.to_csv("measurements.csv", index=False)
```

### Batch Processing

```python
from rbcodes.GUIs.abstools.batch_process_abstools import batch_process_abstools_files, plot_ew_comparison

# Process all JSON files in a directory
json_pattern = "data/Spectrum_Analysis_*.json"
measurements_df = batch_process_abstools_files(json_pattern)

# Create comparison plot of specific ions across spectra
ions_to_compare = ['OVI 1031.93', 'OVI 1037.62', 'CIV 1548.20', 'CIV 1550.78']
fig, ax = plot_ew_comparison(measurements_df, ions_to_compare)
fig.savefig("ew_comparison.png", dpi=300)
```

### Extracting Raw Data

```python
from rbcodes.GUIs.abstools.extract_abstools_data import load_abstools_analysis

# Load the full analysis data
file_path = "Spectrum_Analysis_z_0.348.json"
ions = load_abstools_analysis(file_path)

# Access flux, wavelength, and other data for an ion
ion_name = list(ions.keys())[0]  # First ion
wavelength = ions[ion_name]['wave']
normalized_flux = ions[ion_name]['flux'] / ions[ion_name]['cont']
velocity = ions[ion_name]['vel']

# Now you can work with the raw data in your own code
import matplotlib.pyplot as plt
plt.figure()
plt.step(velocity, normalized_flux, 'k', where='mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Normalized Flux')
plt.show()
```

## Architecture Overview

AbsTools follows a modular architecture with components communicating via a signal-slot mechanism:

1. **AbsToolsLauncher**: User interface for configuring and starting analysis sessions
2. **Absorber**: Creates ion dictionaries from spectrum data
3. **Metal_Plot / MainWindow**: Main application window and controller
4. **PageManager**: Handles page (tab) creation and initialization
5. **Plotting**: Manages all visualization aspects
6. **EventHandler**: Processes user interactions (mouse/keyboard)
7. **EquivalentWidth**: Performs line measurements and calculations
8. **Utility Modules**: Support JSON handling, data extraction, batch processing

The signal-slot communication system ensures robust coordination between components while preventing segmentation faults that could occur with direct coupling.

### Architecture Diagram

The following diagram illustrates the relationship between all components of AbsTools:

```mermaid
flowchart TB
    subgraph "Entry Points"
        launcher[AbsTools Launcher GUI]
        api[Direct API Usage]
        loadfile[Load Saved File]
    end
    
    subgraph "Core Components"
        absorber[Absorber Class]
        metal_plot[Metal_Plot / MainWindow]
        
        subgraph "GUI Components"
            page_manager[Page Manager]
            save_page[Save Dialog]
            help_window[Help Window]
        end
        
        subgraph "Processing"
            plotting[Plotting Engine]
            event_handler[Event Handler]
            equivalent_width[Equivalent Width Calculator]
            signals[Signal-Slot Communication]
        end
        
        subgraph "Data Handling"
            json_utils[JSON Utilities]
            cleanup[Resource Cleanup]
        end
    end
    
    subgraph "Utilities"
        extract_measurements[Extract Measurements]
        batch_processor[Batch Processor]
        extract_data[Extract Raw Data]
    end
    
    subgraph "File I/O"
        save_json[Save JSON]
        save_pickle[Save Pickle]
        save_pdf[Save PDF]
        save_table[Save Data Table]
    end
    
    launcher --> absorber
    api --> absorber
    loadfile --> metal_plot
    
    absorber --> metal_plot
    
    metal_plot --> page_manager
    metal_plot --> save_page
    metal_plot --> help_window
    metal_plot <--> signals
    
    signals <--> plotting
    signals <--> event_handler
    signals <--> equivalent_width
    
    metal_plot --> plotting
    metal_plot --> event_handler
    metal_plot --> equivalent_width
    
    save_page --> save_json
    save_page --> save_pickle
    save_page --> save_pdf
    save_page --> save_table
    
    metal_plot --> cleanup
    json_utils <--> save_json
    
    save_json -.-> extract_measurements
    save_json -.-> batch_processor
    save_json -.-> extract_data
    save_pickle -.-> extract_measurements
    
    classDef core fill:#f9f,stroke:#333,stroke-width:2px;
    classDef io fill:#bbf,stroke:#333,stroke-width:1px;
    classDef utils fill:#bfb,stroke:#333,stroke-width:1px;
    classDef entry fill:#fbb,stroke:#333,stroke-width:1px;
    
    class absorber,metal_plot,signals core;
    class save_json,save_pickle,save_pdf,save_table io;
    class extract_measurements,batch_processor,extract_data utils;
    class launcher,api,loadfile entry;
```


### User Workflow Diagram

This diagram illustrates the typical user workflow when working with AbsTools:

![AbsTools User Workflow](https://github.com/rongmon/rbcodes/blob/master/docs/GUIs/AbsTools/images/abstools_workflow.svg)

## Troubleshooting

### Common Issues

1. **Missing Libraries**: Ensure all dependencies are properly installed. The error messages will indicate which libraries are missing.

2. **JSON Support Unavailable**: If you see "JSON utilities not found" messages, make sure the `json_utils.py` module is accessible within the rbcodes package.

3. **Spectrum Loading Errors**: Verify that your spectrum file is in a supported format (FITS, ASCII, HDF5). Check that the file contains flux, wavelength, and error arrays.

4. **Lines Outside Observed Range**: If you specify lines that fall outside the wavelength range of your spectrum, AbsTools will warn you during initialization.

5. **Launcher Script Path**: If the launcher can't find required modules, ensure the rbcodes package is in your Python path.

### Error Recovery

The updated AbsTools includes improved error handling:

- Signals communicate errors to a central handler
- Error messages are displayed in status bars and dialog boxes
- Critical errors are caught without causing application crashes
- Resource cleanup on exit prevents segmentation faults

## Acknowledgments

AbsTools was originally developed by Sean Clark and Rongmon Bordoloi [2021]
Major refactoring and update with enhanced features by Rongmon Bordoloi [2025]