# MultispecViewer

[Back to Main Page](../../main_readme.md)

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
- Use velocity plots (vStack) for detailed line analysis
- Save and load line identifications and absorber systems

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
- Action buttons for loading, saving, and displaying identified lines

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

### Velocity Plot Analysis with vStack

The vStack feature allows you to analyze absorption lines in velocity space:

| Key | Action |
|-----|--------|
| `v` | Launch vStack with default velocity window (±1000 km/s) |
| `V` | Launch vStack with custom velocity window (prompts for limits) |
| `>` | Navigate to next page of lines in vStack |
| `<` | Navigate to previous page of lines in vStack |
| `w` | Toggle line detection status in vStack (cycles through Non-Detection, Detection, Blended-detection, Low-Confidence) |
| `S` | Save line identifications and return to main display |
| `Y` | Manually set y-limits for the current vStack panel |

### Action Buttons

The new version includes several action buttons:

| Button | Action |
|--------|--------|
| Load | Load previously saved line identifications and absorber systems |
| Save | Save current line identifications and absorber systems |
| Show | Toggle visibility of all identified lines |
| List | Display a table of all identified lines |

## Complete Workflow Examples

### Identifying and Cataloging an Absorber

Here's a typical workflow for identifying and cataloging an absorber:

1. **Load spectral data**:
   - Click "Select FITS Files" and choose your FITS files
   - The spectra will be displayed in the main window

2. **Initial redshift identification**:
   - Right-click on a spectral feature
   - Select a possible line from the dialog to calculate the implied redshift
   - The redshift will be applied to the redshift input widget

   [Insert screenshot of right-click identification here]

3. **Apply redshift and check for additional lines**:
   - Verify the selected line list in the redshift widget
   - Choose a color for the lines
   - Click "Submit" to mark all lines at that redshift
   - Visually inspect the marked lines to confirm the identification

   [Insert screenshot of applied redshift with marked lines here]

4. **Detailed analysis with vStack**:
   - Position your cursor on the spectrum and press `v` for default velocity windows or `V` to set custom velocity limits
   - Navigate through pages of lines using `<` and `>` keys
   - Mark each line as Detection, Non-Detection, Blended, or Low-Confidence using the `w` key
   - Press `S` to save all marked lines and return to the main display

   [Insert screenshot of vStack interface here]

5. **Catalog the absorber system**:
   - Click "Catalog" to add the current redshift system to the absorber manager
   - The absorber will appear in the left panel with its redshift, line list, and color
   - Use the checkbox to toggle visibility of the system

   [Insert screenshot of absorber manager with cataloged system here]

6. **Save your work**:
   - Click the "Save" button to save your line identifications and absorber systems
   - Choose a filename to save the data
   - Two files will be created: a .txt file with line identifications and a .csv file with absorber systems

### Identifying Multiple Systems

To identify multiple absorber systems:

1. Follow the workflow above to identify your first system
2. After cataloging the first system, clear the current redshift markings by pressing `R`
3. Identify a feature from a different system using right-click or direct redshift input
4. Apply a different color to distinguish the second system from the first
5. Confirm the identification and add it to the catalog
6. Repeat for additional systems

[Insert screenshot showing multiple absorber systems in different colors here]

### Managing Line Identifications

After identifying multiple systems:

1. Click the "List" button to view all identified lines in a tabular format
2. Use the "Show" button to toggle the visibility of all identified lines
3. Use the absorber manager checkboxes to control which systems are displayed
4. Click the "Save" button to save all identifications

[Insert screenshot of line list dialog here]

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
- Action buttons for common tasks (Load, Save, Show, List)
- Better visualization of multiple absorber systems

[Insert new version screenshot here]

## Tips and Best Practices

### Keyboard Navigation

- Use the keyboard navigation keys (`x`, `X`, `t`, `b`) for precise control of the display
- Press `r` to reset to the original view if you get lost
- Use `R` to clear line identifications without changing the view

### Redshift Identification

- Start by identifying strong, unambiguous features like MgII, CIV, or Lyman-alpha
- Use right-click to see all possible line identifications at a given wavelength
- After applying a redshift, check for other expected lines at the same redshift
- If uncertain about a feature, use vStack (`v` or `V`) for detailed analysis

### vStack Analysis

- Use `V` to specify custom velocity limits for better visualizing line profiles
- Mark line detections systematically using the `w` key
- Don't forget to press `S` to save your work before returning to the main display
- Use `Y` to adjust y-limits for better visualization of line profiles

### Multiple Absorber Systems

- Use different colors for different absorber systems
- Use the absorber manager to toggle visibility of systems when the display gets crowded
- Save your work frequently using the "Save" button
- For complex spectra with many systems, focus on one redshift at a time

### File Management

- Use consistent file naming conventions for your saved data
- The program creates two files when saving:
  - A .txt file containing line identifications
  - A .csv file with the same base name containing absorber systems
- Both files are needed for a complete restoration of your work

## Development

MultispecViewer is part of the rbcodes package for spectroscopic analysis. It utilizes:
- PyQt5 for the graphical interface
- Matplotlib for plotting
- linetools for spectral data handling
- astropy for astronomical calculations

## Troubleshooting

### Common Issues

- **No files selected**: Ensure you've clicked "Select FITS Files" and chosen valid FITS files
- **No error spectrum**: If no error spectrum is found, the program will assume 5% of flux values
- **No lines visible after applying redshift**: Check that the selected line list contains lines within your wavelength range
- **vStack not launching**: Ensure you have a redshift and line list selected before pressing `v` or `V`

### Error Messages

Watch the message box at the bottom of the screen for helpful error messages and status updates.