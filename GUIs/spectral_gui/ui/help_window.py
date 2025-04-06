#!/usr/bin/env python
"""
Help Window - Documentation for the spectrum viewer application.
"""
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QScrollArea, QTextEdit


HELP_TEXT = """
# Spectrum Viewer Help

## Main Interface Components

The spectrum viewer consists of three main areas:
- **Left Panel**: Absorber manager for organizing identified absorption systems
- **Center**: Spectrum plot canvas
- **Bottom Bar**: Active values and controls

## Keyboard Shortcuts

### Navigation
- `r`: Reset/clear the axes and replot the spectrum
- `R`: Keep the spectrum active but remove all lines/text from the canvas
- `x`: Set left x limit (xmin)
- `X`: Set right x limit (xmax)
- `o`: Zoom out of x range
- `]`: Shift canvas to the right
- `[`: Shift canvas to the left
- `Y`: Manually input y-axis limits

### Display Adjustments
- `t`: Restrict the maximum y value to mouse position
- `b`: Restrict the minimum y value to mouse position
- `S`: Smooth the spectrum (increase smoothing)
- `U`: Unsmooth spectrum (decrease smoothing)

### Measurements
- `E`: Two `E` keystrokes will compute rest frame equivalent width between marked points
- `G`: Three `G` keystrokes to fit a Gaussian profile

### Line Identification
- `j` or Right-click: Open line selection dialog at cursor position
- `v`: Open stack view for transition identification
- `V`: Open stack view with custom velocity limits

### Doublet/Multiplet Marking
- `M`: Mark MgII doublet (2796, 2803)
- `C`: Mark CIV doublet (1548, 1550)
- `F`: Mark FeII triplet (2600, 2586, 2382)
- `6`: Mark OVI doublet (1031, 1037)
- `4`: Mark SiIV doublet (1393, 1402)
- `8`: Mark NeVIII doublet (778, 770)
- `2`: Mark Lyb/Lya pair
- `1`: Mark Lya/Lyb pair

### Other
- `H`: Show this help window
- `q`: Quit application

## Absorber Manager

The absorber manager in the left panel allows you to:
- View all identified absorbers
- Add new absorbers by entering a redshift
- Select line lists for each absorber
- Choose display colors
- Plot, hide, or remove absorbers

### Working with Absorbers

1. **Adding Absorbers**: Enter a redshift in the last empty row
2. **Plotting Lines**: Click the "Plot" button to display absorption lines
3. **Hiding Lines**: Click "Hide" to temporarily hide lines without removing them
4. **Removing Absorbers**: Click "Remove" to delete an absorber

## Spectrum Measurements

### Equivalent Width
To measure equivalent width:
1. Press `E` at one edge of the absorption feature
2. Press `E` at the other edge
3. The EW will be calculated and displayed

### Gaussian Fitting
To fit a Gaussian profile:
1. Press `G` at the left continuum point
2. Press `G` at the center of the feature
3. Press `G` at the right continuum point
4. The Gaussian parameters will be displayed

## Stack View

The stack view (`v` key) helps identify multiple transitions at the same redshift:
- Navigate pages with `<` and `>` keys
- Toggle detection status with `w` key
- Press `S` to save and return to main view

## Saving and Loading

- Click "Save" to save the current list of absorbers
- Click "Load" to load a previously saved catalog
- Catalogs are saved as CSV files
- Line lists are saved as separate text files

## Identified Lines

Click "Plot Identified Lines" to show all detected lines across all absorbers.
This provides a comprehensive view of all identified features.
"""


class HelpWindow(QWidget):
    """
    A window displaying help documentation for the application.
    
    Parameters
    ----------
    parent : QWidget, optional
        Parent widget, by default None
    """
    
    def __init__(self, parent=None):
        """Initialize the help window."""
        super(HelpWindow, self).__init__(parent)
        
        # Set window properties
        self.setWindowTitle("Spectrum Viewer Help")
        self.resize(800, 600)
        
        # Create layout
        layout = QVBoxLayout(self)
        
        # Create text editor for displaying help
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setMarkdown(HELP_TEXT)
        
        # Add to layout
        layout.addWidget(text_edit)
        
        # Set layout
        self.setLayout(layout)