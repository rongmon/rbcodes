# rbcodes/GUIs/multispecviewer/utils.py

def read_line_options():
    """
    Read available line list options from the configuration file.
    
    Returns:
        list: List of line list option names. Returns default options if the
              configuration file cannot be read.
    """
    try:
        # Use pkg_resources to find the configuration file
        from pkg_resources import resource_filename
        config_file = resource_filename('rbcodes.GUIs.multispecviewer', 'line_options.conf')
        
        # Default options in case file is not found
        default_options = ["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong", "AGN"]
        
        # Try to read the file
        with open(config_file, 'r') as f:
            # Read lines and strip whitespace
            options = [line.strip() for line in f.readlines()]
            
            # Filter out empty lines and comments (lines starting with #)
            options = [opt for opt in options if opt and not opt.startswith('#')]
            
            # If no options were read, use defaults
            if not options:
                print("Warning: No line options found in configuration file, using defaults.")
                return default_options
                
            return options
            
    except Exception as e:
        # If there's any error reading the file, use defaults
        print(f"Error reading line options configuration: {str(e)}")
        print("Using default line options instead.")
        return ["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong", "AGN"]


# Add to rbcodes/GUIs/multispecviewer/utils.py

# Add this function to rbcodes/GUIs/multispecviewer/utils.py

def show_help_dialog(parent=None):
    """
    Display a help dialog with keyboard shortcuts and usage information.
    
    Parameters:
    -----------
    parent : QWidget, optional
        Parent widget for the dialog
    """
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QTextEdit, 
                                QPushButton, QHBoxLayout, QLabel)
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QFont
    import os
    
    # Help text
    help_text = """MULTISPECVIEWER HELP

OVERVIEW:
---------
MultispecViewer is a tool for visualizing and analyzing multiple spectroscopic datasets simultaneously.

INTERFACE:
----------
- Left Widget: Absorber manager for adding, tracking, and plotting absorber systems
- Main Canvas: Main spectral display area with keyboard shortcuts (click here to enable keyboard events)
- Redshift Widget: Shows active redshift for line identification; use Catalog button to add to manager

BASIC NAVIGATION:
----------------
- r: Reset view to original state and clear all lines
- R: Keep spectra but remove all lines/text from canvas
- x: Set left limit (xmin) to cursor position
- X: Set right limit (xmax) to cursor position
- [: Shift view left
- ]: Shift view right
- o: Zoom out (x-axis)
- Y: Manually input y-limits range

DISPLAY CONTROLS:
----------------
- t: Set maximum y-limit to cursor position
- b: Set minimum y-limit to cursor position
- S: Increase smoothing (convolution kernel size)
- U: Decrease smoothing (convolution kernel size)

ANALYSIS TOOLS:
--------------
- v: Open vStack GUI for identifying transitions (press S to exit)
- V: Same as 'v' but allows manual velocity axis selection
- Right-Click: Shows a menu of possible spectral lines at cursor position for identification

QUICK LINE IDENTIFICATION:
-------------------------
- M: Mark MgII doublet (2796, 2803)
- C: Mark CIV doublet (1548, 1550)
- F: Mark FeII multiplet (2600, 2586, 2382)
- 6: Mark OVI doublet (1031, 1037)
- 4: Mark SiIV doublet (1393, 1402)
- 8: Mark NeVIII doublet (778, 770)
- 2: Mark Lyb/Lya
- 1: Mark Lya/Lyb

VSTACK CONTROLS:
---------------
- >: Next page
- <: Previous page
- w: Toggle transition flag between:
   - Detection
   - Blended-detection
   - Low-Confidence Detection
   - Non-Detection
- S: Save transition list and return to main view

FILE MANAGEMENT:
---------------
- Save: Save data in different formats:
   - JSON: Combined data with absorbers, line lists, and metadata
   - CSV: Absorber data only (.csv)
   - TXT: Line list data only (.txt)
- Load: Load previously saved data:
   - Select a JSON file for complete data
   - Select a CSV or TXT file for individual components
   - Choose to append or overwrite existing data

GETTING HELP:
------------
- H: Show this help window
- For complete documentation, see the multispec.md file
"""
    
    # Create dialog
    dialog = QDialog(parent)
    dialog.setWindowTitle("MultispecViewer Help")
    dialog.setMinimumSize(700, 500)
    
    # Create layout
    layout = QVBoxLayout(dialog)
    
    # Create text edit for help content
    text_edit = QTextEdit()
    text_edit.setReadOnly(True)
    text_edit.setPlainText(help_text)
    
    # Set monospace font for better formatting
    font = QFont("Courier New", 10)
    text_edit.setFont(font)
    
    # Add to layout
    layout.addWidget(text_edit)
    
    # Add a label for README link
    readme_label = QLabel("For more detailed documentation:")
    layout.addWidget(readme_label)
    
    # Add buttons for documentation access
    button_layout = QHBoxLayout()

    # Add GitHub button - always available
    github_link = "https://github.com/rongmon/rbcodes/blob/master/docs/GUIs/multispec/multispec.md"
    github_button = QPushButton("View README on GitHub")

    def open_github():
        import webbrowser
        webbrowser.open(github_link)

    github_button.clicked.connect(open_github)
    button_layout.addWidget(github_button)

    # Try to find local README file as well
    try:
        import os
        from pkg_resources import resource_filename
        
        # Check possible local locations
        potential_paths = [
            resource_filename('rbcodes', 'docs/GUIs/multispec/multispec.md'),
            resource_filename('rbcodes.GUIs.multispecviewer', 'multispec.md'),
            resource_filename('rbcodes', 'GUIs/multispecviewer/multispec.md')
        ]
        
        # Find first existing path
        readme_path = None
        for path in potential_paths:
            if os.path.exists(path):
                readme_path = path
                break
        
        # If found locally, add button to open it
        if readme_path:
            local_button = QPushButton("Open Local README")
            
            def open_local_readme():
                import sys
                import subprocess
                
                # Different methods to open file based on OS
                if sys.platform.startswith('darwin'):  # macOS
                    subprocess.call(('open', readme_path))
                elif sys.platform.startswith('win'):   # Windows
                    os.startfile(readme_path)
                else:  # Linux and other Unix-like
                    subprocess.call(('xdg-open', readme_path))
            
            local_button.clicked.connect(open_local_readme)
            button_layout.addWidget(local_button)
            
    except Exception as e:
        print(f"Could not locate local README.md: {str(e)}")

    # Add close button
    close_button = QPushButton("Close")
    close_button.clicked.connect(dialog.accept)
    button_layout.addWidget(close_button)

    # Add button layout to main layout
    layout.addLayout(button_layout)
    
    # Apply dark theme if parent has it
    if parent and hasattr(parent, 'palette'):
        dialog.setPalette(parent.palette())
        
        # Additional styling for dark theme
        dialog.setStyleSheet("""
            QDialog {
                background-color: #353535;
                color: #F2F2F7;
            }
            QTextEdit {
                background-color: #252525;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 6px;
                padding: 8px;
            }
            QLabel {
                color: #F2F2F7;
            }
            QPushButton {
                background-color: #474747;
                color: #F2F2F2;
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050;
            }
            QPushButton:pressed {
                background-color: #2A2A2A;
            }
        """)
    
    # Show dialog
    dialog.exec_()