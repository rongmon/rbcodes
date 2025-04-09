#!/usr/bin/env python
"""
AbsTools Launcher - Entry point for the Absorption Line Analysis Toolbox.

This script provides a convenient entry point for launching the AbsTools
wrapper GUI, which in turn facilitates using the absorption line analysis tool.

Usage:
    python abstools_launcher.py
"""

import sys
import os
import traceback

def run_gui():
    """Run the AbsTools Launcher GUI."""
    try:
        from PyQt5.QtWidgets import QApplication
        from PyQt5.QtGui import QPalette, QColor
        from PyQt5.QtCore import Qt
        
        # Create the application
        app = QApplication(sys.argv)
        
        # Set the application style to Fusion (modern look)
        app.setStyle("Fusion")
        
        # Import the launcher class
        try:
            # Try importing directly from the package
            from rbcodes.gui.abstools.AbsToolsLauncher import AbsToolsLauncher
        except ImportError:
            # If that fails, try to import assuming this script is in the same directory
            current_dir = os.path.dirname(os.path.abspath(__file__))
            sys.path.insert(0, os.path.dirname(os.path.dirname(current_dir)))
            from rbcodes.gui.abstools.AbsToolsLauncher import AbsToolsLauncher
        
        # Create and show the main window
        window = AbsToolsLauncher()
        
        # Start the event loop
        sys.exit(app.exec_())
        
    except ImportError as e:
        print(f"Error: Missing required dependency: {str(e)}")
        print("Please make sure PyQt5 and all dependencies are installed.")
        sys.exit(1)
    except Exception as e:
        print(f"Error launching application: {str(e)}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    run_gui()