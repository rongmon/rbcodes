#!/usr/bin/env python
"""
Simplified AbsTools Launcher - Entry point for the Absorption Line Analysis Toolbox.

This script provides a convenient entry point for launching the AbsTools
wrapper GUI, with improved process isolation.

Usage:
    python simplified_abstools_launcher.py
"""

import sys
import os
import traceback
import subprocess

def run_as_separate_process():
    """
    Launch AbsTools in a completely separate process.
    This function should be called from absorption_manager.
    """
    try:
        # Get the path to this script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Create a detached subprocess
        if sys.platform == 'win32':
            # Windows-specific approach
            CREATE_NEW_PROCESS_GROUP = 0x00000200
            DETACHED_PROCESS = 0x00000008
            process = subprocess.Popen(
                [sys.executable, __file__, "--new-process"],
                creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
                close_fds=True
            )
        else:
            # Unix-specific approach
            process = subprocess.Popen(
                [sys.executable, __file__, "--new-process"],
                start_new_session=True,
                close_fds=True
            )
        
        print(f"AbsTools launched as process ID: {process.pid}")
        return True
    
    except Exception as e:
        print(f"Error launching AbsTools: {e}")
        traceback.print_exc()
        return False

def run_gui():
    """Run the AbsTools Launcher GUI."""
    try:
        from PyQt5.QtWidgets import QApplication
        
        # Create the application
        app = QApplication(sys.argv)
        app.setStyle("Fusion")
        
        # Import the launcher class
        try:
            from AbsToolsLauncher import AbsToolsLauncher
        except ImportError:
            # Try to import from full package path
            try:
                from GUIs.abstools.AbsToolsLauncher import AbsToolsLauncher
            except ImportError:
                # If that fails, try to import assuming this script is in the same directory
                current_dir = os.path.dirname(os.path.abspath(__file__))
                sys.path.insert(0, os.path.dirname(os.path.dirname(current_dir)))
                from GUIs.abstools.AbsToolsLauncher import AbsToolsLauncher
        
        # Create and show the main window
        window = AbsToolsLauncher()
        window.show()
        
        # Start the event loop
        sys.exit(app.exec_())
        
    except Exception as e:
        print(f"Error launching application: {str(e)}")
        with open("abstool_launcher_error.log", "w") as f:
            f.write(f"Error: {str(e)}\n")
            traceback.print_exc(file=f)
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    # Check if we're being called with the --new-process flag
    if len(sys.argv) > 1 and sys.argv[1] == "--new-process":
        # Run the GUI directly
        run_gui()
    else:
        # Check if we're being imported or run directly
        is_imported = __name__ != "__main__"
        
        if is_imported:
            # If imported, provide the run_as_separate_process function
            # The calling code should use this function
            pass
        else:
            # If run directly (without --new-process), launch as separate process
            run_as_separate_process()