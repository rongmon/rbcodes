# launch_batch.py
"""
Launch script for rb_spec batch processing GUI.

This script can be used to start the batch processing application
from the command line or imported to start it programmatically.

Usage:
    python launch_batch.py [config_file]
    
Examples:
    python launch_batch.py
    python launch_batch.py my_batch_config.json
"""

import sys
import os

def launch():
    """Launch the batch processing GUI."""
    try:
        # Try to import from the installed package
        from rbcodes.GUIs.specgui.batch.batch_main import main
        return main()
    except ImportError:
        try:
            # Try to import from local development version
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from batch_main import main
            return main()
        except ImportError as e:
            print(f"Error: Could not import batch processing modules: {e}")
            print("Make sure rbcodes is properly installed or you're running from the correct directory.")
            return 1

if __name__ == "__main__":
    sys.exit(launch())