# rbcodes/GUIs/launch_specgui.py
"""
Unified launcher for rb_spec GUI applications.

This script replaces both the old launch_specgui.py and launch_batch.py files.
It can launch either the single-spectrum GUI or the batch processing GUI
depending on the command-line arguments provided.

Usage:
    python launch_specgui.py [filename] [-t filetype]    # Single spectrum mode (default)
    python launch_specgui.py -b [config_file]            # Batch processing mode
    python launch_specgui.py -v                          # Show version
    python launch_specgui.py --help                      # Show help

Examples:
    python launch_specgui.py spectrum.fits                    # Load single spectrum
    python launch_specgui.py spectrum.fits -t linetools       # Load with specific filetype
    python launch_specgui.py analysis.json                    # Load saved analysis
    python launch_specgui.py -b                               # Launch batch processor
    python launch_specgui.py -b batch_config.json             # Load batch configuration
    python launch_specgui.py -b batch_items.csv               # Import batch items from CSV
"""

import sys
import os
import argparse
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QTimer

# Version information
__version__ = "1.0.3"
__author__ = "Rongmon Bordoloi"
__description__ = "rb_spec GUI - Absorption Line Analysis Tool"


def show_version():
    """Show version information."""
    print(f"{__description__}")
    print(f"Version: {__version__}")
    print(f"Author: {__author__}")
    print()
    print("Components:")
    print("  - Single Spectrum GUI: Interactive analysis of individual spectra")
    print("  - Batch Processing GUI: Automated analysis of multiple spectra")


def launch_single_gui(args):
    """Launch the single-spectrum GUI."""
    try:
        from rbcodes.GUIs.specgui.main import RbSpecGUI
    except ImportError:
        print("Error: Could not import single-spectrum GUI modules.")
        print("Make sure rbcodes is properly installed.")
        return 1
    
    # Create application and main window
    app = QApplication(sys.argv)
    app.setApplicationName("rb_spec GUI")
    app.setApplicationVersion(__version__)
    app.setOrganizationName("rbcodes")
    
    window = RbSpecGUI()
    window.show()
    
    # Load file if specified
    if args.filename and os.path.exists(args.filename):
        def load_file():
            filetype = args.filetype if args.filetype else None
            
            # Check if JSON file for special handling
            if args.filename.lower().endswith('.json'):
                success, has_redshift, has_transition = window.controller.load_from_json(args.filename)
                window.input_panel.spectrum_loaded.emit(success, has_redshift, has_transition)
            else:
                success, _, _ = window.controller.load_from_file(args.filename, filetype)
                window.input_panel.spectrum_loaded.emit(success, False, False)
            
            # Update the file path display
            window.input_panel.file_path.setText(args.filename)
        
        # Schedule file loading after GUI is displayed
        QTimer.singleShot(100, load_file)
    elif args.filename:
        print(f"Warning: File '{args.filename}' does not exist.")
    
    return app.exec_()


def launch_batch_gui(args):
    """Launch the batch processing GUI."""
    try:
        from rbcodes.GUIs.specgui.batch.batch_main import BatchSpecGUI
    except ImportError:
        print("Error: Could not import batch processing GUI modules.")
        print("Make sure rbcodes is properly installed.")
        return 1
    
    # Create application and main window
    app = QApplication(sys.argv)
    app.setApplicationName("rb_spec Batch Processing")
    app.setApplicationVersion(__version__)
    app.setOrganizationName("rbcodes")
    
    window = BatchSpecGUI()
    window.show()
    
    # Handle config/CSV file if specified
    if args.config and os.path.exists(args.config):
        def load_file():
            # Check file extension to determine import method
            if args.config.lower().endswith('.csv'):
                # Import from CSV
                print(f"Importing batch items from CSV: {args.config}")
                success, message = window.controller.master_table.import_template_csv(args.config)
                if success:
                    window.config_panel.refresh_table()
                    window.config_panel.configuration_changed.emit()
                    window.update_status(f"Imported batch items from {args.config}")
                    print(f"Success: {message}")
                else:
                    window.update_status(f"Failed to import from {args.config}: {message}")
                    print(f"Error: {message}")
            else:
                # Load as JSON configuration
                print(f"Loading configuration: {args.config}")
                success = window.controller.load_batch_configuration(args.config)
                if success:
                    window.config_panel.refresh_table()
                    window.current_config_file = args.config
                    window.has_unsaved_changes = False
                    window.update_window_title()
                    window.update_tab_states()
                    window.update_status(f"Loaded configuration from {args.config}")
                else:
                    window.update_status(f"Failed to load configuration from {args.config}")
        
        # Schedule file loading after GUI is fully initialized
        QTimer.singleShot(100, load_file)
    elif args.config:
        print(f"Warning: File '{args.config}' does not exist.")
    
    return app.exec_()


def main():
    """Main entry point for the unified launcher."""
    # Create argument parser
    parser = argparse.ArgumentParser(
        description=f'{__description__} v{__version__}',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s spectrum.fits                    Load single spectrum
  %(prog)s spectrum.fits -t linetools       Load with specific filetype  
  %(prog)s analysis.json                    Load saved analysis
  %(prog)s -b                               Launch batch processor
  %(prog)s -b batch_config.json             Load batch configuration
  %(prog)s -b batch_items.csv               Import batch items from CSV
  %(prog)s -v                               Show version information

For more information, visit: https://github.com/rongmon/rbcodes
        """
    )
    
    # Version argument
    parser.add_argument('-v', '--version', action='store_true',
                       help='Show version information and exit')
    
    # Batch mode argument
    parser.add_argument('-b', '--batch', action='store_true',
                       help='Launch batch processing GUI')
    
    # Positional argument (filename for single mode, config/csv for batch mode)
    parser.add_argument('filename', nargs='?', 
                       help='Spectrum file (single mode) or config/CSV file (batch mode)')
    
    # Single mode specific arguments
    parser.add_argument('-t', '--filetype', 
                       help='File type for single mode (e.g., fits, ascii, linetools)')
    
    # Hidden argument for batch config (for backward compatibility)
    parser.add_argument('--config', help=argparse.SUPPRESS)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Handle version request
    if args.version:
        show_version()
        return 0
    
    # Determine mode and validate arguments
    if args.batch:
        # Batch mode
        if args.filetype:
            parser.error("--filetype is not used in batch mode")
        
        # In batch mode, filename becomes config file
        args.config = args.filename
        args.filename = None
        
        print(f"Launching rb_spec Batch Processing GUI v{__version__}")
        if args.config:
            if args.config.lower().endswith('.csv'):
                print(f"Will import batch items from CSV: {args.config}")
            else:
                print(f"Loading configuration: {args.config}")
        
        return launch_batch_gui(args)
    
    else:
        # Single spectrum mode (default)
        print(f"Launching rb_spec GUI v{__version__}")
        if args.filename:
            print(f"Loading spectrum: {args.filename}")
            if args.filetype:
                print(f"File type: {args.filetype}")
        
        return launch_single_gui(args)


def launch():
    """Entry point for setup.py console scripts."""
    return main()


if __name__ == "__main__":
    sys.exit(main())
