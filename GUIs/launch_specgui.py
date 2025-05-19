# launch_specgui.py
import sys
import os
import argparse
from PyQt5.QtWidgets import QApplication
from rbcodes.GUIs.specgui.main import RbSpecGUI
from rbcodes.GUIs.specgui.batch.batch_main import BatchSpecGUI

def launch():
    """
    Launch the rb_spec GUI application with optional file loading.
    
    Command-line options:
    - Positional argument: filename to load
    - -t/--filetype: Specify file type (e.g., 'fits', 'ascii', 'linetools')
    - -b/--batch: Launch in batch processing mode
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Launch rb_spec GUI')
    parser.add_argument('filename', nargs='?', help='Spectrum file to load on startup')
    parser.add_argument('-t', '--filetype', help='File type (e.g., fits, ascii, linetools)')
    parser.add_argument('-b', '--batch', action='store_true', help='Launch in batch processing mode')
    
    args = parser.parse_args()
    
    # Create application
    app = QApplication(sys.argv)
    
    # Launch appropriate GUI based on mode
    if args.batch:
        # Launch batch processing GUI
        window = BatchSpecGUI()
        window.show()
        
        # If a file was specified and it's a batch configuration file
        if args.filename and os.path.exists(args.filename) and args.filename.lower().endswith('.json'):
            # Use Qt's timer to load the configuration after the UI is fully initialized
            from PyQt5.QtCore import QTimer
            
            def load_config():
                # Load batch configuration
                success = window.controller.load_batch_configuration(args.filename)
                if success:
                    # Update UI with loaded configuration
                    window.config_panel.update_ui_with_config(window.controller.batch_items)
            
            # Schedule configuration loading after GUI is displayed
            QTimer.singleShot(100, load_config)
    else:
        # Launch regular GUI
        window = RbSpecGUI()
        window.show()
        
        # Load file if specified
        if args.filename and os.path.exists(args.filename):
            # Use Qt's timer to load the file after the UI is fully initialized
            from PyQt5.QtCore import QTimer
            
            
            def load_file():
                # If filetype is not specified and it's a FITS file, use 'linetools' as default
                if args.filetype is None and args.filename.lower().endswith('.fits'):
                    filetype = 'linetools'
                else:
                    filetype = args.filetype
            
                # If the file is JSON, properly update the UI state through the input panel
                if args.filename.lower().endswith('.json'):
                    success, has_redshift, has_transition = window.controller.load_from_json(args.filename)
                    window.input_panel.spectrum_loaded.emit(success, has_redshift, has_transition)
                else:
                    success, _, _ = window.controller.load_from_file(args.filename, filetype)
                    window.input_panel.spectrum_loaded.emit(success, False, False)
                
                # Update the file path text field
                window.input_panel.file_path.setText(args.filename)


            # Schedule file loading after GUI is displayed
            QTimer.singleShot(100, load_file)
    
    return app.exec_()

if __name__ == "__main__":
    sys.exit(launch())