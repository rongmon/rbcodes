# launch_specgui.py or launch.py
import sys
import os
import argparse
from PyQt5.QtWidgets import QApplication
from rbcodes.GUIs.specgui.main import RbSpecGUI

def launch():
    """
    Launch the rb_spec GUI application with optional file loading.
    
    Command-line options:
    - Positional argument: filename to load
    - -t/--filetype: Specify file type (e.g., 'fits', 'ascii', 'linetools')
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Launch rb_spec GUI')
    parser.add_argument('filename', nargs='?', help='Spectrum file to load on startup')
    parser.add_argument('-t', '--filetype', help='File type (e.g., fits, ascii, linetools)')
    
    args = parser.parse_args()
    
    # Create application and main window
    app = QApplication(sys.argv)
    window = RbSpecGUI()
    window.show()
    
    # Load file if specified
    if args.filename and os.path.exists(args.filename):
        # Use Qt's timer to load the file after the UI is fully initialized
        from PyQt5.QtCore import QTimer
        
        def load_file():
            filetype = args.filetype if args.filetype else None
            window.controller.load_from_file(args.filename, filetype)
            window.input_panel.file_path.setText(args.filename)
            
            # If the file is JSON, properly update the UI state through the input panel
            if args.filename.lower().endswith('.json'):
                success, has_redshift, has_transition = window.controller.load_from_json(args.filename)
                window.input_panel.spectrum_loaded.emit(success, has_redshift, has_transition)
            else:
                success, _, _ = window.controller.load_from_file(args.filename, filetype)
                window.input_panel.spectrum_loaded.emit(success, False, False)
        
        # Schedule file loading after GUI is displayed
        QTimer.singleShot(100, load_file)
    
    return app.exec_()

if __name__ == "__main__":
    sys.exit(launch())