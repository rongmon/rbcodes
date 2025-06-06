#!/usr/bin/env python
import sys
import os
import argparse
import importlib
from rbcodes.GUIs.multispecviewer.io_manager import IOManager

# Version information
__version__ = "1.2.0"

def display_examples():
    """Display detailed usage examples for the rb_multispec tool."""
    print("MultispecViewer Usage Examples")
    print("=============================")
    print("")
    print("Basic Usage:")
    print("  rb_multispec.py                         # Run with GUI file selection")
    print("  rb_multispec.py file1.fits              # Load a single FITS file")
    print("  rb_multispec.py file1.fits file2.fits   # Load multiple specific FITS files")
    print("")
    print("Using Wildcards:")
    print("  rb_multispec.py *.fits                  # Load all FITS files in current directory")
    print("  rb_multispec.py data/*.fits             # Load all FITS files in 'data' directory")
    print("")
    print("Using Full Paths:")
    print("  rb_multispec.py /path/to/data/file.fits             # Load file with absolute path")
    print("  rb_multispec.py path1/file1.fits path2/file2.fits   # Mix files from different paths")
    print("  rb_multispec.py /path/to/data/*.fits                # Load all FITS with pattern")
    print("")
    print("Classic Mode:")
    print("  rb_multispec.py -c                      # Run classic version with GUI file selection")
    print("  rb_multispec.py -c file.fits            # Run classic version with specified file")
    print("  rb_multispec.py --classic *.fits        # Run classic version with all FITS files")
    print("")
    print("Programmatic Usage (from Python):")
    print("  from rbcodes.GUIs.multispecviewer import rb_multispec")
    print("  from rbcodes.utils.rb_spectrum import rb_spectrum")
    print("")
    print("  # Single spectrum from arrays")
    print("  spec = rb_spectrum.from_tuple((wave, flux, error))")
    print("  window = rb_multispec.from_data(spec)")
    print("")
    print("  # Multiple spectra from 2D arrays")
    print("  spec_list = rb_spectrum.from_arrays(wave_2d, flux_2d, error_2d)")
    print("  combined = rb_spectrum.append(spec_list)")
    print("  window = rb_multispec.from_data(combined)")
    print("")
    print("Tips:")
    print("  - Use quotes for paths with spaces: rb_multispec.py \"My Data/file.fits\"")
    print("  - You can mix options and filenames: rb_multispec.py -c /path/to/file.fits")
    print("  - For batch processing, use shell expansion: rb_multispec.py {A,B,C}/spec.fits")
    return 0

def main():
    """
    Wrapper script for rbcodes multispecviewer tools.
    Allows users to choose between the new version and the classic version.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="rb_multispec - A tool for visualizing spectral data",
        epilog="By default, the script runs the new multispecviewer."
    )
    
    # Add command line arguments
    parser.add_argument(
        "-c", "--classic", 
        action="store_true", 
        help="Run the classic version of multispecviewer"
    )
    
    parser.add_argument(
        "-v", "--version", 
        action="store_true", 
        help="Display version information"
    )
    
    parser.add_argument(
        "-e", "--examples",
        action="store_true",
        help="Display detailed usage examples"
    )
    
    # Add positional argument for filenames
    parser.add_argument(
        "filenames", 
        nargs="*",  # Zero or more arguments
        help="FITS files to load automatically on startup"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check for version flag
    if args.version:
        print(f"rb_multispec - MultispecViewer Wrapper v{__version__}")
        print("rbcodes Spectral Analysis Tools")
        return 0
    
    # Check for examples flag
    if args.examples:
        return display_examples()
    
    # Filter out our custom arguments for both versions, keeping filenames
    filtered_argv = [sys.argv[0]]  # Keep the script name
    
    # Attempt to import and run the appropriate module
    try:
        if args.classic:
            print("Starting classic MultispecViewer...")
            # Import the classic version
            from rbcodes.GUIs.multispecviewer.classic import multispec_classic as multispec
        else:
            print("Starting new MultispecViewer...")
            # Import the new version
            from rbcodes.GUIs.multispecviewer import multispec
        
        # Check if filenames were provided
        filenames = args.filenames
        if filenames:
            print(f"Loading {len(filenames)} FITS file(s)...")
            for filename in filenames:
                if not os.path.exists(filename):
                    print(f"Warning: File not found: {filename}")
        
        # Initialize QApplication
        app = multispec.QApplication(filtered_argv)
        
        # Create MainWindow instance
        window = multispec.MainWindow()
        
        # Show the window
        window.show()
        
        # If filenames were provided, load them after the window is shown
        if filenames:
            valid_paths = [path for path in filenames if os.path.exists(path)]
            if valid_paths:
                # Use a short delay to ensure UI is fully initialized before loading files
                from PyQt5.QtCore import QTimer
                QTimer.singleShot(100, lambda: window.load_fits_files(valid_paths))
        
        # Start the event loop
        sys.exit(app.exec_())
        
    except ImportError as e:
        print(f"Error importing module: {e}")
        print("Make sure rbcodes package is installed and in your Python path.")
        return 1
    except Exception as e:
        print(f"Error running multispecviewer: {e}")
        print(f"Error details: {str(e)}")
        return 1
    
    return 0


def from_data(spectrum_data, show=True):
    """
    Launch MultispecViewer GUI with spectral data from arrays or spectrum objects.
    """
    try:
        # Import required modules
        from PyQt5.QtWidgets import QApplication
        from rbcodes.GUIs.multispecviewer import multispec
        
        # Handle QApplication - reuse existing or create new
        app = QApplication.instance()
        if app is None:
            # No existing QApplication, create one
            app = QApplication(sys.argv)
            need_exec = True
        else:
            # Use existing QApplication
            need_exec = False
            
        # Normalize input to a list of spectrum objects
        spectra_list = _prepare_spectra_list(spectrum_data)  # Fix the typo here
        
        if not spectra_list:
            raise ValueError("No valid spectra provided")
            
        # Create MainWindow instance
        window = multispec.MainWindow()
        
        # Load spectra directly into the window
        window.spectra = spectra_list
        window.canvas.plot_spectra(spectra_list)
        
        # Update file label to show data source
        if hasattr(window, 'file_label'):
            window.file_label.setText(f"{len(spectra_list)} spectra loaded from arrays")
            
        # Show the window if requested
        if show:
            window.show()
            
        # Set focus to canvas for keyboard shortcuts
        window.canvas.setFocus()
        
        # Start event loop if we created the QApplication
        if need_exec:
            # Check if we're in IPython
            try:
                ipython = get_ipython()
                # IPython case - enable Qt integration
                try:
                    ipython.magic('gui qt')
                    print("Enabled Qt integration in IPython - GUI should stay open.")
                except Exception as e:
                    print(f"Could not enable Qt integration: {e}")
                    print("You may need to run '%gui qt' manually.")
            except NameError:
                # Regular Python case
                print("GUI launched. Close window to return to Python prompt.")
                app.exec_()  # Remove sys.exit() here
                print("GUI closed, returning to Python.")
            
        return window
        
    except Exception as e:
        print(f"Error launching MultispecViewer: {str(e)}")
        import traceback
        traceback.print_exc()
        raise


def _prepare_spectra_list(spectrum_data):
    """
    Convert various input formats to a list of spectrum objects compatible with MultispecViewer.
    
    Parameters
    ----------
    spectrum_data : various
        Input spectrum data in various formats
        
    Returns
    -------
    list
        List of spectrum objects (rb_spectrum or XSpectrum1D)
    """
    # Import here to avoid circular imports
    try:
        from rbcodes.utils.rb_spectrum import rb_spectrum, rb_spectrum_collection
    except ImportError:
        rb_spectrum = None
        rb_spectrum_collection = None
        
    try:
        from linetools.spectra.xspectrum1d import XSpectrum1D
    except ImportError:
        XSpectrum1D = None
    
    # Handle different input types
    if rb_spectrum and isinstance(spectrum_data, rb_spectrum):
        # Single rb_spectrum object
        return [spectrum_data]
        
    elif rb_spectrum_collection and isinstance(spectrum_data, rb_spectrum_collection):
        # rb_spectrum_collection object
        return list(spectrum_data.spectra)
        
    elif XSpectrum1D and isinstance(spectrum_data, XSpectrum1D):
        # Single XSpectrum1D object
        return [spectrum_data]
        
    elif isinstance(spectrum_data, (list, tuple)):
        # List or tuple of spectrum objects
        validated_list = []
        for i, spec in enumerate(spectrum_data):
            if rb_spectrum and isinstance(spec, rb_spectrum):
                validated_list.append(spec)
            elif XSpectrum1D and isinstance(spec, XSpectrum1D):
                validated_list.append(spec)
            else:
                print(f"Warning: Skipping invalid spectrum at index {i} (type: {type(spec)})")
                
        return validated_list
        
    else:
        raise ValueError(f"Unsupported spectrum data type: {type(spectrum_data)}. "
                        f"Expected rb_spectrum, rb_spectrum_collection, XSpectrum1D, or list of these.")


def launch_empty(show=True):
    """
    Launch empty MultispecViewer GUI.
    """
    try:
        # Import required modules
        from PyQt5.QtWidgets import QApplication
        from rbcodes.GUIs.multispecviewer import multispec
        
        # Handle QApplication - reuse existing or create new
        app = QApplication.instance()
        if app is None:
            # No existing QApplication, create one
            app = QApplication(sys.argv)
            need_exec = True
        else:
            # Use existing QApplication
            need_exec = False
            
        # Create MainWindow instance
        window = multispec.MainWindow()
        
        # Show the window if requested
        if show:
            window.show()
            
        # Start event loop if we created the QApplication
        if need_exec:
            # Check if we're in IPython
            try:
                ipython = get_ipython()
                # IPython case - enable Qt integration
                try:
                    ipython.magic('gui qt')
                    print("Enabled Qt integration in IPython - GUI should stay open.")
                except Exception as e:
                    print(f"Could not enable Qt integration: {e}")
                    print("You may need to run '%gui qt' manually.")
            except NameError:
                # Regular Python case
                print("GUI launched. Close window to return to Python prompt.")
                app.exec_()  # This will block until GUI closes
                print("GUI closed, returning to Python.")
        
        return window
        
    except Exception as e:
        print(f"Error launching MultispecViewer: {str(e)}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    sys.exit(main())