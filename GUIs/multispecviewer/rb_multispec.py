#!/usr/bin/env python
import sys
import os
import argparse
import importlib

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
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check for version flag
    if args.version:
        print("rb_multispec - MultispecViewer Wrapper")
        print("rbcodes Spectral Analysis Tools")
        return 0
    
    # Filter out our custom arguments for both versions
    filtered_argv = [arg for arg in sys.argv if arg not in ['-c', '--classic', '-v', '--version']]
    sys.argv = filtered_argv
    
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
        
        # Both versions use the same initialization pattern
        app = multispec.QApplication(sys.argv)
        window = multispec.MainWindow()
        window.show()
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

if __name__ == "__main__":
    sys.exit(main())