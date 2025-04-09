#!/usr/bin/env python
"""
Utility script to run a single test module with proper path setup.
This is helpful for debugging specific test modules.

Usage:
    python run_single_test.py test_module_name.py
"""

import sys
import os
import unittest
import importlib.util

def run_single_test(test_file):
    """
    Run a single test module.
    
    Parameters:
    -----------
    test_file : str
        Path to the test file to run
    """
    if not os.path.exists(test_file):
        print(f"Error: Test file '{test_file}' not found.")
        sys.exit(1)
    
    # Get the directory containing this script
    test_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Add the parent directory to the Python path
    parent_dir = os.path.dirname(test_dir)
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)
    
    # Load the test module
    module_name = os.path.splitext(os.path.basename(test_file))[0]
    spec = importlib.util.spec_from_file_location(module_name, test_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    # Run the tests in the module
    print(f"\nRunning tests from: {test_file}\n")
    unittest.main(module)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Error: Please specify a test file to run.")
        print("Usage: python run_single_test.py test_module_name.py")
        sys.exit(1)
    
    run_single_test(sys.argv[1])