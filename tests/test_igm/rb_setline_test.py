"""
Unit tests for rb_setline.py

These tests verify the functionality of the rb_setline module,
ensuring that the improvements maintain backward compatibility
and correctly implement the enhanced features.
"""

import unittest
import numpy as np
from unittest.mock import patch, mock_open
import os
import sys
import tempfile
import logging

# Disable logging during tests to reduce noise
logging.disable(logging.CRITICAL)

# Import the module to test - use direct import without try/except
# since we know where it is in your environment
from rbcodes.IGM.rb_setline import rb_setline, read_line_list

# Sample data for mocking
MOCK_ATOM_DATA = """
col1 col2 col3 col4
H I 1215.67 0.4164 6.265e8
C IV 1548.20 0.1908 2.643e8
C IV 1550.77 0.09522 2.628e8
Si II 1260.42 1.007 2.533e9
Mg II 2796.35 0.6123 2.612e8
Mg II 2803.53 0.3054 2.592e8
"""

MOCK_LLS_DATA = """# LLS list
# wrest ion fval
1215.67 H I 0.4164
1025.72 H I 0.07912
972.54 H I 0.02901
"""

class TestRbSetline(unittest.TestCase):
    """Test cases for rb_setline function."""

    def setUp(self):
        """Set up the test case."""
        # Create a temporary mock file
        self.temp_file = tempfile.NamedTemporaryFile(delete=False)
        self.temp_file.write(MOCK_ATOM_DATA.encode('utf-8'))
        self.temp_file.close()

    def tearDown(self):
        """Clean up after the test case."""
        # Delete the temporary file
        if hasattr(self, 'temp_file') and os.path.exists(self.temp_file.name):
            os.unlink(self.temp_file.name)
        
        # Clear the line list cache - using the actual module name
        if '_LINE_LIST_CACHE' in dir(sys.modules['rbcodes.IGM.rb_setline']):
            sys.modules['rbcodes.IGM.rb_setline']._LINE_LIST_CACHE.clear()

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    @patch('astropy.io.ascii.read')
    def test_closest_match(self, mock_read, mock_resource):
        """Test the 'closest' method of rb_setline."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.temp_file.name
        
        # Mock the astropy read function to return structured data
        mock_read.return_value = {
            'col1': ['H I', 'C IV', 'C IV', 'Si II', 'Mg II', 'Mg II'],
            'col2': [1215.67, 1548.20, 1550.77, 1260.42, 2796.35, 2803.53],
            'col3': [0.4164, 0.1908, 0.09522, 1.007, 0.6123, 0.3054],
            'col4': [6.265e8, 2.643e8, 2.628e8, 2.533e9, 2.612e8, 2.592e8]
        }
        
        # Test finding the closest match to Mg II 2796
        result = rb_setline(2796.3, 'closest')
        self.assertIsInstance(result, dict)
        self.assertAlmostEqual(result['wave'], 2796.35)
        # The actual name format should match what your implementation creates
        self.assertEqual(result['name'], 'Mg II 2796')  # Updated to match the actual format
        self.assertAlmostEqual(result['fval'], 0.6123)
        self.assertAlmostEqual(result['gamma'], 2.612e8)
        
        # Test finding a match for a wavelength between two lines
        result = rb_setline(1549.9, 'closest')
        self.assertIsInstance(result, dict)
        self.assertAlmostEqual(result['wave'], 1550.77)  # Adjusted to match actual implementation
        self.assertEqual(result['name'], 'C IV 1550')

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    @patch('astropy.io.ascii.read')
    def test_exact_match(self, mock_read, mock_resource):
        """Test the 'Exact' method of rb_setline."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.temp_file.name
        
        # Mock the astropy read function to return structured data
        mock_read.return_value = {
            'col1': ['H I', 'C IV', 'C IV', 'Si II', 'Mg II', 'Mg II'],
            'col2': [1215.67, 1548.20, 1550.77, 1260.42, 2796.35, 2803.53],
            'col3': [0.4164, 0.1908, 0.09522, 1.007, 0.6123, 0.3054],
            'col4': [6.265e8, 2.643e8, 2.628e8, 2.533e9, 2.612e8, 2.592e8]
        }
        
        # Test finding an exact match
        result = rb_setline(1215.67, 'Exact')
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['wave']), 1)
        self.assertAlmostEqual(result['wave'][0], 1215.67)
        self.assertEqual(result['name'][0], 'H I 1215')
        
        # Test finding an exact match with a small tolerance
        result = rb_setline(1215.671, 'Exact')
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['wave']), 1)
        self.assertAlmostEqual(result['wave'][0], 1215.67)
        
        # Test no match found
        result = rb_setline(1500.0, 'Exact')
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['wave']), 0)

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    @patch('astropy.io.ascii.read')
    def test_name_match(self, mock_read, mock_resource):
        """Test the 'Name' method of rb_setline."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.temp_file.name
        
        # Mock the astropy read function to return structured data
        mock_read.return_value = {
            'col1': ['H I', 'C IV', 'C IV', 'Si II', 'Mg II', 'Mg II'],
            'col2': [1215.67, 1548.20, 1550.77, 1260.42, 2796.35, 2803.53],
            'col3': [0.4164, 0.1908, 0.09522, 1.007, 0.6123, 0.3054],
            'col4': [6.265e8, 2.643e8, 2.628e8, 2.533e9, 2.612e8, 2.592e8]
        }
        
        # Test finding a match by name
        result = rb_setline(0, 'Name', target_name='H I 1215')
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['wave']), 1)
        self.assertAlmostEqual(result['wave'][0], 1215.67)
        self.assertEqual(result['name'][0], 'H I 1215')
        
        # Test finding multiple matches by name (for C IV)
        result = rb_setline(0, 'Name', target_name='C IV 1548')
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['wave']), 1)
        self.assertAlmostEqual(result['wave'][0], 1548.20)
        
        # Test no match found
        result = rb_setline(0, 'Name', target_name='Fe II 2600')
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['wave']), 0)

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    def test_invalid_method(self, mock_resource):
        """Test that an invalid method raises a ValueError."""
        mock_resource.return_value = self.temp_file.name
        
        with self.assertRaises(ValueError):
            rb_setline(1215.67, 'invalid_method')

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    def test_name_without_target(self, mock_resource):
        """Test that using 'Name' method without target_name raises a ValueError."""
        mock_resource.return_value = self.temp_file.name
        
        with self.assertRaises(ValueError):
            rb_setline(0, 'Name')

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    @patch('astropy.io.ascii.read')
    def test_different_line_list(self, mock_read, mock_resource):
        """Test using a different line list."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.temp_file.name
        
        # Mock the astropy read function to return structured data
        mock_read.return_value = {
            'wrest': [1215.67, 1025.72, 972.54],
            'name': ['H I', 'H I', 'H I'],
            'transition': ['Ly-a', 'Ly-b', 'Ly-g'],
            'ID': [1, 2, 3]
        }
        
        # Test with the LBG line list
        result = rb_setline(1215.67, 'Exact', linelist='LBG')
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result['wave']), 1)
        self.assertAlmostEqual(result['wave'][0], 1215.67)
        self.assertEqual(result['name'][0], 'H I Ly-a')

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    def test_file_not_found(self, mock_resource):
        """Test handling of a missing file."""
        # Mock the resource_filename to point to a non-existent file
        mock_resource.return_value = '/path/to/nonexistent/file.dat'
        
        # Ensure that a FileNotFoundError is raised
        with self.assertRaises(FileNotFoundError):
            rb_setline(1215.67, 'closest')

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    @patch('astropy.io.ascii.read')
    def test_line_list_caching(self, mock_read, mock_resource):
        """Test that line lists are properly cached."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.temp_file.name
        
        # Mock the astropy read function to return structured data
        mock_read.return_value = {
            'col1': ['H I', 'C IV', 'C IV'],
            'col2': [1215.67, 1548.20, 1550.77],
            'col3': [0.4164, 0.1908, 0.09522],
            'col4': [6.265e8, 2.643e8, 2.628e8]
        }
        
        # First call should read the file
        rb_setline(1215.67, 'Exact')
        self.assertEqual(mock_read.call_count, 1)
        
        # Second call should use the cached data
        rb_setline(1548.20, 'Exact')
        self.assertEqual(mock_read.call_count, 1)  # Still just one call
        
        # Call with a different line list should read again
        mock_read.return_value = {
            'wrest': [1215.67, 1025.72, 972.54],
            'name': ['H I', 'H I', 'H I'],
            'transition': ['Ly-a', 'Ly-b', 'Ly-g'],
            'ID': [1, 2, 3]
        }
        rb_setline(1215.67, 'Exact', linelist='LBG')
        self.assertEqual(mock_read.call_count, 2)  # Now two calls


class TestReadLineList(unittest.TestCase):
    """Test cases for read_line_list function."""

    def setUp(self):
        """Set up the test case."""
        # Create temporary files for different line list formats
        self.atom_file = tempfile.NamedTemporaryFile(delete=False)
        self.atom_file.write(MOCK_ATOM_DATA.encode('utf-8'))
        self.atom_file.close()
        
        self.lls_file = tempfile.NamedTemporaryFile(delete=False)
        self.lls_file.write(MOCK_LLS_DATA.encode('utf-8'))
        self.lls_file.close()

    def tearDown(self):
        """Clean up after the test case."""
        # Delete temporary files
        if hasattr(self, 'atom_file') and os.path.exists(self.atom_file.name):
            os.unlink(self.atom_file.name)
        if hasattr(self, 'lls_file') and os.path.exists(self.lls_file.name):
            os.unlink(self.lls_file.name)
        
        # Clear the line list cache - using the actual module name
        if '_LINE_LIST_CACHE' in dir(sys.modules['rbcodes.IGM.rb_setline']):
            sys.modules['rbcodes.IGM.rb_setline']._LINE_LIST_CACHE.clear()

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    @patch('astropy.io.ascii.read')
    def test_read_atom_line_list(self, mock_read, mock_resource):
        """Test reading the atom line list."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.atom_file.name
        
        # Mock the astropy read function to return structured data
        mock_read.return_value = {
            'col1': ['H I', 'C IV', 'C IV'],
            'col2': [1215.67, 1548.20, 1550.77],
            'col3': [0.4164, 0.1908, 0.09522],
            'col4': [6.265e8, 2.643e8, 2.628e8]
        }
        
        # Read the line list
        data = read_line_list('atom')
        
        # Verify the data is correctly structured
        self.assertEqual(len(data), 3)
        self.assertEqual(data[0]['ion'], 'H I 1215')
        self.assertAlmostEqual(data[0]['wrest'], 1215.67)
        self.assertAlmostEqual(data[0]['fval'], 0.4164)
        self.assertAlmostEqual(data[0]['gamma'], 6.265e8)

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    @patch('astropy.io.ascii.read')
    def test_read_lbg_line_list(self, mock_read, mock_resource):
        """Test reading the LBG line list."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.atom_file.name
        
        # Mock the astropy read function to return structured data
        mock_read.return_value = {
            'wrest': [1215.67, 1025.72],
            'name': ['H I', 'H I'],
            'transition': ['Ly-a', 'Ly-b'],
            'ID': [1, 2]
        }
        
        # Read the line list
        data = read_line_list('LBG')
        
        # Verify the data is correctly structured
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0]['ion'], 'H I Ly-a')
        self.assertAlmostEqual(data[0]['wrest'], 1215.67)
        self.assertAlmostEqual(data[0]['fval'], 1.0)

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    def test_read_generic_line_list(self, mock_resource):
        """Test reading a generic format line list."""
        # Mock the resource_filename to point to our temporary file
        mock_resource.return_value = self.lls_file.name
        
        # Read the line list (this will use the actual file)
        data = read_line_list('LLS')
        
        # Verify the data is correctly structured
        self.assertEqual(len(data), 3)
        self.assertEqual(data[0]['ion'], 'H I')  # Updated to match actual implementation
        self.assertAlmostEqual(data[0]['wrest'], 1215.67)
        
    @patch('rbcodes.IGM.rb_setline.resource_filename')
    def test_invalid_line_list(self, mock_resource):
        """Test that an invalid line list label raises a ValueError."""
        mock_resource.return_value = self.atom_file.name
        
        with self.assertRaises(ValueError):
            read_line_list('invalid_label')

    @patch('rbcodes.IGM.rb_setline.resource_filename')
    def test_file_not_found(self, mock_resource):
        """Test handling of a missing file."""
        # Mock the resource_filename to point to a non-existent file
        mock_resource.return_value = '/path/to/nonexistent/file.dat'
        
        # Ensure that a FileNotFoundError is raised
        with self.assertRaises(FileNotFoundError):
            read_line_list('atom')


if __name__ == '__main__':
    unittest.main()