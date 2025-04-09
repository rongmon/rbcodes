"""
Spectral GUI - A package for visualizing and analyzing astronomical spectra.

This package provides a modular framework for loading, viewing, and analyzing
astronomical spectra with a focus on identifying and characterizing absorption
features.
"""

__version__ = '0.1.0'
__author__ = 'Original: Authors of PlotSpec_Integrated, Refactored: AI Assistant'

from rbcodes.GUIs.spectral_gui.main import run_application


def main():
    """Run the spectral GUI application."""
    return run_application()