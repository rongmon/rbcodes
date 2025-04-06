#!/usr/bin/env python
"""
Setup script for installing the spectral_gui package.
"""
from setuptools import setup, find_packages

setup(
    name="rbcodes-spectral-gui",
    version="0.1.0",
    description="GUI for visualizing and analyzing astronomical spectra",
    author="Original: Authors of PlotSpec_Integrated, Refactored: AI Assistant",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "astropy",
        "scipy",
        "PyQt5",
        "pandas",
        "linetools"
    ],
    entry_points={
        'console_scripts': [
            'spectral_gui=rbcodes.GUIs.spectral_gui:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.6",
)