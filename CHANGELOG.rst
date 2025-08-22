=========
Changelog
=========

Version 2.0.0 (2025-08-21)
=======================

Breaking Changes
--------------
* Removed deprecated rb_plot_spec.py GUI tool (superseded by multispecviewer)
* Removed lmfit dependency and associated Gaussian fitting in deprecated GUI
* All PYTHONPATH-based installation methods removed in favor of pip install
* Fixed ASCII file reading in rb_specgui_new.py to properly handle Table objects

Dependencies
-----------
* Added tqdm>=4.65.0
* Added photutils>=1.0
* Removed lmfit
* All dependencies now have explicit minimum versions

Version 0.1
===========

- Initial release with basic functionality
- Support for spectra visualization and analysis
- Various GUI tools for astronomical data processing