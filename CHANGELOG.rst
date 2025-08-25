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

Version 1.0.0
===========
* Major Release with significant enhancements
* New Multispecviewer GUI
* New launch_specgui GUI [Wrapper for rb_spec] as well as its batch mode implmenetation
* Mature rb_spec 1d spectrum analysis toolkit with a lot of bells and whistles
* Large volume of documentation included
* New io modules included
* New LLSFitter GUI included
* New zgui inlcuded for NIRCam data visualizationa nd redshift estimation [optimzied for the JWST EIGER survey data]

Version 0.1
===========

- Initial release with basic functionality
- Support for spectra visualization and analysis
- Various GUI tools for astronomical data processing