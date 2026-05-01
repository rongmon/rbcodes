=========
Changelog
=========

Version 2.2.0 (2026-05-01)
==========================

Compatibility
-------------
* pandas 2.x compatibility: replaced all deprecated ``DataFrame.append()`` calls
  with ``pd.concat()`` across 8 GUI files (19 call sites in total).
  Affected: ``PlotSpec_Integrated``, ``AbsorberManager``, ``vStack``,
  ``prep_guidb``, ``linelist_selection``, ``tableview_pandas``,
  ``spec_fit_gauss2d``, ``menu_toolbars``.

Bug Fixes
---------
* ``rb_multispec``: fixed "Error showing lines" crash caused by custom color names
  (e.g. ``sky_blue``, ``vermillion``) being passed directly to matplotlib in
  ``toggle_show_identified_lines``; colors are now resolved through ``rb_set_color()``.
* ``rb_zgui -x``: fixed ``ModuleNotFoundError: No module named 'utils'`` in
  ``gui_io_xspec.py`` caused by a bare import that broke under package installation;
  replaced with a fully-qualified ``rbcodes.GUIs.zgui.utils`` import.

Version 2.1.2 (2026-04-23)
==========================

Enhancements
------------
* launch_specgui batch mode (Review tab): Added manual "Mark as Failed" button.
  - New ``✗ Mark as Failed`` button in the navigation row lets users flag corrupt,
    noisy, or otherwise unusable spectra without re-running processing.
  - Sets a distinct ``manually_failed`` status (separate from auto-generated ``error``
    status) so manual rejections are clearly distinguishable in exports and logs.
  - Button label toggles to ``↩ Reset to Pending`` when the item is already
    manually failed, allowing the mark to be undone in one click.
  - Manually-failed items appear in the ``Failed`` filter view and are excluded from
    the "Ready Items Only" processing selector.
  - Status bar and details line render ``Manually Failed`` in dark-red italic for
    quick visual identification.

Version 2.1.1 (2026-04-23)
==========================

Bug Fixes
---------
* launch_specgui batch mode: Fixed incorrect ion/filename assignment when saving completed items
  alongside failed items.
  - When "Export Only Completed Items" was checked, the JSON export used positions in the
    filtered list as master table row indices, causing rb_spec objects to be looked up at
    the wrong rows (e.g., a MgII spectrum saved as FeII).
  - Fix: actual master table row indices are now tracked during filtering and passed through
    to the export function, ensuring each completed item is matched to its correct rb_spec object.

Version 2.1.0 (2026-02-19)
==========================

Enhancements
------------
* LLSFitter GUI: Added "Restore Session" feature to Load Results
  - Can now restore filepath, redshift, continuum regions, and all fit parameters from saved JSON
  - Automatically loads spectrum if file exists, shows warning if file moved/missing
  - Backward compatible with older JSON files (uses defaults if fit_parameters not present)
* LLSFitter GUI: Save Results now includes full fit_parameters section
  - Stores initial guesses (C0, C1, logNHI)
  - Stores parameter bounds
  - Stores MCMC settings (walkers, steps, burnin_frac)
  - Stores sigma_clip and plot options

Version 2.0.0 (2025-08-21)
==========================

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