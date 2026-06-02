=========
Changelog
=========

Version 2.4.0 (2026-06-01)
==========================

New Features
------------
* ``rb_align`` — new astrometry alignment module: interactive and batch WCS alignment
  of IFU datacubes and 2D images via chi2-matched source pairs and ``astropy`` WCS fitting.
  Supports interactive, batch, DAO, Gaia, knots, cross-correlation, and auto strategies.
  ``box`` parameter is in arcsec (default 0.1"), converted per-frame via WCS pixel scale.
  Integrated into ``rb_ifuview`` via **Analysis > WCS Alignment** (``Ctrl+Shift+A``).

* ``rb_zfind`` — semi-automated redshift finder: PCA, template, and picket-fence chi2
  search against a spectrum. Curated linelist presets (``zfind_em``, ``zfind_galaxy``,
  ``zfind_igm``, ``zfind_qso``, ``zfind_stellar``). Integrates with ``rb_zgui`` via
  **Find z** button. Standalone CLI: ``rb_zfind spectrum.fits``.

Bug Fixes
---------
* ``rb_align`` interactive window: orphaned cyan edit-mode circle remained visible
  after right-click delete or double-click toggle-off. Fixed by removing matplotlib
  artists explicitly before clearing the tracking dict.

Documentation
-------------
* Added ``docs/GUIs/rb_align/rb_align.md`` and ``docs/GUIs/zfind/rb_zfind.md``.
* ``rb_ifuview`` docs updated with WCS Alignment workflow, keyboard shortcut, and menu entry.

Version 2.3.0 (2026-05-07)
==========================

Enhancements
------------
* ``interactive_continuum_fit``: Added cursor-based and dialog-based navigation shortcuts,
  consistent with ``multispecviewer`` keybindings.

  - ``x`` / ``X``: set left / right x-limit to cursor position.
  - ``t`` / ``b``: set top / bottom y-limit of the panel under the cursor
    (``b`` retains its "add exact spline anchor" behaviour in spline mode).
  - ``[`` / ``]``: pan left / right by 50 % of the current x-axis range.
  - ``W``: open a dialog to type an exact x-range (e.g. ``4100, 4300``);
    pre-filled with the current limits.
  - ``Y``: open a dialog to type an exact y-range for the panel under the cursor;
    pre-filled with the current limits.
  - Overhauled HTML help file: added quick-reference card, grouped shortcut
    table with category headers, and a dedicated Navigation section documenting
    all new keys.

* ``multispecviewer`` updated to v1.5.0.

  - ``a`` / ``A`` keys: autoscale y-axis to the flux range visible in the current
    x-window (``a`` = all panels, ``A`` = panel under cursor).
  - Always-visible spectral-coordinates label in the toolbar showing cursor
    wavelength (Å) and velocity offset (km/s) from the nearest line in the
    active line list at the current redshift.
  - ``L`` key: toggle line-label text on/off without removing tick marks.
  - Help toolbar button (``?``) opens the help dialog without needing a
    keyboard shortcut.
  - Replaced plain-text help popup with a ``QTableWidget``-based dialog
    organised into five tabs: Navigation, Display, Quick Line ID, vStack,
    and Overview.

* ``launch_specgui`` updated to v1.0.6: minor under-the-hood cleanup to
  prevent spurious exception errors in batch mode.

Bug Fixes
---------
* ``interactive_continuum_fit``: zoom state was reset to the initial domain
  whenever a mask was added, removed, or modified via click.  Root cause: six
  call sites captured ``xlim`` before redrawing but passed it to
  ``update_plots()`` without the ``input_xrange`` keyword, causing the method
  to fall back to ``self.domain`` (only updated by the ``+``/``-`` keyboard
  zoom, not the matplotlib toolbar).  All six sites now pass
  ``input_xrange=xlim`` correctly.
  Affected methods: ``handle_add_mask``, ``handle_remove_mask``,
  ``add_exact_spline_point``, ``add_median_spline_point``,
  ``remove_closest_spline_point``, ``clear_spline_points``.

* ``multispecviewer``: JSON load now respects absorber visibility — systems
  whose checkboxes were saved as unchecked were being replotted on load.
  Fixed by wrapping ``_populate_row`` calls with ``blockSignals(True/False)``.
* ``multispecviewer``: ``handle_convert_clicked`` no longer treats multi-value
  return tuples as a bool (always truthy); tuples are now unpacked and the
  error field is checked explicitly.
* ``multispecviewer``: removed double render in ``handle_redshift_submission``
  (``set_redshift_data`` already calls ``plot_redshift_lines`` internally).
* ``multispecviewer``: submitting linelist="None" no longer triggers a
  validation warning; it now passes through and correctly clears plotted lines.
* ``multispecviewer``: fixed cumulative-drift merging bug in
  ``reconcile_linelists`` — clustering now always compares against the first
  element of the current cluster, not the last.
* ``multispecviewer``: simplified vStack canvas replacement to use
  ``main_window.right_layout`` directly, removing fragile attribute-walking
  fallback.

Removed
-------
* ``multispecviewer``: ~30 debug ``print()`` statements from ``vStack.py``.
* ``multispecviewer``: ~299 lines of dead commented-out ``display_line_list()``
  implementation from ``multispec.py``.
* ``multispecviewer``: duplicate ``import matplotlib.pyplot as plt`` in
  ``multispec.py``.

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