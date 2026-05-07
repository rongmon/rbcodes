# Changelog
All notable changes to MultispecViewer will be documented in this file.

## [1.5.0] - 2026-05-07

### Added
- **Y-axis autoscale to visible x-range** (`a`/`A` keys): `a` rescales all panels, `A` rescales only the panel under the cursor. Scales to the flux range within the current zoom window, like IRAF's `y` key.
- **Spectral coordinates status bar**: Always-visible toolbar label showing cursor wavelength (λ in Å) and velocity offset (Δv in km/s) from the nearest line in the active line list at the current redshift. Works with all line lists; computes at z=0 when line list is "None".
- **Line label toggle** (`L` key): Hide/show all text labels on identified lines without removing the vertical tick marks. Useful when many absorbers overlap.
- **Help toolbar button** (`?`): Opens the help dialog directly from the toolbar without needing to know the keyboard shortcut.
- **Tabbed help dialog**: Replaced plain-text help popup with a `QTableWidget`-based dialog organized into five tabs — Navigation, Display, Quick Line ID, vStack, and Overview.

### Fixed
- **JSON load plots all absorbers regardless of visibility** (`AbsorberManager`): Loading a saved JSON file no longer plots absorber systems whose checkboxes were saved as unchecked. Root cause was `_populate_row` firing `cellChanged` → `plot_absorber_lines` during programmatic table population. Fixed by wrapping `_populate_row` calls with `table.blockSignals(True/False)` in both `add_absorber` and `populate_table_from_df`.
- **`handle_convert_clicked` ignores return value** (`multispec.py`): `convert_json_to_text()` returns a 4-tuple and `convert_text_to_json()` returns a 2-tuple; both were assigned to `success` and treated as a bool (always truthy). Fixed by unpacking tuples and checking `error is None`.
- **Double render in `handle_redshift_submission`** (`multispec.py`): `set_redshift_data()` already calls `plot_redshift_lines()` internally; the redundant second call was removed.
- **"None" linelist blocked by `RedshiftInputWidget` validation**: Submitting with linelist="None" now passes through to the canvas (which correctly clears lines) instead of showing a warning.
- **`reconcile_linelists` chaining bug** (`utils.py`): Velocity clustering compared each new line to the last item in the cluster rather than the anchor. Lines that drifted cumulatively could be merged far beyond the threshold. Fixed by always comparing against `current_cluster[0]`.
- **vStack canvas replacement** (`vStack.py`): Simplified canvas swap to use `main_window.right_layout` directly, removing fragile attribute-walking fallback logic.

### Removed
- ~30 debug `print()` statements from `vStack.py` (key press, page turn, layout events)
- ~299 lines of dead commented-out old `display_line_list()` implementation from `multispec.py`
- Duplicate `import matplotlib.pyplot as plt` in `multispec.py`

## [1.3.1] - 2026-04-03

### Fixed
- **vStack linelist selection:** vStack ('v' key) now uses the currently selected linelist from the GUI dropdown instead of requiring Submit to be clicked first
- **Linelist propagation to Absorber:** Fixed bug where Absorber class always used 'atom' linelist regardless of selection. Now correctly passes the selected linelist (e.g., LLS_EUV) through the entire chain
- **Gamma key error:** Fixed KeyError when using non-atom linelists that don't include gamma values (now defaults to 0.0)

## [1.3.0] - 2026-03-12

### Added
- **New linelists in IGM/lines:**
  - `hi.lst`: HI Lyman series lines (30 lines) with f-values and gamma values
  - `euv.lst`: EUV lines (816 lines) with f-values and gamma values
  - `lls_euv.lst`: Merged LLS + EUV lines (567 lines) without duplicates, covering 1.35-5897 Å
  - All generated from linetools LineList data

### Fixed
- **Dark theme color improvements:**
  - Changed default linelist colors from 'white' to 'sky_blue' for better dark theme visibility
  - Updated color exclusion list to filter out jarring colors (white, cream, light_gray, black)
  - Replaced hardcoded color cycling list with dark theme friendly palette
  - Applied consistent color filtering across AbsorberManager, RedshiftInputWidget, and multispec canvas

## [1.2.2] - 2025-12-31
### Fixed
- Right-click line identification dialog now uses the linelist selected in the main GUI instead of hardcoding to 'LLS'
- Improved linelist synchronization between redshift widget and right-click dialog

## [1.2.1] - 2025-11-26
### Added
- **Metadata support for linelists:**
  - Display metadata (target, author, comment) in message box when loading JSON files
  - Multi-line formatted metadata display with proper spacing for readability
  - Pre-fill save dialog with previously loaded metadata for seamless workflow
  - Automatic metadata persistence across save/load cycles
  - Only non-empty metadata fields are displayed

### Enhanced
- **Message box improvements:**
  - Messages now properly append instead of replace previous content
  - Multiple load messages and metadata all visible simultaneously
  - Better workflow visibility with sequential message display

## [1.2.0] - 2025-06-06
### Added
- rb_spectrum: Lightweight spectrum reader/writer as alternative to linetools.XSpectrum1D
- Automatic fallback mechanism: tries rb_spectrum first, falls back to linetools on failure
- Enhanced FITS format support including SDSS binary tables
- Safe exit handling with `_read_failed` flag for robust error recovery
- **New programmatic API for launching GUI with arrays:**
  - `rb_spectrum.from_arrays()`: Create multiple spectra from 2D numpy arrays
  - `rb_spectrum.append()`: Combine multiple rb_spectrum objects into collections
  - `rb_multispec.from_data()`: Launch GUI programmatically with spectrum objects
  - `rb_multispec.launch_empty()`: Launch empty GUI programmatically
- **Enhanced SDSS support:**
  - Binary table parser now handles `loglam` wavelength columns
  - Automatic conversion from log wavelength to linear wavelength
  - Support for `ivar` (inverse variance) error columns with proper conversion
- **Improved IPython/Jupyter integration:**
  - Automatic Qt event loop management for IPython environments
  - Fallback to `%gui qt` magic command when needed
  - Better handling of QApplication conflicts

### Changed
- IO operations now use rb_spectrum by default with linetools fallback
- Improved binary FITS table detection (fixes SDSS spectrum loading)
- Enhanced error handling in spectrum loading pipeline
- **Event loop management:** Removed problematic `sys.exit()` calls that caused Python to exit
- **Documentation:** Updated examples to include programmatic usage patterns

### Enhanced
- **Workflow support:** Complete pipeline from arrays → spectrum objects → GUI display
- **Multi-format compatibility:** Seamless integration between rb_spectrum and XSpectrum1D
- **Error robustness:** Better handling of failed spectrum reads with graceful fallbacks

## [1.1.1] - 2025-05-28
### Fixed
- Improved vStack panel layout with proper axis labeling
- X-axis labels now only appear on bottom-most panels in each column
- Y-axis labels only appear on left-most panels
- Eliminated overlapping velocity labels between panels

## [1.1.0] - 2025-05-28
### Added
- Sortable line list table - click column headers to sort data
- Robust deletion system that works after sorting operations

### Fixed
- Table row deletion now works correctly after sorting
- Improved handling of duplicate line identifications

## [1.0.0] - 2025-05-01
### Added
- Initial stable release
- Multiple FITS spectrum display and analysis
- Absorber cataloging and redshift management
- Line identification tools (right-click dialog, keyboard shortcuts, vStack)
- Dark theme interface
- Data persistence (JSON, CSV, TXT formats)
- Comprehensive keyboard navigation
- Help system and documentation
