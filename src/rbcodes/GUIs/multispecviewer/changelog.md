# Changelog
All notable changes to MultispecViewer will be documented in this file.

## [Unreleased]

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