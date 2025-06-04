# Changelog
All notable changes to MultispecViewer will be documented in this file.

## [Unreleased]

## [1.2.0] - 2025-01-XX
### Added
- rb_spectrum: Lightweight spectrum reader/writer as alternative to linetools.XSpectrum1D
- Automatic fallback mechanism: tries rb_spectrum first, falls back to linetools on failure
- Enhanced FITS format support including SDSS binary tables
- Safe exit handling with `_read_failed` flag for robust error recovery

### Changed
- IO operations now use rb_spectrum by default with linetools fallback
- Improved binary FITS table detection (fixes SDSS spectrum loading)
- Enhanced error handling in spectrum loading pipeline


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