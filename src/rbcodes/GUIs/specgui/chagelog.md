# Changelog - rb_spec GUI

## [1.0.2] - 2025-11-25
### Added
- **Batch Mode - Use Existing Continuum**: New feature to use continuum data embedded in spectrum files
  - Added "Use Existing Continuum" option to Processing panel continuum method dropdown
  - Added "Use Existing Continuum" button in Review panel for per-item application
  - Graceful fallback to flat continuum if no embedded continuum is found
  - Seamless integration with existing batch processing workflow

## [1.0.1] - 2025-10-20
### Changed
   - Bug fix: interactive equivalent width range correctly applied.



## [1.0.0] - 2025-05-28

### Added
- **Unified Launcher**: Single `launch_specgui.py` command for both modes
  - `python launch_specgui.py` - Single spectrum analysis
  - `python launch_specgui.py -b` - Batch processing
  - `python launch_specgui.py -v` - Version information

### Single Spectrum GUI
- Interactive workflow: Input → Redshift → Transition → Continuum → Measurement → Output
- Support for FITS, ASCII, JSON, and linetools file formats
- Interactive continuum fitting with point-and-click masking
- Equivalent width and column density measurements
- Save/load complete analysis sessions as JSON
- Export plots in PNG, PDF, SVG formats

### Batch Processing GUI
- Configure multiple absorption systems for automated analysis
- Import batch items from CSV templates
- Interactive review and editing of results
- Export results to CSV or individual JSON files
- Smart parameter validation and error handling

### Changed
- **UI Improvements**: Better button layout and sizing
  - Centered and compact load file button in input panel
  - Horizontal alignment for continuum fitting buttons
  - Improved button sizing for cleaner interface
- **Version Management**: Centralized version handling from launcher
  - About dialogs now dynamically reference launcher version
  - Consistent version display across all GUI components

### Removed
- Old separate `launch_batch.py` file (functionality moved to unified launcher)
- Standalone argparse handling from main.py and batch_main.py (moved to unified launcher)

---

## Usage Examples
```bash
# Single spectrum analysis
python launch_specgui.py spectrum.fits
python launch_specgui.py analysis.json

# Batch processing  
python launch_specgui.py -b
python launch_specgui.py -b batch_config.json

# Version and help
python launch_specgui.py -v
python launch_specgui.py --help
```

## Dependencies
- Python 3.7+, PyQt5, matplotlib, numpy, scipy, astropy, rbcodes