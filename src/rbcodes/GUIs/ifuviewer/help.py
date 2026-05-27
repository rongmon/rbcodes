"""
IFU Viewer — Help dialog.

Phase 13.
"""


def show_help_dialog(parent=None):
    """Open the tabbed help dialog for the IFU Viewer."""
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout,
                                 QTabWidget, QTextEdit, QTableWidget,
                                 QTableWidgetItem, QHeaderView,
                                 QPushButton, QLabel)
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QFont, QFontDatabase

    # ------------------------------------------------------------------
    # Content
    # ------------------------------------------------------------------

    overview_text = """\
IFU Viewer — QFitsView-style IFU datacube viewer for rbcodes
=============================================================

INTERFACE
  Left sidebar      Load and switch between FITS cubes and 2-D images.
                    Right-click a cube → assign another cube as its variance.
  Image panel       2-D display: whitelight, narrowband, channel slice,
                    continuum-subtracted, moment map, or SNR map.
  Spectrum panel    Real-time 1-D spectrum under cursor; locked extractions
                    shown as colored overlays.
  Image controls    Colormap, scale (linear/log/sqrt/square), normalization
                    (ZScale, percentile, manual), right-click drag for
                    contrast/bias adjustment.
  Aperture controls Select extraction mode and parameters (always visible
                    when a cube is loaded).
  Extract strip     Manage locked spectra: save, send to rb_multispec,
                    delete, or clear all.
  Channel slider    Scroll through wavelength slices (Channel mode).

COMMAND-LINE USAGE
  rb_ifuview                              # launch empty
  rb_ifuview cube.fits                    # load one cube
  rb_ifuview cube1.fits cube2.fits        # load multiple files
  rb_ifuview --help                       # print CLI usage

LOADING FILES
  File > Open FITS…  (Ctrl+O)  Load one or more FITS cubes or 2-D images.
  Multiple files appear in the sidebar; click to switch active dataset.
  Right-click a cube in the sidebar → "Use as variance for…" to assign
  a separate variance cube.

COORDINATE DISPLAY
  Pixel x/y and RA/Dec (decimal degrees by default) are shown in the
  status bar as you move the cursor over the image.
  View > Coordinates: Sexagesimal  toggles RA/Dec to hh:mm:ss / dd:mm:ss.

FITS HEADER
  View > Show FITS Header…  opens a searchable monospace header viewer.
  Use the extension dropdown for multi-extension files.
"""

    image_modes_text = """\
IMAGE MODES
===========

Four modes are selected with the buttons above the channel slider.

WHITELIGHT
  Mean-collapsed image over all wavelengths (default on load).

CHANNEL
  Displays a single wavelength slice.  Drag the slider or use the
  step buttons [◄] [►] to scroll through channels.

NARROWBAND
  Integrated image over a user-defined wavelength window.
  1. Click the Narrowband button.
  2. Drag a green span on the spectrum panel to define the window.
     Or use Analysis > Band Range Dialog… (Ctrl+B) for exact values.
  The image updates immediately when the span is released.

CONTINUUM-SUBTRACTED (Cont-sub)
  Narrowband image minus a linear continuum baseline.
  1. Click the Cont-sub button.
  2. Draw the on-band window (green span).
  3. Press 'c' in the spectrum panel to activate continuum mode,
     then draw a red continuum window on each side of the line.
     Press '1' and '2' for a second continuum window if needed.
  4. Press Esc to return to neutral mode.
  The image updates after each span is released.

BAND RANGE DIALOG  (Ctrl+B)
  Set on-band and continuum window values numerically.
  Pre-fills from the current SpanSelector extents.
  Clicking OK also switches to the appropriate image mode automatically.

MOMENT MAPS  (Ctrl+M)
  See the Moment Maps tab for the full workflow.
"""

    extraction_text = """\
SPECTRAL EXTRACTION
===================

MODE SELECTION
  Choose extraction mode in the Aperture Controls bar:

  Single pixel       Click a spaxel → spectrum of that pixel locked.
  Rectangle          Left-drag a rectangle on the image → spectrum of
                     all spaxels within the box.
  Circular           Click a spaxel → circular aperture of given radius.
                     A dashed circle follows the cursor as a preview.
  Circular-Annular   Same as Circular, with annular background subtraction.
                     BG inner / BG outer radii are set in the controls bar.

WEIGHTING
  None               Simple sum / mean / median (Method dropdown).
  Var-weighted       Each spaxel weighted by 1/σ² per channel.
                     Requires a variance cube (see sidebar right-click).
  Optimal (Data)     Spatial-profile weights from the collapsed whitelight
                     image.  Works without variance (no error spectrum).
  Optimal (Gaussian) 2-D Gaussian fit to the whitelight for smoother weights.
                     Falls back to Data mode if the fit fails.

EXTRACT STRIP
  After extraction the spectrum appears in the spectrum panel and the
  Extract strip at the bottom:

  Combo box          Select the active extraction.  Clicking one highlights
                     its aperture on the image (thicker line, others dimmed).
                     Clicking the aperture marker on the image also selects it.
  Save…              Save the selected spectrum to FITS or ASCII.
  → rb_multispec     Send all extracted spectra to the rb_multispec viewer.
  Delete  [Del]      Remove the selected extraction and its image marker.
  Clear all          Remove all extractions.

SAVING SPECTRA
  FITS output includes FLUX, WAVELENGTH, and ERROR extensions.
  RA/Dec provenance (from cube WCS) is stored as header keywords.
  If no error array was available, 5% of flux is used as error — a
  UserWarning is printed to the terminal in that case.
"""

    moment_maps_text = """\
MOMENT MAPS  (Ctrl+M)
=====================

Moment maps collapse the cube along the wavelength axis weighted by flux.

  M0  Integrated flux: Σ F(λ) · Δλ  [flux × Å]
  M1  Flux-weighted velocity centroid  [km/s]
  M2  Flux-weighted velocity dispersion  [km/s]

WORKFLOW
  1. Set the image mode to Narrowband and drag a green span on the
     spectrum to bracket the emission line of interest.
  2. Open the dialog: Analysis > Moment Map… or Ctrl+M.
  3. Select moment order (M0 / M1 / M2).
  4. Click [Use current on-band] to fill the wavelength window from
     the current SpanSelector, or enter values manually.
  5. For M1 / M2: enter the rest wavelength λ_rest in Angstroms
     (the wavelength that should map to v=0).
  6. Optional: enable continuum subtraction and set blue/red continuum
     windows.  A linear baseline is fit and subtracted before the moment
     is computed.
  7. Optional: enable SNR mask.  Spaxels below the threshold are set to
     NaN.  Noise can be estimated from a variance cube, a sky region
     drawn on the image, or from continuum windows.
  8. Click OK.  The moment map replaces the current image.

COLORMAPS
  M0 → gray (flux morphology)
  M1 → RdBu_r (red = receding, blue = approaching)
  M2 → plasma (positive-only dispersion)

RESTORING
  Click any image mode button (Whitelight etc.) to return to normal.
  The moment map image is saved per-dataset — switching to another
  cube and back restores it.  The dialog remembers the last parameters
  used for each dataset independently.

SAVE
  File > Save Current Frame… (Ctrl+S) saves the displayed moment map
  as a 2-D FITS file with spatial WCS header and IFUVMODE keyword.
"""

    crop_text = """\
CROP & ANALYSIS
===============

CROPPING A CUBE
  1. Left-drag a rectangle on the image to define the spatial region.
  2. Either:
       • Ctrl+K  — crop immediately
       • Right-click on image → Crop to selection
  The cropped cube is loaded as a new dataset in the sidebar with
  "_crop" appended to the name.  The original cube is unchanged.

  File > Save Subcube…  saves the currently displayed cube (or any
  dataset) to a new 3-D FITS file with updated WCS header.

SKY REGION
  A sky region is used for SNR estimation in moment maps.
  1. Left-drag a rectangle over an empty sky area.
  2. Right-click → Set as sky region  (or Analysis > Set Sky Region).
  A cyan rectangle marks the region.  It is saved per-dataset and
  restored when you switch back.
  Right-click → Clear sky region  removes it.

REGION OVERLAYS (2-D images)
  When a 2-D FITS image is active, imported ds9 regions are drawn as
  shape overlays (no spectrum extraction).
  Right-click → Clear region overlays  or the Clear all button in
  the extract strip removes all overlays.
  Overlays are saved per-dataset and restored on switch.

DS9 INTEGRATION
  See the ds9 Bridge tab for full details on importing regions from ds9,
  cross-instrument region import (e.g. HST → IFU), and sending images.
"""

    ds9_text = """\
DS9 BRIDGE  (optional)
======================

All ds9 functionality is optional.  The GUI works fully without it.

REQUIREMENTS
  pip install pyds9
  ds9 must be running before connecting.

CONNECTING
  Toolbar: [Connect ds9]
    • If pyds9 is not installed → message box with install instructions.
    • If no ds9 running → message box prompting you to start ds9 first.
    • On success → button becomes [Disconnect]; status bar shows "ds9: connected".

SENDING IMAGES
  Toolbar: [→ Send image]  adds the current image panel frame to the
  ds9 frame queue panel (docked below the sidebar).
  [Send all] in the queue panel pushes all queued frames to ds9.
  [Match WCS] locks all ds9 frames to the same WCS.

IMPORTING REGIONS
  Toolbar: [← Import regions]  or  Analysis > Import regions from ds9.
  Choose "Selected region(s) only" or "All regions" in the dialog.

  On a datacube:   each region extracts a spectrum and draws the shape.
  On a 2-D image:  shapes are drawn as overlays only (no extraction).

CROSS-INSTRUMENT IMPORT  (e.g. HST region → IFU)
  1. Open both images in ds9 (different frames).
  2. Draw the region on the HST image in fk5 coordinates.
  3. In IFU Viewer: make the IFU cube the active dataset in the sidebar.
  4. Import regions from ds9.
  The region sky coordinates are converted to IFU pixel coordinates
  via the cube WCS — the region will be placed at the correct sky
  position if the fields overlap.

LOADING REGION FILES  (no ds9 required)
  File > Load Region File…  loads a ds9 .reg file from disk.
  Works for both fk5 and physical coordinate regions.

SAVING APERTURES AS REGION FILE
  File > Save Apertures as Region File…  exports all current extraction
  apertures back to a .reg file in fk5 coordinates (when WCS available).
"""

    # Shortcut table rows: (key, action)
    spectrum_shortcuts = [
        ("x",         "Set cursor x → new xmin"),
        ("X",         "Set cursor x → new xmax"),
        ("t",         "Set cursor y → new ymax"),
        ("b",         "Set cursor y → new ymin"),
        ("r",         "Reset to full wavelength range"),
        ("a",         "Autoscale Y to visible data"),
        ("[",         "Pan left 10%"),
        ("]",         "Pan right 10%"),
        ("-",         "Zoom out X 20%"),
        ("=",         "Zoom in X 20%"),
        ("l",         "Toggle legend"),
        ("c",         "Activate on-band window (Cont-sub mode; press twice to clear)"),
        ("1",         "Activate continuum window 1 (Cont-sub mode; press twice to clear)"),
        ("2",         "Activate continuum window 2 (Cont-sub mode; press twice to clear)"),
        ("Esc",       "Return to neutral mode (Cont-sub mode)"),
    ]

    app_shortcuts = [
        ("Ctrl+O",          "Open FITS file(s)"),
        ("Ctrl+S",          "Save current image frame as FITS"),
        ("Ctrl+M",          "Open Moment Map dialog"),
        ("Ctrl+K",          "Crop to drawn rectangle"),
        ("Ctrl+B",          "Open Band Range dialog"),
        ("Ctrl+R / Ctrl+W", "Open Spectrum Range & Scale dialog"),
        ("Del / Backspace", "Delete selected extraction"),
        ("Ctrl+Q",          "Quit"),
        ("F1",              "Open this help dialog"),
    ]

    # ------------------------------------------------------------------
    # Build dialog
    # ------------------------------------------------------------------

    dialog = QDialog(parent)
    dialog.setWindowTitle("IFU Viewer — Help")
    dialog.setMinimumSize(820, 580)

    layout = QVBoxLayout(dialog)
    tabs = QTabWidget()

    # Fixed-width font for plain text tabs
    fixed_font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
    fixed_font.setPointSize(10)

    def make_text_tab(text):
        w = QTextEdit()
        w.setReadOnly(True)
        w.setPlainText(text)
        w.setFont(fixed_font)
        return w

    def make_table_tab(rows):
        table = QTableWidget(len(rows), 2)
        table.setHorizontalHeaderLabels(["Key / Shortcut", "Action"])
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        table.setSelectionMode(QTableWidget.NoSelection)
        table.setShowGrid(False)
        table.setAlternatingRowColors(True)
        bold = QFont(fixed_font)
        bold.setBold(True)
        bold.setPointSize(10)
        for i, (key, action) in enumerate(rows):
            key_item = QTableWidgetItem(key)
            key_item.setTextAlignment(Qt.AlignCenter)
            key_item.setFont(bold)
            table.setItem(i, 0, key_item)
            table.setItem(i, 1, QTableWidgetItem(action))
        return table

    # Tab: Overview
    tabs.addTab(make_text_tab(overview_text),     "Overview")
    # Tab: Image Modes
    tabs.addTab(make_text_tab(image_modes_text),  "Image Modes")
    # Tab: Extraction
    tabs.addTab(make_text_tab(extraction_text),   "Extraction")
    # Tab: Moment Maps
    tabs.addTab(make_text_tab(moment_maps_text),  "Moment Maps")
    # Tab: Crop & Analysis
    tabs.addTab(make_text_tab(crop_text),         "Crop & Analysis")
    # Tab: ds9 Bridge
    tabs.addTab(make_text_tab(ds9_text),          "ds9 Bridge")
    # Tab: Keyboard Shortcuts (split into two sub-sections via a combined table)
    combined = (
        [("— Spectrum panel (cursor must be in panel) —", "")] +
        spectrum_shortcuts +
        [("— Application-level —", "")] +
        app_shortcuts
    )
    tabs.addTab(make_table_tab(combined),         "Keyboard Shortcuts")

    layout.addWidget(tabs)

    # Bottom button row
    btn_layout = QHBoxLayout()
    btn_layout.addStretch()
    close_btn = QPushButton("Close")
    close_btn.setFixedWidth(80)
    close_btn.clicked.connect(dialog.accept)
    btn_layout.addWidget(close_btn)
    layout.addLayout(btn_layout)

    # Dark theme styling (matches the GUI palette)
    dialog.setStyleSheet("""
        QDialog { background-color: #1e1e2e; color: #cdd6f4; }
        QTabWidget::pane {
            border: 1px solid #45475a;
            background-color: #181825;
        }
        QTabBar::tab {
            background-color: #313244; color: #cdd6f4;
            padding: 6px 14px;
            border: 1px solid #45475a;
            border-bottom: none;
            border-radius: 4px 4px 0 0;
        }
        QTabBar::tab:selected { background-color: #181825; color: #cba6f7; }
        QTabBar::tab:hover    { background-color: #45475a; }
        QTableWidget {
            background-color: #181825; color: #cdd6f4;
            alternate-background-color: #1e1e2e;
            gridline-color: #313244;
            border: none;
        }
        QHeaderView::section {
            background-color: #313244; color: #cdd6f4;
            padding: 4px; border: 1px solid #45475a;
        }
        QTextEdit {
            background-color: #181825; color: #cdd6f4;
            border: 1px solid #45475a;
            border-radius: 4px;
            padding: 8px;
        }
        QPushButton {
            background-color: #313244; color: #cdd6f4;
            border: 1px solid #45475a;
            border-radius: 4px;
            padding: 5px 14px;
        }
        QPushButton:hover   { background-color: #45475a; }
        QPushButton:pressed { background-color: #181825; }
    """)

    if parent is not None:
        dialog.setPalette(parent.palette())

    dialog.exec_()
