"""
help.py — Help dialog for rb_zfind.
"""


def show_help_dialog(parent=None):
    """Open the tabbed help dialog for the Redshift Finder."""
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout,
                                 QTabWidget, QTextEdit, QTableWidget,
                                 QTableWidgetItem, QHeaderView,
                                 QPushButton)
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QFont, QFontDatabase

    # ------------------------------------------------------------------
    # Content
    # ------------------------------------------------------------------

    overview_text = """\
rb_zfind — Semi-automated Redshift Finder
==========================================

PURPOSE
  Find the redshift of a galaxy or QSO from its 1-D spectrum using one
  of three chi² scanning methods: PCA, Template, or Picket Fence.

QUICK START
  1. Select a Method (PCA recommended for most cases).
  2. Set z min / z max to bracket your expected redshift.
  3. Click Run Search.
  4. Click anywhere on the chi² plot to inspect a candidate z.
     — The spectrum panel updates immediately with line markers.
  5. Select a row in the results table to jump to that solution.
  6. Click Accept to return z to the parent GUI (rb_zgui or caller).

INTERFACE
  Chi² panel    Top row.  Click anywhere to set z.  Scroll to zoom.
                Vertical dashed lines mark the top candidate solutions.
  Spectrum      Bottom row.  Shows raw spectrum (and smoothed overlay if
                smoothing > 1).  Template/PCA best-fit model shown in
                green.  Line markers drawn at the selected z.
  Results table Lists the top candidate redshifts, sorted by chi²/dof.

PROGRAMMATIC USE
  from rbcodes.GUIs.zfind.main import launch_zfind
  launch_zfind('spectrum.fits')
  launch_zfind(spec, linelist='zfind_galaxy')

  # headless (no GUI):
  from rbcodes.GUIs.zfind.engine import pca_search
  result = pca_search(spec, template_set='galaxy', z_min=0.0, z_max=1.6)
  print(result.best())
"""

    methods_text = """\
METHODS
=======

PCA  (recommended default)
  Fits a linear combination of DESI redrock PCA eigenvectors at each
  trial z.  Covers the full spectral shape — continuum slope, D4000
  break, Balmer break, emission lines.  Equivalent to what redrock
  does but without the GPU backend.

  PCA sets:
    galaxy   — star-forming and passive galaxies, z = 0.0–1.6
               DESI eigenvectors, rest 1228–11 000 Å.  Best for
               low-z galaxies with detectable continuum.
    qso_loz  — low-z QSOs / broad-line AGN, z = 0.0–2.5
               rest 1329–9634 Å.  Use for Type-1 AGN, Seyferts,
               and low-z quasars with broad emission lines.
    qso_hiz  — high-z QSOs and emission-line galaxies, z = 1.0–6.0
               rest 444–4499 Å (Lyα forest, CIV, CIII], [OIII]).
               Best choice for z > 2 objects; use for high-z
               emission-line galaxies when template method falls short.
    Best of all — tries all three sets, picks lowest chi²/dof

  Data norm / Model norm:
    raw (recommended) — preserves D4000, Balmer break, continuum colour;
                        eigenvector coefficients absorb flux calibration
                        and throughput differences naturally.
    normalize         — divides by continuum; destroys broad continuum
                        features; only use for pure emission-line spectra.

  Resolution — convolves eigenvectors to match instrument LSF before
  the chi² scan.  DESI PCA eigenvectors are high-res; degrading them
  to your instrument resolution improves the fit on resolved lines.

──────────────────────────────────────────────────────────────────────

TEMPLATE
  Chi² scan against bundled MARZ spectral templates.  Good when you
  know the approximate galaxy type.

  Templates:
    EarlyType        — passive/elliptical (D4000, CaII, MgI b)
    Intermediate     — intermediate-type (some emission + absorption)
    LateTypeEmission — star-forming spiral (strong [OII], Hα); also use
                       for high-z emission-line galaxies (z ~ 0–7)
    Composite        — starburst + AGN composite
    QSO              — broad-line QSO
    Best of all      — tries all 5 templates, picks lowest chi²/dof

  Data norm / Model norm:
    raw (recommended) — preserves D4000 and other continuum features;
                        scalar amplitude A absorbs overall flux level.
                        NOTE: A is a single scalar, so large continuum
                        slope differences (reddening, throughput) will
                        degrade the fit.  Use subtract if that is an issue.
    subtract          — removes polynomial continuum; safer if throughput
                        slope differs significantly from the template.
    normalize         — loses broad continuum features.

──────────────────────────────────────────────────────────────────────

PICKET FENCE
  Weighted sum of matched-filter windows centred on each line in the
  chosen linelist.  Scores every trial z without template convolution.
  Fast and interpretable; best for bright emission-line galaxies or
  when you know which lines to expect.

  Modes:
    Direct         — evaluates the linelist window at every trial z
                     and returns a smooth score curve.
    Detect+Match   — first detects peaks in the spectrum, then matches
                     them to the linelist.  Prominence σ controls the
                     peak-detection threshold.

  Linelists:
    zfind_em       — galaxy nebular emission (z ~ 0–7)
    zfind_galaxy   — emission + stellar absorption combined
    zfind_stellar  — stellar absorption lines (MgII, CaII, NaI)
    zfind_igm      — IGM/CGM intervening absorbers on QSO sightlines
    zfind_qso      — QSO/AGN broad emission lines
    (plus full rb_setline lists for power users)

  Data norm:
    subtract (recommended) — continuum-subtracted; lines stand out
                             above zero.
    normalize               — fractional deviation from continuum.
    raw                     — not recommended for picket fence.
"""

    params_text = """\
PARAMETERS
==========

COMMON PARAMETERS
  Method      PCA / Template / Picket Fence (see Methods tab).
  z min/max   Redshift search range.  Auto-filled when you switch method
              or template to a sensible default for that template.
  n steps     Number of trial redshifts in the scan (default 5000).
              Increase for a finer grid (slower).
  Overlay     Linelist drawn on the spectrum at the selected z.
              Defaults to zfind_galaxy (emission + absorption markers).
  λ min/max   Mask observed-frame pixels outside this range (ivar→0).
              Leave blank to use the full spectrum.

PCA / TEMPLATE PARAMETERS
  Data norm   How to preprocess the observed spectrum:
                raw       no continuum removal (recommended; see Methods)
                normalize flux / continuum
                subtract  flux − continuum
  Model norm  How to preprocess the template / eigenvectors:
                raw       as loaded (recommended)
                normalize L2-normalise (PCA) or divide by pseudo-cont
                subtract  mean-center (PCA) or minus pseudo-cont
  Smooth (pix) Boxcar pre-smoothing on the observed flux before chi²;
               1 = off.  Helps on very noisy spectra.
  Resolution   Convolve the template / eigenvectors to the instrument LSF:
                R           resolving power  λ/Δλ  (e.g. R=2000)
                FWHM (Å)    LSF FWHM in Angstroms
                FWHM (km/s) velocity resolution
                FWHM (pix)  FWHM in detector pixels
               Leave value at 0 to skip convolution.

PICKET FENCE PARAMETERS
  Linelist     Spectral lines used for the picket-fence scan.
  Resolution   Same as above; sets the matched-filter window half-width
               to 1.5× the resolution element.
  Win (pix)    Fallback window half-width in pixels when Resolution = 0.
               Rule of thumb: 1–2× the resolution element.
  Smooth (pix) Gaussian pre-smoothing FWHM in pixels (improves SNR on
               broad or noisy lines).  1 = off.
  Data norm    subtract (recommended) | normalize | raw.
  Mode         Direct: score curve at every z.
               Detect+Match: peak-finding + linelist matching.
  Prominence σ Peak detection threshold for Detect+Match mode.
"""

    # Keyboard shortcut rows
    shortcuts = [
        ("x",       "Set cursor x → new xmin"),
        ("X",       "Set cursor x → new xmax"),
        ("t",       "Set cursor y → new ymax"),
        ("b",       "Set cursor y → new ymin"),
        ("r",       "Reset to full range"),
        ("o",       "Zoom out x-axis 1.5×"),
        ("[",       "Pan left one window width"),
        ("]",       "Pan right one window width"),
        ("a",       "Autoscale Y to the visible X range"),
        ("S",       "Increase display smoothing (spectrum panel only)"),
        ("U",       "Decrease display smoothing (spectrum panel only)"),
        ("scroll",  "Zoom x-axis (in/out)"),
        ("F1",      "Open this help dialog"),
        ("",        ""),
        ("Click chi² panel",  "Set selected z to cursor position"),
        ("Click results row", "Jump to that candidate z"),
    ]

    # ------------------------------------------------------------------
    # Build dialog
    # ------------------------------------------------------------------

    dialog = QDialog(parent)
    dialog.setWindowTitle("rb_zfind — Help")
    dialog.setMinimumSize(820, 560)

    layout = QVBoxLayout(dialog)
    tabs = QTabWidget()

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
        table.setHorizontalHeaderLabels(["Key / Action", "Description"])
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
            if key:
                key_item.setFont(bold)
            table.setItem(i, 0, key_item)
            table.setItem(i, 1, QTableWidgetItem(action))
        return table

    tabs.addTab(make_text_tab(overview_text),  "Overview")
    tabs.addTab(make_text_tab(methods_text),   "Methods")
    tabs.addTab(make_text_tab(params_text),    "Parameters")
    tabs.addTab(make_table_tab(shortcuts),     "Keyboard Shortcuts")

    layout.addWidget(tabs)

    btn_layout = QHBoxLayout()
    btn_layout.addStretch()
    close_btn = QPushButton("Close")
    close_btn.setFixedWidth(80)
    close_btn.clicked.connect(dialog.accept)
    btn_layout.addWidget(close_btn)
    layout.addLayout(btn_layout)

    # Apply dark stylesheet only when the parent dialog is in dark mode;
    # in light mode let Qt use the system palette so the help window matches.
    if getattr(parent, 'dark_theme', False):
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
        dialog.setPalette(parent.palette())

    dialog.exec_()
