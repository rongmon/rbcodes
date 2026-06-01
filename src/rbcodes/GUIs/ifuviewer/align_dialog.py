"""
align_dialog.py — Interactive WCS alignment dialog for rb_ifuview.

Embeds two matplotlib canvases (reference + target) inside a QDialog and
drives wcs_align via GUI clicks, bypassing rb_align's built-in matplotlib
window.  All source matching uses add_pair() from the wcs_align API.

The target image is passed in directly from rb_ifuview's image panel
(_data_raw), so the user sees and aligns exactly what is on screen —
whitelight, narrowband, moment map, or any other display mode.

Works for both IFU cubes (3D) and 2D FITS images as the target.

State machine
-------------
NO_REF   → reference image not loaded yet
IDLE     → ready: left-click on reference panel picks a source
PENDING  → source N pending: left-click on target confirms; right-click cancels

Zoom/pan: use the navigation toolbar above each panel (or scroll wheel).
Clicks are silently ignored while the toolbar is in zoom/pan mode.
"""
from __future__ import annotations

import copy
import numpy as np

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox,
    QFileDialog, QMessageBox, QSizePolicy, QWidget,
)
from PyQt5.QtCore import Qt

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg, NavigationToolbar2QT as NavToolbar)

from rbcodes.rb_align import wcs_align
from rbcodes.rb_align.sources import _refine_centroid, _compute_norm


# ---------------------------------------------------------------------------
# Scroll zoom helper — independent per panel
# ---------------------------------------------------------------------------

def _connect_scroll_zoom(ax, canvas):
    def _on_scroll(event):
        if event.inaxes is not ax:
            return
        factor = 0.85 if event.button == 'up' else 1.15
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        cx, cy = event.xdata, event.ydata
        ax.set_xlim([cx + (x - cx) * factor for x in xlim])
        ax.set_ylim([cy + (y - cy) * factor for y in ylim])
        canvas.draw_idle()
    canvas.mpl_connect('scroll_event', _on_scroll)


# ---------------------------------------------------------------------------
# Strategy options dialog
# ---------------------------------------------------------------------------

class StrategyOptionsDialog(QDialog):
    """
    Small popup for per-strategy source detection parameters.

    Shown fields depend on the strategy:
      cross_corr / dao  — FWHM (px), threshold (σ), max sources
      knots             — threshold (σ), max knots
      gaia              — search radius (deg), max stars
    """

    def __init__(self, strategy: str, current_opts: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Options — {strategy}")
        self.setFixedWidth(320)
        self._strategy = strategy
        self._widgets  = {}

        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        form = QHBoxLayout()
        left  = QVBoxLayout()
        right = QVBoxLayout()
        form.addLayout(left)
        form.addLayout(right)
        layout.addLayout(form)

        def _add_row(label, widget, tooltip=''):
            left.addWidget(QLabel(label))
            if tooltip:
                widget.setToolTip(tooltip)
            right.addWidget(widget)

        if strategy in ('cross_corr', 'dao'):
            w = QDoubleSpinBox()
            w.setRange(0.0, 200.0)
            w.setDecimals(1)
            w.setSuffix(" px")
            w.setValue(current_opts.get('fwhm', 0.0))
            _add_row("Source FWHM:", w,
                     "Expected source FWHM in pixels.\n"
                     "Set to 0 to auto-estimate as box × 0.5.\n"
                     "Increase for seeing-broadened QSOs in KCWI/MUSE.")
            self._widgets['fwhm'] = w

            w = QDoubleSpinBox()
            w.setRange(0.5, 20.0)
            w.setDecimals(1)
            w.setSuffix(" σ")
            w.setValue(current_opts.get('threshold_sigma', 3.0))
            _add_row("Detection threshold:", w,
                     "Minimum source brightness in units of image σ.\n"
                     "Lower = more sources found (but more false positives).")
            self._widgets['threshold_sigma'] = w

            w = QSpinBox()
            w.setRange(1, 500)
            w.setValue(current_opts.get('max_sources', 50))
            _add_row("Max sources:", w,
                     "Maximum number of sources to use from each image.")
            self._widgets['max_sources'] = w

        elif strategy == 'knots':
            w = QDoubleSpinBox()
            w.setRange(0.5, 20.0)
            w.setDecimals(1)
            w.setSuffix(" σ")
            w.setValue(current_opts.get('threshold_sigma', 2.0))
            _add_row("Detection threshold:", w,
                     "Segmentation threshold in units of image σ.\n"
                     "Lower = more knots detected.")
            self._widgets['threshold_sigma'] = w

            w = QSpinBox()
            w.setRange(1, 100)
            w.setValue(current_opts.get('max_knots', 10))
            _add_row("Max knots:", w,
                     "Maximum number of knots/segments to return.")
            self._widgets['max_knots'] = w

        elif strategy == 'gaia':
            w = QDoubleSpinBox()
            w.setRange(0.01, 2.0)
            w.setDecimals(3)
            w.setSuffix(" deg")
            w.setValue(current_opts.get('radius_deg', 0.5))
            _add_row("Search radius:", w,
                     "Cone search radius for Gaia DR3 query.\n"
                     "Capped at 0.5° if auto-estimated from WCS.")
            self._widgets['radius_deg'] = w

            w = QSpinBox()
            w.setRange(1, 500)
            w.setValue(current_opts.get('max_stars', 50))
            _add_row("Max stars:", w,
                     "Maximum number of Gaia stars to match (brightest first).")
            self._widgets['max_stars'] = w

        # OK / Cancel
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        ok  = QPushButton("OK")
        ok.setDefault(True)
        ok.clicked.connect(self.accept)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        btn_row.addWidget(ok)
        btn_row.addWidget(cancel)
        layout.addLayout(btn_row)

    def get_options(self) -> dict:
        return {key: w.value() for key, w in self._widgets.items()}


# ---------------------------------------------------------------------------
# AlignDialog
# ---------------------------------------------------------------------------

class AlignDialog(QDialog):
    """
    Two-panel interactive WCS alignment dialog.

    Parameters
    ----------
    cube : IFUCube or FITSImage
        Active dataset (provides flux, wcs, header for the aligner).
    target_image : np.ndarray
        The 2D image currently displayed in rb_ifuview (ImageCanvas._data_raw).
        This is exactly what the user sees — whitelight, narrowband, moment
        map, etc.  Passed straight into aligner.targets[0].image2d.
    target_label : str
        Human-readable name for the target panel title (e.g. "Narrowband").
    parent : QWidget or None
    """

    def __init__(self, cube, target_image: np.ndarray,
                 target_label: str = '', parent=None):
        super().__init__(parent)
        self.setWindowTitle("Align to Reference — rb_align")
        self.setMinimumSize(1100, 720)

        self._cube         = cube
        self._tgt_img      = np.asarray(target_image, dtype=np.float64)
        self._target_label = target_label
        self._aligner      = None
        self._state        = 'NO_REF'
        self._box          = 5
        self._ref_stretch  = 'zscale'
        self._tgt_stretch  = 'zscale'

        # Per-strategy detection options (populated by StrategyOptionsDialog)
        self._strategy_opts = {
            'cross_corr': dict(fwhm=0.0, threshold_sigma=3.0, max_sources=50),
            'dao':        dict(fwhm=0.0, threshold_sigma=3.0, max_sources=50),
            'knots':      dict(threshold_sigma=2.0, max_knots=10),
            'gaia':       dict(radius_deg=0.5,  max_stars=50),
            'interactive': {},
            'batch':       {},
        }

        # Pending click state
        self._pending_ref_xy   = None
        self._pending_tgt_pred = None

        # Matplotlib artists — one per confirmed pair
        self._ref_crosses = []
        self._ref_labels  = []
        self._tgt_crosses = []
        self._tgt_labels  = []

        # Pending artists (cleared on cancel or confirmation)
        self._pending_ref_artist = None
        self._pending_tgt_circle = None


        # Results — set after align()
        self.corrected_wcs = None
        self.aligner       = None

        self._build_ui()
        # Display target image immediately (stretch combos default to 'zscale')
        self._show_image(self._tgt_ax, self._tgt_canvas,
                         self._tgt_img, self._tgt_stretch)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(4)

        # ---- Controls row ----
        ctrl = QHBoxLayout()

        # Reference loader
        self._ref_btn = QPushButton("Load Reference…")
        self._ref_btn.setToolTip(
            "Load a reference FITS file (2D image or IFU cube).")
        self._ref_btn.clicked.connect(self._on_load_reference)
        ctrl.addWidget(QLabel("Reference:"))
        ctrl.addWidget(self._ref_btn)

        ctrl.addSpacing(16)

        # Reference stretch (independent)
        ctrl.addWidget(QLabel("Ref stretch:"))
        self._ref_stretch_combo = QComboBox()
        self._ref_stretch_combo.addItems(
            ['zscale', '99.5%', '99%', '98%', '95%', 'minmax'])
        self._ref_stretch_combo.setToolTip(
            "Display normalization for the reference panel only.")
        self._ref_stretch_combo.setFixedWidth(80)
        self._ref_stretch_combo.currentTextChanged.connect(self._on_ref_stretch_changed)
        ctrl.addWidget(self._ref_stretch_combo)

        ctrl.addSpacing(8)

        # Target stretch (independent)
        ctrl.addWidget(QLabel("Tgt stretch:"))
        self._tgt_stretch_combo = QComboBox()
        self._tgt_stretch_combo.addItems(
            ['zscale', '99.5%', '99%', '98%', '95%', 'minmax'])
        self._tgt_stretch_combo.setToolTip(
            "Display normalization for the target panel only.")
        self._tgt_stretch_combo.setFixedWidth(80)
        self._tgt_stretch_combo.currentTextChanged.connect(self._on_tgt_stretch_changed)
        ctrl.addWidget(self._tgt_stretch_combo)

        ctrl.addSpacing(16)

        # Centroid box
        ctrl.addWidget(QLabel("Box:"))
        self._box_spin = QSpinBox()
        self._box_spin.setRange(1, 100)
        self._box_spin.setValue(self._box)
        self._box_spin.setSuffix(" px")
        self._box_spin.setFixedWidth(68)
        self._box_spin.setToolTip(
            "Half-width in pixels of the centroid refinement window.\n"
            "Increase for faint/diffuse sources, decrease for crowded fields.")
        self._box_spin.valueChanged.connect(lambda v: setattr(self, '_box', v))
        ctrl.addWidget(self._box_spin)

        ctrl.addStretch()
        layout.addLayout(ctrl)

        # ---- Reference preprocess row (hidden until a cube is loaded as ref) ----
        self._ref_prep_row = QHBoxLayout()
        self._ref_prep_row_widget = QWidget()
        self._ref_prep_row_widget.setLayout(self._ref_prep_row)

        self._ref_prep_row.addWidget(QLabel("Ref preprocess:"))
        self._ref_prep_combo = QComboBox()
        self._ref_prep_combo.addItems(['Whitelight', 'Narrowband'])
        self._ref_prep_combo.setToolTip(
            "How to collapse the reference IFU cube to a 2D alignment image.")
        self._ref_prep_combo.currentIndexChanged.connect(self._on_ref_prep_changed)
        self._ref_prep_row.addWidget(self._ref_prep_combo)

        self._ref_prep_row.addWidget(QLabel("Collapse:"))
        self._ref_collapse_combo = QComboBox()
        self._ref_collapse_combo.addItems(['mean', 'median', 'sum'])
        self._ref_collapse_combo.setToolTip(
            "Statistic used to combine flux along the wavelength axis.")
        self._ref_collapse_combo.setFixedWidth(72)
        self._ref_prep_row.addWidget(self._ref_collapse_combo)

        self._ref_wl_min_lbl = QLabel("λ min:")
        self._ref_wl_min_spin = QDoubleSpinBox()
        self._ref_wl_min_spin.setRange(0, 99999)
        self._ref_wl_min_spin.setSuffix(" Å")
        self._ref_wl_min_spin.setDecimals(1)
        self._ref_wl_min_spin.setFixedWidth(95)

        self._ref_wl_max_lbl = QLabel("λ max:")
        self._ref_wl_max_spin = QDoubleSpinBox()
        self._ref_wl_max_spin.setRange(0, 99999)
        self._ref_wl_max_spin.setSuffix(" Å")
        self._ref_wl_max_spin.setDecimals(1)
        self._ref_wl_max_spin.setFixedWidth(95)

        self._ref_wl_apply_btn = QPushButton("Apply")
        self._ref_wl_apply_btn.setToolTip(
            "Recompute reference image with this λ range (clears all pairs).")
        self._ref_wl_apply_btn.clicked.connect(self._on_ref_repreprocess)

        for w in (self._ref_wl_min_lbl, self._ref_wl_min_spin,
                  self._ref_wl_max_lbl, self._ref_wl_max_spin):
            self._ref_prep_row.addWidget(w)

        # Apply is always visible (needed for whitelight + collapse method too)
        self._ref_prep_row.addWidget(self._ref_wl_apply_btn)
        self._ref_prep_row.addStretch()
        layout.addWidget(self._ref_prep_row_widget)
        self._ref_prep_row_widget.setVisible(False)   # shown only for cube refs
        self._set_ref_narrowband_visible(False)

        # ---- Strategy row ----
        strat_row = QHBoxLayout()

        strat_row.addWidget(QLabel("Strategy:"))
        self._strategy_combo = QComboBox()
        self._strategy_combo.addItems([
            'interactive',
            'cross_corr',
            'dao',
            'knots',
            'gaia',
            'batch',
        ])
        self._strategy_combo.setToolTip(
            "interactive — click pairs manually (always works)\n"
            "cross_corr  — cross-correlation shift; best for same-instrument "
            "IFU↔IFU with similar PSF\n"
            "dao         — DAOStarFinder point-source detection (photutils)\n"
            "knots       — segmentation centroids; best for arcs/extended sources\n"
            "gaia        — Gaia DR3 catalog matching (requires astroquery, wide fields)\n"
            "batch       — reproject a saved FITS catalog onto both images"
        )
        self._strategy_combo.currentTextChanged.connect(self._on_strategy_changed)
        strat_row.addWidget(self._strategy_combo)

        # Catalog file picker — visible only for 'batch'
        self._catalog_lbl = QLabel("Catalog:")
        self._catalog_path_btn = QPushButton("Load catalog…")
        self._catalog_path_btn.setToolTip(
            "Load a FITS catalog previously saved with save_catalog=... "
            "from an interactive rb_align session.")
        self._catalog_path_btn.clicked.connect(self._on_load_catalog)
        self._catalog_path_label = QLabel("")   # shows filename after load
        self._catalog_path_label.setStyleSheet("color: #a6e3a1;")
        strat_row.addWidget(self._catalog_lbl)
        strat_row.addWidget(self._catalog_path_btn)
        strat_row.addWidget(self._catalog_path_label)

        strat_row.addSpacing(8)

        self._strat_opts_btn = QPushButton("Options…")
        self._strat_opts_btn.setToolTip(
            "Set detection parameters for the selected strategy\n"
            "(FWHM, threshold, max sources, etc.).")
        self._strat_opts_btn.clicked.connect(self._on_strategy_options)
        strat_row.addWidget(self._strat_opts_btn)

        strat_row.addSpacing(8)

        self._find_btn = QPushButton("Find Sources")
        self._find_btn.setToolTip(
            "Run the selected strategy to find and match sources automatically.\n"
            "Results are shown as markers on both panels.\n"
            "You can then add, remove, or adjust pairs manually before clicking Align.")
        self._find_btn.clicked.connect(self._on_find_sources)
        self._find_btn.setEnabled(False)   # enabled once reference is loaded
        strat_row.addWidget(self._find_btn)

        strat_row.addStretch()

        strat_row_widget = QWidget()
        strat_row_widget.setLayout(strat_row)
        layout.addWidget(strat_row_widget)

        # State: track catalog path for batch strategy
        self._catalog_path = None
        self._set_catalog_visible(False)

        # ---- Two-panel area ----
        panels = QHBoxLayout()
        panels.setSpacing(4)

        # Reference panel
        ref_box = QGroupBox("Reference")
        ref_vbox = QVBoxLayout(ref_box)
        ref_vbox.setContentsMargins(4, 2, 4, 4)
        ref_vbox.setSpacing(2)
        self._ref_fig    = Figure(figsize=(5, 5))
        self._ref_ax     = self._ref_fig.add_axes([0.1, 0.08, 0.86, 0.88])
        self._ref_canvas = FigureCanvasQTAgg(self._ref_fig)
        self._ref_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._ref_toolbar = NavToolbar(self._ref_canvas, self)
        self._ref_canvas.mpl_connect('button_press_event', self._on_ref_click)
        self._ref_canvas.mpl_connect('button_press_event', self._on_right_click_delete)
        _connect_scroll_zoom(self._ref_ax, self._ref_canvas)
        ref_vbox.addWidget(self._ref_toolbar)
        ref_vbox.addWidget(self._ref_canvas)
        panels.addWidget(ref_box)

        # Target panel
        tgt_title = f"Target — {self._target_label}" if self._target_label else "Target"
        name = getattr(self._cube, 'name', '')
        if name:
            tgt_title += f"  ({name})"
        self._tgt_box_widget = QGroupBox(tgt_title)
        tgt_vbox = QVBoxLayout(self._tgt_box_widget)
        tgt_vbox.setContentsMargins(4, 2, 4, 4)
        tgt_vbox.setSpacing(2)
        self._tgt_fig    = Figure(figsize=(5, 5))
        self._tgt_ax     = self._tgt_fig.add_axes([0.1, 0.08, 0.86, 0.88])
        self._tgt_canvas = FigureCanvasQTAgg(self._tgt_fig)
        self._tgt_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._tgt_toolbar = NavToolbar(self._tgt_canvas, self)
        self._tgt_canvas.mpl_connect('button_press_event', self._on_tgt_click)
        self._tgt_canvas.mpl_connect('button_press_event', self._on_right_click_delete)
        _connect_scroll_zoom(self._tgt_ax, self._tgt_canvas)
        tgt_vbox.addWidget(self._tgt_toolbar)
        tgt_vbox.addWidget(self._tgt_canvas)
        panels.addWidget(self._tgt_box_widget)

        layout.addLayout(panels, stretch=1)

        # ---- Status / info ----
        self._state_label = QLabel("Load a reference image to begin.")
        self._state_label.setAlignment(Qt.AlignLeft)
        layout.addWidget(self._state_label)

        self._info_label = QLabel("")
        self._info_label.setAlignment(Qt.AlignLeft)
        layout.addWidget(self._info_label)

        # ---- Button row ----
        btn_row = QHBoxLayout()

        self._undo_btn = QPushButton("Undo")
        self._undo_btn.setToolTip("Remove the last confirmed source pair.")
        self._undo_btn.clicked.connect(self._on_undo)
        self._undo_btn.setEnabled(False)

        self._clear_btn = QPushButton("Clear all")
        self._clear_btn.setToolTip("Remove all source pairs and start over.")
        self._clear_btn.clicked.connect(self._on_clear)
        self._clear_btn.setEnabled(False)

        self._align_btn = QPushButton("Align")
        self._align_btn.setToolTip(
            "Fit WCS correction from stored pairs.\n"
            "1 pair  → shift only\n"
            "2 pairs → shift + rotation\n"
            "3+ pairs → full affine (shift, rotation, scale)")
        self._align_btn.clicked.connect(self._on_align)
        self._align_btn.setEnabled(False)

        self._apply_btn = QPushButton("Apply WCS")
        self._apply_btn.setToolTip(
            "Update the active dataset's WCS in memory.\n"
            "Run Align first.")
        self._apply_btn.clicked.connect(self._on_apply)
        self._apply_btn.setEnabled(False)

        self._write_btn = QPushButton("Write corrected FITS…")
        self._write_btn.setToolTip(
            "Save a FITS file with the corrected WCS patched into the header.\n"
            "Original file is never modified.")
        self._write_btn.clicked.connect(self._on_write)
        self._write_btn.setEnabled(False)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)

        btn_row.addWidget(self._undo_btn)
        btn_row.addWidget(self._clear_btn)
        btn_row.addStretch()
        btn_row.addWidget(self._align_btn)
        btn_row.addWidget(self._apply_btn)
        btn_row.addWidget(self._write_btn)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

    # ------------------------------------------------------------------
    # Load reference
    # ------------------------------------------------------------------

    def _on_load_reference(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Reference FITS",
            filter="FITS files (*.fits *.fit *.fts);;All files (*)")
        if not path:
            return

        try:
            from rbcodes.rb_align.io import _load_frame_from_file
            from rbcodes.rb_align.preprocess import _preprocess_frame
            ref_frame = _load_frame_from_file(path, auto_detect=True)
            _preprocess_frame(ref_frame, method='whitelight')
        except Exception as e:
            QMessageBox.critical(self, "Load failed", str(e))
            return

        self._ref_frame = ref_frame
        self._ref_btn.setText(path.split('/')[-1])

        # Show/hide reference preprocess row depending on whether it's a cube
        is_cube = (ref_frame.input_type == 'ifu' and ref_frame.data.ndim == 3)
        self._ref_prep_row_widget.setVisible(is_cube)
        if is_cube and ref_frame.wave is not None:
            self._ref_wl_min_spin.setValue(float(ref_frame.wave[0]))
            self._ref_wl_max_spin.setValue(float(ref_frame.wave[-1]))
        self._ref_prep_combo.setCurrentIndex(0)   # reset to Whitelight
        self._set_ref_narrowband_visible(False)

        self._init_aligner()
        self._show_image(self._ref_ax, self._ref_canvas,
                         ref_frame.image2d, self._ref_stretch)

        self._clear_markers()
        self._state = 'IDLE'
        self._update_state_label()
        self._undo_btn.setEnabled(False)
        self._clear_btn.setEnabled(False)
        self._align_btn.setEnabled(False)
        self._apply_btn.setEnabled(False)
        self._write_btn.setEnabled(False)
        self._info_label.setText("")
        # Enable Find Sources now that a reference is loaded
        strategy = self._strategy_combo.currentText()
        self._find_btn.setEnabled(strategy != 'interactive')

    def _init_aligner(self):
        """Build a fresh wcs_align instance from current reference + target."""
        cube       = self._cube
        ref        = self._ref_frame
        tgt_arr    = cube.flux
        input_type = 'ifu' if tgt_arr.ndim == 3 else 'image'

        self._aligner = wcs_align.from_data(
            reference=(ref.data, ref.wcs, ref.header),
            targets=[(tgt_arr, cube.wcs, cube.header)],
            input_type=input_type,
        )
        # Restore reference image2d (from_data builds fresh Frames)
        self._aligner.reference.image2d = ref.image2d

        # Patch wave (from_data discards it)
        if hasattr(cube, 'wave') and cube.wave is not None:
            self._aligner.targets[0].wave = cube.wave

        # Use the exact image the user is looking at in rb_ifuview
        self._aligner.targets[0].image2d = self._tgt_img

        self._aligner.on_pair_callback = self._on_pair_added

    # ------------------------------------------------------------------
    # Reference preprocess (cube refs only)
    # ------------------------------------------------------------------

    def _set_ref_narrowband_visible(self, visible: bool):
        for w in (self._ref_wl_min_lbl, self._ref_wl_min_spin,
                  self._ref_wl_max_lbl, self._ref_wl_max_spin):
            w.setVisible(visible)

    def _on_ref_prep_changed(self, idx):
        self._set_ref_narrowband_visible(idx == 1)

    def _on_ref_repreprocess(self):
        """Recompute reference image with current settings; clears all pairs."""
        if not hasattr(self, '_ref_frame'):
            return
        ref = self._ref_frame
        if ref.data.ndim != 3:
            return

        method   = self._ref_prep_combo.currentText()
        collapse = self._ref_collapse_combo.currentText()
        try:
            from rbcodes.rb_align.preprocess import _preprocess_frame
            if method == 'Narrowband':
                wmin = self._ref_wl_min_spin.value()
                wmax = self._ref_wl_max_spin.value()
                _preprocess_frame(ref, method='narrowband',
                                  wl_range=[wmin, wmax], collapse=collapse)
            else:
                _preprocess_frame(ref, method='whitelight', collapse=collapse)
        except Exception as e:
            QMessageBox.critical(self, "Preprocess failed", str(e))
            return

        if self._aligner is not None:
            self._aligner.reference.image2d = ref.image2d

        self._show_image(self._ref_ax, self._ref_canvas,
                         ref.image2d, self._ref_stretch)
        self._on_clear()

    # ------------------------------------------------------------------
    # Strategy controls
    # ------------------------------------------------------------------

    def _set_catalog_visible(self, visible: bool):
        for w in (self._catalog_lbl, self._catalog_path_btn,
                  self._catalog_path_label):
            w.setVisible(visible)

    def _on_strategy_changed(self, strategy: str):
        is_interactive = (strategy == 'interactive')
        is_batch = (strategy == 'batch')
        self._set_catalog_visible(is_batch)
        # Options button: no settings for interactive or batch
        self._strat_opts_btn.setEnabled(strategy not in ('interactive', 'batch'))
        # Find Sources makes no sense for interactive — clicks drive that
        ref_loaded = hasattr(self, '_ref_frame')
        self._find_btn.setEnabled(not is_interactive and ref_loaded)

    def _on_strategy_options(self):
        strategy = self._strategy_combo.currentText()
        opts = self._strategy_opts.get(strategy, {})
        dlg = StrategyOptionsDialog(strategy, opts, parent=self)
        if dlg.exec_():
            self._strategy_opts[strategy] = dlg.get_options()

    def _on_load_catalog(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Pairs Catalog",
            filter="FITS files (*.fits *.fit);;All files (*)")
        if path:
            self._catalog_path = path
            self._catalog_path_label.setText(path.split('/')[-1])

    def _on_find_sources(self):
        """Run the selected non-interactive strategy and display results."""
        if self._aligner is None:
            QMessageBox.information(self, "No reference",
                                    "Load a reference image first.")
            return

        strategy = self._strategy_combo.currentText()
        if strategy == 'interactive':
            return   # handled by click state machine

        # Batch requires a catalog
        if strategy == 'batch' and not self._catalog_path:
            QMessageBox.information(self, "No catalog",
                                    "Load a catalog file first.")
            return

        # Clear existing pairs before running automated strategy
        self._aligner.clear_pairs()
        self._clear_markers()

        opts = self._strategy_opts.get(strategy, {})
        try:
            if strategy == 'batch':
                self._aligner.find_sources(
                    strategy='batch',
                    catalog=self._catalog_path,
                    box=self._box)
            else:
                self._aligner.find_sources(
                    strategy=strategy,
                    box=self._box,
                    **opts)
        except Exception as e:
            QMessageBox.critical(self, f"Strategy '{strategy}' failed", str(e))
            return

        n = len(self._aligner.pairs)
        if n == 0:
            QMessageBox.warning(self, "No sources found",
                                f"Strategy '{strategy}' returned 0 pairs.\n"
                                "Try a different strategy or add pairs manually.")
            return

        # Display found pairs on both panels
        self._redraw_all_markers()

        self._undo_btn.setEnabled(True)
        self._clear_btn.setEnabled(True)
        self._align_btn.setEnabled(True)
        self._update_info_label()
        self._state_label.setText(
            f"Found {n} pair(s) via '{strategy}'.  "
            "You can still click to add more or Undo to remove the last one.")

    # ------------------------------------------------------------------
    # Stretch
    # ------------------------------------------------------------------

    def _on_ref_stretch_changed(self, text):
        self._ref_stretch = text
        if hasattr(self, '_ref_frame'):
            self._show_image(self._ref_ax, self._ref_canvas,
                             self._ref_frame.image2d, self._ref_stretch)
        self._redraw_all_markers()

    def _on_tgt_stretch_changed(self, text):
        self._tgt_stretch = text
        self._show_image(self._tgt_ax, self._tgt_canvas,
                         self._tgt_img, self._tgt_stretch)
        self._redraw_all_markers()

    # ------------------------------------------------------------------
    # Image display
    # ------------------------------------------------------------------

    def _show_image(self, ax, canvas, img2d, stretch='zscale'):
        """
        Display img2d on ax, preserving the current zoom window.

        stretch is applied only to this panel; the other panel is unaffected.
        After imshow, autoscale is locked so markers never expand the view.
        """
        if img2d is None:
            return

        ny, nx = img2d.shape
        xlim_saved = ax.get_xlim()
        ylim_saved = ax.get_ylim()
        zoomed = (bool(ax.images)
                  and xlim_saved != (0.0, 1.0)
                  and xlim_saved != (-0.5, nx - 0.5))

        ax.cla()
        norm = _compute_norm(img2d, stretch)
        ax.imshow(img2d, origin='lower', norm=norm, cmap='gray',
                  interpolation='nearest', aspect='auto')
        ax.set_xlabel('x (px)', fontsize=8)
        ax.set_ylabel('y (px)', fontsize=8)

        ax.autoscale(False)

        if zoomed:
            ax.set_xlim(xlim_saved)
            ax.set_ylim(ylim_saved)

        canvas.draw_idle()

    # ------------------------------------------------------------------
    # Toolbar mode guard
    # ------------------------------------------------------------------

    def _toolbar_active(self, toolbar) -> bool:
        return toolbar.mode != ''

    # ------------------------------------------------------------------
    # Click handlers (state machine)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Right-click to delete nearest pair (works on either panel)
    # ------------------------------------------------------------------

    def _on_right_click_delete(self, event):
        """Right-click near a marker on either panel → delete nearest pair."""
        if event.button != 3:
            return
        if self._toolbar_active(self._ref_toolbar) or \
           self._toolbar_active(self._tgt_toolbar):
            return
        if self._state == 'PENDING':
            return   # right-click in PENDING is handled by _on_tgt_click
        if self._aligner is None or not self._aligner.pairs:
            return
        if event.xdata is None or event.ydata is None:
            return

        cx, cy = event.xdata, event.ydata
        if event.inaxes is self._ref_ax:
            dists = [np.sqrt((p[0]-cx)**2 + (p[1]-cy)**2)
                     for p in self._aligner.pairs]
        elif event.inaxes is self._tgt_ax:
            dists = [np.sqrt((p[4]-cx)**2 + (p[5]-cy)**2)
                     for p in self._aligner.pairs]
        else:
            return

        self._aligner.remove_pair(int(np.argmin(dists)))
        self._redraw_all_markers()
        n = len(self._aligner.pairs)
        self._undo_btn.setEnabled(n > 0)
        self._clear_btn.setEnabled(n > 0)
        self._align_btn.setEnabled(n > 0)
        self._update_info_label()

    def _on_ref_click(self, event):
        if self._toolbar_active(self._ref_toolbar):
            return
        if event.inaxes is not self._ref_ax:
            return
        if event.button != 1:
            return
        if self._state != 'IDLE':
            return
        if self._aligner is None:
            return

        ref_img = self._aligner.reference.image2d
        cx, cy = _refine_centroid(ref_img, event.xdata, event.ydata, self._box)

        art, = self._ref_ax.plot(cx, cy, 'r+', markersize=14, mew=2,
                                  linestyle='none', zorder=10)
        self._pending_ref_artist = art
        self._ref_canvas.draw_idle()

        # Predict position on target via WCS chain
        tx_pred, ty_pred = None, None
        try:
            sky = self._aligner.reference.wcs.all_pix2world([[cx, cy]], 0)
            ra, dec = float(sky[0, 0]), float(sky[0, 1])
            pix = self._aligner.targets[0].wcs.all_world2pix([[ra, dec]], 0)
            tx_pred, ty_pred = float(pix[0, 0]), float(pix[0, 1])
        except Exception:
            pass

        if tx_pred is not None:
            ny, nx = self._tgt_img.shape
            if -nx < tx_pred < 2 * nx and -ny < ty_pred < 2 * ny:
                from matplotlib.patches import Circle
                circ = Circle((tx_pred, ty_pred), radius=self._box,
                              edgecolor='cyan', facecolor='none',
                              linestyle='--', linewidth=1.5, zorder=9)
                self._tgt_ax.add_patch(circ)
                self._pending_tgt_circle = circ
                self._tgt_canvas.draw_idle()

        self._pending_ref_xy   = (cx, cy)
        self._pending_tgt_pred = (tx_pred, ty_pred)
        self._state = 'PENDING'
        self._update_state_label()

    def _on_tgt_click(self, event):
        if self._toolbar_active(self._tgt_toolbar):
            return
        if event.inaxes is not self._tgt_ax:
            return
        if self._state != 'PENDING':
            return
        if self._aligner is None:
            return

        if event.button == 3:
            self._cancel_pending()
            return
        if event.button != 1:
            return

        tgt_img = self._aligner.targets[0].image2d
        cx, cy = _refine_centroid(tgt_img, event.xdata, event.ydata, self._box)

        ref_x, ref_y = self._pending_ref_xy
        self._aligner.add_pair(ref_x, ref_y, cx, cy)
        # _on_pair_added draws permanent markers

        self._state = 'IDLE'
        self._update_state_label()

    def _cancel_pending(self):
        for attr, canvas in (('_pending_ref_artist', self._ref_canvas),
                              ('_pending_tgt_circle', self._tgt_canvas)):
            art = getattr(self, attr)
            if art is not None:
                try:
                    art.remove()
                except ValueError:
                    pass
                setattr(self, attr, None)
                canvas.draw_idle()

        self._pending_ref_xy   = None
        self._pending_tgt_pred = None
        self._state = 'IDLE'
        self._update_state_label()

    # ------------------------------------------------------------------
    # Pair callback
    # ------------------------------------------------------------------

    def _on_pair_added(self, ref_x, ref_y, tgt_x, tgt_y, ra, dec):
        n = len(self._aligner.pairs)

        for attr, canvas in (('_pending_ref_artist', self._ref_canvas),
                              ('_pending_tgt_circle', self._tgt_canvas)):
            art = getattr(self, attr)
            if art is not None:
                try:
                    art.remove()
                except ValueError:
                    pass
                setattr(self, attr, None)

        art_r, = self._ref_ax.plot(ref_x, ref_y, 'r+', markersize=14, mew=2,
                                    linestyle='none', zorder=10)
        lbl_r = self._ref_ax.text(ref_x + 3, ref_y + 3, str(n),
                                   color='red', fontsize=8, zorder=11)
        self._ref_crosses.append(art_r)
        self._ref_labels.append(lbl_r)
        self._ref_canvas.draw_idle()

        art_t, = self._tgt_ax.plot(tgt_x, tgt_y, 'r+', markersize=14, mew=2,
                                    linestyle='none', zorder=10)
        lbl_t = self._tgt_ax.text(tgt_x + 3, tgt_y + 3, str(n),
                                   color='red', fontsize=8, zorder=11)
        self._tgt_crosses.append(art_t)
        self._tgt_labels.append(lbl_t)
        self._tgt_canvas.draw_idle()

        self._undo_btn.setEnabled(True)
        self._clear_btn.setEnabled(True)
        self._align_btn.setEnabled(True)
        self._update_info_label()

    # ------------------------------------------------------------------
    # Undo / Clear
    # ------------------------------------------------------------------

    def _on_undo(self):
        if self._aligner is None or not self._aligner.pairs:
            return

        self._aligner.remove_last_pair()

        for lst, canvas in (
            (self._ref_crosses, self._ref_canvas),
            (self._ref_labels,  self._ref_canvas),
            (self._tgt_crosses, self._tgt_canvas),
            (self._tgt_labels,  self._tgt_canvas),
        ):
            if lst:
                try:
                    lst.pop().remove()
                except ValueError:
                    pass

        self._ref_canvas.draw_idle()
        self._tgt_canvas.draw_idle()

        n = len(self._aligner.pairs)
        self._undo_btn.setEnabled(n > 0)
        self._clear_btn.setEnabled(n > 0)
        self._align_btn.setEnabled(n > 0)
        self._update_info_label()

        if self._state == 'PENDING':
            self._cancel_pending()

    def _on_clear(self):
        if self._aligner is not None:
            self._aligner.clear_pairs()
        self._clear_markers()
        self._undo_btn.setEnabled(False)
        self._clear_btn.setEnabled(False)
        self._align_btn.setEnabled(False)
        self._info_label.setText("")
        if self._state == 'PENDING':
            self._state = 'IDLE'
            self._update_state_label()

    def _clear_markers(self):
        for lst, canvas in (
            (self._ref_crosses, self._ref_canvas),
            (self._ref_labels,  self._ref_canvas),
            (self._tgt_crosses, self._tgt_canvas),
            (self._tgt_labels,  self._tgt_canvas),
        ):
            for art in lst:
                try:
                    art.remove()
                except ValueError:
                    pass
            lst.clear()

        self._ref_canvas.draw_idle()
        self._tgt_canvas.draw_idle()

        for attr, canvas in (('_pending_ref_artist', self._ref_canvas),
                              ('_pending_tgt_circle', self._tgt_canvas)):
            art = getattr(self, attr)
            if art is not None:
                try:
                    art.remove()
                except ValueError:
                    pass
                setattr(self, attr, None)

        self._pending_ref_xy   = None
        self._pending_tgt_pred = None

    def _redraw_all_markers(self):
        """Redraw confirmed-pair markers after stretch change."""
        if self._aligner is None:
            return
        for art_list in (self._ref_crosses, self._ref_labels,
                          self._tgt_crosses, self._tgt_labels):
            for art in art_list:
                try:
                    art.remove()
                except ValueError:
                    pass
            art_list.clear()

        for i, (rx, ry, ra, dec, tx, ty) in enumerate(self._aligner.pairs, 1):
            art_r, = self._ref_ax.plot(rx, ry, 'r+', markersize=14, mew=2,
                                        linestyle='none', zorder=10)
            lbl_r = self._ref_ax.text(rx + 3, ry + 3, str(i),
                                       color='red', fontsize=8, zorder=11)
            self._ref_crosses.append(art_r)
            self._ref_labels.append(lbl_r)

            art_t, = self._tgt_ax.plot(tx, ty, 'r+', markersize=14, mew=2,
                                        linestyle='none', zorder=10)
            lbl_t = self._tgt_ax.text(tx + 3, ty + 3, str(i),
                                       color='red', fontsize=8, zorder=11)
            self._tgt_crosses.append(art_t)
            self._tgt_labels.append(lbl_t)

        self._ref_canvas.draw_idle()
        self._tgt_canvas.draw_idle()

    # ------------------------------------------------------------------
    # Align
    # ------------------------------------------------------------------

    def _on_align(self):
        if self._aligner is None or not self._aligner.pairs:
            return

        try:
            self._aligner.align()
        except Exception as e:
            QMessageBox.critical(self, "Alignment failed", str(e))
            return

        new_wcs  = self._aligner._new_wcs[0]
        fit_type = self._aligner.fit_type[0]
        n_used   = self._aligner.n_sources_used[0]
        rms      = self._aligner.rms_residuals[0]

        if new_wcs is None:
            QMessageBox.warning(self, "Alignment failed",
                                "No corrected WCS produced — check pairs.")
            return

        self.corrected_wcs = new_wcs
        self.aligner       = self._aligner

        self._info_label.setText(
            f"Fit: {fit_type}   |   Sources used: {n_used}   "
            f"|   RMS: {rms:.3f}\"")
        self._apply_btn.setEnabled(True)
        self._write_btn.setEnabled(True)

    # ------------------------------------------------------------------
    # Apply / Write
    # ------------------------------------------------------------------

    def _on_apply(self):
        self.accept()

    def _on_write(self):
        if self._aligner is None or self.corrected_wcs is None:
            return

        cube    = self._cube
        default = ''
        for attr in ('path', 'name'):
            val = getattr(cube, attr, None)
            if val:
                default = str(val).replace('.fits', '_wcsfix.fits')
                break

        out_path, _ = QFileDialog.getSaveFileName(
            self, "Save corrected FITS", default,
            filter="FITS files (*.fits *.fit);;All files (*)")
        if not out_path:
            return

        try:
            from rbcodes.rb_align.io import _update_3d_header
            from astropy.io import fits as afits

            frame = self._aligner.targets[0]
            hdr   = copy.deepcopy(frame.header)

            if frame.input_type == 'ifu':
                _update_3d_header(hdr, self.corrected_wcs)
            else:
                for key, val in self.corrected_wcs.to_header().items():
                    hdr[key] = val

            hdul = afits.HDUList([afits.PrimaryHDU(data=frame.data, header=hdr)])
            if hasattr(cube, 'var') and cube.var is not None:
                hdul.append(afits.ImageHDU(data=cube.var, name='VAR'))
            hdul.writeto(out_path, overwrite=True)
        except Exception as e:
            QMessageBox.critical(self, "Write failed", str(e))
            return

        QMessageBox.information(self, "Saved",
                                f"Corrected FITS written to:\n{out_path}")

    # ------------------------------------------------------------------
    # Status / info labels
    # ------------------------------------------------------------------

    def _update_state_label(self):
        n_pending = (len(self._aligner.pairs) + 1 if self._aligner else 1)
        msgs = {
            'NO_REF':  "Load a reference image to begin.",
            'IDLE':    ("Left-click on the reference panel to select a source.  "
                        "Right-click near a marker on either panel to delete that pair.  "
                        "Use the toolbar or scroll wheel to zoom/pan each panel."),
            'PENDING': (f"Source {n_pending} pending — "
                        f"left-click on the target panel to confirm.  "
                        f"Right-click to cancel."),
        }
        self._state_label.setText(msgs.get(self._state, ""))

    def _update_info_label(self):
        n = len(self._aligner.pairs) if self._aligner else 0
        if n == 0:
            self._info_label.setText("")
            return
        fit_hint = {1: "shift only", 2: "shift + rotation"}.get(
            n, "full affine")
        self._info_label.setText(
            f"{n} pair(s) stored  →  fit type when aligned: {fit_hint}.")
