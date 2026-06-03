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
    QFileDialog, QMessageBox, QSizePolicy, QWidget, QLineEdit,
)
from PyQt5.QtCore import Qt

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg, NavigationToolbar2QT as NavToolbar)

from rbcodes.rb_align import wcs_align
from rbcodes.rb_align.sources import (_refine_centroid, _compute_norm,
                                       _box_to_pixels)


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
                     "Minimum source brightness in units of background σ.\n"
                     "Lower = more sources found (but more false positives).")
            self._widgets['threshold_sigma'] = w

            w = QSpinBox()
            w.setRange(1, 500)
            w.setValue(current_opts.get('max_sources', 50))
            _add_row("Max sources:", w,
                     "Maximum number of sources to use from each image.")
            self._widgets['max_sources'] = w

            w = QComboBox()
            w.addItems(['sigma_clip', 'mad'])
            w.setCurrentText(current_opts.get('bg_method', 'sigma_clip'))
            _add_row("BG estimator:", w,
                     "sigma_clip — iterative sigma-clipped std (default).\n"
                     "mad        — median absolute deviation × 1.4826.\n"
                     "Use 'mad' for IFU data with correlated noise or\n"
                     "when the source fraction is high (e.g. QSO fields).")
            self._widgets['bg_method'] = w

        elif strategy == 'knots':
            w = QDoubleSpinBox()
            w.setRange(0.5, 20.0)
            w.setDecimals(1)
            w.setSuffix(" σ")
            w.setValue(current_opts.get('threshold_sigma', 2.0))
            _add_row("Detection threshold:", w,
                     "Segmentation threshold in units of background σ.\n"
                     "Lower = more knots detected.")
            self._widgets['threshold_sigma'] = w

            w = QSpinBox()
            w.setRange(1, 100)
            w.setValue(current_opts.get('max_knots', 10))
            _add_row("Max knots:", w,
                     "Maximum number of knots/segments to return.")
            self._widgets['max_knots'] = w

            w = QComboBox()
            w.addItems(['sigma_clip', 'mad'])
            w.setCurrentText(current_opts.get('bg_method', 'sigma_clip'))
            _add_row("BG estimator:", w,
                     "sigma_clip — iterative sigma-clipped std (default).\n"
                     "mad        — median absolute deviation × 1.4826.")
            self._widgets['bg_method'] = w

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

        elif strategy == 'arc_ncc':
            w = QDoubleSpinBox()
            w.setRange(0.5, 60.0)
            w.setDecimals(1)
            w.setSuffix(" \"")
            w.setValue(current_opts.get('search_radius', 5.0))
            _add_row("Search radius:", w,
                     "Half-width in arcsec to search around the WCS-predicted\n"
                     "target position.  Increase if the WCS offset is large.")
            self._widgets['search_radius'] = w

            from PyQt5.QtWidgets import QCheckBox
            w = QCheckBox()
            w.setChecked(current_opts.get('use_LoG', False))
            _add_row("Use LoG filter:", w,
                     "Apply a Laplacian-of-Gaussian filter before NCC.\n"
                     "Suppresses seeing halos; useful when the two datasets\n"
                     "have noticeably different seeing.")
            self._widgets['use_LoG'] = w

            w = QDoubleSpinBox()
            w.setRange(0.0, 20.0)
            w.setDecimals(1)
            w.setSuffix(" px")
            w.setValue(current_opts.get('LoG_sigma', 0.0))
            _add_row("LoG sigma:", w,
                     "Gaussian scale for LoG in pixels.\n"
                     "0 = auto (1.5 px).  Match to expected knot size.")
            self._widgets['LoG_sigma'] = w

            w = QComboBox()
            w.addItems(['mad', 'sigma_clip'])
            w.setCurrentText(current_opts.get('bg_method', 'mad'))
            _add_row("BG estimator:", w,
                     "mad        — median absolute deviation (default, robust).\n"
                     "sigma_clip — iterative sigma-clipped std.")
            self._widgets['bg_method'] = w

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
        from PyQt5.QtWidgets import QCheckBox
        result = {}
        for key, w in self._widgets.items():
            if isinstance(w, QComboBox):
                result[key] = w.currentText()
            elif isinstance(w, QCheckBox):
                result[key] = w.isChecked()
            else:
                result[key] = w.value()
        return result


# ---------------------------------------------------------------------------
# Canvas key filter — catches key events even when the canvas has focus
# ---------------------------------------------------------------------------

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

    # Class-level: persists across dialog instances within the same app session
    _session_ref_path: str = None
    _session_catalog_path: str = None

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
        self._box          = 0.1
        self._ref_stretch  = 'zscale'
        self._ref_scale    = 'linear'
        self._ref_vmin     = None     # used only when stretch == 'manual'
        self._ref_vmax     = None
        self._tgt_stretch  = 'zscale'
        self._tgt_scale    = 'linear'
        self._tgt_vmin     = None
        self._tgt_vmax     = None
        self._last_ref_path = None   # for Reload button

        # Per-strategy detection options (populated by StrategyOptionsDialog)
        self._strategy_opts = {
            'cross_corr': dict(fwhm=0.0, threshold_sigma=3.0, max_sources=50,
                               bg_method='sigma_clip'),
            'dao':        dict(fwhm=0.0, threshold_sigma=3.0, max_sources=50,
                               bg_method='sigma_clip'),
            'knots':      dict(threshold_sigma=2.0, max_knots=10,
                               bg_method='sigma_clip'),
            'gaia':       dict(radius_deg=0.5, max_stars=50),
            'arc_ncc':    dict(search_radius=5.0, use_LoG=False,
                               LoG_sigma=0.0, bg_method='mad'),
            'interactive': {},
            'batch':       {},
        }

        # arc_ncc drag state
        self._arc_ncc_press_xy  = None   # (x, y) where drag started on reference
        self._arc_ncc_rect      = None   # rubber-band Rectangle on reference panel
        self._arc_ncc_tgt_rect  = None   # result box on target panel
        self._arc_ncc_tgt_cross = None   # result cross on target panel

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

        # Edit mode: {pair_idx: (circle_patch, label_text)}
        # Pairs in edit mode show a cyan circle on the target instead of red cross
        self._edit_mode_circles: dict = {}


        # Results — set after align()
        self.corrected_wcs = None
        self.aligner       = None

        self._build_ui()
        # Display target image immediately
        self._refresh_tgt_image()

        # Auto-load reference from previous session if available
        if AlignDialog._session_ref_path:
            self._load_reference_from_path(AlignDialog._session_ref_path)

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

        self._reload_ref_btn = QPushButton("Reload")
        self._reload_ref_btn.setToolTip(
            "Reload the same reference file without browsing.\n"
            "Clears all pairs and resets the aligner.")
        self._reload_ref_btn.clicked.connect(self._on_reload_reference)
        self._reload_ref_btn.setEnabled(False)
        ctrl.addWidget(self._reload_ref_btn)

        ctrl.addSpacing(16)

        # ---- Reference display controls ----
        ctrl.addWidget(QLabel("Ref:"))

        self._ref_stretch_combo = QComboBox()
        self._ref_stretch_combo.addItems(
            ['zscale', '99.5%', '99%', '98%', '95%', 'minmax', 'manual'])
        self._ref_stretch_combo.setToolTip(
            "Normalization for the reference panel.\n"
            "'manual' reveals vmin/vmax fields for exact control.")
        self._ref_stretch_combo.setFixedWidth(72)
        self._ref_stretch_combo.currentTextChanged.connect(self._on_ref_stretch_changed)
        ctrl.addWidget(self._ref_stretch_combo)

        self._ref_scale_combo = QComboBox()
        self._ref_scale_combo.addItems(['linear', 'log', 'sqrt'])
        self._ref_scale_combo.setToolTip(
            "Transfer function for the reference panel.\n"
            "log: compress bright features; sqrt: intermediate.")
        self._ref_scale_combo.setFixedWidth(62)
        self._ref_scale_combo.currentTextChanged.connect(self._on_ref_scale_changed)
        ctrl.addWidget(self._ref_scale_combo)

        self._ref_vmin_edit = QLineEdit()
        self._ref_vmin_edit.setPlaceholderText("vmin")
        self._ref_vmin_edit.setFixedWidth(68)
        self._ref_vmin_edit.setToolTip("Lower display limit (any float or sci notation, e.g. 1e-18)")
        self._ref_vmin_edit.setVisible(False)
        self._ref_vmin_edit.editingFinished.connect(self._refresh_ref_image)
        ctrl.addWidget(self._ref_vmin_edit)

        self._ref_vmax_edit = QLineEdit()
        self._ref_vmax_edit.setPlaceholderText("vmax")
        self._ref_vmax_edit.setFixedWidth(68)
        self._ref_vmax_edit.setToolTip("Upper display limit (any float or sci notation)")
        self._ref_vmax_edit.setVisible(False)
        self._ref_vmax_edit.editingFinished.connect(self._refresh_ref_image)
        ctrl.addWidget(self._ref_vmax_edit)

        ctrl.addSpacing(12)

        # ---- Target display controls ----
        ctrl.addWidget(QLabel("Tgt:"))

        self._tgt_stretch_combo = QComboBox()
        self._tgt_stretch_combo.addItems(
            ['zscale', '99.5%', '99%', '98%', '95%', 'minmax', 'manual'])
        self._tgt_stretch_combo.setToolTip(
            "Normalization for the target panel.\n"
            "'manual' reveals vmin/vmax fields for exact control.")
        self._tgt_stretch_combo.setFixedWidth(72)
        self._tgt_stretch_combo.currentTextChanged.connect(self._on_tgt_stretch_changed)
        ctrl.addWidget(self._tgt_stretch_combo)

        self._tgt_scale_combo = QComboBox()
        self._tgt_scale_combo.addItems(['linear', 'log', 'sqrt'])
        self._tgt_scale_combo.setToolTip(
            "Transfer function for the target panel.\n"
            "log: compress bright features; sqrt: intermediate.")
        self._tgt_scale_combo.setFixedWidth(62)
        self._tgt_scale_combo.currentTextChanged.connect(self._on_tgt_scale_changed)
        ctrl.addWidget(self._tgt_scale_combo)

        self._tgt_vmin_edit = QLineEdit()
        self._tgt_vmin_edit.setPlaceholderText("vmin")
        self._tgt_vmin_edit.setFixedWidth(68)
        self._tgt_vmin_edit.setToolTip("Lower display limit (any float or sci notation)")
        self._tgt_vmin_edit.setVisible(False)
        self._tgt_vmin_edit.editingFinished.connect(self._refresh_tgt_image)
        ctrl.addWidget(self._tgt_vmin_edit)

        self._tgt_vmax_edit = QLineEdit()
        self._tgt_vmax_edit.setPlaceholderText("vmax")
        self._tgt_vmax_edit.setFixedWidth(68)
        self._tgt_vmax_edit.setToolTip("Upper display limit (any float or sci notation)")
        self._tgt_vmax_edit.setVisible(False)
        self._tgt_vmax_edit.editingFinished.connect(self._refresh_tgt_image)
        ctrl.addWidget(self._tgt_vmax_edit)

        ctrl.addSpacing(16)

        # Centroid box
        ctrl.addWidget(QLabel("Box:"))
        self._box_spin = QDoubleSpinBox()
        self._box_spin.setRange(0.01, 30.0)
        self._box_spin.setDecimals(2)
        self._box_spin.setSingleStep(0.05)
        self._box_spin.setValue(self._box)
        self._box_spin.setSuffix(" \"")
        self._box_spin.setFixedWidth(80)
        self._box_spin.setToolTip(
            "Half-width in arcsec of the centroid refinement window.\n"
            "Converted to pixels internally using the WCS pixel scale.\n"
            "Typical values: 0.05–0.5\" for HST, 0.3–1.0\" for KCWI/MUSE.")
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
            'arc_ncc',
            'cross_corr',
            'dao',
            'knots',
            'gaia',
            'batch',
        ])
        self._strategy_combo.setToolTip(
            "interactive — click pairs manually (always works)\n"
            "arc_ncc     — NCC template matching on a user-drawn box;\n"
            "              best for extended arcs/lensed sources with different seeing\n"
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
        self._catalog_path = AlignDialog._session_catalog_path  # restore last used
        self._set_catalog_visible(False)
        if self._catalog_path:
            self._catalog_path_label.setText(self._catalog_path.split('/')[-1])

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
        self._ref_canvas.mpl_connect('button_press_event',   self._on_ref_click)
        self._ref_canvas.mpl_connect('button_press_event',   self._on_right_click_delete)
        self._ref_canvas.mpl_connect('button_press_event',   self._on_double_click)
        self._ref_canvas.mpl_connect('button_release_event', self._on_ref_release)
        self._ref_canvas.mpl_connect('motion_notify_event',  self._on_ref_motion)
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
        self._tgt_canvas.mpl_connect('button_press_event', self._on_double_click)
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

        self._save_catalog_btn = QPushButton("Save Catalog…")
        self._save_catalog_btn.setToolTip(
            "Save current source pairs to a FITS catalog.\n"
            "Reload later with the 'batch' strategy to skip re-detection.")
        self._save_catalog_btn.clicked.connect(self._on_save_catalog)
        self._save_catalog_btn.setEnabled(False)

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

        self._apply_to_btn = QPushButton("Apply WCS to file…")
        self._apply_to_btn.setToolTip(
            "Apply the same WCS correction to another FITS file\n"
            "(e.g. a variance cube, mask, or second exposure).\n"
            "Run Align first.")
        self._apply_to_btn.clicked.connect(self._on_apply_to)
        self._apply_to_btn.setEnabled(False)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)

        btn_row.addWidget(self._undo_btn)
        btn_row.addWidget(self._clear_btn)
        btn_row.addWidget(self._save_catalog_btn)
        btn_row.addStretch()
        btn_row.addWidget(self._align_btn)
        btn_row.addWidget(self._apply_btn)
        btn_row.addWidget(self._write_btn)
        btn_row.addWidget(self._apply_to_btn)
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
        self._load_reference_from_path(path)

    def _on_reload_reference(self):
        """Reload the last-used reference path without a file dialog."""
        if self._last_ref_path:
            self._load_reference_from_path(self._last_ref_path)

    def _load_reference_from_path(self, path: str):
        """Shared logic for initial load and reload."""
        try:
            from rbcodes.rb_align.io import _load_frame_from_file
            from rbcodes.rb_align.preprocess import _preprocess_frame
            ref_frame = _load_frame_from_file(path, auto_detect=True)
            _preprocess_frame(ref_frame, method='whitelight')
        except Exception as e:
            QMessageBox.critical(self, "Load failed", str(e))
            return

        self._ref_frame     = ref_frame
        self._last_ref_path = path
        AlignDialog._session_ref_path = path   # persist across dialog instances
        self._ref_btn.setText(path.split('/')[-1])
        self._reload_ref_btn.setEnabled(True)

        # Show/hide reference preprocess row depending on whether it's a cube
        is_cube = (ref_frame.input_type == 'ifu' and ref_frame.data.ndim == 3)
        self._ref_prep_row_widget.setVisible(is_cube)
        if is_cube and ref_frame.wave is not None:
            self._ref_wl_min_spin.setValue(float(ref_frame.wave[0]))
            self._ref_wl_max_spin.setValue(float(ref_frame.wave[-1]))
        self._ref_prep_combo.setCurrentIndex(0)   # reset to Whitelight
        self._set_ref_narrowband_visible(False)

        self._init_aligner()
        # Reset manual fields when loading a new reference so they get re-initialised
        self._ref_vmin_edit.clear()
        self._ref_vmax_edit.clear()
        self._show_image(self._ref_ax, self._ref_canvas,
                         ref_frame.image2d,
                         stretch=self._ref_stretch, scale=self._ref_scale,
                         reset_zoom=True)

        self._clear_markers()
        self._state = 'IDLE'
        self._update_state_label()
        self._undo_btn.setEnabled(False)
        self._clear_btn.setEnabled(False)
        self._save_catalog_btn.setEnabled(False)
        self._align_btn.setEnabled(False)
        self._apply_btn.setEnabled(False)
        self._write_btn.setEnabled(False)
        self._apply_to_btn.setEnabled(False)
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

        self._ref_vmin_edit.clear()
        self._ref_vmax_edit.clear()
        self._show_image(self._ref_ax, self._ref_canvas,
                         ref.image2d,
                         stretch=self._ref_stretch, scale=self._ref_scale,
                         reset_zoom=False)   # same image size, preserve zoom
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
        is_batch       = (strategy == 'batch')
        is_arc_ncc     = (strategy == 'arc_ncc')
        self._set_catalog_visible(is_batch)
        # Options button: no settings for interactive or batch
        self._strat_opts_btn.setEnabled(strategy not in ('interactive', 'batch'))
        ref_loaded = hasattr(self, '_ref_frame')
        self._find_btn.setEnabled(not is_interactive and ref_loaded)
        # arc_ncc uses drag-to-draw on the reference panel
        if is_arc_ncc:
            self._find_btn.setText("Draw NCC box…")
        else:
            self._find_btn.setText("Find Sources")
        # Exit arc_ncc draw mode if user switches strategy
        if not is_arc_ncc and self._state == 'ARC_NCC':
            self._exit_arc_ncc_mode()
            self._state = 'IDLE'
            self._update_state_label()

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
            AlignDialog._session_catalog_path = path   # persist across dialog instances
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

        # arc_ncc: enter draw mode — user drags a box on the reference panel
        if strategy == 'arc_ncc':
            self._exit_arc_ncc_mode()
            self._state = 'ARC_NCC'
            self._update_state_label()
            return

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
                    box=self._box,
                    headless=True)
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
        self._save_catalog_btn.setEnabled(True)
        self._align_btn.setEnabled(True)
        self._update_info_label()
        self._state_label.setText(
            f"Found {n} pair(s) via '{strategy}'.  "
            "You can still click to add more or Undo to remove the last one.")

    # ------------------------------------------------------------------
    # Display controls (stretch, scale, manual vmin/vmax)
    # ------------------------------------------------------------------

    def _on_ref_stretch_changed(self, text):
        self._ref_stretch = text
        is_manual = (text == 'manual')
        self._ref_vmin_edit.setVisible(is_manual)
        self._ref_vmax_edit.setVisible(is_manual)
        if is_manual:
            self._init_manual_fields(
                self._ref_vmin_edit, self._ref_vmax_edit,
                getattr(self, '_ref_frame', None) and self._ref_frame.image2d,
                self._ref_scale)
        self._refresh_ref_image()

    def _on_tgt_stretch_changed(self, text):
        self._tgt_stretch = text
        is_manual = (text == 'manual')
        self._tgt_vmin_edit.setVisible(is_manual)
        self._tgt_vmax_edit.setVisible(is_manual)
        if is_manual:
            self._init_manual_fields(
                self._tgt_vmin_edit, self._tgt_vmax_edit,
                self._tgt_img, self._tgt_scale)
        self._refresh_tgt_image()

    def _on_ref_scale_changed(self, text):
        self._ref_scale = text
        self._refresh_ref_image()

    def _on_tgt_scale_changed(self, text):
        self._tgt_scale = text
        self._refresh_tgt_image()

    def _init_manual_fields(self, vmin_edit, vmax_edit, img2d, scale):
        """Pre-fill vmin/vmax fields with current zscale limits when first entering manual mode."""
        if img2d is None or not isinstance(img2d, np.ndarray):
            return
        if vmin_edit.text() and vmax_edit.text():
            return   # already set by user — don't overwrite
        try:
            norm = _compute_norm(img2d, 'zscale', scale=scale)
            lo = norm.vmin if norm.vmin is not None else float(np.nanmin(img2d))
            hi = norm.vmax if norm.vmax is not None else float(np.nanmax(img2d))
            vmin_edit.setText(f"{lo:.4g}")
            vmax_edit.setText(f"{hi:.4g}")
        except Exception:
            pass

    def _parse_manual(self, vmin_edit, vmax_edit):
        """Return (vmin, vmax) floats from the edit fields, or (None, None) on error."""
        try:
            lo = float(vmin_edit.text())
            hi = float(vmax_edit.text())
            return lo, hi
        except (ValueError, AttributeError):
            return None, None

    def _refresh_ref_image(self):
        """Redisplay the reference panel with current stretch/scale/vmin/vmax."""
        if not hasattr(self, '_ref_frame') or self._ref_frame is None:
            return
        img = self._ref_frame.image2d
        if img is None:
            return
        vmin, vmax = self._parse_manual(self._ref_vmin_edit, self._ref_vmax_edit) \
            if self._ref_stretch == 'manual' else (None, None)
        self._show_image(self._ref_ax, self._ref_canvas, img,
                         stretch=self._ref_stretch, scale=self._ref_scale,
                         vmin=vmin, vmax=vmax)
        self._redraw_all_markers()

    def _refresh_tgt_image(self):
        """Redisplay the target panel with current stretch/scale/vmin/vmax."""
        vmin, vmax = self._parse_manual(self._tgt_vmin_edit, self._tgt_vmax_edit) \
            if self._tgt_stretch == 'manual' else (None, None)
        self._show_image(self._tgt_ax, self._tgt_canvas, self._tgt_img,
                         stretch=self._tgt_stretch, scale=self._tgt_scale,
                         vmin=vmin, vmax=vmax)
        self._redraw_all_markers()

    # ------------------------------------------------------------------
    # Image display
    # ------------------------------------------------------------------

    def _show_image(self, ax, canvas, img2d, stretch='zscale', scale='linear',
                    vmin=None, vmax=None, reset_zoom: bool = False):
        """
        Display img2d on ax, preserving the current zoom window.

        Parameters
        ----------
        stretch : str
            Normalization: 'zscale' | percentile string | 'minmax' | 'manual'
        scale : str
            Transfer function: 'linear' | 'log' | 'sqrt'
        vmin, vmax : float or None
            Manual limits — used only when stretch == 'manual'.
        reset_zoom : bool
            True → discard saved zoom (use when loading a new reference).
        """
        if img2d is None:
            return

        ny, nx = img2d.shape
        xlim_saved = ax.get_xlim()
        ylim_saved = ax.get_ylim()
        zoomed = (not reset_zoom
                  and bool(ax.images)
                  and xlim_saved != (0.0, 1.0)
                  and xlim_saved != (-0.5, nx - 0.5)
                  and abs(ylim_saved[1] - ylim_saved[0]) < 2 * ny)

        ax.cla()
        norm = _compute_norm(img2d, stretch=stretch, scale=scale,
                             vmin=vmin, vmax=vmax)
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

    # ------------------------------------------------------------------
    # Hover tracking and 'e' key edit mode
    # ------------------------------------------------------------------

    def _on_double_click(self, event):
        """Double-click on either panel → toggle edit mode for nearest pair."""
        if not event.dblclick or event.button != 1:
            return
        if self._toolbar_active(self._ref_toolbar) or \
           self._toolbar_active(self._tgt_toolbar):
            return
        if self._state != 'IDLE':
            return
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
        self._toggle_edit_pair(int(np.argmin(dists)))

    def _toggle_edit_pair(self, pair_idx: int):
        """Enter or exit edit mode for a specific pair index."""
        if pair_idx in self._edit_mode_circles:
            # Toggle off: remove the artist BEFORE deleting from dict,
            # otherwise _redraw_all_markers can't find it to clean up
            val = self._edit_mode_circles.pop(pair_idx)
            if isinstance(val, tuple):
                for a in val:
                    try:
                        a.remove()
                    except Exception:
                        pass
        else:
            # Enter edit mode (placeholder rebuilt in _redraw_all_markers)
            self._edit_mode_circles[pair_idx] = True
        self._redraw_all_markers()

    def _clear_edit_mode(self):
        """Remove all edit-mode circles and reset the dict."""
        for val in self._edit_mode_circles.values():
            if isinstance(val, tuple):
                circ, lbl = val
                for a in (circ, lbl):
                    try:
                        a.remove()
                    except Exception:
                        pass
        self._edit_mode_circles.clear()

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
        # Indices shifted after deletion — remove artists then clear dict
        self._clear_edit_mode()
        self._redraw_all_markers()
        n = len(self._aligner.pairs)
        self._undo_btn.setEnabled(n > 0)
        self._clear_btn.setEnabled(n > 0)
        self._save_catalog_btn.setEnabled(n > 0)
        self._align_btn.setEnabled(n > 0)
        self._update_info_label()

    # ------------------------------------------------------------------
    # arc_ncc drag helpers
    # ------------------------------------------------------------------

    def _exit_arc_ncc_mode(self):
        """Remove any arc_ncc rubber-band and result artists."""
        for attr, canvas in (
            ('_arc_ncc_rect',      self._ref_canvas),
            ('_arc_ncc_tgt_rect',  self._tgt_canvas),
            ('_arc_ncc_tgt_cross', self._tgt_canvas),
        ):
            art = getattr(self, attr)
            if art is not None:
                try:
                    art.remove()
                except Exception:
                    pass
                setattr(self, attr, None)
        self._arc_ncc_press_xy = None
        self._ref_canvas.draw_idle()
        self._tgt_canvas.draw_idle()

    def _on_ref_motion(self, event):
        """Update rubber-band rectangle while dragging in ARC_NCC mode."""
        if self._state != 'ARC_NCC':
            return
        if self._arc_ncc_press_xy is None:
            return
        if event.inaxes is not self._ref_ax or event.xdata is None:
            return
        x0, y0 = self._arc_ncc_press_xy
        x1, y1 = event.xdata, event.ydata
        if self._arc_ncc_rect is not None:
            try:
                self._arc_ncc_rect.remove()
            except Exception:
                pass
        import matplotlib.patches as mpatches
        bx0, bx1 = sorted([x0, x1])
        by0, by1 = sorted([y0, y1])
        self._arc_ncc_rect = mpatches.Rectangle(
            (bx0, by0), bx1 - bx0, by1 - by0,
            edgecolor='cyan', facecolor='none', linewidth=1.5,
            linestyle='--', zorder=9)
        self._ref_ax.add_patch(self._arc_ncc_rect)
        self._ref_canvas.draw_idle()

    def _on_ref_release(self, event):
        """On mouse release in ARC_NCC mode: run NCC with the drawn box."""
        if self._state != 'ARC_NCC':
            return
        if self._arc_ncc_press_xy is None:
            return
        if event.button != 1:
            return
        if event.xdata is None or event.ydata is None:
            return

        x0, y0 = self._arc_ncc_press_xy
        x1, y1 = event.xdata, event.ydata
        self._arc_ncc_press_xy = None

        bx0, bx1 = sorted([int(round(x0)), int(round(x1))])
        by0, by1 = sorted([int(round(y0)), int(round(y1))])

        if bx1 - bx0 < 3 or by1 - by0 < 3:
            self._state_label.setText("Box too small — try again.")
            return

        # Remove old result artists (keep the ref box — already drawn as rubber-band)
        for attr, canvas in (
            ('_arc_ncc_tgt_rect',  self._tgt_canvas),
            ('_arc_ncc_tgt_cross', self._tgt_canvas),
        ):
            art = getattr(self, attr)
            if art is not None:
                try:
                    art.remove()
                except Exception:
                    pass
                setattr(self, attr, None)

        opts = self._strategy_opts.get('arc_ncc', {})
        try:
            from rbcodes.rb_align.sources import _arc_ncc
            _arc_ncc(self._aligner,
                     box_ref=(bx0, by0, bx1, by1),
                     search_radius=opts.get('search_radius', 5.0),
                     use_LoG=opts.get('use_LoG', False),
                     LoG_sigma=opts.get('LoG_sigma', 0.0),
                     bg_method=opts.get('bg_method', 'mad'))
        except Exception as e:
            self._state_label.setText(f"arc_ncc failed: {e}")
            return

        if not self._aligner.pairs:
            self._state_label.setText("NCC produced no result — try a different region.")
            return

        _, _, _, _, tgt_cx, tgt_cy = self._aligner.pairs[-1]
        box_w = bx1 - bx0
        box_h = by1 - by0

        import matplotlib.patches as mpatches
        tgt_rect = mpatches.Rectangle(
            (tgt_cx - box_w / 2.0, tgt_cy - box_h / 2.0), box_w, box_h,
            edgecolor='cyan', facecolor='none', linewidth=1.5,
            linestyle='--', zorder=9)
        self._tgt_ax.add_patch(tgt_rect)
        self._arc_ncc_tgt_rect = tgt_rect

        cross, = self._tgt_ax.plot(tgt_cx, tgt_cy, 'r+',
                                    markersize=14, mew=2, zorder=10)
        self._arc_ncc_tgt_cross = cross

        from rbcodes.rb_align.sources import _estimate_pixel_scale_arcsec
        tgt_scale = _estimate_pixel_scale_arcsec(
            self._aligner.targets[0], self._aligner.targets[0].image2d)
        try:
            ra, dec = self._aligner.pairs[-1][2], self._aligner.pairs[-1][3]
            pred = self._aligner.targets[0].wcs.all_world2pix([[ra, dec]], 0)
            shift_dx = tgt_cx - float(pred[0, 0])
            shift_dy = tgt_cy - float(pred[0, 1])
        except Exception:
            shift_dx = shift_dy = 0.0
        shift_as = np.sqrt((shift_dx * tgt_scale)**2 + (shift_dy * tgt_scale)**2)

        self._tgt_canvas.draw_idle()
        self._undo_btn.setEnabled(True)
        self._clear_btn.setEnabled(True)
        self._save_catalog_btn.setEnabled(True)
        self._align_btn.setEnabled(True)
        self._state_label.setText(
            f"arc_ncc: shift = ({shift_dx:+.2f}, {shift_dy:+.2f}) px"
            f" = {shift_as:.2f}\"  |  Draw another box to add a second pair,"
            f" or click Align.")

    def _on_ref_click(self, event):
        if self._toolbar_active(self._ref_toolbar):
            return
        if event.inaxes is not self._ref_ax:
            return
        if event.button != 1:
            return
        if event.dblclick:
            return   # handled by _on_double_click
        # In ARC_NCC mode, record drag start; actual work done in _on_ref_release
        if self._state == 'ARC_NCC':
            if event.xdata is not None and event.ydata is not None:
                self._arc_ncc_press_xy = (event.xdata, event.ydata)
            return
        if self._state != 'IDLE':
            return
        if self._aligner is None:
            return

        ref_frame = self._aligner.reference
        ref_img   = ref_frame.image2d
        ref_box_px = _box_to_pixels(self._box, ref_frame, ref_img)
        cx, cy = _refine_centroid(ref_img, event.xdata, event.ydata, ref_box_px)

        art, = self._ref_ax.plot(cx, cy, 'r+', markersize=14, mew=2,
                                  linestyle='none', zorder=10)
        self._pending_ref_artist = art
        self._ref_canvas.draw_idle()

        # Predict position on target via WCS chain
        tx_pred, ty_pred = None, None
        try:
            sky = ref_frame.wcs.all_pix2world([[cx, cy]], 0)
            ra, dec = float(sky[0, 0]), float(sky[0, 1])
            pix = self._aligner.targets[0].wcs.all_world2pix([[ra, dec]], 0)
            tx_pred, ty_pred = float(pix[0, 0]), float(pix[0, 1])
        except Exception:
            pass

        if tx_pred is not None:
            ny, nx = self._tgt_img.shape
            if -nx < tx_pred < 2 * nx and -ny < ty_pred < 2 * ny:
                from matplotlib.patches import Circle
                tgt_frame  = self._aligner.targets[0]
                tgt_box_px = _box_to_pixels(self._box, tgt_frame, self._tgt_img)
                circ = Circle((tx_pred, ty_pred), radius=tgt_box_px,
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
        if self._aligner is None:
            return

        if event.button == 3:
            self._cancel_pending()
            return
        if event.button != 1:
            return
        if event.dblclick:
            return   # handled by _on_double_click
        if event.xdata is None or event.ydata is None:
            return

        # ---- Edit mode: re-place the nearest edit circle --------------------
        if self._state == 'IDLE' and self._edit_mode_circles:
            cx, cy = event.xdata, event.ydata
            # Find nearest edit circle by its stored center
            dists = {}
            for idx, val in self._edit_mode_circles.items():
                if isinstance(val, tuple):
                    circ, _ = val
                    ccx, ccy = circ.center
                else:
                    # placeholder — use current pair position
                    ccx, ccy = self._aligner.pairs[idx][4], self._aligner.pairs[idx][5]
                dists[idx] = np.sqrt((ccx - cx)**2 + (ccy - cy)**2)
            nearest_idx = min(dists, key=dists.get)
            tgt_frame  = self._aligner.targets[0]
            tgt_img    = tgt_frame.image2d
            tgt_box_px = _box_to_pixels(self._box, tgt_frame, tgt_img)
            new_tx, new_ty = _refine_centroid(tgt_img, cx, cy, tgt_box_px)
            rx, ry, ra, dec, _, _ = self._aligner.pairs[nearest_idx]
            self._aligner.pairs[nearest_idx] = (rx, ry, ra, dec, new_tx, new_ty)
            # Remove artist before popping from dict
            val = self._edit_mode_circles.pop(nearest_idx)
            if isinstance(val, tuple):
                for a in val:
                    try:
                        a.remove()
                    except Exception:
                        pass
            self._redraw_all_markers()
            self._update_info_label()
            return

        if self._state != 'PENDING':
            return

        tgt_frame  = self._aligner.targets[0]
        tgt_img    = tgt_frame.image2d
        tgt_box_px = _box_to_pixels(self._box, tgt_frame, tgt_img)
        cx, cy = _refine_centroid(tgt_img, event.xdata, event.ydata, tgt_box_px)

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
        self._save_catalog_btn.setEnabled(True)
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
        self._save_catalog_btn.setEnabled(n > 0)
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
        self._save_catalog_btn.setEnabled(False)
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
                if art is not None:
                    try:
                        art.remove()
                    except ValueError:
                        pass
            lst.clear()

        # Clear edit-mode circles
        self._clear_edit_mode()

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

        # Clean up arc_ncc artists
        self._exit_arc_ncc_mode()
        if self._state == 'ARC_NCC':
            self._state = 'IDLE'
            self._update_state_label()

    def _redraw_all_markers(self):
        """Redraw confirmed-pair markers.
        Pairs in edit mode show a numbered cyan circle on the target instead of red cross."""
        if self._aligner is None:
            return

        from matplotlib.patches import Circle

        # Remove old confirmed artists
        for art_list in (self._ref_crosses, self._ref_labels,
                          self._tgt_crosses, self._tgt_labels):
            for art in art_list:
                if art is not None:
                    try:
                        art.remove()
                    except Exception:
                        pass
            art_list.clear()

        # Remove old edit circles (will be re-created below)
        for val in self._edit_mode_circles.values():
            if isinstance(val, tuple):
                for a in val:
                    try:
                        a.remove()
                    except Exception:
                        pass

        # Compute target box in pixels once for circle radii
        tgt_frame  = self._aligner.targets[0]
        tgt_box_px = _box_to_pixels(self._box, tgt_frame,
                                    tgt_frame.image2d if tgt_frame.image2d is not None
                                    else self._tgt_img)

        for i, (rx, ry, ra, dec, tx, ty) in enumerate(self._aligner.pairs):
            n = i + 1   # 1-based label
            art_r, = self._ref_ax.plot(rx, ry, 'r+', markersize=14, mew=2,
                                        linestyle='none', zorder=10)
            lbl_r = self._ref_ax.text(rx + 3, ry + 3, str(n),
                                       color='red', fontsize=8, zorder=11)
            self._ref_crosses.append(art_r)
            self._ref_labels.append(lbl_r)

            if i in self._edit_mode_circles:
                # Draw cyan circle instead of red cross on target
                circ = Circle((tx, ty), radius=tgt_box_px,
                              edgecolor='cyan', facecolor='none',
                              linestyle='--', linewidth=1.5, zorder=9)
                self._tgt_ax.add_patch(circ)
                lbl_c = self._tgt_ax.text(tx + tgt_box_px + 2, ty + tgt_box_px + 2,
                                           str(n), color='cyan', fontsize=8, zorder=10)
                self._edit_mode_circles[i] = (circ, lbl_c)
                self._tgt_crosses.append(None)
                self._tgt_labels.append(None)
            else:
                art_t, = self._tgt_ax.plot(tx, ty, 'r+', markersize=14, mew=2,
                                            linestyle='none', zorder=10)
                lbl_t = self._tgt_ax.text(tx + 3, ty + 3, str(n),
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
        self._apply_to_btn.setEnabled(True)

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

    def _on_save_catalog(self):
        """Save current source pairs to a FITS catalog for later batch reuse."""
        if self._aligner is None or not self._aligner.pairs:
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Save Source Pairs Catalog", "sources.fits",
            filter="FITS files (*.fits *.fit);;All files (*)")
        if not path:
            return

        try:
            from rbcodes.rb_align.sources import _save_pairs_catalog
            _save_pairs_catalog(self._aligner.pairs, path)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))
            return

        QMessageBox.information(self, "Saved",
                                f"Catalog written to:\n{path}\n\n"
                                "Load it with strategy='batch' to skip re-detection.")

    def _on_apply_to(self):
        """Apply the corrected WCS to one or more external FITS files."""
        if self.corrected_wcs is None:
            return

        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select FITS file(s) to patch",
            filter="FITS files (*.fits *.fit *.fts);;All files (*)")
        if not paths:
            return

        try:
            from rbcodes.rb_align.io import _apply_wcs_to_files
            _apply_wcs_to_files(self.corrected_wcs,
                                paths,
                                [None] * len(paths),
                                suffix='_wcsfix')
        except Exception as e:
            QMessageBox.critical(self, "Apply failed", str(e))
            return

        written = '\n'.join(
            str(p).replace('.fits', '_wcsfix.fits').replace('.fit', '_wcsfix.fit')
            for p in paths)
        QMessageBox.information(self, "Done",
                                f"WCS patched and written to:\n{written}")

    # ------------------------------------------------------------------
    # Status / info labels
    # ------------------------------------------------------------------

    def _update_state_label(self):
        n_pending = (len(self._aligner.pairs) + 1 if self._aligner else 1)
        msgs = {
            'NO_REF':  "Load a reference image to begin.",
            'IDLE':    ("Left-click ref: add pair.  "
                        "Double-click any marker: toggle edit mode (cyan circle on target).  "
                        "Left-click near circle: re-place.  "
                        "Right-click: delete pair."),
            'PENDING': (f"Source {n_pending} pending — "
                        f"left-click on the target panel to confirm.  "
                        f"Right-click to cancel."),
            'ARC_NCC': ("arc_ncc mode — click-drag on the reference panel to draw a box"
                        " around a feature.  NCC runs on release.  "
                        "Draw another box to add a second pair."),
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
