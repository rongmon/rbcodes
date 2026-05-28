"""
IFU Viewer — MainWindow

Phase 3: sidebar + image panel (whitelight / plain 2D image).
"""
import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QStatusBar,
                              QHBoxLayout, QVBoxLayout, QLabel, QAction,
                              QFileDialog, QMessageBox, QSplitter,
                              QDialog, QFormLayout, QLineEdit,
                              QDialogButtonBox, QDoubleSpinBox, QSpinBox,
                              QCheckBox, QComboBox, QPushButton,
                              QRadioButton, QButtonGroup,
                              QGroupBox, QListWidget, QListWidgetItem,
                              QProgressBar, QProgressDialog,
                              QPlainTextEdit, QTextEdit)
from PyQt5.QtCore import Qt, QThread, pyqtSignal as _Signal
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavToolbar

from rbcodes.GUIs.ifuviewer.sidebar import DatasetSidebar
from rbcodes.GUIs.ifuviewer.image_panel import ImageCanvas
from rbcodes.GUIs.ifuviewer.image_controls import ImageControls
from rbcodes.GUIs.ifuviewer.channel_slider import ChannelSlider
from rbcodes.GUIs.ifuviewer.spectrum_panel import SpectrumCanvas
from rbcodes.GUIs.ifuviewer.io.image2d import FITSImage
from rbcodes.GUIs.ifuviewer.processing.cube_collapse import (
    build_whitelight, build_narrowband, build_continuum_sub)
from rbcodes.GUIs.ifuviewer.processing.moment_maps import (
    moment_map, subtract_linear_continuum, compute_snr_map)
from rbcodes.GUIs.ifuviewer.processing.aperture_extract import (
    extract_aperture, extract_variance_weighted,
    extract_with_method, subtract_background,
    make_circular_mask, make_annulus_mask)
from rbcodes.utils.rb_optimal_extract import extract_optimal_weighted
from rbcodes.GUIs.ifuviewer.aperture_controls import ApertureControls


class MainWindow(QMainWindow):
    """Top-level window for the IFU Viewer."""

    def __init__(self):
        super().__init__()
        self._active_cube          = None
        self._status_base          = ""   # dataset portion of status bar text
        self._specgui_windows      = []   # keep references to launched GUI windows
        self._extraction_marker_groups = []  # per-extraction list of image artists
        self._extraction_colors    = []   # parallel color list for highlighting
        self._extraction_marker_info = [] # parallel: how to recreate each marker on restore
        self._image_overlays       = []   # (reg_dict, color) for 2D-image-only overlays
        self._dataset_states       = {}   # id(cube) → saved per-dataset state
        self._sky_mask             = None  # 2-D bool array — sky region for SNR
        self._sky_rect             = None  # (x1, y1, x2, y2) pixel bounds
        self._ds9_bridge           = None  # DS9Bridge, created on connect
        self._ds9_frame_queue      = []   # list of (label, data2d, header, frame_no)

        self._build_ui()
        self._build_menu()

        self.setWindowTitle("IFU Viewer")
        self.resize(1200, 800)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        self.setStatusBar(QStatusBar(self))

        central = QWidget()
        self.setCentralWidget(central)

        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # Left: dataset sidebar + ds9 queue panel
        self._sidebar = DatasetSidebar()
        self._sidebar.dataset_changed.connect(self._on_dataset_changed)
        self._sidebar.use_as_variance.connect(self._on_use_as_variance)

        self._ds9_queue_panel = self._build_ds9_queue_panel()

        left_widget = QWidget()
        left_vbox   = QVBoxLayout(left_widget)
        left_vbox.setContentsMargins(0, 0, 0, 0)
        left_vbox.setSpacing(0)
        left_vbox.addWidget(self._sidebar, stretch=1)
        left_vbox.addWidget(self._ds9_queue_panel)
        root.addWidget(left_widget)

        # Right side — VBoxLayout added directly (no wrapper widget)
        right = QVBoxLayout()
        right.setContentsMargins(0, 0, 0, 0)
        right.setSpacing(0)

        self._image_canvas = ImageCanvas(dpi=100)
        self._image_canvas.setMinimumSize(600, 300)
        self._image_canvas.cursor_info.connect(self._on_cursor_info)
        self._image_canvas.spaxel_hovered.connect(self._on_spaxel_hovered)
        self._image_canvas.spaxel_locked.connect(self._on_spaxel_locked)
        self._image_canvas.aperture_drawn.connect(self._on_aperture_drawn)
        self._image_canvas.aperture_picked.connect(self._on_aperture_picked)
        self._image_canvas.crop_requested.connect(self._on_crop_selection)
        self._image_canvas.sky_region_requested.connect(self._on_set_sky_region)
        self._image_canvas.clear_sky_requested.connect(self._on_clear_sky_region)
        self._image_canvas.clear_overlays_requested.connect(self._on_clear_extractions)

        self._mpl_toolbar    = NavToolbar(self._image_canvas, self)
        self._image_canvas._nav_toolbar = self._mpl_toolbar  # for mode check
        self._image_controls = ImageControls(self._image_canvas)
        self._channel_slider = ChannelSlider()
        self._channel_slider.channel_changed.connect(self._on_channel_changed)
        self._channel_slider.mode_changed.connect(self._on_mode_changed)
        self._channel_slider.method_changed.connect(self._on_method_changed)

        self._spectrum_canvas = SpectrumCanvas()
        self._spectrum_canvas.setMinimumHeight(150)
        self._spectrum_canvas.cursor_info.connect(self._on_cursor_info)
        self._spectrum_canvas.onband_changed.connect(self._on_onband_changed)
        self._spectrum_canvas.contband_changed.connect(self._on_contband_changed)
        self._spectrum_canvas.window_cleared.connect(self._on_window_cleared)
        self._spectrum_canvas.xlim_changed.connect(self._sync_spec_range_strip)
        self._spectrum_canvas.ylim_changed.connect(self._sync_spec_range_strip)

        self._placeholder = QLabel("Load a cube to begin")
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setObjectName("placeholder")

        # Spectrum canvas + compact X/Y range strip below it
        self._spec_range_strip = self._build_spec_range_strip()
        _spec_wrapper = QWidget()
        _spec_vbox = QVBoxLayout(_spec_wrapper)
        _spec_vbox.setContentsMargins(0, 0, 0, 0)
        _spec_vbox.setSpacing(0)
        _spec_vbox.addWidget(self._spectrum_canvas, stretch=1)
        _spec_vbox.addWidget(self._spec_range_strip)

        # Splitter: image (top) ↕ spectrum+controls (bottom) — user-resizable
        self._splitter = QSplitter(Qt.Vertical)
        self._splitter.addWidget(self._image_canvas)
        self._splitter.addWidget(_spec_wrapper)
        self._splitter.setStretchFactor(0, 3)
        self._splitter.setStretchFactor(1, 1)
        self._splitter.setHandleWidth(5)

        # Aperture + extraction controls [Phase 9]
        self._aperture_controls = ApertureControls()
        self._aperture_controls.settings_changed.connect(self._on_aperture_settings_changed)
        self._aperture_controls.extract_requested.connect(self._on_circular_extract)
        self._aperture_controls.mode_changed.connect(self._on_aperture_mode_changed)
        self._extract_strip = self._build_extract_strip()

        # ds9 button strip [Phase 11]
        self._ds9_strip = self._build_ds9_strip()

        # Fixed widgets above/below the splitter
        right.addWidget(self._mpl_toolbar)
        right.addWidget(self._ds9_strip)
        right.addWidget(self._channel_slider)
        right.addWidget(self._placeholder)
        right.addWidget(self._splitter, stretch=1)
        right.addWidget(self._aperture_controls)
        right.addWidget(self._extract_strip)
        right.addWidget(self._image_controls)

        self._splitter.hide()
        self._mpl_toolbar.hide()
        self._image_controls.hide()
        self._channel_slider.hide()
        self._aperture_controls.hide()
        self._extract_strip.hide()

        root.addLayout(right, stretch=1)   # ← layout, not widget — fills correctly

    def _build_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("File")

        open_action = QAction("Open…", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._on_open)
        file_menu.addAction(open_action)

        load_reg_action = QAction("Load Region File…", self)
        load_reg_action.setToolTip(
            "Load a ds9 .reg region file from disk (no ds9 required).\n"
            "Converts regions to extraction masks.")
        load_reg_action.triggered.connect(self._on_load_region_file)
        file_menu.addAction(load_reg_action)

        file_menu.addSeparator()

        save_frame_action = QAction("Save Current Frame…", self)
        save_frame_action.setShortcut("Ctrl+S")
        save_frame_action.setToolTip(
            "Save the currently displayed 2D image (whitelight, narrowband, "
            "moment map, …) as a FITS file with spatial WCS header.")
        save_frame_action.triggered.connect(self._on_save_current_frame)
        file_menu.addAction(save_frame_action)

        save_subcube_action = QAction("Save Subcube…", self)
        save_subcube_action.setToolTip(
            "Save the active cube (or cropped subcube) as a 3D FITS file "
            "with updated WCS header.  Includes VAR extension if present.")
        save_subcube_action.triggered.connect(self._on_save_subcube)
        file_menu.addAction(save_subcube_action)

        save_reg_action = QAction("Save Apertures as Region File…", self)
        save_reg_action.setToolTip(
            "Save all current aperture markers as a ds9 .reg file.")
        save_reg_action.triggered.connect(self._on_save_apertures_as_region)
        file_menu.addAction(save_reg_action)

        save_figure_action = QAction("Save Figure…", self)
        save_figure_action.setShortcut("Ctrl+Shift+S")
        save_figure_action.setToolTip(
            "Save a publication-quality PNG or PDF figure.\n"
            "Cube: image panel + extracted spectra (color-matched).\n"
            "2D image: image + region overlays.\n"
            "Uses the exact colormap, clim, and axis limits currently shown.")
        save_figure_action.triggered.connect(self._on_save_figure)
        file_menu.addAction(save_figure_action)

        file_menu.addSeparator()

        quit_action = QAction("Quit", self)
        quit_action.setShortcut("Ctrl+Q")
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        view_menu = menubar.addMenu("View")

        spec_range_action = QAction("Set Spectrum Range…", self)
        spec_range_action.setShortcut("Ctrl+R")
        spec_range_action.setToolTip(
            "Set X (wavelength) range, Y (flux) range, and Y scale for the spectrum panel")
        spec_range_action.triggered.connect(self._on_set_spec_range)
        view_menu.addAction(spec_range_action)

        # Ctrl+W kept as alias for convenience
        spec_range_w = QAction("Set Spectrum Range (Ctrl+W alias)", self)
        spec_range_w.setShortcut("Ctrl+W")
        spec_range_w.setVisible(False)
        spec_range_w.triggered.connect(self._on_set_spec_range)
        self.addAction(spec_range_w)   # window-level, not in menu

        view_menu.addSeparator()

        self._coord_sex_action = QAction("Coordinates: Sexagesimal (hh:mm:ss)", self)
        self._coord_sex_action.setCheckable(True)
        self._coord_sex_action.setToolTip(
            "Toggle cursor RA/Dec display between decimal degrees and sexagesimal (hh:mm:ss / dd:mm:ss)")
        self._coord_sex_action.toggled.connect(self._on_toggle_coord_fmt)
        view_menu.addAction(self._coord_sex_action)

        view_menu.addSeparator()

        header_action = QAction("Show FITS Header…", self)
        header_action.setToolTip(
            "Show the FITS header of the active dataset.  "
            "If the file has multiple extensions you can choose which one to inspect.")
        header_action.triggered.connect(self._on_show_header)
        view_menu.addAction(header_action)

        view_menu.addSeparator()

        band_range_action = QAction("Set Band Ranges…", self)
        band_range_action.setShortcut("Ctrl+B")
        band_range_action.triggered.connect(self._on_set_band_ranges)
        view_menu.addAction(band_range_action)

        analysis_menu = menubar.addMenu("Analysis")

        moment_action = QAction("Moment Map…", self)
        moment_action.setShortcut("Ctrl+M")
        moment_action.setToolTip(
            "Compute a zeroth, first, or second moment map over a wavelength window.")
        moment_action.triggered.connect(self._on_moment_map)
        analysis_menu.addAction(moment_action)

        analysis_menu.addSeparator()

        crop_sel_action = QAction("Crop to Drawn Rectangle", self)
        crop_sel_action.setShortcut("Ctrl+K")
        crop_sel_action.setToolTip(
            "Crop the active cube/image to the last rectangle drawn on the image.\n"
            "Creates a new dataset in the sidebar — original is preserved.")
        crop_sel_action.triggered.connect(self._on_crop_selection)
        analysis_menu.addAction(crop_sel_action)

        crop_coord_action = QAction("Crop by Coordinates…", self)
        crop_coord_action.setToolTip(
            "Crop by entering pixel (X/Y) or sky (RA/Dec) coordinates.")
        crop_coord_action.triggered.connect(self._on_crop_by_coords)
        analysis_menu.addAction(crop_coord_action)

        analysis_menu.addSeparator()

        sky_action = QAction("Set Drawn Rectangle as Sky Region", self)
        sky_action.setToolTip(
            "Mark the last drawn rectangle as the sky/background region.\n"
            "Used for SNR-based moment map masking.")
        sky_action.triggered.connect(self._on_set_sky_region)
        analysis_menu.addAction(sky_action)

        clear_sky_action = QAction("Clear Sky Region", self)
        clear_sky_action.triggered.connect(self._on_clear_sky_region)
        analysis_menu.addAction(clear_sky_action)

        analysis_menu.addSeparator()

        batch_action = QAction("Batch Extract from Regions…", self)
        batch_action.setToolTip(
            "Extract 1D spectra from all regions in a .reg file or from ds9.")
        batch_action.triggered.connect(self._on_batch_extract)
        analysis_menu.addAction(batch_action)

        help_menu = menubar.addMenu("Help")

        help_action = QAction("Help", self)
        help_action.setShortcut("F1")
        help_action.setToolTip("Open the IFU Viewer help dialog")
        help_action.triggered.connect(self._on_show_help)
        help_menu.addAction(help_action)

    def _build_extract_strip(self):
        """Build the extraction list bar (always shown when cube loaded)."""
        strip = QWidget()
        layout = QHBoxLayout(strip)
        layout.setContentsMargins(6, 3, 6, 3)
        layout.setSpacing(6)

        layout.addWidget(QLabel("Extracted:"))

        self._extract_combo = QComboBox()
        self._extract_combo.setMinimumWidth(220)
        self._extract_combo.setToolTip("Select extracted spectrum  (Del = delete selected)")
        self._extract_combo.currentIndexChanged.connect(self._on_extract_selected)
        layout.addWidget(self._extract_combo)

        self._save_btn = QPushButton("Save…")
        self._save_btn.setFixedWidth(65)
        self._save_btn.setToolTip("Save selected spectrum to FITS/ASCII")
        self._save_btn.clicked.connect(self._on_save_extraction)
        layout.addWidget(self._save_btn)

        self._multispec_btn = QPushButton("→ rb_multispec")
        self._multispec_btn.setFixedWidth(105)
        self._multispec_btn.setToolTip("Send all extracted spectra to rb_multispec viewer")
        self._multispec_btn.clicked.connect(self._on_send_to_multispec)
        layout.addWidget(self._multispec_btn)

        self._del_btn = QPushButton("Delete")
        self._del_btn.setFixedWidth(62)
        self._del_btn.setToolTip("Delete selected extraction and its image marker  [Del]")
        self._del_btn.clicked.connect(self._on_delete_extraction)
        layout.addWidget(self._del_btn)

        self._clear_btn = QPushButton("Clear all")
        self._clear_btn.setFixedWidth(72)
        self._clear_btn.setToolTip("Remove all extracted spectra and markers")
        self._clear_btn.clicked.connect(self._on_clear_extractions)
        layout.addWidget(self._clear_btn)

        layout.addStretch()
        return strip

    def _set_extract_strip_enabled(self, enabled):
        """Enable or disable the interactive extraction strip controls."""
        for w in (self._extract_combo, self._save_btn, self._multispec_btn,
                  self._del_btn, self._clear_btn):
            w.setEnabled(enabled)

    # ------------------------------------------------------------------
    # ds9 UI builders
    # ------------------------------------------------------------------

    def _build_ds9_strip(self):
        """Thin button strip above the image for ds9 operations."""
        strip = QWidget()
        strip.setMaximumHeight(32)
        layout = QHBoxLayout(strip)
        layout.setContentsMargins(4, 2, 4, 2)
        layout.setSpacing(6)

        self._ds9_connect_btn = QPushButton("Connect ds9")
        self._ds9_connect_btn.setFixedWidth(110)
        self._ds9_connect_btn.setToolTip("Connect to a running ds9 instance")
        self._ds9_connect_btn.clicked.connect(self._on_ds9_connect)
        layout.addWidget(self._ds9_connect_btn)

        self._ds9_send_btn = QPushButton("→ Add to ds9")
        self._ds9_send_btn.setFixedWidth(100)
        self._ds9_send_btn.setToolTip("Add the current image to the ds9 frame queue")
        self._ds9_send_btn.setEnabled(False)
        self._ds9_send_btn.clicked.connect(self._on_ds9_add_to_queue)
        layout.addWidget(self._ds9_send_btn)

        self._ds9_import_btn = QPushButton("← Import regions")
        self._ds9_import_btn.setFixedWidth(120)
        self._ds9_import_btn.setToolTip("Import regions from the running ds9 instance")
        self._ds9_import_btn.setEnabled(False)
        self._ds9_import_btn.clicked.connect(self._on_ds9_import_regions)
        layout.addWidget(self._ds9_import_btn)

        self._ds9_push_btn = QPushButton("→ Push regions")
        self._ds9_push_btn.setFixedWidth(110)
        self._ds9_push_btn.setToolTip(
            "Push current apertures to ds9 as regions (adds to ds9, does not replace)")
        self._ds9_push_btn.setEnabled(False)
        self._ds9_push_btn.clicked.connect(self._on_push_regions_to_ds9)
        layout.addWidget(self._ds9_push_btn)

        layout.addStretch()

        self._ds9_status_label = QLabel("ds9: not connected")
        self._ds9_status_label.setStyleSheet("color: #585b70; font-size: 8pt;")
        layout.addWidget(self._ds9_status_label)

        layout.addSpacing(12)
        help_btn = QPushButton("?")
        help_btn.setFixedWidth(26)
        help_btn.setToolTip("Open help  (F1)")
        help_btn.clicked.connect(self._on_show_help)
        layout.addWidget(help_btn)

        return strip

    def _build_ds9_queue_panel(self):
        """Frame queue panel shown below the sidebar when ds9 is connected."""
        panel = QWidget()
        panel.setMaximumHeight(200)
        vbox = QVBoxLayout(panel)
        vbox.setContentsMargins(4, 4, 4, 4)
        vbox.setSpacing(4)

        hdr = QHBoxLayout()
        hdr.addWidget(QLabel("ds9 Frames"))
        hdr.addStretch()
        vbox.addLayout(hdr)

        self._ds9_queue_list = QListWidget()
        self._ds9_queue_list.setMaximumHeight(100)
        self._ds9_queue_list.setToolTip("Images queued to send to ds9")
        vbox.addWidget(self._ds9_queue_list)

        btn_row = QHBoxLayout()
        add_btn = QPushButton("+ Add current")
        add_btn.setFixedWidth(90)
        add_btn.setToolTip("Add the current image to the queue")
        add_btn.clicked.connect(self._on_ds9_add_to_queue)

        send_btn = QPushButton("Send all")
        send_btn.setFixedWidth(70)
        send_btn.setToolTip("Send all queued images to ds9")
        send_btn.clicked.connect(self._on_ds9_send_all)

        match_btn = QPushButton("Match WCS")
        match_btn.setFixedWidth(80)
        match_btn.setToolTip("Lock all ds9 frames to the same WCS")
        match_btn.clicked.connect(self._on_ds9_match_wcs)

        del_btn = QPushButton("×")
        del_btn.setFixedWidth(28)
        del_btn.setToolTip("Remove selected entry from queue")
        del_btn.clicked.connect(self._on_ds9_queue_delete)

        btn_row.addWidget(add_btn)
        btn_row.addWidget(send_btn)
        btn_row.addWidget(match_btn)
        btn_row.addStretch()
        btn_row.addWidget(del_btn)
        vbox.addLayout(btn_row)

        panel.hide()   # hidden until connected
        return panel

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_open(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Open FITS file(s)",
            "",
            "FITS files (*.fits *.fit *.fits.gz *.fit.gz);;All files (*)",
        )
        for p in paths:
            try:
                self._sidebar.add_dataset(p)
            except RuntimeError as exc:
                QMessageBox.warning(self, "Load error", str(exc))

    # ------------------------------------------------------------------
    # Per-dataset state save / restore
    # ------------------------------------------------------------------

    def _save_dataset_state(self):
        """Snapshot current extraction + image state for the active dataset."""
        cube = self._active_cube
        if cube is None:
            return
        locked = []
        lines  = self._spectrum_canvas._locked_lines
        specs  = self._spectrum_canvas.locked_spectra   # list of (label, rb_spec)
        for i, (label, rb_spec) in enumerate(specs):
            color  = (self._extraction_colors[i]
                      if i < len(self._extraction_colors) else '#f9e2af')
            minfo  = (self._extraction_marker_info[i]
                      if i < len(self._extraction_marker_info) else {})
            wave_d = lines[i].get_xdata() if i < len(lines) else None
            flux_d = lines[i].get_ydata() if i < len(lines) else None
            locked.append({'label': label, 'color': color,
                           'wave': wave_d, 'flux': flux_d,
                           'rb_spec': rb_spec, 'marker_info': minfo})
        im           = getattr(self._image_canvas, '_im', None)
        norm         = im.norm     if im is not None else None
        clim         = im.get_clim() if im is not None else (None, None)
        display_state = self._image_controls.get_display_state()
        sc   = self._spectrum_canvas
        # Snapshot MomentMapDialog class-level params so they restore per-dataset
        mm = MomentMapDialog
        moment_params = {
            'order':    mm._last_order,
            'wmin':     mm._last_wmin,
            'wmax':     mm._last_wmax,
            'lrest':    mm._last_lrest,
            'cont_sub': mm._last_cont_sub,
            'bcont':    mm._last_bcont,
            'rcont':    mm._last_rcont,
            'snr':      mm._last_snr,
            'snr_thr':  mm._last_snr_thr,
            'snr_meth': mm._last_snr_meth,
        }
        self._dataset_states[id(cube)] = {
            'image_data':     self._image_canvas._data_raw,
            'image_norm':     norm,
            'image_clim':     clim,
            'display_state':  display_state,
            'slider_mode':    self._channel_slider.current_mode,
            'locked':         locked,
            'color_idx':      sc._color_idx,
            'sky_mask':       self._sky_mask,
            'sky_rect':       self._sky_rect,
            'onband_ext':     sc._onband_ext,
            'cont1_ext':      sc._cont1_ext,
            'cont2_ext':      sc._cont2_ext,
            'spec_mode':      sc._mode,
            'spec_xlim':      sc._ax.get_xlim(),
            'spec_ylim':      sc._ax.get_ylim(),
            'aperture_state': self._aperture_controls.get_state(),
            'moment_params':  moment_params,
            'image_overlays': list(self._image_overlays),
        }

    def _restore_dataset_state(self, cube):
        """Restore a previously saved state for *cube*. Returns True if found."""
        state = self._dataset_states.get(id(cube))
        if not state:
            return False
        # Override image with saved one, then restore full display state
        img = state.get('image_data')
        if img is not None:
            import matplotlib.colors as mcolors
            hdr  = cube.spatial_header() if cube is not None else None
            # Rebuild norm from saved clim so vmin/vmax survive the round-trip
            clim = state.get('image_clim', (None, None))
            ds   = state.get('display_state', {})
            cmap = ds.get('cmap', 'gray')
            norm = None
            if clim[0] is not None and clim[1] is not None:
                norm = mcolors.Normalize(vmin=clim[0], vmax=clim[1])
            self._image_canvas.show_image(img, header=hdr, cmap=cmap, norm=norm)
            # show_image calls ax.cla() which removes all artists; sync the marks list
            self._image_canvas._extraction_marks.clear()
            # Restore scale/norm/cmap controls (without triggering a reset)
            if ds:
                self._image_controls.set_display_state(ds)
        # Restore channel slider button (Narrowband / Cont-sub / etc.)
        slider_mode = state.get('slider_mode', 'Whitelight')
        btns = self._channel_slider._mode_btns
        if slider_mode in btns:
            btns[slider_mode].setChecked(True)
            # Enable slider for band modes (without triggering a new image build)
            self._channel_slider._set_slider_enabled(
                slider_mode in ('Channel', 'Narrowband', 'Cont-sub'))
        # Sky region
        self._sky_mask = state.get('sky_mask')
        self._sky_rect = state.get('sky_rect')
        # Color cycle position
        sc = self._spectrum_canvas
        sc._color_idx = state.get('color_idx', 0)
        # Re-lock spectra and redraw markers
        for entry in state.get('locked', []):
            wave, flux = entry.get('wave'), entry.get('flux')
            if wave is None or flux is None:
                continue
            color   = entry['color']
            label   = entry['label']
            minfo   = entry.get('marker_info', {})
            rb_spec = entry.get('rb_spec')
            sc.lock_spectrum(wave, flux, None, label, rb_spec=rb_spec, color=color)
            new_markers = self._redraw_marker(minfo, color)
            self._extract_combo.addItem(label)
            self._extraction_marker_groups.append(list(new_markers))
            self._extraction_colors.append(color)
            self._extraction_marker_info.append(minfo)
        if self._extract_combo.count() > 0:
            self._extract_combo.setCurrentIndex(0)
            self._set_extract_strip_enabled(True)
        # Restore spectrum band spans
        onband = state.get('onband_ext', (np.nan, np.nan))
        cont1  = state.get('cont1_ext',  (np.nan, np.nan))
        cont2  = state.get('cont2_ext',  (np.nan, np.nan))
        spec_mode = state.get('spec_mode', 'Whitelight')
        sc.set_mode(spec_mode)
        sc.apply_band_ranges(onband, cont1, cont2)
        # Restore y-limits (autoscale if spectra were re-locked)
        ylim = state.get('spec_ylim')
        if ylim is not None:
            sc._ax.set_ylim(ylim)
        if state.get('locked'):
            sc.autoscale_y()
        # Restore spectrum x-limits
        xlim = state.get('spec_xlim')
        if xlim is not None:
            sc._ax.set_xlim(xlim)
        # Restore aperture controls state
        ap_state = state.get('aperture_state')
        if ap_state:
            self._aperture_controls.set_state(ap_state)
        # Redraw sky region overlay if one was set
        sky_rect = state.get('sky_rect')
        if sky_rect is not None:
            self._image_canvas.draw_sky_region(*sky_rect)
        # Redraw 2D-image region overlays
        self._image_overlays = list(state.get('image_overlays', []))
        for reg, color in self._image_overlays:
            self._image_canvas.draw_region_shape(reg, cube.wcs, color)
        if self._image_overlays and isinstance(cube, FITSImage):
            self._clear_btn.setEnabled(True)
        # Restore MomentMapDialog class-level params for this dataset
        mp = state.get('moment_params')
        if mp:
            mm = MomentMapDialog
            mm._last_order    = mp.get('order',    mm._last_order)
            mm._last_wmin     = mp.get('wmin',     mm._last_wmin)
            mm._last_wmax     = mp.get('wmax',     mm._last_wmax)
            mm._last_lrest    = mp.get('lrest',    mm._last_lrest)
            mm._last_cont_sub = mp.get('cont_sub', mm._last_cont_sub)
            mm._last_bcont    = mp.get('bcont',    mm._last_bcont)
            mm._last_rcont    = mp.get('rcont',    mm._last_rcont)
            mm._last_snr      = mp.get('snr',      mm._last_snr)
            mm._last_snr_thr  = mp.get('snr_thr',  mm._last_snr_thr)
            mm._last_snr_meth = mp.get('snr_meth', mm._last_snr_meth)
        self._image_canvas.draw_idle()
        return True

    def _redraw_marker(self, marker_info, color):
        """Recreate image markers from stored params; returns list of artists."""
        if not marker_info:
            return []
        mtype  = marker_info.get('type')
        canvas = self._image_canvas
        cube   = self._active_cube
        n_before = len(canvas._extraction_marks)
        try:
            if mtype == 'spaxel':
                mk, = canvas._ax.plot(marker_info['x'], marker_info['y'], '+',
                                      color=color, ms=9, mew=1.4, zorder=12, picker=5)
                canvas._extraction_marks.append(mk)
            elif mtype == 'circle':
                canvas.draw_aperture_marker(
                    marker_info['cx'], marker_info['cy'], 'Circle',
                    marker_info['radius'],
                    marker_info.get('bg_inner'), marker_info.get('bg_outer'),
                    color=color)
            elif mtype == 'rect':
                from matplotlib.patches import Rectangle as _Rect
                xi0, yi0 = marker_info['xi0'], marker_info['yi0']
                xi1, yi1 = marker_info['xi1'], marker_info['yi1']
                rect = _Rect((xi0 - 0.5, yi0 - 0.5), xi1 - xi0, yi1 - yi0,
                             linewidth=1.2, edgecolor=color, facecolor='none',
                             linestyle='--', zorder=11, picker=5)
                canvas._ax.add_patch(rect)
                canvas._extraction_marks.append(rect)
            elif mtype == 'region':
                reg = marker_info.get('region')
                if reg and cube:
                    canvas.draw_region_shape(reg, cube.wcs, color)
        except Exception as exc:
            print(f"[restore] marker redraw failed ({mtype}): {exc}")
        return canvas._extraction_marks[n_before:]

    # ------------------------------------------------------------------

    def _on_dataset_changed(self, cube):
        # Save current dataset's state before switching
        self._save_dataset_state()
        # Clear extraction UI (spectrum lines, image marks, combo)
        self._on_clear_extractions(_from_switch=True)

        self._active_cube = cube

        if cube is None:
            self.statusBar().clearMessage()
            self._splitter.hide()
            self._mpl_toolbar.hide()
            self._image_controls.hide()
            self._channel_slider.hide()
            self._placeholder.show()
            self._placeholder.setText("Load a cube to begin")
            self._channel_slider.set_cube(None)
            self._spectrum_canvas.clear()
            self._aperture_controls.hide()
            self._extract_strip.hide()
            return

        self._status_base = _status_text(cube)
        self.statusBar().showMessage(self._status_base)
        self._show_default_image(cube)
        # Restore state if returning to a previously viewed dataset
        self._restore_dataset_state(cube)

    def _on_use_as_variance(self, var_cube, target_cube):
        """Right-click sidebar: assign var_cube.flux as variance for target_cube."""
        if var_cube.flux.shape != target_cube.flux.shape:
            QMessageBox.warning(
                self, "Shape mismatch",
                f"Variance cube shape {var_cube.flux.shape} does not match "
                f"flux cube shape {target_cube.flux.shape}.")
            return
        target_cube.var = var_cube.flux
        self.statusBar().showMessage(
            f"Variance for '{target_cube.name}' set from '{var_cube.name}'.", 4000)

    def _on_save_current_frame(self):
        """Save the currently displayed 2D image as a FITS file (Ctrl+S)."""
        from astropy.io import fits as _fits

        data2d = self._image_canvas._data_raw
        if data2d is None:
            self.statusBar().showMessage("No image to save.", 3000)
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Save Current Frame", "",
            "FITS files (*.fits);;All files (*)")
        if not path:
            return

        # Build header: spatial WCS from active cube if available
        cube = self._active_cube
        if cube is not None and not isinstance(cube, FITSImage):
            try:
                hdr = cube.spatial_header()
            except Exception:
                hdr = _fits.Header()
        elif cube is not None:
            try:
                hdr = cube.spatial_header()
            except Exception:
                hdr = _fits.Header()
        else:
            hdr = _fits.Header()

        # Tag what this frame is
        mode = self._channel_slider.current_mode
        hdr['IFUVMODE'] = (mode, 'IFU viewer display mode when saved')

        try:
            _fits.writeto(path, data2d.astype('float32'), header=hdr, overwrite=True)
            self.statusBar().showMessage(f"Saved frame to {path}", 4000)
        except Exception as exc:
            QMessageBox.warning(self, "Save failed", str(exc))

    def _on_save_figure(self):
        """
        Save a publication-quality PNG or PDF figure (Ctrl+Shift+S).

        Cube with locked spectra  → 2-panel: image top, extracted spectra bottom.
        Cube with no extractions  → 1-panel: image only.
        2D image                  → 1-panel: image + region overlays.

        Exact colormap, clim, and axis limits from the current display are used.
        """
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.patches import Circle, Rectangle as MplRect

        cube = self._active_cube
        if cube is None:
            self.statusBar().showMessage("No data loaded.", 3000)
            return

        # ---- ask for path ------------------------------------------------
        default = (cube.name or 'figure') + '.pdf'
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", default,
            "PDF (*.pdf);;PNG (*.png);;All files (*)")
        if not path:
            return

        # ---- gather current display state --------------------------------
        ic   = self._image_canvas
        sc   = self._spectrum_canvas
        im   = ic._im
        if im is None:
            self.statusBar().showMessage("No image to save.", 3000)
            return

        cmap_name  = im.get_cmap().name
        vmin, vmax = im.get_clim()
        img_xlim   = ic._ax.get_xlim()
        img_ylim   = ic._ax.get_ylim()
        wcs        = ic._wcs
        data2d     = ic._data_raw

        is_image   = isinstance(cube, FITSImage)
        locked     = sc._locked_lines          # list of Line2D
        locked_specs = sc.locked_spectra       # list of (label, rb_spec)
        has_spectra  = len(locked) > 0 and not is_image

        onband = sc._onband_ext
        cont1  = sc._cont1_ext
        cont2  = sc._cont2_ext
        spec_xlim = sc._ax.get_xlim()
        spec_ylim = sc._ax.get_ylim()

        # ---- build figure ------------------------------------------------
        if has_spectra:
            fig = plt.figure(figsize=(8, 6), constrained_layout=True)
            gs  = fig.add_gridspec(2, 1, height_ratios=[2.5, 1])
            ax_im   = fig.add_subplot(gs[0], projection=wcs) if wcs else fig.add_subplot(gs[0])
            ax_spec = fig.add_subplot(gs[1])
        else:
            fig_h = 5 if is_image else 4.5
            fig   = plt.figure(figsize=(6, fig_h), constrained_layout=True)
            ax_im = fig.add_subplot(111, projection=wcs) if wcs else fig.add_subplot(111)
            ax_spec = None

        # ---- image panel -------------------------------------------------
        ax_im.imshow(data2d, origin='lower', interpolation='nearest',
                     cmap=cmap_name, vmin=vmin, vmax=vmax, aspect='auto')
        ax_im.set_xlim(img_xlim)
        ax_im.set_ylim(img_ylim)

        # colorbar
        sm = plt.cm.ScalarMappable(
            cmap=cmap_name,
            norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        fig.colorbar(sm, ax=ax_im, pad=0.02, fraction=0.046)

        # axis labels
        if wcs:
            ax_im.set_xlabel('RA', fontsize=9)
            ax_im.set_ylabel('Dec', fontsize=9)
        else:
            ax_im.set_xlabel('x (px)', fontsize=9)
            ax_im.set_ylabel('y (px)', fontsize=9)

        # aperture / region markers
        for i, minfo in enumerate(self._extraction_marker_info):
            color = (self._extraction_colors[i]
                     if i < len(self._extraction_colors) else '#888888')
            mtype = minfo.get('type')
            if mtype == 'spaxel':
                ax_im.plot(minfo['x'], minfo['y'], '+',
                           color=color, ms=9, mew=1.4, zorder=12)
            elif mtype == 'circle':
                cx, cy, r = minfo['cx'], minfo['cy'], minfo['radius']
                ax_im.add_patch(Circle((cx, cy), r,
                                       edgecolor=color, facecolor='none',
                                       lw=1.2, zorder=11))
                ax_im.plot(cx, cy, '+', color=color, ms=6, mew=1.2, zorder=12)
                bg_in  = minfo.get('bg_inner')
                bg_out = minfo.get('bg_outer')
                if bg_in and bg_out:
                    ax_im.add_patch(Circle((cx, cy), bg_in,
                                           edgecolor=color, facecolor='none',
                                           lw=0.8, linestyle='--', zorder=11))
                    ax_im.add_patch(Circle((cx, cy), bg_out,
                                           edgecolor=color, facecolor='none',
                                           lw=0.8, linestyle='--', zorder=11))
            elif mtype == 'rect':
                xi0, yi0 = minfo['xi0'], minfo['yi0']
                xi1, yi1 = minfo['xi1'], minfo['yi1']
                ax_im.add_patch(MplRect(
                    (xi0 - 0.5, yi0 - 0.5), xi1 - xi0, yi1 - yi0,
                    edgecolor=color, facecolor='none',
                    lw=1.2, linestyle='--', zorder=11))

        # 2D-image region overlays
        if is_image:
            for reg, color in self._image_overlays:
                try:
                    pixel_region = reg.to_pixel(wcs) if wcs else None
                    if pixel_region is not None:
                        patch = pixel_region.as_artist(
                            edgecolor=color, facecolor='none', lw=1.2)
                        ax_im.add_patch(patch)
                except Exception:
                    pass

        # title
        mode_label = self._channel_slider.current_mode if not is_image else '2D Image'
        ax_im.set_title(f"{cube.name}  —  {mode_label}", fontsize=9, pad=4)

        # ---- spectrum panel (cube with extractions only) -----------------
        if ax_spec is not None:
            for line, (label, _) in zip(locked, locked_specs):
                wave = line.get_xdata()
                flux = line.get_ydata()
                ax_spec.plot(wave, flux,
                             color=line.get_color(), lw=0.9, alpha=0.9,
                             label=label)

            # band spans
            if not (np.isnan(onband[0]) or np.isnan(onband[1])):
                ax_spec.axvspan(onband[0], onband[1],
                                color='#40a040', alpha=0.15, zorder=0)
            if not (np.isnan(cont1[0]) or np.isnan(cont1[1])):
                ax_spec.axvspan(cont1[0], cont1[1],
                                color='#c04040', alpha=0.15, zorder=0)
            if not (np.isnan(cont2[0]) or np.isnan(cont2[1])):
                ax_spec.axvspan(cont2[0], cont2[1],
                                color='#c08020', alpha=0.15, zorder=0)

            ax_spec.set_xlim(spec_xlim)
            ax_spec.set_ylim(spec_ylim)
            ax_spec.set_xlabel('Wavelength (Å)', fontsize=9)
            ax_spec.set_ylabel('Flux', fontsize=9)
            ax_spec.tick_params(labelsize=8)
            if locked:
                ncol = max(1, min(4, len(locked) // 4 + 1))
                ax_spec.legend(fontsize=6, framealpha=0.8,
                               loc='upper right', ncol=ncol)

        # ---- save --------------------------------------------------------
        try:
            dpi = 150 if path.lower().endswith('.png') else None
            fig.savefig(path, dpi=dpi)
            plt.close(fig)
            self.statusBar().showMessage(f"Figure saved to {path}", 4000)
        except Exception as exc:
            plt.close(fig)
            QMessageBox.warning(self, "Save failed", str(exc))

    def _on_save_subcube(self):
        """Save the active 3D cube (or cropped subcube) as a FITS file."""
        from astropy.io import fits as _fits

        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            QMessageBox.information(self, "No cube",
                                    "Load a cube before saving a subcube.")
            return

        default_name = cube.name + '.fits'
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Subcube", default_name,
            "FITS files (*.fits);;All files (*)")
        if not path:
            return

        try:
            hdr = cube.header.copy() if cube.header is not None else _fits.Header()
            primary = _fits.PrimaryHDU(cube.flux.astype('float32'), header=hdr)
            hdul = _fits.HDUList([primary])
            if cube.var is not None:
                var_hdu = _fits.ImageHDU(cube.var.astype('float32'), name='VAR')
                hdul.append(var_hdu)
            hdul.writeto(path, overwrite=True)
            var_note = ' + VAR' if cube.var is not None else ''
            self.statusBar().showMessage(
                f"Saved subcube{var_note} to {path}  "
                f"({cube.flux.shape[0]} channels, {cube.ny}×{cube.nx} spaxels)", 6000)
        except Exception as exc:
            QMessageBox.warning(self, "Save failed", str(exc))

    # ------------------------------------------------------------------
    # Crop
    # ------------------------------------------------------------------

    def _on_crop_selection(self):
        """Crop active cube/image to the last drawn rectangle (Ctrl+K)."""
        rect = self._image_canvas._last_rect
        if rect is None:
            self.statusBar().showMessage(
                "Draw a rectangle on the image first, then crop.", 4000)
            return
        self._do_crop(*rect)

    def _on_crop_by_coords(self):
        """Open the coordinate-input crop dialog."""
        cube = self._active_cube
        if cube is None:
            self.statusBar().showMessage("Load a dataset first.", 3000)
            return
        ny = cube.ny if hasattr(cube, 'ny') else cube.flux.shape[0]
        nx = cube.nx if hasattr(cube, 'nx') else cube.flux.shape[1]
        wcs = cube.wcs
        dlg = CropDialog(nx, ny, wcs, parent=self)
        if not dlg.exec_():
            return
        x1, y1, x2, y2 = dlg.pixel_bounds()
        self._do_crop(x1, y1, x2, y2)

    def _do_crop(self, x1, y1, x2, y2):
        """Perform the crop and add the result to the sidebar."""
        cube = self._active_cube
        if cube is None:
            return
        ny = cube.flux.shape[-2]
        nx = cube.flux.shape[-1]
        x1 = max(0, x1); y1 = max(0, y1)
        x2 = min(nx, x2); y2 = min(ny, y2)
        if x2 <= x1 or y2 <= y1:
            QMessageBox.warning(self, "Crop failed", "Crop region is empty.")
            return
        try:
            cropped = cube.crop(x1, y1, x2, y2)
        except Exception as exc:
            QMessageBox.warning(self, "Crop failed", str(exc))
            return
        self._sidebar.add_cube_object(cropped)
        self.statusBar().showMessage(
            f"Cropped '{cube.name}' [{x1}:{x2}, {y1}:{y2}] → '{cropped.name}'  "
            f"({x2-x1}×{y2-y1} spaxels)", 6000)

    # ------------------------------------------------------------------
    # Sky region
    # ------------------------------------------------------------------

    def _on_set_sky_region(self):
        """Mark the last drawn rectangle as the sky/background region."""
        rect = self._image_canvas._last_rect
        if rect is None:
            self.statusBar().showMessage(
                "Draw a rectangle on the image first, then set as sky region.", 4000)
            return
        cube = self._active_cube
        if cube is None:
            return
        x1, y1, x2, y2 = rect
        ny = cube.flux.shape[-2]
        nx = cube.flux.shape[-1]
        mask = np.zeros((ny, nx), dtype=bool)
        mask[y1:y2, x1:x2] = True
        self._sky_mask = mask
        self._sky_rect = rect
        self._image_canvas.draw_sky_region(x1, y1, x2, y2)
        self.statusBar().showMessage(
            f"Sky region set: x=[{x1},{x2}) y=[{y1},{y2})  "
            f"({mask.sum()} spaxels)", 5000)

    def _on_clear_sky_region(self):
        """Remove the sky region."""
        self._sky_mask = None
        self._sky_rect = None
        self._image_canvas.clear_sky_region()
        self.statusBar().showMessage("Sky region cleared.", 3000)

    # ------------------------------------------------------------------
    # Aperture highlighting
    # ------------------------------------------------------------------

    def _on_extract_selected(self, idx):
        """Highlight the selected extraction's image marker, dim others."""
        if idx < 0:
            return
        n = len(self._extraction_marker_groups)
        for i, artists in enumerate(self._extraction_marker_groups):
            selected = (i == idx)
            for artist in artists:
                try:
                    artist.set_linewidth(3.5 if selected else 2.0)
                except AttributeError:
                    pass
        self._image_canvas.draw_idle()

    def _on_aperture_picked(self, artist):
        """Click on an aperture marker → select the corresponding extraction."""
        for idx, artists in enumerate(self._extraction_marker_groups):
            if artist in artists:
                self._extract_combo.setCurrentIndex(idx)
                return

    # ------------------------------------------------------------------
    # ds9 slots
    # ------------------------------------------------------------------

    def _on_ds9_connect(self):
        """Connect to (or disconnect from) ds9."""
        from rbcodes.GUIs.ifuviewer.ds9.bridge import DS9Bridge

        if self._ds9_bridge is not None and self._ds9_bridge.available:
            # Disconnect
            self._ds9_bridge.disconnect()
            self._ds9_bridge = None
            self._ds9_connect_btn.setText("Connect ds9")
            self._ds9_status_label.setText("ds9: not connected")
            self._ds9_status_label.setStyleSheet("color: #585b70; font-size: 8pt;")
            self._ds9_send_btn.setEnabled(False)
            self._ds9_import_btn.setEnabled(False)
            self._ds9_push_btn.setEnabled(False)
            self._ds9_queue_panel.hide()
            return

        bridge = DS9Bridge()
        ok, reason = bridge.connect()

        if not ok:
            if reason == 'no_pyds9':
                QMessageBox.warning(
                    self, "pyds9 not found",
                    "pyds9 is not installed.\n\nInstall with:\n    pip install pyds9")
            else:
                QMessageBox.warning(
                    self, "ds9 not running",
                    "No running ds9 instance found.\n\nStart ds9 first, then connect.")
            return

        self._ds9_bridge = bridge
        self._ds9_connect_btn.setText("Disconnect ds9")
        self._ds9_status_label.setText("ds9: connected ●")
        self._ds9_status_label.setStyleSheet(
            "color: #a6e3a1; font-size: 8pt; font-weight: bold;")
        self._ds9_send_btn.setEnabled(True)
        self._ds9_import_btn.setEnabled(True)
        self._ds9_push_btn.setEnabled(True)
        self._ds9_queue_panel.show()
        self.statusBar().showMessage("ds9 connected.", 3000)

    def _on_ds9_add_to_queue(self):
        """Add the current image panel data to the ds9 frame queue."""
        data2d = self._image_canvas._data_raw
        if data2d is None:
            self.statusBar().showMessage("No image to add.", 3000)
            return

        cube = self._active_cube
        header = None
        if cube is not None:
            try:
                header = cube.spatial_header()
            except Exception:
                pass

        # Auto-name from current mode
        mode = getattr(self._channel_slider, 'current_mode', 'Image')
        sc = self._spectrum_canvas
        if mode == 'Narrowband' and np.isfinite(sc._onband_ext[0]):
            w0, w1 = sc._onband_ext
            label = f"NB {w0:.0f}–{w1:.0f} Å"
        elif mode == 'Cont-sub' and np.isfinite(sc._onband_ext[0]):
            w0, w1 = sc._onband_ext
            label = f"Cont-sub {w0:.0f}–{w1:.0f} Å"
        else:
            label = mode

        frame_no = len(self._ds9_frame_queue) + 1
        self._ds9_frame_queue.append((label, data2d.copy(), header, frame_no))
        item = QListWidgetItem(f"[{frame_no}]  {label}")
        self._ds9_queue_list.addItem(item)

    def _on_ds9_queue_delete(self):
        row = self._ds9_queue_list.currentRow()
        if row < 0:
            return
        self._ds9_queue_list.takeItem(row)
        if row < len(self._ds9_frame_queue):
            del self._ds9_frame_queue[row]

    def _on_ds9_send_all(self):
        """Send all queued images to ds9."""
        if not self._ds9_bridge or not self._ds9_bridge.available:
            self.statusBar().showMessage("Not connected to ds9.", 3000)
            return
        if not self._ds9_frame_queue:
            self.statusBar().showMessage("Frame queue is empty.", 3000)
            return

        for label, data2d, header, frame_no in self._ds9_frame_queue:
            self._ds9_bridge.send_image(data2d, header, frame=frame_no)

        self.statusBar().showMessage(
            f"Sent {len(self._ds9_frame_queue)} frame(s) to ds9.", 4000)

    def _on_ds9_match_wcs(self):
        if self._ds9_bridge and self._ds9_bridge.available:
            self._ds9_bridge.match_wcs()
            self.statusBar().showMessage("ds9 frames WCS-matched.", 3000)

    def _on_push_regions_to_ds9(self):
        """Push all current apertures to the live ds9 instance as region overlays."""
        if not self._ds9_bridge or not self._ds9_bridge.available:
            return
        from rbcodes.GUIs.ifuviewer.processing.spatial_mask import regions_to_ds9_text

        cube  = self._active_cube
        wcs   = cube.wcs if cube is not None else None

        # Build aperture list from marker_info
        ap_list = []
        for i, minfo in enumerate(self._extraction_marker_info):
            label = ''
            specs = self._spectrum_canvas.locked_spectra
            if i < len(specs):
                label = specs[i][0]
            mtype = minfo.get('type', 'single')
            if mtype == 'spaxel':
                ap_list.append({'type': 'single', 'label': label,
                                'cx': minfo['x'], 'cy': minfo['y']})
            elif mtype == 'circle':
                ap_list.append({'type': 'circle', 'label': label,
                                'cx': minfo['cx'], 'cy': minfo['cy'],
                                'radius': minfo['radius'],
                                'bg_inner': minfo.get('bg_inner'),
                                'bg_outer': minfo.get('bg_outer')})
            elif mtype == 'annulus':
                ap_list.append({'type': 'annulus', 'label': label,
                                'cx': minfo['cx'], 'cy': minfo['cy'],
                                'bg_inner': minfo.get('bg_inner', minfo.get('radius', 3)),
                                'bg_outer': minfo.get('bg_outer', minfo.get('radius', 3) + 3)})
            elif mtype == 'rect':
                ap_list.append({'type': 'rect', 'label': label,
                                'x1': minfo['xi0'], 'y1': minfo['yi0'],
                                'x2': minfo['xi1'], 'y2': minfo['yi1']})
            # mtype == 'region': skip — these already came from ds9

        if not ap_list:
            QMessageBox.information(self, "Nothing to push",
                                    "No apertures have been defined.")
            return

        reg_text = regions_to_ds9_text(ap_list, wcs=wcs)
        ok = self._ds9_bridge.push_regions(reg_text)
        if ok:
            self.statusBar().showMessage(
                f"Pushed {len(ap_list)} aperture(s) to ds9.", 4000)
        else:
            QMessageBox.warning(self, "Push failed",
                                "Could not push regions to ds9.\n"
                                "Try saving the region file and loading it manually in ds9.")

    def _on_ds9_import_regions(self):
        """Import regions from the live ds9 instance."""
        if not self._ds9_bridge or not self._ds9_bridge.available:
            return
        if self._active_cube is None:
            QMessageBox.information(self, "No data", "Load a cube or image first.")
            return
        dlg = ImportRegionsDialog(source='ds9', parent=self)
        if not dlg.exec_():
            return
        selected_only = dlg.selected_only()
        reg_text = self._ds9_bridge.get_regions(selected_only=selected_only)
        if not reg_text.strip():
            QMessageBox.information(self, "No regions", "No regions found in ds9.")
            return
        self._extract_from_region_text(reg_text)

    # ------------------------------------------------------------------
    # Region file I/O
    # ------------------------------------------------------------------

    def _on_load_region_file(self):
        """Load a ds9 .reg file; extract spectra if cube, overlay only if 2D image."""
        if self._active_cube is None:
            QMessageBox.information(self, "No data", "Load a cube or image first.")
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Region File", "",
            "ds9 region files (*.reg);;All files (*)")
        if not path:
            return
        try:
            with open(path, 'r') as f:
                reg_text = f.read()
        except Exception as exc:
            QMessageBox.warning(self, "Load failed", str(exc))
            return

        dlg = ImportRegionsDialog(source='file', filename=path, parent=self)
        if not dlg.exec_():
            return

        self._extract_from_region_text(reg_text)

    def _extract_from_region_text(self, reg_text, method='sum', weighting='None'):
        """Parse region text; extract spectra if cube, overlay shapes only if 2D image."""
        from rbcodes.GUIs.ifuviewer.processing.spatial_mask import (
            parse_ds9_regions, region_to_mask, region_center_sky)

        cube = self._active_cube
        if cube is None:
            return

        is_image = isinstance(cube, FITSImage)

        regions = parse_ds9_regions(reg_text)
        if not regions:
            QMessageBox.information(self, "No regions",
                                    "No regions found in the file.")
            return

        # Color cycle for image-only overlays (no spectrum panel involvement)
        from rbcodes.GUIs.ifuviewer.spectrum_panel import _EXTRACT_COLORS
        _img_color_idx = 0

        errors = []
        n_ok   = 0
        for i, reg in enumerate(regions):
            try:
                extractable = reg.get('extract', True)

                if is_image or not extractable:
                    # 2D image OR annotation shape — draw shape overlay, no spectrum
                    color = _EXTRACT_COLORS[_img_color_idx % len(_EXTRACT_COLORS)]
                    _img_color_idx += 1
                    self._image_canvas.draw_region_shape(reg, cube.wcs, color)
                    self._image_overlays.append((reg, color))
                    if is_image:
                        self._clear_btn.setEnabled(True)
                    n_ok += 1
                    continue

                # IFU cube with extractable shape — extract spectrum
                mask = region_to_mask(reg, cube.wcs, cube.ny, cube.nx)
                if not mask.any():
                    continue

                flux_arr, err_arr = extract_optimal_weighted(
                    cube.flux, cube.var, mask, method='data')
                name  = reg.get('name') or f"region_{i+1:03d}"
                label = f"[{reg['shape']}] {name}"
                rb_spec = _make_rb_spectrum(cube.wave, flux_arr, err_arr, label)
                ra, dec = region_center_sky(reg, cube.wcs)
                if rb_spec is not None and ra is not None:
                    rb_spec.meta['ra']  = ra
                    rb_spec.meta['dec'] = dec
                    rb_spec.meta['object']   = cube.header.get('OBJECT', '') if cube.header else ''
                    rb_spec.meta['instrume'] = cube.header.get('INSTRUME', '') if cube.header else ''
                color = self._spectrum_canvas.lock_spectrum(
                    cube.wave, flux_arr, err_arr, label, rb_spec=rb_spec)
                new_markers = self._image_canvas.draw_region_shape(
                    reg, cube.wcs, color)
                if not new_markers:
                    yy, xx = np.mgrid[0:cube.ny, 0:cube.nx]
                    cx = int(np.round(np.mean(xx[mask])))
                    cy = int(np.round(np.mean(yy[mask])))
                    mk, = self._image_canvas._ax.plot(
                        cx, cy, '+', color=color, ms=9, mew=1.4, zorder=12, picker=5)
                    self._image_canvas._extraction_marks.append(mk)
                    new_markers = [mk]
                self._register_extraction(label, rb_spec, new_markers, color=color,
                                          marker_info={'type': 'region', 'region': reg})
                n_ok += 1
            except Exception as e:
                import traceback
                errors.append(f"Region {i+1}: {e}")
                traceback.print_exc()

        self._image_canvas.draw_idle()
        msg = f"Imported {n_ok} region(s)."
        if errors:
            msg += f"  {len(errors)} error(s) — see console."
            for e in errors:
                print(f"[region import] {e}")
        self.statusBar().showMessage(msg, 5000)

    def _on_save_apertures_as_region(self):
        """Save current apertures back to a ds9 .reg file."""
        from rbcodes.GUIs.ifuviewer.processing.spatial_mask import regions_to_ds9_text

        if not self._extraction_marker_groups:
            QMessageBox.information(self, "Nothing to save",
                                    "No apertures have been defined.")
            return

        cube  = self._active_cube
        wcs   = cube.wcs if cube is not None else None
        spectra = self._spectrum_canvas.locked_spectra

        ap_list = []
        for i, minfo in enumerate(self._extraction_marker_info):
            label = spectra[i][0] if i < len(spectra) else ''
            mtype = minfo.get('type', 'single')
            if mtype == 'spaxel':
                ap_list.append({'type': 'single', 'label': label,
                                'cx': minfo['x'], 'cy': minfo['y']})
            elif mtype == 'circle':
                ap_list.append({'type': 'circle', 'label': label,
                                'cx': minfo['cx'], 'cy': minfo['cy'],
                                'radius': minfo['radius'],
                                'bg_inner': minfo.get('bg_inner'),
                                'bg_outer': minfo.get('bg_outer')})
            elif mtype == 'annulus':
                ap_list.append({'type': 'annulus', 'label': label,
                                'cx': minfo['cx'], 'cy': minfo['cy'],
                                'bg_inner': minfo.get('bg_inner', minfo.get('radius', 3)),
                                'bg_outer': minfo.get('bg_outer', minfo.get('radius', 3) + 3)})
            elif mtype == 'rect':
                ap_list.append({'type': 'rect', 'label': label,
                                'x1': minfo['xi0'], 'y1': minfo['yi0'],
                                'x2': minfo['xi1'], 'y2': minfo['yi1']})
            # mtype == 'region': skip — these came from ds9, not user-drawn apertures

        reg_text = regions_to_ds9_text(ap_list, wcs=wcs)

        path, _ = QFileDialog.getSaveFileName(
            self, "Save Region File", "apertures.reg",
            "ds9 region files (*.reg);;All files (*)")
        if not path:
            return
        try:
            with open(path, 'w') as f:
                f.write(reg_text)
            self.statusBar().showMessage(f"Saved {len(ap_list)} region(s) to {path}", 4000)
        except Exception as exc:
            QMessageBox.warning(self, "Save failed", str(exc))

    # ------------------------------------------------------------------
    # Batch extraction
    # ------------------------------------------------------------------

    def _on_batch_extract(self):
        """Open the batch extract dialog."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            QMessageBox.information(self, "No cube",
                                    "Load a cube before batch extracting.")
            return

        has_ds9 = bool(self._ds9_bridge and self._ds9_bridge.available)
        dlg = BatchExtractDialog(has_ds9=has_ds9, parent=self)
        if not dlg.exec_():
            return

        # Get settings
        source, reg_text_or_path, method, weighting, \
            load_strip, save_folder, fname_style, prefix, \
            open_multispec, run_bg = dlg.values()

        # Get region text
        reg_text = ''
        if source == 'file':
            try:
                with open(reg_text_or_path, 'r') as f:
                    reg_text = f.read()
            except Exception as exc:
                QMessageBox.warning(self, "Load failed", str(exc))
                return
        elif source == 'ds9':
            selected_only = (reg_text_or_path == 'selected')
            reg_text = self._ds9_bridge.get_regions(selected_only=selected_only)

        from rbcodes.GUIs.ifuviewer.processing.spatial_mask import (
            parse_ds9_regions, region_to_mask, region_center_sky, iau_name)

        regions = parse_ds9_regions(reg_text)
        if not regions:
            QMessageBox.information(self, "No regions",
                                    "No extractable regions found.")
            return

        if load_strip and len(regions) > 10:
            ret = QMessageBox.question(
                self, "Many regions",
                f"Loading {len(regions)} spectra into the strip may be slow.\nContinue?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if ret != QMessageBox.Yes:
                load_strip = False

        if run_bg:
            # Launch background thread
            self._batch_thread = BatchExtractThread(
                cube, regions, method, weighting,
                load_strip, save_folder, fname_style, prefix, open_multispec,
                wcs=cube.wcs)
            self._batch_thread.progress.connect(self._on_batch_progress)
            self._batch_thread.finished.connect(self._on_batch_finished_dict)
            self._batch_thread.error_log.connect(
                lambda msgs: self.statusBar().showMessage(
                    f"Batch errors: {len(msgs)} — see console.", 5000))
            self._batch_thread.start()
            self.statusBar().showMessage(
                f"Batch extraction started: {len(regions)} regions…", 0)
        else:
            # Run synchronously with a progress dialog
            prog = QProgressDialog(
                "Extracting regions…", "Cancel", 0, len(regions), self)
            prog.setWindowTitle("Batch Extract")
            prog.setMinimumDuration(0)
            results = []
            errors  = []
            for i, reg in enumerate(regions):
                if prog.wasCanceled():
                    break
                prog.setValue(i)
                name = reg.get('name') or f"region_{i+1:03d}"
                prog.setLabelText(f"Extracting {name}  ({i+1}/{len(regions)})")
                try:
                    mask = region_to_mask(reg, cube.wcs, cube.ny, cube.nx)
                    if not mask.any():
                        continue
                    flux_arr = extract_with_method(cube.flux, mask, method)
                    _, err_arr = extract_aperture(cube.flux, cube.var, mask)
                    ra, dec = region_center_sky(reg, cube.wcs)
                    results.append((name, flux_arr, err_arr, ra, dec))
                except Exception as e:
                    errors.append(f"{name}: {e}")
            prog.setValue(len(regions))
            self._apply_batch_results(
                results, cube, method, load_strip, save_folder,
                fname_style, prefix, open_multispec, errors)

    def _on_batch_progress(self, current, total, name):
        self.statusBar().showMessage(
            f"Extracting regions… {current}/{total}  ({name})", 0)

    def _on_batch_finished_dict(self, d):
        self._apply_batch_results(
            d['results'], d['cube'], d['method'],
            d['load_strip'], d['save_folder'],
            d['fname_style'], d['prefix'],
            d['open_multispec'], d['errors'])

    def _apply_batch_results(self, results, cube, method, load_strip,
                             save_folder, fname_style, prefix,
                             open_multispec, errors):
        """Load spectra into strip / save to disk after batch extraction."""
        from rbcodes.GUIs.ifuviewer.processing.spatial_mask import iau_name

        saved = 0
        loaded = 0
        rb_specs = []

        for name, flux_arr, err_arr, ra, dec in results:
            label   = name
            rb_spec = _make_rb_spectrum(cube.wave, flux_arr, err_arr, label)
            if rb_spec is not None:
                rb_spec.meta['object']   = cube.header.get('OBJECT', '') if cube.header else ''
                rb_spec.meta['instrume'] = cube.header.get('INSTRUME', '') if cube.header else ''
                if ra is not None:
                    rb_spec.meta['ra']  = ra
                    rb_spec.meta['dec'] = dec

            if load_strip and rb_spec is not None:
                color = self._spectrum_canvas.lock_spectrum(
                    cube.wave, flux_arr, err_arr, label, rb_spec=rb_spec)
                self._register_extraction(label, rb_spec, [], color=color)
                loaded += 1

            if save_folder and rb_spec is not None:
                if fname_style == 'coords' and ra is not None:
                    fname = iau_name(ra, dec, prefix=prefix)
                elif fname_style == 'name':
                    safe = name.replace(' ', '_').replace('/', '_')
                    fname = f"{prefix}{safe}.fits"
                else:
                    fname = f"{prefix}{name}.fits"
                import os
                fpath = os.path.join(save_folder, fname)
                try:
                    rb_spec.write(fpath)
                    saved += 1
                except Exception as e:
                    errors.append(f"Save {fname}: {e}")

            if rb_spec is not None:
                rb_specs.append(rb_spec)

        self._image_canvas.draw_idle()

        msg = f"Batch complete: {len(results)} extracted"
        if loaded:
            msg += f", {loaded} loaded into strip"
        if saved:
            msg += f", {saved} saved"
        if errors:
            msg += f", {len(errors)} error(s)"
            for e in errors:
                print(f"[batch extract] {e}")
        self.statusBar().showMessage(msg, 6000)

        if open_multispec and rb_specs:
            try:
                from rbcodes.GUIs.multispecviewer.rb_multispec import from_data
                win = from_data(rb_specs, show=True)
                win.setWindowTitle(f"rb_multispec — {len(rb_specs)} batch spectra")
                self._specgui_windows.append(win)
            except Exception as exc:
                print(f"[batch extract] rb_multispec launch failed: {exc}")

    def _on_moment_map(self):
        """Open the Moment Map dialog (Ctrl+M)."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            QMessageBox.information(self, "No cube",
                                    "Load a cube before computing moment maps.")
            return

        sc = self._spectrum_canvas
        onband = sc._onband_ext   # (wmin, wmax) or (nan, nan)
        cont1  = sc._cont1_ext
        cont2  = sc._cont2_ext

        dlg = MomentMapDialog(
            wave          = cube.wave,
            onband        = onband,
            has_variance  = cube.var is not None,
            has_sky       = self._sky_mask is not None,
            has_cont      = np.isfinite(cont1[0]) if cont1 else False,
            parent        = self,
        )
        if not dlg.exec_():
            return

        order, wmin, wmax, lambda_rest, subtract_cont, \
            bcont_min, bcont_max, rcont_min, rcont_max, \
            apply_snr, snr_thresh, snr_method = dlg.values()

        try:
            flux = cube.flux
            if subtract_cont:
                flux = subtract_linear_continuum(
                    flux, cube.wave, bcont_min, bcont_max, rcont_min, rcont_max)
            img = moment_map(flux, cube.wave, wmin, wmax, order, lambda_rest)
        except Exception as exc:
            QMessageBox.warning(self, "Moment map failed", str(exc))
            return

        # Apply SNR mask — always derived from M0 (integrated flux).
        # For M1/M2 compute M0 internally; M1/M2 are only meaningful
        # where the line is detected, which is a flux question, not a
        # velocity question.
        if apply_snr:
            try:
                c1 = tuple(cont1) if np.isfinite(cont1[0]) else None
                c2 = tuple(cont2) if np.isfinite(cont2[0]) else None
                m0_for_snr = (img if order == 0
                              else moment_map(flux, cube.wave, wmin, wmax, 0))
                snr = compute_snr_map(
                    m0_for_snr, flux, cube.wave, wmin, wmax,
                    var      = cube.var       if snr_method == 'Variance cube' else None,
                    sky_mask = self._sky_mask if snr_method == 'Sky region'    else None,
                    cont1    = c1             if snr_method == 'Cont windows'  else None,
                    cont2    = c2             if snr_method == 'Cont windows'  else None,
                )
                if snr is not None:
                    img = np.where(snr >= snr_thresh, img, np.nan)
            except Exception as exc:
                QMessageBox.warning(self, "SNR mask failed", str(exc))

        # Choose colormap per moment order
        cmaps = {0: 'gray', 1: 'RdBu_r', 2: 'plasma'}
        cmap  = cmaps.get(order, 'gray')

        header = cube.spatial_header()
        self._image_canvas.show_image(img, header=header, cmap=cmap)
        self._image_controls.set_cmap(cmap)
        self._image_controls.refresh(img)

        cont_tag = (f"  cont [{bcont_min:.0f}–{bcont_max:.0f}, "
                    f"{rcont_min:.0f}–{rcont_max:.0f}] Å") if subtract_cont else ""
        snr_tag  = f"  SNR≥{snr_thresh} ({snr_method})" if apply_snr else ""
        units = {0: 'flux·Å', 1: 'km/s', 2: 'km/s'}
        self.statusBar().showMessage(
            f"Moment {order} map  [{wmin:.1f}–{wmax:.1f} Å]"
            + (f"  λ_center={lambda_rest:.2f} Å" if lambda_rest else "")
            + cont_tag + snr_tag
            + f"  ({units[order]})  — click any mode button to return",
            0)

    def _on_show_help(self):
        """Open the IFU Viewer help dialog (F1)."""
        from rbcodes.GUIs.ifuviewer.help import show_help_dialog
        show_help_dialog(self)

    def _on_toggle_coord_fmt(self, checked):
        """Toggle cursor RA/Dec display between decimal degrees and sexagesimal."""
        self._image_canvas.coord_fmt = 'sex' if checked else 'deg'

    def _on_show_header(self):
        """Open the FITS header viewer for the active dataset."""
        cube = self._active_cube
        if cube is None:
            QMessageBox.information(self, "No dataset", "Load a dataset first.")
            return
        dlg = HeaderDialog(cube, parent=self)
        dlg.exec_()

    def _on_set_spec_range(self):
        """Open the spectrum X/Y range + scale dialog (Ctrl+R / Ctrl+W)."""
        wave = (self._active_cube.wave
                if self._active_cube and not isinstance(self._active_cube, FITSImage)
                else None)
        dlg = SpectrumRangeDialog(self._spectrum_canvas, wave, parent=self)
        dlg.exec_()

    def _on_set_band_ranges(self):
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            QMessageBox.information(self, "No cube",
                                    "Load a cube before setting band ranges.")
            return

        mode = self._channel_slider.current_mode
        sc   = self._spectrum_canvas
        dlg  = BandRangeDialog(
            wave     = cube.wave,
            mode     = mode,
            onband   = sc._onband_ext,
            cont1    = sc._cont1_ext,
            cont2    = sc._cont2_ext,
            parent   = self,
        )
        if not dlg.exec_():
            return

        onband, cont1, cont2 = dlg.values()

        # Push values into spectrum canvas spans and internal state
        sc.apply_band_ranges(onband, cont1, cont2)

        # Switch to the right mode automatically
        if not np.isnan(onband[0]):
            target_mode = 'Cont-sub' if not np.isnan(cont1[0]) else 'Narrowband'
            self._channel_slider._mode_btns[target_mode].setChecked(True)
            self._channel_slider._on_mode_clicked(
                self._channel_slider._mode_btns[target_mode])

        # Recompute image — pick the right path once
        if not np.isnan(onband[0]):
            if not np.isnan(cont1[0]):
                self._recompute_contsub()
            else:
                self._on_onband_changed(*onband)

    def _on_spaxel_locked(self, x, y, flux_single):
        """Left-click on image — extract using current Mode setting."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return

        ny, nx   = cube.ny, cube.nx
        ac       = self._aperture_controls
        mode     = ac.mode

        if mode == 'Rectangle':
            self.statusBar().showMessage("Rectangle mode: drag to define extraction box.", 3000)
            return

        # ---- Single pixel ------------------------------------------------
        if mode == 'Single pixel':
            flux_arr = flux_single
            err_arr  = None
            if cube.var is not None:
                try:
                    err_arr = np.sqrt(cube.var[:, y, x])
                except Exception:
                    pass
            label = f"Spaxel ({x},{y})"
            n_markers_before = len(self._image_canvas._extraction_marks)
            mk, = self._image_canvas._ax.plot(
                x, y, '+', color='#f9e2af', ms=9, mew=1.4, zorder=12)
            self._image_canvas._extraction_marks.append(mk)
            self._image_canvas.draw_idle()
            new_markers = self._image_canvas._extraction_marks[n_markers_before:]
            marker_info = {'type': 'spaxel', 'x': x, 'y': y}

        # ---- Circular / Circular-Annular ---------------------------------
        else:
            radius    = ac.radius
            method    = ac.method
            weighting = ac.extraction_weighting
            src_mask  = make_circular_mask(ny, nx, x, y, radius)
            if src_mask.sum() == 0:
                return

            if weighting == 'Var-weighted' and cube.var is None:
                self.statusBar().showMessage(
                    "Var-weighted: no variance cube loaded — using unweighted extraction.", 4000)
                weighting = 'None'

            try:
                if weighting == 'Var-weighted':
                    flux_arr, err_arr = extract_variance_weighted(
                        cube.flux, cube.var, src_mask)
                elif weighting == 'Optimal (Data)':
                    flux_arr, err_arr = extract_optimal_weighted(
                        cube.flux, cube.var, src_mask, method='data')
                elif weighting == 'Optimal (Gaussian)':
                    flux_arr, err_arr = extract_optimal_weighted(
                        cube.flux, cube.var, src_mask, method='gaussian')
                else:
                    flux_arr = extract_with_method(cube.flux, src_mask, method)
                    _, err_arr = extract_aperture(cube.flux, cube.var, src_mask)
            except Exception as exc:
                QMessageBox.warning(self, "Extraction error", str(exc))
                return

            bg_inner = bg_outer = None
            bg_tag = ""
            if mode == 'Circular-Annular':
                bg_inner = ac.bg_inner
                bg_outer = ac.bg_outer
                bg_mask  = make_annulus_mask(ny, nx, x, y, bg_inner, bg_outer)
                if bg_mask.sum() > 0:
                    bg_method = 'median' if method == 'median' else 'mean'
                    flux_arr  = subtract_background(
                        flux_arr, cube.flux, bg_mask, bg_method)
                    bg_tag = " [bg]"

            wtag  = {"Var-weighted": " [var]",
                     "Optimal (Data)": " [opt]",
                     "Optimal (Gaussian)": " [opt-g]"}.get(weighting, "")
            label = f"r{radius}px ({x},{y}){wtag}{bg_tag}"

            n_markers_before = len(self._image_canvas._extraction_marks)
            self._image_canvas.draw_aperture_marker(x, y, 'Circle', radius,
                                                    bg_inner, bg_outer)
            new_markers = self._image_canvas._extraction_marks[n_markers_before:]
            marker_info = {'type': 'circle', 'cx': x, 'cy': y, 'radius': radius,
                           'bg_inner': bg_inner, 'bg_outer': bg_outer}

        rb_spec = _make_rb_spectrum(cube.wave, flux_arr, err_arr, label)
        color   = self._spectrum_canvas.lock_spectrum(cube.wave, flux_arr, err_arr,
                                                      label, rb_spec=rb_spec)
        _inject_provenance(rb_spec, cube, x, y, mode, ac)
        self._register_extraction(label, rb_spec, new_markers, color=color,
                                  marker_info=marker_info)

    def _on_aperture_settings_changed(self):
        """Radius/bg spinbox changed — update live preview radius on canvas."""
        ac = self._aperture_controls
        ic = self._image_canvas
        ic.preview_radius   = ac.radius
        ic.preview_bg_inner = ac.bg_inner if ac.mode == 'Circular-Annular' else None
        ic.preview_bg_outer = ac.bg_outer if ac.mode == 'Circular-Annular' else None

    def _on_aperture_mode_changed(self, mode):
        """Mode switched — enable/disable cursor-following preview circle."""
        ic = self._image_canvas
        if mode in ('Circular', 'Circular-Annular'):
            ic.preview_mode = 'circle'
            self._on_aperture_settings_changed()
        else:
            ic.preview_mode = None
            ic.clear_preview_circle()

    def _on_circular_extract(self):
        """[Extract] button — no-op (kept for signal connection; click extracts directly)."""
        pass

    def _on_aperture_drawn(self, mask):
        """Left-drag on image — extract over the drawn rectangle."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return
        n = int(mask.sum())
        if n == 0:
            return

        ac       = self._aperture_controls
        # Drag always uses Circular method settings (mode is secondary for drag)
        method    = ac.method if ac.mode != 'Single pixel' else 'sum'
        weighting = ac.extraction_weighting

        if weighting == 'Var-weighted' and cube.var is None:
            self.statusBar().showMessage(
                "Var-weighted: no variance cube loaded — using unweighted extraction.", 4000)
            weighting = 'None'

        try:
            if weighting == 'Var-weighted':
                flux_arr, err_arr = extract_variance_weighted(
                    cube.flux, cube.var, mask)
            elif weighting == 'Optimal (Data)':
                flux_arr, err_arr = extract_optimal_weighted(
                    cube.flux, cube.var, mask, method='data')
            elif weighting == 'Optimal (Gaussian)':
                flux_arr, err_arr = extract_optimal_weighted(
                    cube.flux, cube.var, mask, method='gaussian')
            else:
                flux_arr = extract_with_method(cube.flux, mask, method)
                _, err_arr = extract_aperture(cube.flux, cube.var, mask)
        except Exception as exc:
            QMessageBox.warning(self, "Extraction error", str(exc))
            return

        wtag  = {"Var-weighted": " [var]",
                 "Optimal (Data)": " [opt]",
                 "Optimal (Gaussian)": " [opt-g]"}.get(weighting, "")
        label = f"Box {n}px/{method}{wtag}"

        # The rectangle was added to _extraction_marks by _handle_aperture_drag
        # just before this signal fired — capture whatever was added since
        # the last registered extraction.
        prev_total = sum(len(g) for g in self._extraction_marker_groups)
        new_markers = list(
            self._image_canvas._extraction_marks[prev_total:])

        rb_spec = _make_rb_spectrum(cube.wave, flux_arr, err_arr, label)
        color   = self._spectrum_canvas.lock_spectrum(cube.wave, flux_arr, err_arr,
                                                      label, rb_spec=rb_spec)
        # For rectangle drag, inject provenance using mask centroid
        if cube.wcs is not None:
            try:
                yy, xx = np.mgrid[0:cube.ny, 0:cube.nx]
                cx = float(np.mean(xx[mask])); cy = float(np.mean(yy[mask]))
                _inject_provenance(rb_spec, cube, cx, cy, 'rect', None)
            except Exception:
                pass
        lr = getattr(self._image_canvas, '_last_rect', None)
        marker_info = ({'type': 'rect',
                        'xi0': lr[0], 'yi0': lr[1], 'xi1': lr[2], 'yi1': lr[3]}
                       if lr else {})
        self._register_extraction(label, rb_spec, new_markers, color=color,
                                  marker_info=marker_info)

    def _register_extraction(self, label, rb_spec, markers=None, color=None,
                             marker_info=None):
        """Add entry to extraction combo and enable strip controls."""
        self._extract_combo.addItem(label)
        self._extract_combo.setCurrentIndex(self._extract_combo.count() - 1)
        self._extraction_marker_groups.append(list(markers or []))
        self._extraction_colors.append(color or '#f9e2af')
        self._extraction_marker_info.append(marker_info or {})
        # Recolor markers to match spectrum line color
        if color and markers:
            for artist in markers:
                try:
                    artist.set_edgecolor(color)
                except AttributeError:
                    try:
                        artist.set_color(color)
                    except Exception:
                        pass
            self._image_canvas.draw_idle()
        self._set_extract_strip_enabled(True)

    def _on_save_extraction(self):
        idx = self._extract_combo.currentIndex()
        if idx < 0:
            return
        spectra = self._spectrum_canvas.locked_spectra
        if idx >= len(spectra):
            return
        label, rb_spec = spectra[idx]
        if rb_spec is None:
            QMessageBox.information(self, "Save", "No spectrum object available.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save spectrum", f"{label}.fits",
            "FITS (*.fits);;ASCII (*.txt *.dat);;All files (*)",
        )
        if not path:
            return
        try:
            rb_spec.write(path)
        except Exception as exc:
            QMessageBox.warning(self, "Save failed", str(exc))

    def _on_send_to_multispec(self):
        """Let the user pick which extracted spectra to send to rb_multispec."""
        all_spectra = [(label, rb_spec)
                       for label, rb_spec in self._spectrum_canvas.locked_spectra
                       if rb_spec is not None]
        if not all_spectra:
            QMessageBox.information(self, "rb_multispec",
                                    "No extracted spectra available.")
            return

        dlg = MultispecSendDialog(all_spectra, parent=self)
        if not dlg.exec_():
            return

        chosen = dlg.selected_spectra()
        if not chosen:
            return

        spec_list = [rb_spec for _, rb_spec in chosen]
        no_err = [label for label, rb_spec in chosen
                  if getattr(rb_spec, 'sig', None) is None]
        if no_err:
            self.statusBar().showMessage(
                f"Note: {len(no_err)} spectrum/spectra sent without error array "
                f"({', '.join(no_err[:3])}{'…' if len(no_err) > 3 else ''}) "
                f"— error panel will be blank in rb_multispec", 8000)
        try:
            from rbcodes.GUIs.multispecviewer.rb_multispec import from_data
            win = from_data(spec_list, show=True)
            win.setWindowTitle(f"rb_multispec — {len(spec_list)} spectra")
            self._specgui_windows.append(win)
        except Exception as exc:
            QMessageBox.warning(self, "rb_multispec launch failed", str(exc))

    def _on_delete_extraction(self):
        """Delete the currently selected extraction (spectrum line + image markers)."""
        idx = self._extract_combo.currentIndex()
        if idx < 0:
            return

        # Remove image markers
        if idx < len(self._extraction_marker_groups):
            for artist in self._extraction_marker_groups[idx]:
                try:
                    artist.remove()
                except Exception:
                    pass
            del self._extraction_marker_groups[idx]
            if idx < len(self._extraction_colors):
                del self._extraction_colors[idx]
            if idx < len(self._extraction_marker_info):
                del self._extraction_marker_info[idx]
            self._image_canvas.draw_idle()

        # Remove spectrum line
        self._spectrum_canvas.remove_locked_at(idx)

        # Remove combo entry
        self._extract_combo.removeItem(idx)

        if self._extract_combo.count() == 0:
            self._set_extract_strip_enabled(False)

    def _get_selected_spectrum(self):
        """Return (rb_spectrum, label) for the currently selected extraction."""
        idx = self._extract_combo.currentIndex()
        if idx < 0:
            return None, None
        spectra = self._spectrum_canvas.locked_spectra
        if idx >= len(spectra):
            return None, None
        label, rb_spec = spectra[idx]
        if rb_spec is None:
            QMessageBox.information(self, "No spectrum",
                                    "No spectrum object available for this entry.")
            return None, None
        return rb_spec, label

    def _on_clear_extractions(self, _from_switch=False):
        self._spectrum_canvas.clear_locked()
        self._image_canvas.clear_extraction_marks()
        self._image_canvas.clear_preview_circle()
        self._extract_combo.clear()
        self._extraction_marker_groups.clear()
        self._extraction_colors.clear()
        self._extraction_marker_info.clear()
        self._image_overlays.clear()
        self._set_extract_strip_enabled(False)
        # Persist the cleared state so overlays don't reappear on next switch
        if not _from_switch and self._active_cube is not None:
            self._save_dataset_state()

    def _on_spaxel_hovered(self, x, y, spectrum):
        cube = self._active_cube
        ac   = self._aperture_controls
        if (cube is not None and not isinstance(cube, FITSImage)
                and ac.mode in ('Circular', 'Circular-Annular')):
            # Skip re-extraction if same pixel and same radius as last call
            cache = getattr(self, '_hover_cache', None)
            key   = (x, y, ac.radius, ac.mode)
            if cache is not None and cache[0] == key:
                spectrum = cache[1]
            else:
                try:
                    mask = make_circular_mask(cube.ny, cube.nx, x, y, ac.radius)
                    if mask.sum() > 0:
                        spectrum = extract_with_method(cube.flux, mask, ac.method)
                except Exception:
                    pass
                self._hover_cache = (key, spectrum)
        else:
            self._hover_cache = None
        self._spectrum_canvas.on_spaxel_hovered(x, y, spectrum)

    def _on_cursor_info(self, info):
        if info:
            self.statusBar().showMessage(f"{self._status_base}    |    {info}")
        else:
            self.statusBar().showMessage(self._status_base)

    def _on_channel_changed(self, idx, wavelength):
        """Slider moved — show slice (Channel) or shift band windows (Narrowband/Cont-sub)."""
        if self._active_cube is None or isinstance(self._active_cube, FITSImage):
            return
        mode = self._channel_slider.current_mode
        if mode == 'Channel':
            self._image_canvas.show_slice(idx)
            self._image_controls.refresh(self._image_canvas._data_raw)
        elif mode in ('Narrowband', 'Cont-sub'):
            self._shift_bands_to_channel(idx)

    def _on_mode_changed(self, mode):
        """Image mode button clicked — recompute and display the right image."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return

        # Tell spectrum panel which selectors to activate
        self._spectrum_canvas.set_mode(mode)

        if mode == 'Channel':
            idx = self._channel_slider.current_index
            self._image_canvas.show_slice(idx)
            data2d = self._image_canvas._data_raw
        elif mode == 'Whitelight':
            method = self._channel_slider.current_method
            data2d = build_whitelight(cube.flux, cube.wave, method=method)
            self._image_canvas.update_data(data2d)
        elif mode in ('Narrowband', 'Cont-sub'):
            sc = self._spectrum_canvas
            wmin, wmax = sc._onband_ext
            if not np.isnan(wmin):
                # Span already set — recompute and sync slider to band center
                center = (wmin + wmax) / 2.0
                init_idx = int(np.argmin(np.abs(cube.wave - center)))
                self._channel_slider._slider.blockSignals(True)
                self._channel_slider._slider.setValue(init_idx)
                self._channel_slider._slider.blockSignals(False)
                if mode == 'Narrowband':
                    method = self._channel_slider.current_method
                    data2d = build_narrowband(cube.flux, cube.wave, wmin, wmax, method=method)
                    self._image_canvas.update_data(data2d)
                else:
                    self._recompute_contsub()
                    data2d = self._image_canvas._data_raw
            else:
                # No span yet — show whitelight as placeholder
                method = self._channel_slider.current_method
                data2d = build_whitelight(cube.flux, cube.wave, method=method)
                self._image_canvas.update_data(data2d)
        else:
            return

        self._image_controls.refresh(data2d)

    def _on_onband_changed(self, wmin, wmax):
        """On-band span dragged — recompute narrowband or cont-sub image."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return
        mode   = self._channel_slider.current_mode
        method = self._channel_slider.current_method
        if mode == 'Narrowband':
            data2d = build_narrowband(cube.flux, cube.wave, wmin, wmax,
                                      method=method)
            self._image_canvas.update_data(data2d)
            self._image_controls.refresh(data2d)
        elif mode == 'Cont-sub':
            # Recompute if continuum is already defined
            c1min, c1max = self._spectrum_canvas._cont1_ext
            if not np.isnan(c1min):
                self._recompute_contsub()

        # Sync slider to on-band center so step buttons stay in sync
        if mode in ('Narrowband', 'Cont-sub'):
            center = (wmin + wmax) / 2.0
            idx = int(np.argmin(np.abs(cube.wave - center)))
            self._channel_slider._slider.blockSignals(True)
            self._channel_slider._slider.setValue(idx)
            self._channel_slider._slider.blockSignals(False)

    def _on_contband_changed(self, c1min, c1max, c2min, c2max):
        """Continuum span dragged — recompute cont-sub image."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return
        if self._channel_slider.current_mode == 'Cont-sub':
            self._recompute_contsub()

    def _on_window_cleared(self, name):
        """A band window was cleared — show appropriate fallback image."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return
        sc = self._spectrum_canvas
        wmin, _ = sc._onband_ext

        if name == 'onband' or np.isnan(wmin):
            # No on-band window → whitelight
            method = self._channel_slider.current_method
            data2d = build_whitelight(cube.flux, cube.wave, method=method)
            self._image_canvas.update_data(data2d)
            self._image_controls.refresh(data2d)
        else:
            # A cont window was cleared — recompute with whatever remains
            # _recompute_contsub handles: cont2-only → cont-sub, no cont → narrowband
            self._recompute_contsub()

    def _recompute_contsub(self):
        """Build and display the continuum-subtracted image from current spans."""
        cube   = self._active_cube
        method = self._channel_slider.current_method
        wmin,  wmax  = self._spectrum_canvas._onband_ext
        c1min, c1max = self._spectrum_canvas._cont1_ext
        c2min, c2max = self._spectrum_canvas._cont2_ext

        if np.isnan(wmin):
            return   # no on-band — nothing to compute

        has_c1 = not np.isnan(c1min)
        has_c2 = not np.isnan(c2min)

        if not has_c1 and not has_c2:
            # No continuum windows — fall back to plain narrowband
            data2d = build_narrowband(cube.flux, cube.wave, wmin, wmax, method=method)
            self._image_canvas.update_data(data2d)
            self._image_controls.refresh(data2d)
            return

        # build_continuum_sub requires c1; if only c2 is set, promote it to c1
        if has_c1:
            use_c1min, use_c1max = c1min, c1max
            c2min_arg = c2min if has_c2 else None
            c2max_arg = c2max if has_c2 else None
        else:
            use_c1min, use_c1max = c2min, c2max
            c2min_arg = None
            c2max_arg = None

        try:
            data2d = build_continuum_sub(
                cube.flux, cube.wave,
                wmin, wmax, use_c1min, use_c1max,
                c2min=c2min_arg, c2max=c2max_arg,
                method=method,
            )
        except Exception:
            return

        self._image_canvas.update_data(data2d)
        self._image_controls.refresh(data2d)

    # ------------------------------------------------------------------
    # Band-mode slider: shift on-band (and cont) windows to a new channel
    # ------------------------------------------------------------------

    def _shift_bands_to_channel(self, idx):
        """Slide the on-band window center to wave[idx]; cont windows follow."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return
        sc   = self._spectrum_canvas
        wave = cube.wave

        wmin, wmax = sc._onband_ext
        if np.isnan(wmin):
            return   # no span drawn yet — do nothing

        onband_center = (wmin + wmax) / 2.0
        dlam          = wave[idx] - onband_center
        new_wmin      = wmin + dlam
        new_wmax      = wmax + dlam

        # On-band edge check
        if new_wmin < wave[0] or new_wmax > wave[-1]:
            # Snap slider back to the current on-band center
            prev_idx = int(np.argmin(np.abs(wave - onband_center)))
            self._channel_slider._slider.blockSignals(True)
            self._channel_slider._slider.setValue(prev_idx)
            self._channel_slider._slider.blockSignals(False)
            return

        mode = self._channel_slider.current_mode

        if mode == 'Narrowband':
            sc._onband_ext = (new_wmin, new_wmax)
            if sc._onband_sel is not None:
                try:
                    sc._onband_sel.extents = (float(new_wmin), float(new_wmax))
                except Exception:
                    pass
            sc.draw_idle()
            method = self._channel_slider.current_method
            data2d = build_narrowband(cube.flux, wave, new_wmin, new_wmax, method=method)
            self._image_canvas.update_data(data2d)
            self._image_controls.refresh(data2d)

        elif mode == 'Cont-sub':
            c1min, c1max = sc._cont1_ext
            c2min, c2max = sc._cont2_ext

            new_c1 = (np.nan, np.nan)
            if not np.isnan(c1min):
                nc1min, nc1max = c1min + dlam, c1max + dlam
                if nc1min < wave[0] or nc1max > wave[-1]:
                    self.statusBar().showMessage(
                        "Cont 1 window hit spectral edge — re-draw with mouse to reset.", 4000)
                    prev_idx = int(np.argmin(np.abs(wave - onband_center)))
                    self._channel_slider._slider.blockSignals(True)
                    self._channel_slider._slider.setValue(prev_idx)
                    self._channel_slider._slider.blockSignals(False)
                    return
                new_c1 = (nc1min, nc1max)

            new_c2 = (np.nan, np.nan)
            if not np.isnan(c2min):
                nc2min, nc2max = c2min + dlam, c2max + dlam
                if nc2min < wave[0] or nc2max > wave[-1]:
                    self.statusBar().showMessage(
                        "Cont 2 window hit spectral edge — re-draw with mouse to reset.", 4000)
                    prev_idx = int(np.argmin(np.abs(wave - onband_center)))
                    self._channel_slider._slider.blockSignals(True)
                    self._channel_slider._slider.setValue(prev_idx)
                    self._channel_slider._slider.blockSignals(False)
                    return
                new_c2 = (nc2min, nc2max)

            sc.apply_band_ranges((new_wmin, new_wmax), new_c1, new_c2)
            self._recompute_contsub()

    # ------------------------------------------------------------------
    # Spectrum range strip
    # ------------------------------------------------------------------

    def _build_spec_range_strip(self):
        """Thin strip below the spectrum: Auto X/Y buttons + X/Y range text boxes."""
        strip = QWidget()
        strip.setMaximumHeight(30)
        layout = QHBoxLayout(strip)
        layout.setContentsMargins(4, 1, 4, 1)
        layout.setSpacing(4)

        auto_x = QPushButton("Auto X")
        auto_x.setFixedWidth(56)
        auto_x.setToolTip("Reset wavelength axis to full range")
        auto_x.clicked.connect(self._on_spec_auto_x)
        layout.addWidget(auto_x)

        layout.addWidget(QLabel("λ:"))
        self._spec_xmin = QLineEdit()
        self._spec_xmin.setFixedWidth(72)
        self._spec_xmin.setToolTip("Wavelength min (press Enter to apply)")
        self._spec_xmin.returnPressed.connect(self._on_spec_x_entered)
        layout.addWidget(self._spec_xmin)
        layout.addWidget(QLabel("–"))
        self._spec_xmax = QLineEdit()
        self._spec_xmax.setFixedWidth(72)
        self._spec_xmax.setToolTip("Wavelength max (press Enter to apply)")
        self._spec_xmax.returnPressed.connect(self._on_spec_x_entered)
        layout.addWidget(self._spec_xmax)

        layout.addSpacing(10)

        auto_y = QPushButton("Auto Y")
        auto_y.setFixedWidth(56)
        auto_y.setToolTip("Auto-scale Y to visible data (1st–99th percentile)")
        auto_y.clicked.connect(self._on_spec_auto_y)
        layout.addWidget(auto_y)

        layout.addWidget(QLabel("Y:"))
        self._spec_ymin = QLineEdit()
        self._spec_ymin.setFixedWidth(72)
        self._spec_ymin.setToolTip("Y min (press Enter to apply)")
        self._spec_ymin.returnPressed.connect(self._on_spec_y_entered)
        layout.addWidget(self._spec_ymin)
        layout.addWidget(QLabel("–"))
        self._spec_ymax = QLineEdit()
        self._spec_ymax.setFixedWidth(72)
        self._spec_ymax.setToolTip("Y max (press Enter to apply)")
        self._spec_ymax.returnPressed.connect(self._on_spec_y_entered)
        layout.addWidget(self._spec_ymax)

        layout.addStretch()
        return strip

    def _sync_spec_range_strip(self, *_):
        """Populate X/Y text boxes from current spectrum axes limits."""
        try:
            ax = self._spectrum_canvas._ax
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            self._spec_xmin.setText(f'{xmin:.2f}')
            self._spec_xmax.setText(f'{xmax:.2f}')
            self._spec_ymin.setText(f'{ymin:.4g}')
            self._spec_ymax.setText(f'{ymax:.4g}')
        except Exception:
            pass

    def _on_spec_auto_x(self):
        self._spectrum_canvas.reset_wave_limits()
        self._sync_spec_range_strip()

    def _on_spec_auto_y(self):
        self._spectrum_canvas.autoscale_y()
        self._sync_spec_range_strip()

    def _on_spec_x_entered(self):
        try:
            xmin = float(self._spec_xmin.text())
            xmax = float(self._spec_xmax.text())
        except ValueError:
            return
        if xmin >= xmax:
            return
        self._spectrum_canvas.set_wave_limits(xmin, xmax)

    def _on_spec_y_entered(self):
        try:
            ymin = float(self._spec_ymin.text())
            ymax = float(self._spec_ymax.text())
        except ValueError:
            return
        if ymin >= ymax:
            return
        sc = self._spectrum_canvas
        sc._ax.set_ylim(ymin, ymax)
        sc.draw_idle()

    def _on_method_changed(self, method):
        """Collapse method dropdown changed — recompute current image mode."""
        cube = self._active_cube
        if cube is None or isinstance(cube, FITSImage):
            return
        mode = self._channel_slider.current_mode
        if mode == 'Whitelight':
            data2d = build_whitelight(cube.flux, cube.wave, method=method)
            self._image_canvas.update_data(data2d)
            self._image_controls.refresh(data2d)
        elif mode == 'Narrowband':
            wmin, wmax = self._spectrum_canvas._onband_ext
            if not np.isnan(wmin):
                data2d = build_narrowband(cube.flux, cube.wave, wmin, wmax, method=method)
                self._image_canvas.update_data(data2d)
                self._image_controls.refresh(data2d)
        elif mode == 'Cont-sub':
            self._recompute_contsub()

    def _show_default_image(self, cube):
        """Compute the default display image and hand it to ImageCanvas."""
        is_image = isinstance(cube, FITSImage)

        # Store cube reference in canvas for show_slice
        self._image_canvas._cube = cube

        if is_image:
            data2d = cube.flux
            header = cube.spatial_header()
        else:
            data2d = build_whitelight(cube.flux, cube.wave)
            header = cube.spatial_header()

        self._placeholder.hide()
        self._splitter.show()
        self._mpl_toolbar.show()
        self._image_controls.show()
        self._image_canvas.show_image(data2d, header=header)
        self._image_controls.reset(data2d)

        # Collapsed spectrum
        if not is_image:
            collapsed = np.nanmedian(cube.flux, axis=(1, 2))
            self._spectrum_canvas.show_spectrum(cube.wave, collapsed,
                                                label='Collapsed (median)')
            self._sync_spec_range_strip()
            self._splitter.setSizes([600, 200])
            self._spectrum_canvas.show()
            self._aperture_controls.show()
            self._aperture_controls.set_has_variance(cube.var is not None)
            for w in (self._extract_combo, self._save_btn,
                      self._multispec_btn, self._del_btn, self._clear_btn):
                w.setVisible(True)
            self._extract_strip.show()
            self._set_extract_strip_enabled(False)   # empty state
        else:
            self._spectrum_canvas.hide()
            self._splitter.setSizes([800, 0])
            self._aperture_controls.hide()
            # Show strip with only Clear all enabled (no spectra, but overlays can be cleared)
            for w in (self._extract_combo, self._save_btn,
                      self._multispec_btn, self._del_btn):
                w.setVisible(False)
            self._clear_btn.setVisible(True)
            self._clear_btn.setEnabled(False)
            self._extract_strip.show()

        # Channel slider: only for cubes
        if is_image:
            self._channel_slider.hide()
            self._channel_slider.set_cube(None)
        else:
            self._channel_slider.show()
            self._channel_slider.set_cube(cube)
            self._channel_slider._mode_btns['Whitelight'].setChecked(True)
            self._channel_slider._set_slider_enabled(False)

    # ------------------------------------------------------------------
    # Keyboard shortcuts
    # ------------------------------------------------------------------

    def keyPressEvent(self, event):
        """Delete / Backspace → delete the selected extraction."""
        if event.key() in (Qt.Key_Delete, Qt.Key_Backspace):
            if self._extract_combo.count() > 0:
                self._on_delete_extraction()
                return
        super().keyPressEvent(event)

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    @property
    def active_cube(self):
        return self._active_cube


# ---------------------------------------------------------------------------
# rb_spectrum factory
# ---------------------------------------------------------------------------

def _write_temp_fits(rb_spec):
    """Write rb_spectrum to a temp FITS file; return path or None on failure."""
    import tempfile
    try:
        path = tempfile.mktemp(suffix='.fits')
        rb_spec.write(path)
        return path
    except Exception:
        return None


def _make_rb_spectrum(wave, flux, err, label):
    """Wrap extracted arrays into an rb_spectrum object. Returns None on failure."""
    try:
        from rbcodes.utils.rb_spectrum import rb_spectrum
        return rb_spectrum(
            wavelength=wave,
            flux=flux,
            error=err,
            filename=label,
        )
    except Exception:
        return None


def _inject_provenance(rb_spec, cube, cx, cy, mode, ac):
    """
    Inject RA/Dec and extraction metadata into rb_spec.meta.

    Parameters
    ----------
    rb_spec : rb_spectrum or None
    cube    : IFUCube
    cx, cy  : float  — pixel centre (0-indexed)
    mode    : str    — 'Single pixel', 'Circular', 'Circular-Annular', 'rect'
    ac      : ApertureControls or None
    """
    if rb_spec is None:
        return
    try:
        if cube.wcs is not None:
            xy  = cube.wcs.all_pix2world([[cx, cy]], 0)
            rb_spec.meta['ra']  = float(xy[0][0])
            rb_spec.meta['dec'] = float(xy[0][1])
        rb_spec.meta['extr_x'] = int(round(cx))
        rb_spec.meta['extr_y'] = int(round(cy))
        if ac is not None and mode not in ('Single pixel', 'rect'):
            rb_spec.meta['extr_rad'] = int(ac.radius)
        else:
            rb_spec.meta['extr_rad'] = 1
        if cube.header is not None:
            rb_spec.meta['object']   = cube.header.get('OBJECT', '')
            rb_spec.meta['instrume'] = cube.header.get('INSTRUME', '')
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Status bar helper
# ---------------------------------------------------------------------------

def _status_text(cube):
    if isinstance(cube, FITSImage):
        return f"{cube.name}  |  2D image  |  {cube.ny} × {cube.nx} px"
    w0, w1 = cube.wave_range
    wave_str = f"{w0:.1f} – {w1:.1f} Å" if w0 is not None else "unknown"
    return (f"{cube.name}  |  {cube.__class__.__name__}  "
            f"|  shape {cube.flux.shape}  |  λ {wave_str}")



# ---------------------------------------------------------------------------
# Wavelength range dialog
# ---------------------------------------------------------------------------

class WaveRangeDialog(QDialog):
    """
    Small popup for manually setting the spectrum x-axis range.

    Pre-fills with the current axis limits.  Shows the full cube wavelength
    range as a hint when available.
    """

    def __init__(self, current_xlim, wave=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Set Spectrum Wavelength Range")
        self.setMinimumWidth(320)

        layout = QFormLayout(self)
        layout.setContentsMargins(12, 12, 12, 8)
        layout.setSpacing(8)

        if wave is not None:
            hint = QLabel(f"Full range: {wave[0]:.1f} – {wave[-1]:.1f} Å")
            hint.setStyleSheet("color: #a6adc8; font-size: 8pt;")
            layout.addRow(hint)

        self._wmin = QDoubleSpinBox()
        self._wmin.setDecimals(2)
        self._wmin.setRange(0, 1e7)
        self._wmin.setValue(current_xlim[0])
        self._wmin.setSuffix("  Å")
        self._wmin.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
        layout.addRow("λ min:", self._wmin)

        self._wmax = QDoubleSpinBox()
        self._wmax.setDecimals(2)
        self._wmax.setRange(0, 1e7)
        self._wmax.setValue(current_xlim[1])
        self._wmax.setSuffix("  Å")
        self._wmax.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
        layout.addRow("λ max:", self._wmax)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addRow(btns)

    def values(self):
        wmin = self._wmin.value()
        wmax = self._wmax.value()
        if wmin > wmax:
            wmin, wmax = wmax, wmin
        return wmin, wmax


# ---------------------------------------------------------------------------
# Band range dialog
# ---------------------------------------------------------------------------

def _make_spinbox(lo, hi, value):
    """Create a QDoubleSpinBox pre-configured for wavelength entry."""
    sb = QDoubleSpinBox()
    sb.setDecimals(2)
    sb.setRange(max(0, lo * 0.5), hi * 1.5)
    sb.setValue(value)
    sb.setSuffix("  Å")
    sb.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
    return sb


class BandRangeDialog(QDialog):
    """
    Popup for manually entering on-band and continuum wavelength windows.

    Sections
    --------
    On-band (green)         — always active
    Continuum 1 (red)       — enable checkbox
    Continuum 2 (orange)    — enable checkbox (optional second window)
    """

    def __init__(self, wave, mode='Narrowband',
                 onband=(np.nan, np.nan),
                 cont1=(np.nan, np.nan),
                 cont2=(np.nan, np.nan),
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Set Band Ranges")
        self.setMinimumWidth(380)

        wlo = float(wave[0])
        whi = float(wave[-1])

        layout = QFormLayout(self)
        layout.setContentsMargins(12, 12, 12, 8)
        layout.setSpacing(6)

        # Full-range hint
        hint = QLabel(f"Full range: {wlo:.1f} – {whi:.1f} Å")
        hint.setStyleSheet("color: #a6adc8; font-size: 8pt;")
        layout.addRow(hint)

        # --- On-band ---
        ob_label = QLabel("On-band window (green)")
        ob_label.setStyleSheet("color: #a6e3a1; font-weight: bold; margin-top: 6px;")
        layout.addRow(ob_label)

        ob_lo = onband[0] if not np.isnan(onband[0]) else wlo
        ob_hi = onband[1] if not np.isnan(onband[1]) else whi
        self._ob_min = _make_spinbox(wlo, whi, ob_lo)
        self._ob_max = _make_spinbox(wlo, whi, ob_hi)
        layout.addRow("λ min:", self._ob_min)
        layout.addRow("λ max:", self._ob_max)

        # --- Continuum 1 ---
        c1_label = QLabel("Continuum window 1 (red)")
        c1_label.setStyleSheet("color: #f38ba8; font-weight: bold; margin-top: 6px;")
        layout.addRow(c1_label)

        self._c1_cb = QCheckBox("Enable")
        self._c1_cb.setChecked(not np.isnan(cont1[0]))
        layout.addRow(self._c1_cb)

        c1_lo = cont1[0] if not np.isnan(cont1[0]) else wlo
        c1_hi = cont1[1] if not np.isnan(cont1[1]) else whi
        self._c1_min = _make_spinbox(wlo, whi, c1_lo)
        self._c1_max = _make_spinbox(wlo, whi, c1_hi)
        layout.addRow("λ min:", self._c1_min)
        layout.addRow("λ max:", self._c1_max)

        def _toggle_c1(checked):
            self._c1_min.setEnabled(checked)
            self._c1_max.setEnabled(checked)
        self._c1_cb.toggled.connect(_toggle_c1)
        _toggle_c1(self._c1_cb.isChecked())

        # --- Continuum 2 (optional) ---
        c2_label = QLabel("Continuum window 2 (orange, optional)")
        c2_label.setStyleSheet("color: #fab387; font-weight: bold; margin-top: 6px;")
        layout.addRow(c2_label)

        self._c2_cb = QCheckBox("Enable")
        self._c2_cb.setChecked(not np.isnan(cont2[0]))
        layout.addRow(self._c2_cb)

        c2_lo = cont2[0] if not np.isnan(cont2[0]) else wlo
        c2_hi = cont2[1] if not np.isnan(cont2[1]) else whi
        self._c2_min = _make_spinbox(wlo, whi, c2_lo)
        self._c2_max = _make_spinbox(wlo, whi, c2_hi)
        layout.addRow("λ min:", self._c2_min)
        layout.addRow("λ max:", self._c2_max)

        def _toggle_c2(checked):
            self._c2_min.setEnabled(checked)
            self._c2_max.setEnabled(checked)
        self._c2_cb.toggled.connect(_toggle_c2)
        _toggle_c2(self._c2_cb.isChecked())

        # Buttons
        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addRow(btns)

    # ------------------------------------------------------------------

    def values(self):
        """Return (onband, cont1, cont2) tuples.  Disabled windows → (nan, nan)."""
        def _sorted(lo, hi):
            return (lo, hi) if lo <= hi else (hi, lo)

        onband = _sorted(self._ob_min.value(), self._ob_max.value())

        if self._c1_cb.isChecked():
            cont1 = _sorted(self._c1_min.value(), self._c1_max.value())
        else:
            cont1 = (np.nan, np.nan)

        if self._c2_cb.isChecked():
            cont2 = _sorted(self._c2_min.value(), self._c2_max.value())
        else:
            cont2 = (np.nan, np.nan)

        return onband, cont1, cont2


# ---------------------------------------------------------------------------
# rb_multispec selection dialog
# ---------------------------------------------------------------------------

class MultispecSendDialog(QDialog):
    """
    Checkbox list of extracted spectra — user picks which to send to rb_multispec.
    All items checked by default.
    """

    def __init__(self, spectra, parent=None):
        """
        Parameters
        ----------
        spectra : list of (label, rb_spectrum)
        """
        super().__init__(parent)
        self.setWindowTitle("Send to rb_multispec")
        self.setMinimumWidth(320)

        self._spectra = spectra

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 8)
        layout.setSpacing(6)

        layout.addWidget(QLabel("Select spectra to send:"))

        self._list = QListWidget()
        self._list.setAlternatingRowColors(True)
        for label, _ in spectra:
            item = QListWidgetItem(label)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)
            self._list.addItem(item)
        layout.addWidget(self._list)

        # Select-all / none row
        btn_row = QHBoxLayout()
        all_btn  = QPushButton("Select all")
        none_btn = QPushButton("Select none")
        all_btn.setFixedWidth(90)
        none_btn.setFixedWidth(90)
        all_btn.clicked.connect(lambda: self._set_all(Qt.Checked))
        none_btn.clicked.connect(lambda: self._set_all(Qt.Unchecked))
        btn_row.addWidget(all_btn)
        btn_row.addWidget(none_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def _set_all(self, state):
        for i in range(self._list.count()):
            self._list.item(i).setCheckState(state)

    def selected_spectra(self):
        """Return list of (label, rb_spectrum) for checked items."""
        result = []
        for i in range(self._list.count()):
            if self._list.item(i).checkState() == Qt.Checked:
                result.append(self._spectra[i])
        return result


# ---------------------------------------------------------------------------
# Spectrum range / scale dialog
# ---------------------------------------------------------------------------

class SpectrumRangeDialog(QDialog):
    """
    Pop-up to set X range (wavelength), Y range (flux), and Y scale for the
    spectrum panel.  Changes apply live on [Apply]; cancelled on [Cancel].

    Keyboard shortcut: Ctrl+R  (also Ctrl+W alias).
    """

    def __init__(self, spectrum_canvas, wave=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Spectrum Range & Scale")
        self.setMinimumWidth(360)

        self._sc   = spectrum_canvas
        self._wave = wave

        # Snapshot so Cancel can restore
        self._orig_xlim   = spectrum_canvas._ax.get_xlim()
        self._orig_ylim   = spectrum_canvas._ax.get_ylim()
        self._orig_yscale = spectrum_canvas._ax.get_yscale()

        layout = QFormLayout(self)
        layout.setContentsMargins(12, 12, 12, 8)
        layout.setSpacing(6)

        # ---- Wavelength hint ----
        if wave is not None:
            hint = QLabel(f"Full range: {wave[0]:.1f} – {wave[-1]:.1f} Å")
            hint.setStyleSheet("color: #a6adc8; font-size: 8pt;")
            layout.addRow(hint)

        # ---- X range ----
        x_label = QLabel("X range (wavelength)")
        x_label.setStyleSheet("font-weight: bold; margin-top: 4px;")
        layout.addRow(x_label)

        lo, hi = self._orig_xlim
        self._xmin = _dblspin(lo, suffix=" Å")
        self._xmax = _dblspin(hi, suffix=" Å")
        layout.addRow("λ min:", self._xmin)
        layout.addRow("λ max:", self._xmax)

        auto_x = QPushButton("Auto X")
        auto_x.setFixedWidth(72)
        auto_x.clicked.connect(self._auto_x)
        layout.addRow("", auto_x)

        # ---- Y range ----
        y_label = QLabel("Y range (flux)")
        y_label.setStyleSheet("font-weight: bold; margin-top: 4px;")
        layout.addRow(y_label)

        ylo, yhi = self._orig_ylim
        # Allow wide negative range for continuum-subtracted cubes
        self._ymin = _dblspin(ylo, lo=-1e15, hi=1e15, suffix="")
        self._ymax = _dblspin(yhi, lo=-1e15, hi=1e15, suffix="")
        layout.addRow("Y min:", self._ymin)
        layout.addRow("Y max:", self._ymax)

        auto_y = QPushButton("Auto Y")
        auto_y.setFixedWidth(72)
        auto_y.clicked.connect(self._auto_y)
        layout.addRow("", auto_y)

        # ---- Y scale ----
        scale_label = QLabel("Y scale")
        scale_label.setStyleSheet("font-weight: bold; margin-top: 4px;")
        layout.addRow(scale_label)

        self._linear_rb = QRadioButton("Linear")
        self._log_rb    = QRadioButton("Log")
        self._scale_grp = QButtonGroup(self)
        self._scale_grp.addButton(self._linear_rb)
        self._scale_grp.addButton(self._log_rb)

        if self._orig_yscale == 'log':
            self._log_rb.setChecked(True)
        else:
            self._linear_rb.setChecked(True)

        scale_row = QHBoxLayout()
        scale_row.addWidget(self._linear_rb)
        scale_row.addWidget(self._log_rb)
        scale_row.addStretch()
        layout.addRow("", scale_row)

        # ---- Buttons ----
        apply_btn = QPushButton("Apply")
        apply_btn.setToolTip("Apply without closing")
        apply_btn.clicked.connect(self._apply)

        reset_btn = QPushButton("Reset")
        reset_btn.setToolTip("Restore full wavelength range, auto Y, linear scale")
        reset_btn.clicked.connect(self._reset)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.addButton(apply_btn, QDialogButtonBox.ActionRole)
        btns.addButton(reset_btn, QDialogButtonBox.ResetRole)
        btns.accepted.connect(self._on_ok)
        btns.rejected.connect(self._on_cancel)
        layout.addRow(btns)

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    def _apply(self):
        xmin = self._xmin.value(); xmax = self._xmax.value()
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        ymin = self._ymin.value(); ymax = self._ymax.value()
        if ymin > ymax:
            ymin, ymax = ymax, ymin

        scale = 'log' if self._log_rb.isChecked() else 'linear'
        self._sc.set_yscale(scale)
        self._sc.set_wave_limits(xmin, xmax)
        self._sc._ax.set_ylim(ymin, ymax)
        self._sc.draw_idle()

    def _auto_x(self):
        if self._wave is not None:
            self._xmin.setValue(float(self._wave[0]))
            self._xmax.setValue(float(self._wave[-1]))

    def _auto_y(self):
        self._sc.autoscale_y()
        ylo, yhi = self._sc._ax.get_ylim()
        self._ymin.setValue(ylo)
        self._ymax.setValue(yhi)

    def _reset(self):
        self._linear_rb.setChecked(True)
        self._sc.set_yscale('linear')
        self._sc.reset_wave_limits()
        self._sc.autoscale_y()
        if self._wave is not None:
            self._xmin.setValue(float(self._wave[0]))
            self._xmax.setValue(float(self._wave[-1]))
        ylo, yhi = self._sc._ax.get_ylim()
        self._ymin.setValue(ylo)
        self._ymax.setValue(yhi)

    def _on_ok(self):
        self._apply()
        self.accept()

    def _on_cancel(self):
        # Restore original state
        self._sc.set_yscale(self._orig_yscale)
        self._sc._ax.set_xlim(*self._orig_xlim)
        self._sc._ax.set_ylim(*self._orig_ylim)
        self._sc.draw_idle()
        self.reject()


def _dblspin(value, lo=0.0, hi=1e7, suffix=""):
    sb = QDoubleSpinBox()
    sb.setDecimals(2)
    sb.setRange(lo, hi)
    sb.setValue(value)
    if suffix:
        sb.setSuffix(f"  {suffix}")
    sb.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
    return sb


# ---------------------------------------------------------------------------
# Stylesheet
# ---------------------------------------------------------------------------

_QSS = """
QMainWindow, QWidget {
    background-color: #1e1e2e;
    color: #cdd6f4;
}

QLabel {
    font-size: 9pt;
    color: #cdd6f4;
}

QLabel#placeholder {
    font-size: 14pt;
    color: #585b70;
}

QListWidget {
    background-color: #181825;
    color: #cdd6f4;
    border: 1px solid #45475a;
    font-size: 9pt;
    alternate-background-color: #1e1e2e;
}

QListWidget::item:selected {
    background-color: #89b4fa;
    color: #1e1e2e;
}

QListWidget::item:hover {
    background-color: #313244;
}

QPushButton {
    font-size: 9pt;
    padding: 4px 10px;
    border: 1px solid #45475a;
    border-radius: 4px;
    background: #313244;
    color: #cdd6f4;
    min-height: 22px;
}

QPushButton:hover {
    background: #45475a;
    border-color: #89b4fa;
}

QPushButton:pressed {
    background: #89b4fa;
    color: #1e1e2e;
}

QMenuBar {
    background-color: #181825;
    color: #cdd6f4;
}

QMenuBar::item:selected {
    background-color: #313244;
}

QMenu {
    background-color: #181825;
    color: #cdd6f4;
    border: 1px solid #45475a;
}

QMenu::item:selected {
    background-color: #89b4fa;
    color: #1e1e2e;
}

QStatusBar {
    background-color: #181825;
    color: #a6adc8;
    font-size: 8pt;
}

QScrollBar:vertical {
    background: #181825;
    width: 10px;
}

QScrollBar::handle:vertical {
    background: #45475a;
    border-radius: 5px;
    min-height: 20px;
}

QDialog {
    background-color: #1e1e2e;
    color: #cdd6f4;
}

QFormLayout QLabel {
    color: #cdd6f4;
    font-size: 9pt;
}

QDoubleSpinBox {
    background-color: #313244;
    color: #cdd6f4;
    border: 1px solid #45475a;
    border-radius: 3px;
    padding: 3px 6px;
    font-size: 9pt;
}

QDoubleSpinBox:focus { border-color: #89b4fa; }

QDialogButtonBox QPushButton {
    min-width: 70px;
}

QCheckBox {
    color: #cdd6f4;
    font-size: 9pt;
    spacing: 4px;
}

QCheckBox::indicator {
    width: 14px;
    height: 14px;
    border: 1px solid #45475a;
    border-radius: 3px;
    background: #313244;
}

QCheckBox::indicator:checked {
    background: #89b4fa;
    border-color: #89b4fa;
}

QSplitter::handle {
    background: #45475a;
}

QSplitter::handle:vertical {
    height: 5px;
}

QSplitter::handle:hover {
    background: #89b4fa;
}

QSlider::groove:horizontal {
    background: #45475a;
    height: 4px;
    border-radius: 2px;
}

QSlider::handle:horizontal {
    background: #89b4fa;
    width: 12px;
    height: 12px;
    margin: -4px 0;
    border-radius: 6px;
}

QSlider::handle:horizontal:hover {
    background: #cdd6f4;
}

QPushButton[checkable="true"] {
    background: #313244;
    color: #cdd6f4;
}

QPushButton[checkable="true"]:checked {
    background: #89b4fa;
    color: #1e1e2e;
    border-color: #89b4fa;
}

QComboBox {
    background-color: #313244;
    color: #cdd6f4;
    border: 1px solid #45475a;
    border-radius: 3px;
    padding: 2px 6px;
    font-size: 9pt;
}

QComboBox:hover { border-color: #89b4fa; }

QComboBox QAbstractItemView {
    background-color: #181825;
    color: #cdd6f4;
    selection-background-color: #89b4fa;
    selection-color: #1e1e2e;
}

QLineEdit {
    background-color: #313244;
    color: #cdd6f4;
    border: 1px solid #45475a;
    border-radius: 3px;
    padding: 2px 4px;
    font-size: 9pt;
}

QLineEdit:focus { border-color: #89b4fa; }

NavigationToolbar2QT, QToolBar {
    background-color: #cdd6f4;
    border-bottom: 1px solid #45475a;
    spacing: 2px;
    padding: 2px;
}

NavigationToolbar2QT QToolButton, QToolBar QToolButton {
    background: transparent;
    border: 1px solid transparent;
    border-radius: 3px;
    padding: 2px;
    color: #1e1e2e;
}

NavigationToolbar2QT QToolButton:hover, QToolBar QToolButton:hover {
    background: #a6b0d8;
    border-color: #89b4fa;
}

NavigationToolbar2QT QToolButton:checked, QToolBar QToolButton:checked {
    background: #89b4fa;
    border-color: #89b4fa;
}

NavigationToolbar2QT QLabel, QToolBar QLabel {
    color: #1e1e2e;
    font-size: 8pt;
}
"""


# ---------------------------------------------------------------------------
# Batch extraction thread
# ---------------------------------------------------------------------------

class BatchExtractThread(QThread):
    """
    Background thread for batch-extracting spectra from a list of regions.

    Signals
    -------
    progress(int, int, str)       — (current, total, region_name)
    finished(list, ...)           — (results, cube_ref, method, load_strip,
                                      save_folder, fname_style, prefix,
                                      open_multispec, errors)
    error_log(list[str])          — accumulated per-region error messages
    """
    progress  = _Signal(int, int, str)
    finished  = _Signal(object)   # emits dict with all results
    error_log = _Signal(list)

    def __init__(self, cube, regions, method, weighting,
                 load_strip, save_folder, fname_style, prefix, open_multispec,
                 wcs=None, parent=None):
        super().__init__(parent)
        self._cube          = cube
        self._regions       = regions
        self._method        = method
        self._weighting     = weighting
        self._load_strip    = load_strip
        self._save_folder   = save_folder
        self._fname_style   = fname_style
        self._prefix        = prefix
        self._open_multispec = open_multispec
        self._wcs           = wcs
        self._cancel        = False

    def cancel(self):
        self._cancel = True

    def run(self):
        from rbcodes.GUIs.ifuviewer.processing.spatial_mask import (
            region_to_mask, region_center_sky)
        from rbcodes.GUIs.ifuviewer.processing.aperture_extract import (
            extract_aperture, extract_with_method)

        cube    = self._cube
        results = []
        errors  = []
        total   = len(self._regions)

        for i, reg in enumerate(self._regions):
            if self._cancel:
                break
            name = reg.get('name') or f"region_{i+1:03d}"
            self.progress.emit(i + 1, total, name)
            try:
                mask = region_to_mask(reg, self._wcs, cube.ny, cube.nx)
                if not mask.any():
                    continue
                flux_arr = extract_with_method(cube.flux, mask, self._method, var=cube.var)
                _, err_arr = extract_aperture(cube.flux, cube.var, mask)
                ra, dec = region_center_sky(reg, self._wcs)
                results.append((name, flux_arr, err_arr, ra, dec))
            except Exception as e:
                errors.append(f"{name}: {e}")

        if errors:
            self.error_log.emit(errors)

        self.finished.emit({
            'results':       results,
            'cube':          cube,
            'method':        self._method,
            'load_strip':    self._load_strip,
            'save_folder':   self._save_folder or '',
            'fname_style':   self._fname_style,
            'prefix':        self._prefix,
            'open_multispec': self._open_multispec,
            'errors':        errors,
        })


# ---------------------------------------------------------------------------
# Import regions dialog
# ---------------------------------------------------------------------------

class ImportRegionsDialog(QDialog):
    """Simple dialog for choosing which ds9 regions to import."""

    def __init__(self, source='ds9', filename='', parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import Regions")
        self.setMinimumWidth(300)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 8)

        if source == 'ds9':
            layout.addWidget(QLabel("Import regions from ds9:"))
            self._sel_rb = QRadioButton("Selected region(s) only")
            self._all_rb = QRadioButton("All regions")
            self._all_rb.setChecked(True)
            layout.addWidget(self._sel_rb)
            layout.addWidget(self._all_rb)
        else:
            layout.addWidget(QLabel(f"Load regions from:\n{filename}"))
            self._sel_rb = None
            self._all_rb = None

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def selected_only(self):
        return self._sel_rb is not None and self._sel_rb.isChecked()


# ---------------------------------------------------------------------------
# Batch extract dialog
# ---------------------------------------------------------------------------

class BatchExtractDialog(QDialog):
    """
    Dialog for configuring a batch spectral extraction from region file or ds9.
    """

    def __init__(self, has_ds9=False, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Batch Extract from Regions")
        self.setMinimumWidth(400)
        self._has_ds9 = has_ds9
        self._reg_path = ''

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 8)
        layout.setSpacing(8)

        # --- Source ---
        src_grp = QGroupBox("Source")
        src_lay = QVBoxLayout(src_grp)

        self._src_file_rb = QRadioButton("Region file")
        self._src_file_rb.setChecked(True)
        self._src_file_rb.toggled.connect(self._on_source_changed)
        src_lay.addWidget(self._src_file_rb)

        file_row = QHBoxLayout()
        self._file_edit = QLineEdit()
        self._file_edit.setPlaceholderText("path/to/regions.reg")
        browse_btn = QPushButton("Browse…")
        browse_btn.setFixedWidth(70)
        browse_btn.clicked.connect(self._browse_reg)
        file_row.addWidget(self._file_edit)
        file_row.addWidget(browse_btn)
        src_lay.addLayout(file_row)

        if has_ds9:
            self._src_ds9_sel_rb = QRadioButton("Live ds9 — selected region(s)")
            self._src_ds9_all_rb = QRadioButton("Live ds9 — all regions")
            self._src_ds9_sel_rb.toggled.connect(self._on_source_changed)
            self._src_ds9_all_rb.toggled.connect(self._on_source_changed)
            src_lay.addWidget(self._src_ds9_sel_rb)
            src_lay.addWidget(self._src_ds9_all_rb)

        layout.addWidget(src_grp)

        # --- Extraction method ---
        meth_row = QHBoxLayout()
        meth_row.addWidget(QLabel("Method:"))
        self._method_box = QComboBox()
        self._method_box.addItems(['sum', 'mean', 'median'])
        self._method_box.setFixedWidth(90)
        meth_row.addWidget(self._method_box)
        meth_row.addSpacing(16)
        meth_row.addWidget(QLabel("Weighting:"))
        self._weight_box = QComboBox()
        self._weight_box.addItems(['None', 'Var-weighted'])
        self._weight_box.setFixedWidth(110)
        meth_row.addWidget(self._weight_box)
        meth_row.addStretch()
        layout.addLayout(meth_row)

        # --- Output ---
        out_grp = QGroupBox("Output")
        out_lay = QVBoxLayout(out_grp)

        self._load_strip_cb = QCheckBox("Load into extract strip")
        self._load_strip_cb.setChecked(True)
        out_lay.addWidget(self._load_strip_cb)

        save_row = QHBoxLayout()
        self._save_cb = QCheckBox("Save to folder:")
        self._save_cb.setChecked(False)
        self._save_cb.toggled.connect(self._on_save_toggled)
        save_row.addWidget(self._save_cb)
        self._save_edit = QLineEdit()
        self._save_edit.setPlaceholderText("output folder")
        self._save_edit.setEnabled(False)
        save_row.addWidget(self._save_edit)
        save_browse = QPushButton("…")
        save_browse.setFixedWidth(28)
        save_browse.clicked.connect(self._browse_save_dir)
        save_row.addWidget(save_browse)
        out_lay.addLayout(save_row)

        fname_lay = QHBoxLayout()
        fname_lay.addWidget(QLabel("Filename:"))
        self._fname_coords_rb = QRadioButton("IAU coords")
        self._fname_idx_rb    = QRadioButton("Region index")
        self._fname_name_rb   = QRadioButton("Region name")
        self._fname_coords_rb.setChecked(True)
        fname_lay.addWidget(self._fname_coords_rb)
        fname_lay.addWidget(self._fname_idx_rb)
        fname_lay.addWidget(self._fname_name_rb)
        out_lay.addLayout(fname_lay)

        pfx_row = QHBoxLayout()
        pfx_row.addWidget(QLabel("Prefix:"))
        self._prefix_edit = QLineEdit("spec_")
        self._prefix_edit.setFixedWidth(90)
        pfx_row.addWidget(self._prefix_edit)
        pfx_row.addStretch()
        out_lay.addLayout(pfx_row)

        self._multispec_cb = QCheckBox("Open in rb_multispec when done")
        self._multispec_cb.setChecked(False)
        out_lay.addWidget(self._multispec_cb)

        layout.addWidget(out_grp)

        # --- Background ---
        self._bg_cb = QCheckBox("Run in background (non-blocking)")
        self._bg_cb.setChecked(True)
        layout.addWidget(self._bg_cb)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def _browse_reg(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Region File", "",
            "ds9 region files (*.reg);;All files (*)")
        if path:
            self._file_edit.setText(path)
            self._reg_path = path

    def _browse_save_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Select output folder")
        if d:
            self._save_edit.setText(d)

    def _on_source_changed(self):
        is_file = self._src_file_rb.isChecked()
        self._file_edit.setEnabled(is_file)

    def _on_save_toggled(self, checked):
        self._save_edit.setEnabled(checked)

    def values(self):
        """Return tuple of all settings."""
        if self._src_file_rb.isChecked():
            source = 'file'
            reg_ref = self._file_edit.text().strip()
        elif self._has_ds9 and self._src_ds9_sel_rb.isChecked():
            source  = 'ds9'
            reg_ref = 'selected'
        else:
            source  = 'ds9'
            reg_ref = 'all'

        method   = self._method_box.currentText()
        weighting = self._weight_box.currentText()
        load_strip = self._load_strip_cb.isChecked()
        save_folder = self._save_edit.text().strip() if self._save_cb.isChecked() else ''

        if self._fname_coords_rb.isChecked():
            fname_style = 'coords'
        elif self._fname_name_rb.isChecked():
            fname_style = 'name'
        else:
            fname_style = 'index'

        prefix = self._prefix_edit.text().strip() or 'spec_'
        open_multispec = self._multispec_cb.isChecked()
        run_bg = self._bg_cb.isChecked()

        return (source, reg_ref, method, weighting,
                load_strip, save_folder, fname_style, prefix,
                open_multispec, run_bg)


def launch_viewer(files=None):
    """
    Start the IFU Viewer GUI.

    Parameters
    ----------
    files : str or list of str or None
        FITS file(s) to load on startup.
    """
    app = QApplication.instance() or QApplication(sys.argv)
    win = MainWindow()
    win.setAttribute(Qt.WA_DeleteOnClose)
    if files:
        if isinstance(files, str):
            files = [files]
        for f in files:
            try:
                win._sidebar.add_dataset(f)
            except Exception as exc:
                print(f"[rb_ifuview] could not load '{f}': {exc}")
    app.setStyleSheet(_QSS)
    win.show()
    sys.exit(app.exec_())


class CropDialog(QDialog):
    """
    Dialog for entering crop coordinates as pixel (X/Y) or sky (RA/Dec).

    pixel_bounds() returns (x1, y1, x2, y2) in 0-indexed pixel coords
    suitable for numpy slicing.
    """

    def __init__(self, nx, ny, wcs=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Crop by Coordinates")
        self.setMinimumWidth(360)
        self._nx  = nx
        self._ny  = ny
        self._wcs = wcs

        layout = QVBoxLayout(self)
        form   = QFormLayout()
        layout.addLayout(form)

        # Input type selector
        self._mode_box = QComboBox()
        self._mode_box.addItem("Pixel (X / Y)")
        if wcs is not None:
            self._mode_box.addItem("Sky (RA / Dec)")
        self._mode_box.currentIndexChanged.connect(self._on_mode_changed)
        form.addRow("Input type:", self._mode_box)

        # ---- Pixel inputs ----
        self._px_widget = QWidget()
        px_form = QFormLayout(self._px_widget)
        px_form.setContentsMargins(0, 0, 0, 0)

        self._x1_spin = _int_spinbox(0, nx - 1, 0)
        self._x2_spin = _int_spinbox(1, nx,     nx)
        self._y1_spin = _int_spinbox(0, ny - 1, 0)
        self._y2_spin = _int_spinbox(1, ny,     ny)

        px_form.addRow("X  (start – end):", _hpair(self._x1_spin, self._x2_spin))
        px_form.addRow("Y  (start – end):", _hpair(self._y1_spin, self._y2_spin))
        form.addRow(self._px_widget)

        # ---- RA/Dec inputs ----
        self._rd_widget = QWidget()
        rd_form = QFormLayout(self._rd_widget)
        rd_form.setContentsMargins(0, 0, 0, 0)

        self._ra_spin  = QDoubleSpinBox(); self._ra_spin.setRange(0, 360); self._ra_spin.setDecimals(6); self._ra_spin.setSuffix(" deg"); self._ra_spin.setFixedWidth(140)
        self._dec_spin = QDoubleSpinBox(); self._dec_spin.setRange(-90, 90); self._dec_spin.setDecimals(6); self._dec_spin.setSuffix(" deg"); self._dec_spin.setFixedWidth(140)
        self._w_spin   = QDoubleSpinBox(); self._w_spin.setRange(0.1, 36000); self._w_spin.setDecimals(1); self._w_spin.setSuffix(" arcsec"); self._w_spin.setValue(30.0); self._w_spin.setFixedWidth(120)
        self._h_spin   = QDoubleSpinBox(); self._h_spin.setRange(0.1, 36000); self._h_spin.setDecimals(1); self._h_spin.setSuffix(" arcsec"); self._h_spin.setValue(30.0); self._h_spin.setFixedWidth(120)

        # Pre-fill RA/Dec from WCS centre pixel
        if wcs is not None:
            try:
                from astropy.wcs.utils import proj_plane_pixel_scales
                cx, cy = wcs.all_pix2world(nx / 2, ny / 2, 0)
                self._ra_spin.setValue(float(cx))
                self._dec_spin.setValue(float(cy))
                scales = proj_plane_pixel_scales(wcs) * 3600  # arcsec/pix
                self._w_spin.setValue(round(nx * scales[0], 1))
                self._h_spin.setValue(round(ny * scales[1], 1))
            except Exception:
                pass

        rd_form.addRow("RA center:",  self._ra_spin)
        rd_form.addRow("Dec center:", self._dec_spin)
        rd_form.addRow("Width:",      self._w_spin)
        rd_form.addRow("Height:",     self._h_spin)
        form.addRow(self._rd_widget)
        self._rd_widget.setVisible(False)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _on_mode_changed(self, idx):
        self._px_widget.setVisible(idx == 0)
        self._rd_widget.setVisible(idx == 1)

    def pixel_bounds(self):
        """Return (x1, y1, x2, y2) in 0-indexed pixels."""
        if self._mode_box.currentIndex() == 0 or self._wcs is None:
            x1 = self._x1_spin.value()
            x2 = self._x2_spin.value()
            y1 = self._y1_spin.value()
            y2 = self._y2_spin.value()
        else:
            # RA/Dec → pixel
            try:
                from astropy.wcs.utils import proj_plane_pixel_scales
                ra   = self._ra_spin.value()
                dec  = self._dec_spin.value()
                w_as = self._w_spin.value()
                h_as = self._h_spin.value()
                cx, cy = self._wcs.all_pix2world(ra, dec, 0)
                # pixel scale
                scales = proj_plane_pixel_scales(self._wcs) * 3600  # arcsec/pix
                hw = w_as / (2 * scales[0])
                hh = h_as / (2 * scales[1])
                # all_world2pix wants (ra, dec) → (x, y)
                px, py = self._wcs.all_world2pix(ra, dec, 0)
                x1, x2 = int(px - hw), int(px + hw)
                y1, y2 = int(py - hh), int(py + hh)
            except Exception as exc:
                QMessageBox.warning(self, "Coordinate error", str(exc))
                return 0, 0, self._nx, self._ny
        return int(x1), int(y1), int(x2), int(y2)


def _int_spinbox(lo, hi, value):
    sb = QSpinBox()
    sb.setRange(lo, hi)
    sb.setValue(value)
    sb.setFixedWidth(80)
    return sb


def _hpair(w1, w2):
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(4)
    lay.addWidget(w1)
    lay.addWidget(QLabel("–"))
    lay.addWidget(w2)
    return w


def _wave_spinbox(wave, value):
    """QDoubleSpinBox pre-configured for wavelength entry."""
    sb = QDoubleSpinBox()
    sb.setRange(float(wave[0]), float(wave[-1]))
    sb.setDecimals(2)
    sb.setSuffix(" Å")
    sb.setValue(float(np.clip(value, wave[0], wave[-1])))
    sb.setFixedWidth(110)
    return sb


def _range_widget(spin_min, spin_max):
    """Horizontal widget: [spin_min] – [spin_max]."""
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(4)
    lay.addWidget(spin_min)
    lay.addWidget(QLabel("–"))
    lay.addWidget(spin_max)
    return w


class MomentMapDialog(QDialog):
    """
    Dialog for computing moment maps (M0, M1, M2).

    Ctrl+M shortcut from MainWindow.
    Class attributes persist last-used values across dialog opens.
    """

    _last_order    = 0
    _last_wmin     = None
    _last_wmax     = None
    _last_lrest    = None
    _last_cont_sub = False
    _last_bcont    = (None, None)
    _last_rcont    = (None, None)
    _last_snr      = False
    _last_snr_thr  = 3.0
    _last_snr_meth = 'Sky region'

    def __init__(self, wave, onband=None,
                 has_variance=False, has_sky=False, has_cont=False,
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Moment Map")
        self.setMinimumWidth(380)

        self._wave        = wave
        self._onband      = onband
        self._has_variance = has_variance
        self._has_sky      = has_sky
        self._has_cont     = has_cont

        layout = QVBoxLayout(self)
        form   = QFormLayout()
        layout.addLayout(form)

        # ---- Moment order ----
        order_widget = QWidget()
        order_layout = QHBoxLayout(order_widget)
        order_layout.setContentsMargins(0, 0, 0, 0)
        self._order_group = QButtonGroup(self)
        for i, label in enumerate(["M0  (flux)", "M1  (velocity)", "M2  (dispersion)"]):
            rb = QRadioButton(label)
            if i == self._last_order:
                rb.setChecked(True)
            self._order_group.addButton(rb, i)
            order_layout.addWidget(rb)
        self._order_group.buttonClicked.connect(self._on_order_changed)
        form.addRow("Moment:", order_widget)

        # ---- Wavelength window ----
        # Priority: last used → current on-band → full wave range
        if self._last_wmin is not None:
            wmin_def = self._last_wmin
            wmax_def = self._last_wmax
        elif onband is not None and np.isfinite(onband[0]):
            wmin_def, wmax_def = float(onband[0]), float(onband[1])
        else:
            wmin_def = float(wave[0])
            wmax_def = float(wave[-1])

        self._wmin_spin = QDoubleSpinBox()
        self._wmin_spin.setRange(float(wave[0]), float(wave[-1]))
        self._wmin_spin.setDecimals(2)
        self._wmin_spin.setSuffix(" Å")
        self._wmin_spin.setValue(wmin_def)
        self._wmin_spin.setFixedWidth(120)

        self._wmax_spin = QDoubleSpinBox()
        self._wmax_spin.setRange(float(wave[0]), float(wave[-1]))
        self._wmax_spin.setDecimals(2)
        self._wmax_spin.setSuffix(" Å")
        self._wmax_spin.setValue(wmax_def)
        self._wmax_spin.setFixedWidth(120)

        wave_widget = QWidget()
        wave_layout = QHBoxLayout(wave_widget)
        wave_layout.setContentsMargins(0, 0, 0, 0)
        wave_layout.addWidget(self._wmin_spin)
        wave_layout.addWidget(QLabel("–"))
        wave_layout.addWidget(self._wmax_spin)

        use_onband_btn = QPushButton("Use current on-band")
        use_onband_btn.setFixedWidth(160)
        use_onband_btn.clicked.connect(self._fill_from_onband)
        wave_layout.addWidget(use_onband_btn)
        form.addRow("λ window:", wave_widget)

        # ---- Line center (M1 / M2 only) ----
        # This is the observed wavelength of the line that maps to v=0.
        # Must be in the same frame as the cube wavelength axis (observed frame).
        self._lrest_spin = QDoubleSpinBox()
        self._lrest_spin.setRange(100.0, 100000.0)
        self._lrest_spin.setDecimals(2)
        self._lrest_spin.setSuffix(" Å")
        self._lrest_spin.setToolTip(
            "Observed wavelength of the line = v 0 km/s reference.\n"
            "Must be in the observed frame (same as cube wavelength axis).\n"
            "Pre-filled with the center of the on-band window.")
        lrest_def = self._last_lrest if self._last_lrest is not None else (wmin_def + wmax_def) / 2.0
        self._lrest_spin.setValue(lrest_def)
        self._lrest_spin.setFixedWidth(120)
        self._lrest_label = QLabel("Line center (v=0):")
        form.addRow(self._lrest_label, self._lrest_spin)

        self._on_order_changed()   # set initial visibility

        # ---- Continuum subtraction ----
        self._cont_cb = QCheckBox("Subtract linear continuum")
        self._cont_cb.setToolTip(
            "Fit a per-spaxel linear baseline anchored in the blue and red\n"
            "continuum windows and subtract before computing moments.\n"
            "Leave unchecked for pure emission lines with no continuum.")
        self._cont_cb.setChecked(self._last_cont_sub)
        self._cont_cb.toggled.connect(self._on_cont_toggled)
        form.addRow("", self._cont_cb)

        # Blue continuum window
        wspan = wmax_def - wmin_def
        bc_lo = self._last_bcont[0] if self._last_bcont[0] is not None else max(float(wave[0]), wmin_def - wspan)
        bc_hi = self._last_bcont[1] if self._last_bcont[1] is not None else wmin_def - 2.0
        rc_lo = self._last_rcont[0] if self._last_rcont[0] is not None else wmax_def + 2.0
        rc_hi = self._last_rcont[1] if self._last_rcont[1] is not None else min(float(wave[-1]), wmax_def + wspan)

        self._bcont_min = _wave_spinbox(wave, bc_lo)
        self._bcont_max = _wave_spinbox(wave, bc_hi)
        self._rcont_min = _wave_spinbox(wave, rc_lo)
        self._rcont_max = _wave_spinbox(wave, rc_hi)

        self._bcont_widget = _range_widget(self._bcont_min, self._bcont_max)
        self._rcont_widget = _range_widget(self._rcont_min, self._rcont_max)
        self._bcont_label  = QLabel("Blue cont window:")
        self._rcont_label  = QLabel("Red cont window:")
        form.addRow(self._bcont_label, self._bcont_widget)
        form.addRow(self._rcont_label, self._rcont_widget)

        self._on_cont_toggled(self._last_cont_sub)   # set initial visibility

        # ---- SNR mask ----
        self._snr_cb = QCheckBox("Apply SNR mask")
        self._snr_cb.setToolTip(
            "Mask spaxels below the SNR threshold to NaN before displaying.\n"
            "SNR is computed from the moment-0 map divided by the noise estimate.")
        self._snr_cb.setChecked(self._last_snr)
        self._snr_cb.toggled.connect(self._on_snr_toggled)
        form.addRow("", self._snr_cb)

        snr_thr_widget = QWidget()
        snr_thr_layout = QHBoxLayout(snr_thr_widget)
        snr_thr_layout.setContentsMargins(0, 0, 0, 0)
        self._snr_thr_spin = QDoubleSpinBox()
        self._snr_thr_spin.setRange(0.1, 1000.0)
        self._snr_thr_spin.setDecimals(1)
        self._snr_thr_spin.setValue(self._last_snr_thr)
        self._snr_thr_spin.setFixedWidth(80)
        snr_thr_layout.addWidget(self._snr_thr_spin)

        self._snr_meth_box = QComboBox()
        snr_methods = []
        if has_variance: snr_methods.append('Variance cube')
        if has_sky:      snr_methods.append('Sky region')
        if has_cont:     snr_methods.append('Cont windows')
        if not snr_methods:
            snr_methods = ['(none available)']
            self._snr_cb.setEnabled(False)
            self._snr_cb.setToolTip(
                "No noise estimate available.\n"
                "Load a variance cube, set a sky region, or define continuum windows.")
        self._snr_meth_box.addItems(snr_methods)
        if self._last_snr_meth in snr_methods:
            self._snr_meth_box.setCurrentText(self._last_snr_meth)
        self._snr_meth_box.setFixedWidth(130)
        snr_thr_layout.addWidget(QLabel("σ  Method:"))
        snr_thr_layout.addWidget(self._snr_meth_box)
        snr_thr_layout.addStretch()

        self._snr_thr_label  = QLabel("Threshold:")
        self._snr_thr_widget = snr_thr_widget
        form.addRow(self._snr_thr_label, self._snr_thr_widget)

        self._on_snr_toggled(self._last_snr)

        # ---- Buttons ----
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    # ------------------------------------------------------------------

    def _on_order_changed(self, *_):
        needs_rest = self._order_group.checkedId() in (1, 2)
        self._lrest_label.setVisible(needs_rest)
        self._lrest_spin.setVisible(needs_rest)

    def _on_cont_toggled(self, checked):
        for w in (self._bcont_label, self._bcont_widget,
                  self._rcont_label, self._rcont_widget):
            w.setVisible(checked)

    def _on_snr_toggled(self, checked):
        for w in (self._snr_thr_label, self._snr_thr_widget):
            w.setVisible(checked)

    def _fill_from_onband(self):
        if self._onband is not None and np.isfinite(self._onband[0]):
            self._wmin_spin.setValue(float(self._onband[0]))
            self._wmax_spin.setValue(float(self._onband[1]))

    def values(self):
        """
        Return (order, wmin, wmax, lambda_rest_or_None,
                subtract_cont, bcont_min, bcont_max, rcont_min, rcont_max)
        and save to class for next open.
        """
        order  = self._order_group.checkedId()
        wmin   = self._wmin_spin.value()
        wmax   = self._wmax_spin.value()
        if wmin > wmax:
            wmin, wmax = wmax, wmin
        lrest       = self._lrest_spin.value() if order in (1, 2) else None
        subtract    = self._cont_cb.isChecked()
        bcont_min   = self._bcont_min.value()
        bcont_max   = self._bcont_max.value()
        rcont_min   = self._rcont_min.value()
        rcont_max   = self._rcont_max.value()
        if bcont_min > bcont_max:
            bcont_min, bcont_max = bcont_max, bcont_min
        if rcont_min > rcont_max:
            rcont_min, rcont_max = rcont_max, rcont_min
        # Persist for next open
        apply_snr  = self._snr_cb.isChecked() and self._snr_cb.isEnabled()
        snr_thresh = self._snr_thr_spin.value()
        snr_method = self._snr_meth_box.currentText()

        MomentMapDialog._last_order    = order
        MomentMapDialog._last_wmin     = wmin
        MomentMapDialog._last_wmax     = wmax
        MomentMapDialog._last_lrest    = lrest
        MomentMapDialog._last_cont_sub = subtract
        MomentMapDialog._last_bcont    = (bcont_min, bcont_max)
        MomentMapDialog._last_rcont    = (rcont_min, rcont_max)
        MomentMapDialog._last_snr      = apply_snr
        MomentMapDialog._last_snr_thr  = snr_thresh
        MomentMapDialog._last_snr_meth = snr_method
        return (order, wmin, wmax, lrest,
                subtract, bcont_min, bcont_max, rcont_min, rcont_max,
                apply_snr, snr_thresh, snr_method)


class HeaderDialog(QDialog):
    """
    Show FITS header(s) for the active dataset.

    If the source file is available and has multiple extensions, a dropdown
    lets the user switch between them.  For cropped cubes (no file path) or
    when the file is unavailable, the in-memory header is shown directly.
    """

    def __init__(self, cube, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"FITS Header — {cube.name}")
        self.resize(700, 560)

        self._headers = self._collect_headers(cube)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # Extension selector (hidden when only one header)
        self._ext_combo = QComboBox()
        for label, _ in self._headers:
            self._ext_combo.addItem(label)
        if len(self._headers) > 1:
            self._ext_combo.currentIndexChanged.connect(self._on_ext_changed)
            layout.addWidget(self._ext_combo)

        # Header text area
        self._text = QPlainTextEdit()
        self._text.setReadOnly(True)
        self._text.setFont(_mono_font())
        self._text.setLineWrapMode(QPlainTextEdit.NoWrap)
        layout.addWidget(self._text)

        # Search bar
        search_row = QHBoxLayout()
        self._search_box   = QLineEdit()
        self._search_box.setPlaceholderText("Search header…")
        self._search_box.returnPressed.connect(self._find_next)
        self._match_label  = QLabel("")
        prev_btn = QPushButton("Prev")
        next_btn = QPushButton("Next")
        prev_btn.setFixedWidth(48)
        next_btn.setFixedWidth(48)
        prev_btn.clicked.connect(self._find_prev)
        next_btn.clicked.connect(self._find_next)
        self._search_box.textChanged.connect(self._on_search_changed)
        search_row.addWidget(self._search_box)
        search_row.addWidget(prev_btn)
        search_row.addWidget(next_btn)
        search_row.addWidget(self._match_label)
        layout.addLayout(search_row)

        # Buttons
        btn_row = QHBoxLayout()
        copy_btn  = QPushButton("Copy to Clipboard")
        close_btn = QPushButton("Close")
        copy_btn.clicked.connect(self._copy)
        close_btn.clicked.connect(self.accept)
        btn_row.addWidget(copy_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

        self._show_header(0)

    # ------------------------------------------------------------------

    @staticmethod
    def _collect_headers(cube):
        """
        Return list of (label, header) pairs.

        Tries to open cube.path first; falls back to cube.header.
        """
        from astropy.io import fits as _fits
        import os

        headers = []
        path = getattr(cube, 'path', '')
        if path and os.path.isfile(path):
            try:
                with _fits.open(path, memmap=False) as hdul:
                    for i, hdu in enumerate(hdul):
                        name = hdu.name or 'UNNAMED'
                        label = f"[{i}]  {name}"
                        if hasattr(hdu, 'data') and hdu.data is not None:
                            shape = getattr(hdu.data, 'shape', '')
                            if shape:
                                label += f"  {shape}"
                        headers.append((label, hdu.header.copy()))
                if headers:
                    return headers
            except Exception:
                pass

        # Fallback: in-memory header only
        hdr = getattr(cube, 'header', None)
        label = f"[0]  {cube.name}"
        if hdr is None:
            from astropy.io import fits as _fits
            hdr = _fits.Header()
        return [(label, hdr)]

    def _show_header(self, idx):
        _, hdr = self._headers[idx]
        lines = [f"{str(card)}" for card in hdr.cards]
        self._text.setPlainText("\n".join(lines))

    def _on_ext_changed(self, idx):
        self._show_header(idx)
        self._on_search_changed(self._search_box.text())

    def _on_search_changed(self, text):
        """Re-run search from top whenever the query changes."""
        self._match_label.setText("")
        if not text:
            # Clear any existing highlights by resetting extra selections
            self._text.setExtraSelections([])
            return
        self._highlight_all(text)
        # Move cursor to first match
        self._text.moveCursor(self._text.textCursor().Start)
        self._find_next()

    def _highlight_all(self, text):
        """Highlight all occurrences with a yellow background."""
        from PyQt5.QtGui import QTextCharFormat, QColor, QTextCursor
        selections = []
        fmt = QTextCharFormat()
        fmt.setBackground(QColor('#f9e2af'))   # warm yellow (Catppuccin latte)
        fmt.setForeground(QColor('#1e1e2e'))
        cursor = self._text.document().find(text)
        count = 0
        while not cursor.isNull():
            sel = QTextEdit.ExtraSelection()
            sel.cursor = cursor
            sel.format  = fmt
            selections.append(sel)
            count += 1
            cursor = self._text.document().find(text, cursor)
        self._text.setExtraSelections(selections)
        self._match_label.setText(f"{count} match{'es' if count != 1 else ''}" if count else "Not found")

    def _find_next(self):
        text = self._search_box.text()
        if not text:
            return
        found = self._text.find(text)
        if not found:   # wrap around
            self._text.moveCursor(self._text.textCursor().Start)
            self._text.find(text)

    def _find_prev(self):
        from PyQt5.QtGui import QTextDocument
        text = self._search_box.text()
        if not text:
            return
        found = self._text.find(text, QTextDocument.FindBackward)
        if not found:   # wrap around
            self._text.moveCursor(self._text.textCursor().End)
            self._text.find(text, QTextDocument.FindBackward)

    def _copy(self):
        from PyQt5.QtWidgets import QApplication
        QApplication.clipboard().setText(self._text.toPlainText())


def _mono_font():
    from PyQt5.QtGui import QFont, QFontDatabase
    f = QFontDatabase.systemFont(QFontDatabase.FixedFont)
    f.setPointSize(10)
    return f


if __name__ == "__main__":
    from rbcodes.GUIs.ifuviewer.rb_ifuview import main
    main()
