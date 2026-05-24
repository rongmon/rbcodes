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
                              QGroupBox, QListWidget, QListWidgetItem)
from PyQt5.QtCore import Qt
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
        self._sky_mask             = None  # 2-D bool array — sky region for SNR
        self._sky_rect             = None  # (x1, y1, x2, y2) pixel bounds

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

        # Left: dataset sidebar
        self._sidebar = DatasetSidebar()
        self._sidebar.dataset_changed.connect(self._on_dataset_changed)
        self._sidebar.use_as_variance.connect(self._on_use_as_variance)
        root.addWidget(self._sidebar)

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
        self._extract_strip = self._build_extract_strip()

        # Fixed widgets above/below the splitter
        right.addWidget(self._mpl_toolbar)
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

        file_menu.addSeparator()

        save_frame_action = QAction("Save Current Frame…", self)
        save_frame_action.setShortcut("Ctrl+S")
        save_frame_action.setToolTip(
            "Save the currently displayed 2D image (whitelight, narrowband, "
            "moment map, …) as a FITS file with spatial WCS header.")
        save_frame_action.triggered.connect(self._on_save_current_frame)
        file_menu.addAction(save_frame_action)

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

    def _on_dataset_changed(self, cube):
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
            self._on_clear_extractions()
            return

        self._status_base = _status_text(cube)
        self.statusBar().showMessage(self._status_base)
        self._show_default_image(cube)

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
        # Warn if sky region overlaps crop — it won't be valid anymore
        if self._sky_rect is not None:
            sx1, sy1, sx2, sy2 = self._sky_rect
            if not (sx1 >= x1 and sx2 <= x2 and sy1 >= y1 and sy2 <= y2):
                QMessageBox.information(
                    self, "Sky region",
                    "The sky region extends outside the crop area and has been cleared.")
                self._sky_mask = None
                self._sky_rect = None
                self._image_canvas.clear_sky_region()
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

        # ---- Circular / Circular-Annular ---------------------------------
        else:
            radius   = ac.radius
            method   = ac.method
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

        rb_spec = _make_rb_spectrum(cube.wave, flux_arr, err_arr, label)
        self._spectrum_canvas.lock_spectrum(cube.wave, flux_arr, err_arr,
                                            label, rb_spec=rb_spec)
        self._register_extraction(label, rb_spec, new_markers)

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
        self._spectrum_canvas.lock_spectrum(cube.wave, flux_arr, err_arr,
                                            label, rb_spec=rb_spec)
        self._register_extraction(label, rb_spec, new_markers)

    def _register_extraction(self, label, rb_spec, markers=None):
        """Add entry to extraction combo and enable strip controls."""
        self._extract_combo.addItem(label)
        self._extract_combo.setCurrentIndex(self._extract_combo.count() - 1)
        self._extraction_marker_groups.append(list(markers or []))
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

    def _on_clear_extractions(self):
        self._spectrum_canvas.clear_locked()
        self._image_canvas.clear_extraction_marks()
        self._extract_combo.clear()
        self._extraction_marker_groups.clear()
        self._set_extract_strip_enabled(False)

    def _on_spaxel_hovered(self, x, y, spectrum):
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

    def _sync_spec_range_strip(self):
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
            self._extract_strip.show()
            self._set_extract_strip_enabled(False)   # empty state
        else:
            self._spectrum_canvas.hide()
            self._splitter.setSizes([800, 0])
            self._aperture_controls.hide()
            self._extract_strip.hide()

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


def launch_viewer(fitsfile=None):
    app = QApplication.instance() or QApplication(sys.argv)
    win = MainWindow()
    win.setAttribute(Qt.WA_DeleteOnClose)
    if fitsfile:
        try:
            win._sidebar.add_dataset(fitsfile)
        except RuntimeError as exc:
            print(f"[ifuviewer] {exc}")
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


def main():
    import argparse
    parser = argparse.ArgumentParser(description="IFU Viewer")
    parser.add_argument("fitsfile", nargs="?", default=None,
                        help="FITS cube or image to open on startup")
    args = parser.parse_args()
    launch_viewer(fitsfile=args.fitsfile)


if __name__ == "__main__":
    main()
