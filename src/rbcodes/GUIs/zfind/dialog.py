"""
dialog.py — PyQt5 popup dialog for rb_zfind.

Figure layout:
  Row 0 — chi2 / significance vs z   (click anywhere to set z)
  Row 1 — spectrum + linelist overlay at the selected z

Spectrum panel keystrokes  (cursor must be over the spectrum panel):
  x / X   — set left / right x limit to cursor wavelength
  t / b   — set top / bottom y limit to cursor flux
  r       — reset both axes to full range
  o       — zoom out x-axis (1.5×)
  [ / ]   — pan left / right one window width
  S / U   — increase / decrease boxcar smoothing
  a       — autoscale y to the visible x range

Signals
-------
accepted_z         : list  [z, z_err]              → zgui widget_z._on_estZ_changed
accepted_absorbers : list  [{zabs, name, label}]   → rb_multispec AbsorberManager
"""

import numpy as np

try:
    from PyQt5.QtWidgets import (
        QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
        QLabel, QLineEdit, QPushButton, QComboBox,
        QRadioButton, QGroupBox, QButtonGroup,
        QTableWidget, QTableWidgetItem,
        QStackedWidget, QApplication, QSizePolicy,
        QHeaderView,
    )
    from PyQt5.QtCore import pyqtSignal, Qt
    from PyQt5.QtGui import QDoubleValidator, QIntValidator

    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
    import matplotlib.gridspec as gridspec
    from astropy.convolution import convolve, Box1DKernel

    # Curated presets first, then full rb_setline lists for power users
    _LINELIST_NAMES = [
        # -- built-in curated (recommended) --
        'zfind_em', 'zfind_abs',
        # -- full rb_setline lists --
        'Gal_Em', 'Gal_Abs', 'Gal_long', 'Gal',
        'LLS', 'LLS Small', 'DLA', 'LBG',
        'AGN', 'Eiger_Strong',
        'HI', 'HI_recomb', 'HI_recomb_light',
        'EUV', 'LLS_EUV', 'atom',
    ]
    _DEFAULT_LINELIST = {'emission': 'zfind_em', 'absorption': 'zfind_abs'}
    _SOL_COLORS = ['#e74c3c', '#e67e22', '#27ae60', '#8e44ad', '#2980b9',
                   '#c0392b', '#d35400', '#16a085', '#7f8c8d', '#2c3e50']


    def _linelist_to_df(linelist_name):
        """Return a DataFrame for *linelist_name*.

        Curated presets ('zfind_em', 'zfind_abs') are served from the
        built-in registry in linelists.py.  All other names are passed
        through to rb_setline.read_line_list.
        """
        from rbcodes.GUIs.zfind.linelists import CURATED_NAMES, get_curated_df
        if linelist_name in CURATED_NAMES:
            return get_curated_df(linelist_name)

        import pandas as pd
        from rbcodes.IGM.rb_setline import read_line_list
        raw = read_line_list(linelist_name)
        if not raw:
            raise ValueError(f'Empty linelist: {linelist_name!r}')
        if any(map(str.isdigit, raw[0]['ion'])):
            rows = [{'wave': li['wrest'], 'name': li['ion']} for li in raw]
        else:
            rows = [{'wave': li['wrest'],
                     'name': li['ion'] + ' ' + str(round(li['wrest']))}
                    for li in raw]
        df = pd.DataFrame(rows)
        df.attrs['name'] = linelist_name
        return df


    class ZFindDialog(QDialog):
        """
        Self-contained redshift / absorber finder popup.

        Parameters
        ----------
        spec            : rb_spectrum or None
        mode            : 'emission' or 'absorption'
        default_linelist: str or None
        z_qso           : float or None  — caps z_max in absorption mode
        parent          : QWidget or None
        """

        accepted_z         = pyqtSignal(list)
        accepted_absorbers = pyqtSignal(list)

        def __init__(self, spec=None, mode='emission', default_linelist=None,
                     z_qso=None, parent=None):
            super().__init__(parent)
            self.spec             = spec
            self.mode             = mode
            self.default_linelist = default_linelist
            self.z_qso            = z_qso

            self._result          = None
            self._df              = None       # linelist after last Run
            self._selected_em_idx = 0
            self._selected_z      = None       # z currently shown in spectrum

            # smoothing state
            self._smooth_scale    = 1
            self._raw_flux        = None       # unsmoothed flux array
            self._raw_err         = None       # unsmoothed error array

            self.setWindowTitle('Redshift Estimator' if mode == 'emission'
                                else 'Absorber Finder')
            self.resize(960, 880)
            self.setMinimumSize(780, 680)

            self._build_ui()
            self._setup_connections()
            self._set_initial_state()

        # ------------------------------------------------------------------
        # UI
        # ------------------------------------------------------------------

        def _build_ui(self):
            outer = QVBoxLayout(self)
            outer.setSpacing(5)

            # ---- Mode ----
            mode_box = QGroupBox('Mode')
            mode_row = QHBoxLayout(mode_box)
            self._rb_emission   = QRadioButton('Emission')
            self._rb_absorption = QRadioButton('Absorption')
            self._mode_group    = QButtonGroup(self)
            self._mode_group.addButton(self._rb_emission,   0)
            self._mode_group.addButton(self._rb_absorption, 1)
            (self._rb_emission if self.mode == 'emission'
             else self._rb_absorption).setChecked(True)
            mode_row.addWidget(self._rb_emission)
            mode_row.addWidget(self._rb_absorption)
            mode_row.addStretch()
            outer.addWidget(mode_box)

            # ---- Parameters ----
            param_box  = QGroupBox('Parameters')
            param_grid = QGridLayout(param_box)

            param_grid.addWidget(QLabel('Linelist:'), 0, 0)
            self._linelist_combo = QComboBox()
            self._linelist_combo.setMinimumWidth(120)
            self._linelist_combo.addItems(_LINELIST_NAMES)
            param_grid.addWidget(self._linelist_combo, 0, 1)

            param_grid.addWidget(QLabel('z min:'), 0, 2)
            self._z_min = QLineEdit('0.0')
            self._z_min.setFixedWidth(65)
            self._z_min.setValidator(QDoubleValidator(0.0, 30.0, 6))
            param_grid.addWidget(self._z_min, 0, 3)

            param_grid.addWidget(QLabel('z max:'), 0, 4)
            self._z_max = QLineEdit('6.0')
            self._z_max.setFixedWidth(65)
            self._z_max.setValidator(QDoubleValidator(0.0, 30.0, 6))
            param_grid.addWidget(self._z_max, 0, 5)

            param_grid.addWidget(QLabel('n steps:'), 0, 6)
            self._n_steps = QLineEdit('5000')
            self._n_steps.setFixedWidth(65)
            self._n_steps.setValidator(QIntValidator(100, 200000))
            param_grid.addWidget(self._n_steps, 0, 7)

            param_grid.addWidget(QLabel('FWHM (Å):'), 1, 0)
            self._fwhm = QLineEdit('0.0')
            self._fwhm.setFixedWidth(65)
            self._fwhm.setValidator(QDoubleValidator(0.0, 1000.0, 3))
            param_grid.addWidget(self._fwhm, 1, 1)

            param_grid.addWidget(QLabel('Win (pix):'), 1, 2)
            self._win_pix = QLineEdit('5')
            self._win_pix.setFixedWidth(45)
            self._win_pix.setValidator(QIntValidator(1, 200))
            self._win_pix.setToolTip(
                'Half-window in pixels around each line for chi2 evaluation')
            param_grid.addWidget(self._win_pix, 1, 3)

            param_grid.addWidget(QLabel('Smooth (pix):'), 1, 4)
            self._smooth_pix = QLineEdit('3')
            self._smooth_pix.setFixedWidth(45)
            self._smooth_pix.setValidator(QIntValidator(1, 99))
            self._smooth_pix.setToolTip(
                'Boxcar kernel width (pixels) applied to flux before chi2 scan. '
                'Set to 1 to disable.')
            param_grid.addWidget(self._smooth_pix, 1, 5)

            param_grid.addWidget(QLabel('Method:'), 1, 6)
            self._method_combo = QComboBox()
            self._method_combo.addItems(['Line Search', 'Template', 'PCA'])
            param_grid.addWidget(self._method_combo, 1, 7)

            # Template selector — shown only when Template method is active
            from rbcodes.GUIs.zfind.engine import TEMPLATE_NAMES, PCA_NAMES
            self._template_combo = QComboBox()
            self._template_combo.addItems(['Best of all'] + TEMPLATE_NAMES)
            self._template_combo.setVisible(False)
            self._template_combo.setToolTip(
                '"Best of all" tries every MARZ template and picks the lowest chi2/dof')
            param_grid.addWidget(self._template_combo, 1, 8)

            # PCA selector — shown only when PCA method is active
            self._pca_combo = QComboBox()
            _pca_display = ['Best of all'] + PCA_NAMES
            self._pca_combo.addItems(_pca_display)
            self._pca_combo.setVisible(False)
            self._pca_combo.setToolTip(
                'DESI redrock PCA eigenvectors.\n'
                '"Best of all" runs galaxy + qso and picks lowest chi2/dof.\n'
                'Download first: python templates/download_desi_pca.py')
            param_grid.addWidget(self._pca_combo, 1, 8)

            self._run_btn = QPushButton('Run Search')
            self._run_btn.setMinimumWidth(100)
            param_grid.addWidget(self._run_btn, 1, 9)

            # Row 2: overlay linelist | wave range | data norm | model norm
            param_grid.addWidget(QLabel('Overlay lines:'), 2, 0)
            self._overlay_combo = QComboBox()
            self._overlay_combo.setMinimumWidth(120)
            self._overlay_combo.addItems(_LINELIST_NAMES)
            self._overlay_combo.setToolTip(
                'Linelist drawn on the spectrum panel at the selected z.\n'
                'In Template/PCA mode this is the only linelist shown.')
            param_grid.addWidget(self._overlay_combo, 2, 1)

            param_grid.addWidget(QLabel('λ min (Å):'), 2, 2)
            self._wave_min = QLineEdit('')
            self._wave_min.setFixedWidth(65)
            self._wave_min.setValidator(QDoubleValidator(0.0, 1e6, 2))
            self._wave_min.setPlaceholderText('auto')
            self._wave_min.setToolTip(
                'Observed-frame wavelength lower limit for chi2 scan.\n'
                'Pixels below this are masked (ivar→0). Leave blank for full range.')
            param_grid.addWidget(self._wave_min, 2, 3)

            param_grid.addWidget(QLabel('λ max (Å):'), 2, 4)
            self._wave_max = QLineEdit('')
            self._wave_max.setFixedWidth(65)
            self._wave_max.setValidator(QDoubleValidator(0.0, 1e6, 2))
            self._wave_max.setPlaceholderText('auto')
            self._wave_max.setToolTip(
                'Observed-frame wavelength upper limit for chi2 scan.\n'
                'Pixels above this are masked (ivar→0). Leave blank for full range.')
            param_grid.addWidget(self._wave_max, 2, 5)

            param_grid.addWidget(QLabel('Data norm:'), 2, 6)
            self._data_norm_combo = QComboBox()
            self._data_norm_combo.addItems(['normalize', 'subtract', 'raw'])
            self._data_norm_combo.setToolTip(
                'How to preprocess the observed spectrum before chi2 comparison:\n'
                '  normalize  flux / continuum\n'
                '  subtract   flux − continuum\n'
                '  raw        flux as-is (no continuum removal)')
            param_grid.addWidget(self._data_norm_combo, 2, 7)

            param_grid.addWidget(QLabel('Model norm:'), 2, 8)
            self._model_norm_combo = QComboBox()
            self._model_norm_combo.addItems(['normalize', 'subtract', 'raw'])
            self._model_norm_combo.setToolTip(
                'How to preprocess the template / PCA eigenvectors:\n'
                '  normalize  template / pseudo-continuum  (MARZ default)\n'
                '             PCA: L2-normalise each eigenvector\n'
                '  subtract   template − pseudo-continuum\n'
                '             PCA: mean-center each eigenvector\n'
                '  raw        template as-is\n'
                'N/A for Line Search (greyed out).')
            param_grid.addWidget(self._model_norm_combo, 2, 9)

            outer.addWidget(param_box)

            # ---- Two-panel figure ----
            self._fig = Figure(figsize=(9, 6.5))
            gs = gridspec.GridSpec(2, 1, figure=self._fig,
                                   height_ratios=[2, 3], hspace=0.42,
                                   top=0.95, bottom=0.08,
                                   left=0.08, right=0.97)
            self._ax      = self._fig.add_subplot(gs[0])
            self._ax_spec = self._fig.add_subplot(gs[1])

            self._canvas = FigureCanvasQTAgg(self._fig)
            self._canvas.setSizePolicy(QSizePolicy.Expanding,
                                       QSizePolicy.Expanding)
            self._canvas.setMinimumHeight(360)
            # Canvas must accept focus so key events reach it
            self._canvas.setFocusPolicy(Qt.StrongFocus)
            outer.addWidget(self._canvas, stretch=4)

            # ---- Results (stacked) ----
            outer.addWidget(QLabel('Results:'))
            self._stack = QStackedWidget()

            self._em_table = QTableWidget(0, 4)
            self._em_table.setHorizontalHeaderLabels(
                ['z', 'z_err', 'chi2/dof', 'n_features'])
            self._em_table.setSelectionBehavior(QTableWidget.SelectRows)
            self._em_table.setSelectionMode(QTableWidget.SingleSelection)
            self._em_table.setEditTriggers(QTableWidget.NoEditTriggers)
            self._em_table.horizontalHeader().setSectionResizeMode(
                QHeaderView.Stretch)
            self._em_table.setMaximumHeight(120)
            self._stack.addWidget(self._em_table)

            self._abs_table = QTableWidget(0, 5)
            self._abs_table.setHorizontalHeaderLabels(
                ['Accept', 'z', 'Significance', 'n_lines', 'Lines matched'])
            self._abs_table.setSelectionBehavior(QTableWidget.SelectRows)
            self._abs_table.setSelectionMode(QTableWidget.SingleSelection)
            self._abs_table.setEditTriggers(QTableWidget.NoEditTriggers)
            hdr = self._abs_table.horizontalHeader()
            hdr.setSectionResizeMode(0, QHeaderView.ResizeToContents)
            hdr.setSectionResizeMode(4, QHeaderView.Stretch)
            self._abs_table.setMaximumHeight(150)
            self._stack.addWidget(self._abs_table)

            outer.addWidget(self._stack, stretch=1)

            # ---- Status ----
            self._status = QLabel('Ready.   '
                                  'Keys (hover either panel): '
                                  'x/X=xlim  t/b=ylim  r=reset  '
                                  'o=zoom-out  [/]=pan  a=autoscale-y  '
                                  'S/U=smooth (spectrum only)  scroll=zoom')
            self._status.setWordWrap(True)
            outer.addWidget(self._status)

            # ---- Buttons ----
            btn_row = QHBoxLayout()
            self._accept_btn = QPushButton('Accept')
            self._accept_btn.setEnabled(False)
            self._cancel_btn = QPushButton('Cancel')
            btn_row.addStretch()
            btn_row.addWidget(self._accept_btn)
            btn_row.addWidget(self._cancel_btn)
            outer.addLayout(btn_row)

        # ------------------------------------------------------------------
        # Connections
        # ------------------------------------------------------------------

        def _setup_connections(self):
            self._run_btn.clicked.connect(self._on_run)
            self._accept_btn.clicked.connect(self._on_accept)
            self._cancel_btn.clicked.connect(self.reject)
            self._rb_emission.toggled.connect(self._on_mode_toggled)
            self._method_combo.currentIndexChanged.connect(
                self._on_method_changed)
            self._em_table.itemSelectionChanged.connect(
                self._on_em_selection_changed)
            self._abs_table.itemSelectionChanged.connect(
                self._on_abs_selection_changed)
            self._canvas.mpl_connect('button_press_event',  self._on_plot_click)
            self._canvas.mpl_connect('key_press_event',     self._on_key_press)
            self._canvas.mpl_connect('scroll_event',        self._on_scroll)

        # ------------------------------------------------------------------
        # Init
        # ------------------------------------------------------------------

        def _set_initial_state(self):
            target = (self.default_linelist
                      if self.default_linelist in _LINELIST_NAMES
                      else _DEFAULT_LINELIST.get(self.mode, 'Gal_Em'))
            idx = self._linelist_combo.findText(target)
            if idx >= 0:
                self._linelist_combo.setCurrentIndex(idx)
            # Sync overlay combo to same default
            oidx = self._overlay_combo.findText(target)
            if oidx >= 0:
                self._overlay_combo.setCurrentIndex(oidx)

            if self.z_qso is not None and self.mode == 'absorption':
                self._z_max.setText(f'{self.z_qso:.4f}')

            self._update_stack()
            self._reset_plots()

            if self.spec is None:
                self._status.setText(
                    'No spectrum loaded. Pass a spectrum to launch_zfind().')
                self._run_btn.setEnabled(False)
            else:
                wave = self.spec.wavelength.to('AA').value
                self._status.setText(
                    f'Spectrum loaded: {len(wave):,} px, '
                    f'{wave.min():.0f}–{wave.max():.0f} Å.   '
                    'Keys (hover spectrum): x/X t/b r o [/] S/U a')
                self._draw_spectrum_only()

        # ------------------------------------------------------------------
        # Method switch  (Line Search ↔ Template ↔ PCA)
        # ------------------------------------------------------------------

        def _on_method_changed(self, index):
            method = self._method_combo.currentText()
            is_template = (method == 'Template')
            is_pca      = (method == 'PCA')
            is_line     = not is_template and not is_pca

            self._linelist_combo.setVisible(is_line)
            self._template_combo.setVisible(is_template)
            self._pca_combo.setVisible(is_pca)

            # model_norm is not meaningful for line search
            self._model_norm_combo.setEnabled(not is_line)

            if is_template:
                self._data_norm_combo.setCurrentText('normalize')
                self._model_norm_combo.setCurrentText('normalize')
                from rbcodes.GUIs.zfind.engine import _TEMPLATE_Z_RANGE
                tname = self._template_combo.currentText()
                if tname != 'Best of all' and tname in _TEMPLATE_Z_RANGE:
                    zlo, zhi = _TEMPLATE_Z_RANGE[tname]
                    self._z_min.setText(str(zlo))
                    self._z_max.setText(str(zhi))
            elif is_pca:
                self._data_norm_combo.setCurrentText('normalize')
                self._model_norm_combo.setCurrentText('normalize')
                from rbcodes.GUIs.zfind.engine import _PCA_Z_RANGE
                pname = self._pca_combo.currentText()
                if pname != 'Best of all' and pname in _PCA_Z_RANGE:
                    zlo, zhi = _PCA_Z_RANGE[pname]
                    self._z_min.setText(str(zlo))
                    self._z_max.setText(str(zhi))
            else:   # Line Search
                self._data_norm_combo.setCurrentText('subtract')
                self._model_norm_combo.setCurrentText('subtract')  # greyed out

        # ------------------------------------------------------------------
        # Mode switch
        # ------------------------------------------------------------------

        def _on_mode_toggled(self, checked):
            if not checked:
                return
            self.mode = ('emission' if self._rb_emission.isChecked()
                         else 'absorption')
            self._result = self._df = self._selected_z = None
            self._smooth_scale = 1
            self._accept_btn.setEnabled(False)
            self._em_table.setRowCount(0)
            self._abs_table.setRowCount(0)
            default = _DEFAULT_LINELIST.get(self.mode, 'Gal_Em')
            idx = self._linelist_combo.findText(default)
            if idx >= 0:
                self._linelist_combo.setCurrentIndex(idx)
            oidx = self._overlay_combo.findText(default)
            if oidx >= 0:
                self._overlay_combo.setCurrentIndex(oidx)
            if self.mode == 'absorption' and self.z_qso is not None:
                self._z_max.setText(f'{self.z_qso:.4f}')
            self._update_stack()
            self._reset_plots()
            if self.spec is not None:
                self._draw_spectrum_only()

        def _update_stack(self):
            self._stack.setCurrentIndex(0 if self.mode == 'emission' else 1)

        # ------------------------------------------------------------------
        # Run
        # ------------------------------------------------------------------

        def _on_run(self):
            if self.spec is None:
                self._status.setText('No spectrum loaded.')
                return
            try:
                z_min       = float(self._z_min.text())
                z_max       = float(self._z_max.text())
                n_steps     = int(self._n_steps.text())
                fwhm        = float(self._fwhm.text())
                win_pix     = int(self._win_pix.text())
                smooth_pix  = int(self._smooth_pix.text())
            except ValueError:
                self._status.setText('Invalid parameter values.')
                return
            if z_min >= z_max:
                self._status.setText('z_min must be < z_max.')
                return

            method      = self._method_combo.currentText()
            data_norm   = self._data_norm_combo.currentText()
            model_norm  = self._model_norm_combo.currentText()
            wave_min_s  = self._wave_min.text().strip()
            wave_max_s  = self._wave_max.text().strip()
            wave_min    = float(wave_min_s) if wave_min_s else None
            wave_max    = float(wave_max_s) if wave_max_s else None

            # Always load the overlay linelist so line markers are drawn
            # regardless of whether Line Search or Template mode is active.
            overlay_name = self._overlay_combo.currentText()
            try:
                overlay_df = _linelist_to_df(overlay_name)
            except Exception as exc:
                self._status.setText(
                    f'Cannot load overlay linelist "{overlay_name}": {exc}')
                return

            self._status.setText('Running…')
            self._run_btn.setEnabled(False)
            QApplication.processEvents()

            try:
                if method == 'Template':
                    from rbcodes.GUIs.zfind.engine import (
                        template_search, multi_template_search)
                    tname = self._template_combo.currentText()
                    if tname == 'Best of all':
                        result = multi_template_search(
                            self.spec,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=True,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                        )
                    else:
                        result = template_search(
                            self.spec, template_name=tname,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=True,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                        )
                elif method == 'PCA':
                    from rbcodes.GUIs.zfind.engine import (
                        pca_search, multi_pca_search)
                    pname = self._pca_combo.currentText()
                    if pname == 'Best of all':
                        result = multi_pca_search(
                            self.spec,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=True,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                        )
                    else:
                        result = pca_search(
                            self.spec, template_set=pname,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=True,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                        )
                else:   # Line Search
                    linelist_name = self._linelist_combo.currentText()
                    try:
                        search_df = _linelist_to_df(linelist_name)
                    except Exception as exc:
                        self._status.setText(
                            f'Cannot load linelist "{linelist_name}": {exc}')
                        self._run_btn.setEnabled(True)
                        return
                    from rbcodes.GUIs.zfind.engine import line_search
                    result = line_search(
                        self.spec, search_df,
                        z_min=z_min, z_max=z_max, n_steps=n_steps,
                        mode=self.mode, fwhm_ang=fwhm,
                        fit_continuum=True,
                        window_pixels=win_pix,
                        smooth_pixels=smooth_pix,
                        data_norm=data_norm,
                        wave_min=wave_min, wave_max=wave_max,
                    )
            except Exception as exc:
                self._status.setText(f'Search failed: {exc}')
                self._run_btn.setEnabled(True)
                return

            self._run_btn.setEnabled(True)
            self._result       = result
            self._df           = overlay_df   # always the overlay linelist
            self._selected_em_idx = 0
            self._smooth_scale    = 1

            if result.solutions:
                self._selected_z = result.solutions[0].z
            elif hasattr(result, 'candidates') and result.candidates:
                self._selected_z = result.candidates[0].z
            else:
                self._selected_z = None

            self._update_chi2_plot(highlight_z=self._selected_z)
            self._update_spectrum_plot(z=self._selected_z)
            self._update_table()
            self._accept_btn.setEnabled(True)

            warn_str = '; '.join(result.warnings) if result.warnings else ''
            best = result.best() if hasattr(result, 'best') else None
            if best:
                method_tag = best.method if best.method else ''
                msg = (f'Done.  Best z = {best.z:.5f} ± {best.z_err:.5f}  '
                       f'chi2={best.chi2_dof:.3g}  [{method_tag}]')
            elif hasattr(result, 'candidates') and result.candidates:
                msg = f'Done.  {len(result.candidates)} absorber candidates.'
            else:
                msg = 'Done. No solution found.'
            if warn_str:
                msg += f'   ⚠ {warn_str}'
            self._status.setText(msg)

        # ------------------------------------------------------------------
        # Chi2 plot
        # ------------------------------------------------------------------

        def _reset_plots(self):
            self._ax.clear()
            self._ax.set_xlabel('Redshift  z')
            self._ax.set_ylabel('chi2' if self.mode == 'emission'
                                else 'Significance')
            self._ax.set_title('Run search to see results')
            self._ax_spec.clear()
            self._ax_spec.set_xlabel('Wavelength  (Å)')
            self._ax_spec.set_ylabel('Flux')
            self._canvas.draw()

        def _update_chi2_plot(self, highlight_z=None):
            self._ax.clear()
            if self._result is None:
                self._canvas.draw()
                return

            z = self._result.z_array
            curves = self._result.chi2_curves
            method = self._method_combo.currentText()

            # Curve colours — single curve gets steelblue; multiple get palette
            _curve_colors = ['steelblue', '#e67e22', '#27ae60', '#8e44ad',
                             '#e74c3c', '#16a085', '#7f8c8d']

            if hasattr(self._result, 'candidates'):
                # Absorption result
                self._ax.plot(z, self._result.significance_curve,
                              color='steelblue', lw=0.8, alpha=0.85)
                self._ax.set_ylabel('Significance')
                items = self._result.candidates[:10]
            else:
                # Emission / template result — draw all chi2 curves
                for ci, curve in enumerate(curves):
                    col = _curve_colors[ci % len(_curve_colors)]
                    lw  = 0.9 if len(curves) > 1 else 0.8
                    alpha = 0.7 if len(curves) > 1 else 0.85
                    self._ax.plot(z, curve['chi2'], color=col, lw=lw,
                                  alpha=alpha, label=curve['label'])
                self._ax.set_ylabel('chi2 / dof')
                items = self._result.solutions

            for i, item in enumerate(items):
                self._ax.axvline(item.z,
                                 color=_SOL_COLORS[i % len(_SOL_COLORS)],
                                 lw=1.6 if i == 0 else 0.9,
                                 ls='--', alpha=0.8,
                                 label=f'z={item.z:.4f}')

            if highlight_z is not None:
                self._ax.axvline(highlight_z, color='black', lw=2.0,
                                 ls='-', alpha=0.35)

            self._ax.set_xlabel('Redshift  z')
            # Build title tag
            if method == 'Template':
                tag = self._template_combo.currentText()
            elif method == 'PCA':
                tag = f'PCA:{self._pca_combo.currentText()}'
            else:
                tag = self._linelist_combo.currentText()
            self._ax.set_title(
                f'{self.mode.capitalize()} — {tag}   '
                f'(click anywhere to set z)', fontsize=9)
            if items or len(curves) > 1:
                self._ax.legend(fontsize=7, loc='best',
                                ncol=min(5, max(len(items), len(curves))))

        # ------------------------------------------------------------------
        # Spectrum panel
        # ------------------------------------------------------------------

        def _draw_spectrum_only(self):
            self._smooth_scale = 1
            import astropy.units as u
            wave = self.spec.wavelength.to(u.AA).value
            flux = self.spec.flux.value
            self._raw_flux = flux.copy()
            self._raw_err  = (self.spec.sig.value.copy()
                              if self.spec.sig_is_set else None)
            self._spec_wave_full = wave.copy()
            self._update_spectrum_plot(z=None, flux_override=flux)

        def _update_spectrum_plot(self, z=None, flux_override=None, _no_draw=False):
            ax = self._ax_spec
            ax.clear()

            if self.spec is None:
                self._canvas.draw()
                return

            import astropy.units as u
            wave = self.spec.wavelength.to(u.AA).value

            # Use smoothed flux if available, raw otherwise
            if flux_override is not None:
                flux = flux_override
            elif self._raw_flux is not None:
                if self._smooth_scale > 1:
                    flux = convolve(self._raw_flux,
                                    Box1DKernel(self._smooth_scale))
                else:
                    flux = self._raw_flux
            else:
                flux = self.spec.flux.value
                self._raw_flux = flux.copy()

            ax.plot(wave, flux, color='#2c3e50', lw=0.7, alpha=0.9)

            if self._raw_err is not None:
                err = (convolve(self._raw_err, Box1DKernel(self._smooth_scale))
                       if self._smooth_scale > 1 else self._raw_err)
                ax.fill_between(wave, flux - err, flux + err,
                                color='#bdc3c7', alpha=0.35, lw=0)

            if self.spec.co_is_set:
                ax.plot(wave, self.spec.co.value,
                        color='#e74c3c', lw=0.8, ls='--', alpha=0.7)

            ax.set_xlabel('Wavelength  (Å)', fontsize=8)
            ax.set_ylabel('Flux', fontsize=8)
            ax.tick_params(labelsize=7)

            finite = flux[np.isfinite(flux)]
            if len(finite) > 10:
                ylo = np.percentile(finite, 1)
                yhi = np.percentile(finite, 99)
                pad = 0.15 * (yhi - ylo) if yhi > ylo else 0.1 * abs(yhi)
                ax.set_ylim(ylo - pad, yhi + pad)

            # --- line overlay ---
            if z is not None and self._df is not None:
                wmin, wmax = wave.min(), wave.max()
                obs   = self._df['wave'].values * (1.0 + z)
                names = self._df['name'].values
                mask  = (obs > wmin) & (obs < wmax)

                ylo_ax, yhi_ax = ax.get_ylim()
                label_y = yhi_ax - 0.04 * (yhi_ax - ylo_ax)

                for ow, nm in zip(obs[mask], names[mask]):
                    ax.axvline(ow, color='#e74c3c', lw=0.9,
                               ls='--', alpha=0.75)
                    short = nm.split()[0] if ' ' in nm else nm
                    ax.text(ow, label_y, short, rotation=90,
                            fontsize=6, color='#e74c3c',
                            va='top', ha='right', alpha=0.85)

                title = (f'Spectrum — z = {z:.5f}  '
                         f'({int(mask.sum())} lines in range)'
                         + (f'   smooth={self._smooth_scale}'
                            if self._smooth_scale > 1 else ''))
            else:
                title = ('Spectrum'
                         + (f'   smooth={self._smooth_scale}'
                            if self._smooth_scale > 1 else '')
                         + '   (run search then click chi2 plot to set z)')

            ax.set_title(title, fontsize=9)
            if not _no_draw:
                self._canvas.draw()

        # ------------------------------------------------------------------
        # Keystroke navigation (cursor must be over the canvas)
        # ------------------------------------------------------------------

        def _on_key_press(self, event):
            """
            Navigation keys — work on whichever panel the cursor is hovering over.
            S/U (smoothing) only apply to the spectrum panel.

            x/X  set left/right x limit    t/b  set top/bottom y limit
            r    reset axes                 o    zoom out 1.5×
            [/]  pan left/right             a    autoscale y to visible x
            S/U  increase/decrease display smoothing  (spectrum panel only)
            """
            key = event.key
            if key is None:
                return

            # Which axes is the cursor in?  Fall back to spectrum panel.
            ax = event.inaxes if event.inaxes in (self._ax, self._ax_spec) \
                 else self._ax_spec
            x, y = event.xdata, event.ydata

            # --- reset (r) — no cursor position needed ---
            if key == 'r':
                if ax is self._ax_spec:
                    if self._raw_flux is not None and hasattr(self, '_spec_wave_full'):
                        wave = self._spec_wave_full
                        ax.set_xlim(wave.min(), wave.max())
                        finite = self._raw_flux[np.isfinite(self._raw_flux)]
                        if len(finite) > 10:
                            ylo = np.percentile(finite, 1)
                            yhi = np.percentile(finite, 99)
                            pad = 0.15 * (yhi - ylo)
                            ax.set_ylim(ylo - pad, yhi + pad)
                else:  # chi2 axes
                    ax.autoscale()
                self._canvas.draw()
                return

            # S/U smoothing — spectrum panel only, no cursor coords needed
            if key == 'S':
                self._smooth_scale = max(1, self._smooth_scale + 2)
                if self._smooth_scale % 2 == 0:
                    self._smooth_scale += 1
                xlim = self._ax_spec.get_xlim()
                ylim = self._ax_spec.get_ylim()
                self._update_spectrum_plot(z=self._selected_z, _no_draw=True)
                self._ax_spec.set_xlim(xlim)
                self._ax_spec.set_ylim(ylim)
                self._canvas.draw()
                return
            if key == 'U':
                self._smooth_scale = max(1, self._smooth_scale - 2)
                if self._smooth_scale > 1 and self._smooth_scale % 2 == 0:
                    self._smooth_scale -= 1
                xlim = self._ax_spec.get_xlim()
                ylim = self._ax_spec.get_ylim()
                self._update_spectrum_plot(z=self._selected_z, _no_draw=True)
                self._ax_spec.set_xlim(xlim)
                self._ax_spec.set_ylim(ylim)
                self._canvas.draw()
                return

            # All remaining keys need cursor coordinates
            if x is None or y is None:
                return

            if key == 'x':
                ax.set_xlim(x, ax.get_xlim()[1])
            elif key == 'X':
                ax.set_xlim(ax.get_xlim()[0], x)
            elif key == 't':
                ax.set_ylim(ax.get_ylim()[0], y)
            elif key == 'b':
                ax.set_ylim(y, ax.get_ylim()[1])
            elif key == 'o':
                xlim = ax.get_xlim()
                cen  = 0.5 * (xlim[0] + xlim[1])
                hw   = 0.75 * (xlim[1] - xlim[0])
                ax.set_xlim(cen - hw, cen + hw)
            elif key == '[':
                xlim = ax.get_xlim()
                dx   = xlim[1] - xlim[0]
                ax.set_xlim(xlim[0] - dx, xlim[1] - dx)
            elif key == ']':
                xlim = ax.get_xlim()
                dx   = xlim[1] - xlim[0]
                ax.set_xlim(xlim[0] + dx, xlim[1] + dx)
            elif key == 'a':
                # autoscale y to visible x range
                xlim = ax.get_xlim()
                if ax is self._ax_spec:
                    if self._raw_flux is not None and hasattr(self, '_spec_wave_full'):
                        wave = self._spec_wave_full
                        flux = (convolve(self._raw_flux,
                                         Box1DKernel(self._smooth_scale))
                                if self._smooth_scale > 1 else self._raw_flux)
                        mask = (wave >= xlim[0]) & (wave <= xlim[1])
                        vis  = flux[mask][np.isfinite(flux[mask])]
                        if len(vis) > 0:
                            ylo, yhi = vis.min(), vis.max()
                            pad = 0.1 * (yhi - ylo) if yhi != ylo else 0.1 * abs(yhi)
                            ax.set_ylim(ylo - pad, yhi + pad)
                else:  # chi2 axes — autoscale y to visible z range
                    if self._result is not None:
                        z_arr = self._result.z_array
                        curve = (self._result.chi2_curves[0]['chi2']
                                 if self.mode == 'emission'
                                 else self._result.significance_curve)
                        mask = (z_arr >= xlim[0]) & (z_arr <= xlim[1])
                        vis  = curve[mask][np.isfinite(curve[mask])]
                        if len(vis) > 0:
                            ylo, yhi = vis.min(), vis.max()
                            pad = 0.1 * (yhi - ylo) if yhi != ylo else 0.1 * abs(yhi)
                            ax.set_ylim(ylo - pad, yhi + pad)
            else:
                return   # unhandled key — skip redraw

            self._canvas.draw()

        def _on_scroll(self, event):
            """Mouse-wheel zoom on x-axis of whichever panel the cursor is in."""
            if event.inaxes not in (self._ax, self._ax_spec):
                return
            ax   = event.inaxes
            xc   = event.xdata
            xlim = ax.get_xlim()
            hw   = 0.5 * (xlim[1] - xlim[0])
            factor = 0.85 if event.button == 'up' else 1.0 / 0.85
            ax.set_xlim(xc - hw * factor, xc + hw * factor)
            self._canvas.draw()

        # ------------------------------------------------------------------
        # Plot click
        # ------------------------------------------------------------------

        def _on_plot_click(self, event):
            if event.inaxes == self._ax and event.xdata is not None:
                # Any click on chi2/sig axes → immediately use that z
                z = float(event.xdata)
                self._selected_z = z
                self._update_chi2_plot(highlight_z=z)
                # Preserve whatever zoom/pan the user has on the spectrum panel
                xlim = self._ax_spec.get_xlim()
                ylim = self._ax_spec.get_ylim()
                self._update_spectrum_plot(z=z, _no_draw=True)
                self._ax_spec.set_xlim(xlim)
                self._ax_spec.set_ylim(ylim)
                self._canvas.draw()
                return

                # Try to select nearest table row (within Δz=0.05)
                if self.mode == 'emission' and self._result and self._result.solutions:
                    zvals = np.array([s.z for s in self._result.solutions])
                    idx   = int(np.argmin(np.abs(zvals - z)))
                    if abs(zvals[idx] - z) < 0.05:
                        self._em_table.blockSignals(True)
                        self._em_table.selectRow(idx)
                        self._em_table.blockSignals(False)
                elif self.mode == 'absorption' and self._result and self._result.candidates:
                    n     = self._abs_table.rowCount()
                    zvals = np.array([self._result.candidates[i].z
                                      for i in range(n)])
                    idx   = int(np.argmin(np.abs(zvals - z)))
                    if abs(zvals[idx] - z) < 0.05:
                        self._abs_table.blockSignals(True)
                        self._abs_table.selectRow(idx)
                        self._abs_table.blockSignals(False)

        # ------------------------------------------------------------------
        # Table population
        # ------------------------------------------------------------------

        def _update_table(self):
            if self._result is None:
                return
            if self.mode == 'emission':
                self._em_table.setRowCount(0)
                for sol in self._result.solutions:
                    row = self._em_table.rowCount()
                    self._em_table.insertRow(row)
                    self._em_table.setItem(row, 0,
                        QTableWidgetItem(f'{sol.z:.6f}'))
                    z_err_str = (f'{sol.z_err:.6f}'
                                 if np.isfinite(sol.z_err) else 'nan')
                    self._em_table.setItem(row, 1,
                        QTableWidgetItem(z_err_str))
                    self._em_table.setItem(row, 2,
                        QTableWidgetItem(f'{sol.chi2_dof:.4g}'))
                    self._em_table.setItem(row, 3,
                        QTableWidgetItem(str(sol.n_features)))
                if self._em_table.rowCount() > 0:
                    self._em_table.blockSignals(True)
                    self._em_table.selectRow(0)
                    self._em_table.blockSignals(False)
            else:
                self._abs_table.setRowCount(0)
                for i, cand in enumerate(self._result.candidates[:20]):
                    row = self._abs_table.rowCount()
                    self._abs_table.insertRow(row)
                    chk = QTableWidgetItem()
                    chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
                    chk.setCheckState(Qt.Checked if cand.significance > 3.0
                                      else Qt.Unchecked)
                    self._abs_table.setItem(row, 0, chk)
                    self._abs_table.setItem(row, 1,
                        QTableWidgetItem(f'{cand.z:.6f}'))
                    self._abs_table.setItem(row, 2,
                        QTableWidgetItem(f'{cand.significance:.2f}'))
                    self._abs_table.setItem(row, 3,
                        QTableWidgetItem(str(cand.n_lines)))
                    lines_str = ', '.join(cand.lines_matched[:4])
                    if len(cand.lines_matched) > 4:
                        lines_str += '…'
                    self._abs_table.setItem(row, 4,
                        QTableWidgetItem(lines_str))

        # ------------------------------------------------------------------
        # Table selection callbacks
        # ------------------------------------------------------------------

        def _on_em_selection_changed(self):
            idx = self._em_table.currentRow()
            if idx < 0 or self._result is None:
                return
            self._selected_em_idx = idx
            if idx < len(self._result.solutions):
                z = self._result.solutions[idx].z
                self._selected_z = z
                self._update_chi2_plot(highlight_z=z)
                xlim = self._ax_spec.get_xlim()
                ylim = self._ax_spec.get_ylim()
                self._update_spectrum_plot(z=z, _no_draw=True)
                self._ax_spec.set_xlim(xlim)
                self._ax_spec.set_ylim(ylim)
                self._canvas.draw()

        def _on_abs_selection_changed(self):
            idx = self._abs_table.currentRow()
            if idx < 0 or self._result is None:
                return
            if idx < len(self._result.candidates):
                z = self._result.candidates[idx].z
                self._selected_z = z
                self._update_chi2_plot(highlight_z=z)
                xlim = self._ax_spec.get_xlim()
                ylim = self._ax_spec.get_ylim()
                self._update_spectrum_plot(z=z, _no_draw=True)
                self._ax_spec.set_xlim(xlim)
                self._ax_spec.set_ylim(ylim)
                self._canvas.draw()

        # ------------------------------------------------------------------
        # Accept / Cancel
        # ------------------------------------------------------------------

        def _on_accept(self):
            if self._result is None:
                return
            if self.mode == 'emission':
                from rbcodes.GUIs.zfind.adapters import zfind_to_zgui_z
                payload = zfind_to_zgui_z(self._result,
                                          idx=self._selected_em_idx)
                self.accepted_z.emit(payload)
                # Also print for standalone / CLI use
                print(f'rb_zfind accepted z = {payload[0]:.6f}  '
                      f'z_err = {payload[1]:.6f}')
            else:
                from rbcodes.GUIs.zfind.adapters import absorbers_to_multispec
                accepted = [
                    row for row in range(self._abs_table.rowCount())
                    if (self._abs_table.item(row, 0) is not None and
                        self._abs_table.item(row, 0).checkState()
                        == Qt.Checked)
                ]
                payload = absorbers_to_multispec(self._result, accepted)
                self.accepted_absorbers.emit(payload)
                print('rb_zfind accepted absorbers:')
                for ab in payload:
                    print(f"  {ab['label']}")
            self.accept()

except ImportError:
    pass
