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
        QGroupBox,
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
        'zfind_em', 'zfind_galaxy', 'zfind_stellar', 'zfind_igm', 'zfind_qso',
        # -- full rb_setline lists --
        'Gal_Em', 'Gal_Abs', 'Gal_long', 'Gal',
        'LLS', 'LLS Small', 'DLA', 'LBG',
        'AGN', 'Eiger_Strong',
        'HI', 'HI_recomb', 'HI_recomb_light',
        'EUV', 'LLS_EUV', 'atom',
    ]
    _DEFAULT_LINELIST = 'zfind_em'
    # Solution marker colours — visible on both light and dark backgrounds
    # (avoid near-black or near-white entries)
    _SOL_COLORS = ['#e74c3c', '#e67e22', '#27ae60', '#8e44ad', '#2980b9',
                   '#c0392b', '#d35400', '#16a085', '#e91e63', '#00bcd4']


    def _linelist_to_df(linelist_name):
        """Return a DataFrame for *linelist_name*.

        Curated presets ('zfind_em', 'zfind_galaxy', etc.) are served from the
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
        Self-contained redshift finder popup.

        Parameters
        ----------
        spec            : rb_spectrum or None
        default_linelist: str or None
        parent          : QWidget or None
        """

        accepted_z         = pyqtSignal(list)
        accepted_absorbers = pyqtSignal(list)

        def __init__(self, spec=None, default_linelist=None, parent=None,
                     dark_theme=False,
                     # legacy kwargs silently ignored
                     mode=None, z_qso=None):
            super().__init__(parent)
            self.spec             = spec
            self.default_linelist = default_linelist
            self.dark_theme       = dark_theme

            self._result          = None
            self._df              = None       # linelist after last Run
            self._selected_em_idx = 0
            self._selected_z      = None       # z currently shown in spectrum

            # smoothing state
            self._smooth_scale    = 1
            self._raw_flux        = None       # unsmoothed flux array
            self._raw_err         = None       # unsmoothed error array

            # last-run metadata for spectrum overlay
            self._last_run_method        = None   # 'Picket Fence' | 'Template' | 'PCA'
            self._last_run_template_name = None   # template name for overlay
            self._last_run_pca_set       = None   # PCA set name for overlay
            self._last_run_smooth_pix    = 1      # engine smooth used in last run
            self._last_run_data_norm     = 'raw'  # data_norm used in last run
            self._last_run_model_norm    = 'raw'  # model_norm used in last run

            self.setWindowTitle('Redshift Finder')
            self.resize(960, 860)
            self.setMinimumSize(780, 660)

            self._build_ui()
            if dark_theme:
                self._apply_dark_theme()
            self._setup_connections()
            self._set_initial_state()

        # ------------------------------------------------------------------
        # UI
        # ------------------------------------------------------------------

        def _build_ui(self):
            outer = QVBoxLayout(self)
            outer.setSpacing(5)

            # ---- Parameters ----
            param_box  = QGroupBox('Parameters')
            param_vbox = QVBoxLayout(param_box)
            param_vbox.setSpacing(4)

            # Common row 1: Method | z min | z max | n steps | Run btn
            row1 = QHBoxLayout()
            row1.addWidget(QLabel('Method:'))
            self._method_combo = QComboBox()
            self._method_combo.addItems(['PCA', 'Template', 'Picket Fence'])
            self._method_combo.setMinimumWidth(100)
            row1.addWidget(self._method_combo)
            row1.addSpacing(12)
            row1.addWidget(QLabel('z min:'))
            self._z_min = QLineEdit('0.0')
            self._z_min.setFixedWidth(65)
            self._z_min.setValidator(QDoubleValidator(0.0, 30.0, 6))
            row1.addWidget(self._z_min)
            row1.addWidget(QLabel('z max:'))
            self._z_max = QLineEdit('6.0')
            self._z_max.setFixedWidth(65)
            self._z_max.setValidator(QDoubleValidator(0.0, 30.0, 6))
            row1.addWidget(self._z_max)
            self._z_reset_btn = QPushButton('↺')
            self._z_reset_btn.setFixedWidth(26)
            self._z_reset_btn.setToolTip(
                'Reset z min/max to the default range for the selected template or PCA set.')
            row1.addWidget(self._z_reset_btn)
            row1.addWidget(QLabel('n steps:'))
            self._n_steps = QLineEdit('5000')
            self._n_steps.setFixedWidth(65)
            self._n_steps.setValidator(QIntValidator(100, 200000))
            row1.addWidget(self._n_steps)
            row1.addStretch()
            self._run_btn = QPushButton('Run Search')
            self._run_btn.setMinimumWidth(100)
            row1.addWidget(self._run_btn)
            param_vbox.addLayout(row1)

            # Common row 2: Overlay lines | λ min | λ max
            row2 = QHBoxLayout()
            row2.addWidget(QLabel('Overlay lines:'))
            self._overlay_combo = QComboBox()
            self._overlay_combo.setMinimumWidth(120)
            self._overlay_combo.addItems(_LINELIST_NAMES)
            self._overlay_combo.setToolTip(
                'Linelist drawn on the spectrum panel at the selected z.\n'
                'In Template/PCA mode this is the only linelist shown.')
            row2.addWidget(self._overlay_combo)
            row2.addSpacing(12)
            row2.addWidget(QLabel('λ min (Å):'))
            self._wave_min = QLineEdit('')
            self._wave_min.setFixedWidth(65)
            self._wave_min.setValidator(QDoubleValidator(0.0, 1e6, 2))
            self._wave_min.setPlaceholderText('auto')
            self._wave_min.setToolTip(
                'Observed-frame wavelength lower limit for chi2 scan.\n'
                'Pixels below this are masked (ivar→0). Leave blank for full range.')
            row2.addWidget(self._wave_min)
            row2.addWidget(QLabel('λ max (Å):'))
            self._wave_max = QLineEdit('')
            self._wave_max.setFixedWidth(65)
            self._wave_max.setValidator(QDoubleValidator(0.0, 1e6, 2))
            self._wave_max.setPlaceholderText('auto')
            self._wave_max.setToolTip(
                'Observed-frame wavelength upper limit for chi2 scan.\n'
                'Pixels above this are masked (ivar→0). Leave blank for full range.')
            row2.addWidget(self._wave_max)
            row2.addStretch()
            param_vbox.addLayout(row2)

            # Common row 3: Resolution (shared across all methods) | Smooth (pix)
            row3 = QHBoxLayout()
            row3.addWidget(QLabel('Resolution:'))
            self._res_type = QComboBox()
            self._res_type.addItems(['R', 'FWHM (Å)', 'FWHM (km/s)', 'FWHM (pix)'])
            self._res_type.setToolTip(
                'Instrument resolution — shared by all three methods.\n'
                'R:           resolving power λ/Δλ  (e.g. R=1000 for LRIS, R=4000 for DEIMOS)\n'
                'FWHM (Å):    LSF FWHM in Angstroms  (e.g. 3.5 Å for a typical spectrograph)\n'
                'FWHM (km/s): velocity resolution  (e.g. 150 km/s)\n'
                'FWHM (pix):  LSF FWHM in detector pixels\n\n'
                'PCA / Template: convolves model to match instrument LSF before chi² scan.\n'
                'Picket Fence:   sets matched-filter window half-width to 1.5× this element.\n'
                'Leave value at 0 to skip convolution / use Win (pix) fallback.')
            row3.addWidget(self._res_type)
            self._res_val = QLineEdit('0.0')
            self._res_val.setFixedWidth(65)
            self._res_val.setValidator(QDoubleValidator(0.0, 1e6, 3))
            self._res_val.setToolTip('Resolution value in the units shown. 0 = no convolution.')
            row3.addWidget(self._res_val)
            row3.addSpacing(16)
            row3.addWidget(QLabel('Smooth (pix):'))
            self._smooth_pix = QLineEdit('3')
            self._smooth_pix.setFixedWidth(45)
            self._smooth_pix.setValidator(QIntValidator(1, 99))
            self._smooth_pix.setToolTip(
                'Boxcar pre-smoothing applied to the spectrum before the chi² scan.\n'
                'Width in pixels — shared by all three methods.\n'
                'Rule of thumb: ~1 resolution element.  Set to 1 to disable.')
            row3.addWidget(self._smooth_pix)
            row3.addStretch()
            param_vbox.addLayout(row3)

            # ---- Method-specific stacked panels ----
            self._method_stack = QStackedWidget()

            # -- Page 0: Picket Fence --
            pf_box  = QGroupBox('Picket Fence options')
            pf_grid = QGridLayout(pf_box)
            pf_grid.setSpacing(4)

            pf_grid.addWidget(QLabel('Linelist:'), 0, 0)
            self._pf_linelist_combo = QComboBox()
            self._pf_linelist_combo.setMinimumWidth(120)
            self._pf_linelist_combo.addItems(_LINELIST_NAMES)
            pf_grid.addWidget(self._pf_linelist_combo, 0, 1)

            pf_grid.addWidget(QLabel('Win (pix):'), 0, 2)
            self._pf_win_pix = QLineEdit('5')
            self._pf_win_pix.setFixedWidth(45)
            self._pf_win_pix.setValidator(QIntValidator(1, 200))
            self._pf_win_pix.setToolTip(
                'Matched-filter half-window in pixels — used only when Resolution = 0.\n'
                'The window is centred on each predicted line position.\n'
                'Rule of thumb: ~1–2× the resolution element in pixels.')
            pf_grid.addWidget(self._pf_win_pix, 0, 3)

            pf_grid.addWidget(QLabel('Data norm:'), 0, 4)
            self._pf_data_norm = QComboBox()
            self._pf_data_norm.addItems(['subtract', 'normalize', 'raw'])
            self._pf_data_norm.setToolTip(
                'How to preprocess the observed spectrum:\n'
                '  subtract   flux − continuum  (recommended)\n'
                '  normalize  flux / continuum − 1\n'
                '  raw        flux as-is')
            pf_grid.addWidget(self._pf_data_norm, 0, 5)

            pf_grid.addWidget(QLabel('Mode:'), 0, 6)
            self._pf_mode = QComboBox()
            self._pf_mode.addItems(['Direct', 'Detect+Match'])
            self._pf_mode.setToolTip(
                'Direct: score linelist windows at every trial z (Mode A).\n'
                'Detect+Match: first find spectral peaks, then match to linelist (Mode B).')
            pf_grid.addWidget(self._pf_mode, 0, 7)

            self._pf_prom_label = QLabel('Prominence σ:')
            self._pf_prom_label.setVisible(False)
            pf_grid.addWidget(self._pf_prom_label, 0, 8)
            self._pf_prom_val = QLineEdit('3.0')
            self._pf_prom_val.setFixedWidth(55)
            self._pf_prom_val.setValidator(QDoubleValidator(0.5, 100.0, 2))
            self._pf_prom_val.setToolTip(
                'Peak detection threshold in units of local noise σ (Mode B only).')
            self._pf_prom_val.setVisible(False)
            pf_grid.addWidget(self._pf_prom_val, 0, 9)

            self._method_stack.addWidget(pf_box)

            # -- Page 1: Template --
            from rbcodes.GUIs.zfind.engine import TEMPLATE_NAMES, PCA_NAMES
            tmpl_box  = QGroupBox('Template options')
            tmpl_grid = QGridLayout(tmpl_box)
            tmpl_grid.setSpacing(4)

            tmpl_grid.addWidget(QLabel('Template:'), 0, 0)
            self._tmpl_combo = QComboBox()
            self._tmpl_combo.addItems(['Best of all'] + TEMPLATE_NAMES)
            self._tmpl_combo.setToolTip(
                '"Best of all" tries every MARZ template and picks the lowest chi2/dof')
            _tmpl_tips = {
                'Best of all':      'Tries all 5 templates, picks lowest chi²/dof.',
                'EarlyType':        'Passive/elliptical galaxy. D4000 break, CaII H&K, '
                                    'MgI b absorption. Best for z = 0–1.5.',
                'Intermediate':     'Intermediate-type galaxy. Mix of emission and '
                                    'absorption. Best for z = 0–1.5.',
                'LateTypeEmission': 'Star-forming spiral. Strong [OII], Hβ, [OIII], Hα '
                                    'emission. Use for z = 0–7 emission-line galaxies.',
                'Composite':        'Starburst + AGN composite. Best for z = 0–1.5.',
                'QSO':              'Broad-line QSO / Type-1 AGN. Best for z = 0–5.5.',
            }
            for _i, _name in enumerate(['Best of all'] + TEMPLATE_NAMES):
                self._tmpl_combo.setItemData(
                    _i, _tmpl_tips.get(_name, ''), Qt.ToolTipRole)
            tmpl_grid.addWidget(self._tmpl_combo, 0, 1)

            tmpl_grid.addWidget(QLabel('Data norm:'), 0, 2)
            self._tmpl_data_norm = QComboBox()
            self._tmpl_data_norm.addItems(['raw', 'normalize', 'subtract'])
            self._tmpl_data_norm.setToolTip(
                'How to preprocess the observed spectrum:\n'
                '  raw        flux as-is (recommended — preserves D4000,\n'
                '             Balmer break and other continuum features)\n'
                '  normalize  flux / continuum (removes continuum shape)\n'
                '  subtract   flux − continuum')
            tmpl_grid.addWidget(self._tmpl_data_norm, 0, 3)

            tmpl_grid.addWidget(QLabel('Model norm:'), 0, 4)
            self._tmpl_model_norm = QComboBox()
            self._tmpl_model_norm.addItems(['raw', 'normalize', 'subtract'])
            self._tmpl_model_norm.setToolTip(
                'How to preprocess the template:\n'
                '  raw        template as-is (recommended)\n'
                '  normalize  template / pseudo-continuum\n'
                '  subtract   template − pseudo-continuum')
            tmpl_grid.addWidget(self._tmpl_model_norm, 0, 5)

            self._method_stack.addWidget(tmpl_box)

            # -- Page 2: PCA --
            pca_box  = QGroupBox('PCA options')
            pca_grid = QGridLayout(pca_box)
            pca_grid.setSpacing(4)

            pca_grid.addWidget(QLabel('PCA set:'), 0, 0)
            self._pca_combo = QComboBox()
            self._pca_combo.addItems(['Best of all'] + PCA_NAMES)
            self._pca_combo.setToolTip(
                'DESI redrock PCA eigenvectors.\n'
                '"Best of all" runs galaxy + qso and picks lowest chi2/dof.\n'
                'Download first: python templates/download_desi_pca.py')
            _pca_tips = {
                'Best of all': 'Tries galaxy + qso_loz + qso_hiz, picks lowest chi²/dof.',
                'galaxy':      'DESI galaxy eigenvectors. Star-forming and passive galaxies.\n'
                               'Rest 1228–11000 Å, z = 0.0–1.6.\n'
                               'Best for low-z galaxies with detectable continuum.',
                'qso_loz':    'DESI low-z QSO eigenvectors. Broad-line AGN, Seyferts,\n'
                               'low-z quasars. Rest 1329–9634 Å, z = 0.0–2.5.',
                'qso_hiz':    'DESI high-z QSO eigenvectors. Lyα, CIV, CIII], [OIII].\n'
                               'Rest 444–4499 Å, z = 1.0–6.0.\n'
                               'Best choice for z > 2 objects and high-z emission-line galaxies.',
            }
            for _i, _name in enumerate(['Best of all'] + PCA_NAMES):
                self._pca_combo.setItemData(
                    _i, _pca_tips.get(_name, ''), Qt.ToolTipRole)
            pca_grid.addWidget(self._pca_combo, 0, 1)

            pca_grid.addWidget(QLabel('Data norm:'), 0, 2)
            self._pca_data_norm = QComboBox()
            self._pca_data_norm.addItems(['raw', 'normalize', 'subtract'])
            self._pca_data_norm.setToolTip(
                'How to preprocess the observed spectrum:\n'
                '  raw        flux as-is (recommended — preserves D4000,\n'
                '             Balmer break; eigenvector coefficients absorb\n'
                '             flux calibration and throughput differences)\n'
                '  normalize  flux / continuum (removes continuum features)\n'
                '  subtract   flux − continuum')
            pca_grid.addWidget(self._pca_data_norm, 0, 3)

            pca_grid.addWidget(QLabel('Model norm:'), 0, 4)
            self._pca_model_norm = QComboBox()
            self._pca_model_norm.addItems(['raw', 'normalize', 'subtract'])
            self._pca_model_norm.setToolTip(
                'How to preprocess the PCA eigenvectors:\n'
                '  raw        eigenvectors as-is (recommended)\n'
                '  normalize  L2-normalise each eigenvector\n'
                '  subtract   mean-center each eigenvector')
            pca_grid.addWidget(self._pca_model_norm, 0, 5)

            self._method_stack.addWidget(pca_box)

            param_vbox.addWidget(self._method_stack)
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

            # ---- Results ----
            outer.addWidget(QLabel('Results:'))
            self._stack = QStackedWidget()

            self._em_table = QTableWidget(0, 4)
            self._em_table.setHorizontalHeaderLabels(
                ['z', 'z_err', 'score/chi2', 'n_lines'])
            self._em_table.setSelectionBehavior(QTableWidget.SelectRows)
            self._em_table.setSelectionMode(QTableWidget.SingleSelection)
            self._em_table.setEditTriggers(QTableWidget.NoEditTriggers)
            self._em_table.horizontalHeader().setSectionResizeMode(
                QHeaderView.Stretch)
            self._em_table.setMaximumHeight(120)
            self._stack.addWidget(self._em_table)

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
            self._help_btn = QPushButton('?')
            self._help_btn.setFixedWidth(28)
            self._help_btn.setToolTip('Help (F1)')
            self._accept_btn = QPushButton('Accept')
            self._accept_btn.setEnabled(False)
            self._cancel_btn = QPushButton('Cancel')
            btn_row.addWidget(self._help_btn)
            btn_row.addStretch()
            btn_row.addWidget(self._accept_btn)
            btn_row.addWidget(self._cancel_btn)
            outer.addLayout(btn_row)

        # ------------------------------------------------------------------
        # Connections
        # ------------------------------------------------------------------

        def _apply_dark_theme(self):
            """Apply dark palette + matplotlib dark style — matching rb_multispec."""
            from PyQt5.QtGui import QPalette, QColor
            from PyQt5.QtCore import Qt

            # --- Qt palette (same colours as rb_multispec) ---
            pal = QPalette()
            pal.setColor(QPalette.Window,          QColor(53, 53, 53))
            pal.setColor(QPalette.WindowText,      Qt.white)
            pal.setColor(QPalette.Base,            QColor(25, 25, 25))
            pal.setColor(QPalette.AlternateBase,   QColor(53, 53, 53))
            pal.setColor(QPalette.Text,            Qt.white)
            pal.setColor(QPalette.Button,          QColor(53, 53, 53))
            pal.setColor(QPalette.ButtonText,      Qt.white)
            pal.setColor(QPalette.Highlight,       QColor(42, 130, 218))
            pal.setColor(QPalette.HighlightedText, Qt.black)
            pal.setColor(QPalette.Link,            QColor(42, 130, 218))
            pal.setColor(QPalette.BrightText,      Qt.red)
            self.setPalette(pal)

            self.setStyleSheet("""
                QGroupBox { color: #cccccc; border: 1px solid #555; border-radius: 4px;
                             margin-top: 6px; padding-top: 6px; }
                QGroupBox::title { subcontrol-origin: margin; left: 8px; }
                QLabel      { color: #dddddd; }
                QLineEdit   { background: #1a1a1a; color: #ffffff;
                              border: 1px solid #555; border-radius: 3px; padding: 2px; }
                QComboBox   { background: #2a2a2a; color: #ffffff;
                              border: 1px solid #555; border-radius: 3px; padding: 2px; }
                QComboBox QAbstractItemView { background: #2a2a2a; color: #ffffff;
                              selection-background-color: #2a82da; }
                QPushButton { background: #474747; color: #ffffff;
                              border: 1px solid #666; border-radius: 4px; padding: 4px 10px; }
                QPushButton:hover   { background: #555555; }
                QPushButton:pressed { background: #333333; }
                QPushButton:disabled { color: #777; }
                QTableWidget { background: #1a1a1a; color: #ffffff;
                               gridline-color: #444; }
                QHeaderView::section { background: #353535; color: #cccccc;
                                       border: 1px solid #555; padding: 2px; }
            """)

            # --- Matplotlib figure: dark background ---
            _bg   = '#353535'
            _axes = '#252525'
            _fg   = 'white'
            _grid = '#444444'

            self._fig.patch.set_facecolor(_bg)
            for ax in (self._ax, self._ax_spec):
                ax.set_facecolor(_axes)
                ax.tick_params(colors=_fg, labelcolor=_fg)
                ax.xaxis.label.set_color(_fg)
                ax.yaxis.label.set_color(_fg)
                ax.title.set_color(_fg)
                for spine in ax.spines.values():
                    spine.set_edgecolor('#666666')

            self._canvas.draw()

        def _setup_connections(self):
            self._run_btn.clicked.connect(self._on_run)
            self._z_reset_btn.clicked.connect(self._on_z_range_reset)
            self._accept_btn.clicked.connect(self._on_accept)
            self._cancel_btn.clicked.connect(self.reject)
            self._help_btn.clicked.connect(self._on_help)
            self._method_combo.currentIndexChanged.connect(
                self._on_method_changed)
            self._pf_linelist_combo.currentIndexChanged.connect(
                self._on_linelist_changed)
            self._pf_mode.currentIndexChanged.connect(
                self._on_pf_mode_changed)
            self._em_table.itemSelectionChanged.connect(
                self._on_em_selection_changed)
            self._canvas.mpl_connect('button_press_event',  self._on_plot_click)
            self._canvas.mpl_connect('key_press_event',     self._on_key_press)
            self._canvas.mpl_connect('scroll_event',        self._on_scroll)
            # F1 as application-level shortcut (works even without canvas focus)
            from PyQt5.QtWidgets import QShortcut
            from PyQt5.QtGui import QKeySequence
            QShortcut(QKeySequence('F1'), self).activated.connect(self._on_help)

        # ------------------------------------------------------------------
        # Init
        # ------------------------------------------------------------------

        def _set_initial_state(self):
            target = (self.default_linelist
                      if self.default_linelist in _LINELIST_NAMES
                      else _DEFAULT_LINELIST)
            idx = self._pf_linelist_combo.findText(target)
            if idx >= 0:
                self._pf_linelist_combo.setCurrentIndex(idx)
            # Overlay defaults to zfind_galaxy — shows both emission and
            # absorption lines, useful for any galaxy type at a glance.
            # Only follow the caller's linelist choice if one was explicitly given.
            overlay_default = target if self.default_linelist else 'zfind_galaxy'
            oidx = self._overlay_combo.findText(overlay_default)
            if oidx >= 0:
                self._overlay_combo.setCurrentIndex(oidx)

            self._stack.setCurrentIndex(0)        # always show the z-solutions table
            self._method_stack.setCurrentIndex(2) # default to PCA panel (stack page 2)
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
        # Method switch  (Picket Fence ↔ Template ↔ PCA)
        # ------------------------------------------------------------------

        def _on_method_changed(self, index):
            # combo order: PCA=0, Template=1, Picket Fence=2
            # stack order: PF=0, Template=1, PCA=2  (built in that order)
            _STACK_PAGE = [2, 1, 0]
            self._method_stack.setCurrentIndex(_STACK_PAGE[index])
            # z range is NOT auto-filled on method switch — the user's values
            # are preserved.  Use the ↺ button next to z max to reset to the
            # default range for the currently selected template / PCA set.

        # ------------------------------------------------------------------
        # Linelist sync — overlay follows search linelist
        # ------------------------------------------------------------------

        def _on_linelist_changed(self, index):
            text = self._pf_linelist_combo.currentText()
            oidx = self._overlay_combo.findText(text)
            if oidx >= 0:
                self._overlay_combo.setCurrentIndex(oidx)

        def _on_pf_mode_changed(self, index):
            is_detect = (self._pf_mode.currentText() == 'Detect+Match')
            self._pf_prom_label.setVisible(is_detect)
            self._pf_prom_val.setVisible(is_detect)

        def _on_help(self):
            from rbcodes.GUIs.zfind.help import show_help_dialog
            show_help_dialog(parent=self)

        def _on_z_range_reset(self):
            """Reset z min/max to the default range for the active method."""
            method = self._method_combo.currentText()
            if method == 'Template':
                from rbcodes.GUIs.zfind.engine import _TEMPLATE_Z_RANGE
                tname = self._tmpl_combo.currentText()
                if tname != 'Best of all' and tname in _TEMPLATE_Z_RANGE:
                    zlo, zhi = _TEMPLATE_Z_RANGE[tname]
                    self._z_min.setText(str(zlo))
                    self._z_max.setText(str(zhi))
                else:
                    self._z_min.setText('0.0')
                    self._z_max.setText('6.0')
            elif method == 'PCA':
                from rbcodes.GUIs.zfind.engine import _PCA_Z_RANGE
                pname = self._pca_combo.currentText()
                if pname != 'Best of all' and pname in _PCA_Z_RANGE:
                    zlo, zhi = _PCA_Z_RANGE[pname]
                    self._z_min.setText(str(zlo))
                    self._z_max.setText(str(zhi))
                else:
                    self._z_min.setText('0.0')
                    self._z_max.setText('6.0')
            else:  # Picket Fence — no preset, reset to wide default
                self._z_min.setText('0.0')
                self._z_max.setText('6.0')

        # ------------------------------------------------------------------
        # Run
        # ------------------------------------------------------------------

        def _on_run(self):
            if self.spec is None:
                self._status.setText('No spectrum loaded.')
                return
            try:
                z_min   = float(self._z_min.text())
                z_max   = float(self._z_max.text())
                n_steps = int(self._n_steps.text())
            except ValueError:
                self._status.setText('Invalid parameter values.')
                return
            if z_min >= z_max:
                self._status.setText('z_min must be < z_max.')
                return

            method     = self._method_combo.currentText()
            wave_min_s = self._wave_min.text().strip()
            wave_max_s = self._wave_max.text().strip()
            wave_min   = float(wave_min_s) if wave_min_s else None
            wave_max   = float(wave_max_s) if wave_max_s else None

            # Always load the overlay linelist so line markers are drawn
            # regardless of method.
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

            # Read shared instrument parameters (common to all methods)
            try:
                smooth_pix = int(self._smooth_pix.text())
                res_val    = float(self._res_val.text())
            except ValueError:
                self._status.setText('Invalid Resolution or Smooth value.')
                self._run_btn.setEnabled(True)
                return
            res_type = self._res_type.currentText()
            res_kwargs = {}
            if res_val > 0.0:
                if res_type == 'R':
                    res_kwargs['R'] = res_val
                elif res_type == 'FWHM (Å)':
                    res_kwargs['fwhm_ang'] = res_val
                elif res_type == 'FWHM (km/s)':
                    res_kwargs['fwhm_kms'] = res_val
                elif res_type == 'FWHM (pix)':
                    res_kwargs['fwhm_pix'] = res_val

            try:
                if method == 'Template':
                    data_norm  = self._tmpl_data_norm.currentText()
                    model_norm = self._tmpl_model_norm.currentText()
                    from rbcodes.GUIs.zfind.engine import (
                        template_search, multi_template_search)
                    tname = self._tmpl_combo.currentText()
                    need_cont = (data_norm != 'raw')
                    if tname == 'Best of all':
                        result = multi_template_search(
                            self.spec,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=need_cont,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                            template_res_kwargs=res_kwargs,
                        )
                    else:
                        result = template_search(
                            self.spec, template_name=tname,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=need_cont,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                            template_res_kwargs=res_kwargs,
                        )
                elif method == 'PCA':
                    data_norm  = self._pca_data_norm.currentText()
                    model_norm = self._pca_model_norm.currentText()
                    from rbcodes.GUIs.zfind.engine import (
                        pca_search, multi_pca_search)
                    pname = self._pca_combo.currentText()
                    need_cont = (data_norm != 'raw')
                    if pname == 'Best of all':
                        result = multi_pca_search(
                            self.spec,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=need_cont,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                            pca_res_kwargs=res_kwargs,
                        )
                    else:
                        result = pca_search(
                            self.spec, template_set=pname,
                            z_min=z_min, z_max=z_max, n_steps=n_steps,
                            fit_continuum=need_cont,
                            data_norm=data_norm, model_norm=model_norm,
                            smooth_pixels=smooth_pix,
                            wave_min=wave_min, wave_max=wave_max,
                            pca_res_kwargs=res_kwargs,
                        )
                else:   # Picket Fence
                    try:
                        win_pix    = int(self._pf_win_pix.text())
                        prom_sigma = float(self._pf_prom_val.text())
                    except ValueError:
                        self._status.setText('Invalid Picket Fence parameter values.')
                        self._run_btn.setEnabled(True)
                        return
                    data_norm = self._pf_data_norm.currentText()
                    pf_mode   = self._pf_mode.currentText()

                    linelist_name = self._pf_linelist_combo.currentText()
                    try:
                        search_df = _linelist_to_df(linelist_name)
                    except Exception as exc:
                        self._status.setText(
                            f'Cannot load linelist "{linelist_name}": {exc}')
                        self._run_btn.setEnabled(True)
                        return
                    from rbcodes.GUIs.zfind.engine import picket_fence_search
                    result = picket_fence_search(
                        self.spec, search_df,
                        z_min=z_min, z_max=z_max, n_steps=n_steps,
                        smooth_fwhm_pix=(float(smooth_pix) if smooth_pix > 1 else None),
                        window_pixels=win_pix,
                        fit_continuum=True,
                        data_norm=data_norm,
                        wave_min=wave_min, wave_max=wave_max,
                        res_kwargs=res_kwargs,
                        pf_mode='detect_match' if pf_mode == 'Detect+Match' else 'direct',
                        prominence_sigma=prom_sigma,
                    )
            except Exception as exc:
                self._status.setText(f'Search failed: {exc}')
                self._run_btn.setEnabled(True)
                return

            self._run_btn.setEnabled(True)
            self._result          = result
            self._df              = overlay_df   # always the overlay linelist
            self._selected_em_idx = 0

            # store metadata for spectrum overlay and display smooth
            self._last_run_method     = method
            self._last_run_smooth_pix = smooth_pix   # already read from shared widget
            # store norms for correct overlay rendering
            if method == 'Template':
                self._last_run_data_norm  = self._tmpl_data_norm.currentText()
                self._last_run_model_norm = self._tmpl_model_norm.currentText()
            elif method == 'PCA':
                self._last_run_data_norm  = self._pca_data_norm.currentText()
                self._last_run_model_norm = self._pca_model_norm.currentText()
            else:
                self._last_run_data_norm  = 'raw'
                self._last_run_model_norm = 'raw'
            if method == 'Template':
                tname_stored = self._tmpl_combo.currentText()
                # For 'Best of all', extract winning template from result
                if tname_stored == 'Best of all' and result.solutions:
                    best_method = result.solutions[0].method  # e.g. 'Template:EarlyType'
                    if ':' in best_method:
                        tname_stored = best_method.split(':', 1)[1]
                self._last_run_template_name = tname_stored
            elif method == 'PCA':
                pset_stored = self._pca_combo.currentText()
                if pset_stored == 'Best of all' and result.solutions:
                    best_method = result.solutions[0].method  # e.g. 'PCA:galaxy'
                    if ':' in best_method:
                        pset_stored = best_method.split(':', 1)[1]
                self._last_run_pca_set       = pset_stored
                self._last_run_template_name = None
            else:
                self._last_run_template_name = None

            # Initialise display smooth to match the engine smooth so the
            # plot shows what the engine actually used
            self._smooth_scale = self._last_run_smooth_pix

            if result.solutions:
                self._selected_z = result.solutions[0].z
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
            else:
                msg = 'Done. No solution found.'
            if warn_str:
                msg += f'   ⚠ {warn_str}'
            self._status.setText(msg)

        # ------------------------------------------------------------------
        # Chi2 plot
        # ------------------------------------------------------------------

        def _dark_ax(self, ax):
            """Re-apply dark colours to ax after ax.clear() (only in dark_theme mode)."""
            if not self.dark_theme:
                return
            _fg = 'white'
            ax.set_facecolor('#252525')
            ax.tick_params(colors=_fg, labelcolor=_fg)
            ax.xaxis.label.set_color(_fg)
            ax.yaxis.label.set_color(_fg)
            ax.title.set_color(_fg)
            for spine in ax.spines.values():
                spine.set_edgecolor('#666666')

        def _theme_color(self, key):
            """Return a plot colour appropriate for the current theme."""
            _DARK = {
                'spectrum_raw':    '#aacce8',   # light steel-blue
                'spectrum_smooth': '#f5a623',   # amber
                'template':        '#2ecc71',   # bright green
                'chi2_highlight':  '#eeeeee',   # near-white
                'line_marker':     '#ff6b6b',   # salmon-red (brighter than #e74c3c on dark)
                'line_label':      '#ff6b6b',
                'continuum':       '#ff6b6b',
                'error_fill':      '#4a6278',   # muted blue-gray
            }
            _LIGHT = {
                'spectrum_raw':    '#2c3e50',   # dark navy
                'spectrum_smooth': '#d35400',   # burnt orange
                'template':        '#1a7a3f',   # dark green
                'chi2_highlight':  '#333333',   # near-black
                'line_marker':     '#c0392b',   # dark red
                'line_label':      '#c0392b',
                'continuum':       '#c0392b',
                'error_fill':      '#bdc3c7',   # light gray
            }
            return (_DARK if self.dark_theme else _LIGHT).get(key, 'gray')

        def _reset_plots(self):
            self._ax.clear()
            self._dark_ax(self._ax)
            self._ax.set_xlabel('Redshift  z')
            self._ax.set_ylabel('score / chi2')
            self._ax.set_title('Run search to see results')
            self._ax_spec.clear()
            self._dark_ax(self._ax_spec)
            self._ax_spec.set_xlabel('Wavelength  (Å)')
            self._ax_spec.set_ylabel('Flux')
            self._canvas.draw()

        def _update_chi2_plot(self, highlight_z=None):
            self._ax.clear()
            self._dark_ax(self._ax)
            if self._result is None:
                self._canvas.draw()
                return

            z = self._result.z_array
            curves = self._result.chi2_curves
            method = self._method_combo.currentText()

            # Curve colours — single curve gets steelblue; multiple get palette
            _curve_colors = ['steelblue', '#e67e22', '#27ae60', '#8e44ad',
                             '#e74c3c', '#16a085', '#7f8c8d']

            for ci, curve in enumerate(curves):
                col   = _curve_colors[ci % len(_curve_colors)]
                lw    = 0.9 if len(curves) > 1 else 0.8
                alpha = 0.7 if len(curves) > 1 else 0.85
                self._ax.plot(z, curve['chi2'], color=col, lw=lw,
                              alpha=alpha, label=curve['label'])
            self._ax.set_ylabel('score / chi2')
            items = self._result.solutions

            for i, item in enumerate(items):
                self._ax.axvline(item.z,
                                 color=_SOL_COLORS[i % len(_SOL_COLORS)],
                                 lw=1.6 if i == 0 else 0.9,
                                 ls='--', alpha=0.8,
                                 label=f'z={item.z:.4f}')

            if highlight_z is not None:
                self._ax.axvline(highlight_z,
                                 color=self._theme_color('chi2_highlight'),
                                 lw=2.0, ls='-', alpha=0.45)

            self._ax.set_xlabel('Redshift  z')
            # Build title tag
            if method == 'Template':
                tag = self._tmpl_combo.currentText()
            elif method == 'PCA':
                tag = f'PCA:{self._pca_combo.currentText()}'
            else:
                tag = self._pf_linelist_combo.currentText()
            self._ax.set_title(
                f'{tag}   (click anywhere to set z)', fontsize=9)
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
            self._raw_flux = self.spec.flux.value.copy()
            self._raw_err  = (self.spec.sig.value.copy()
                              if self.spec.sig_is_set else None)
            self._spec_wave_full = wave.copy()
            self._update_spectrum_plot(z=None)

        def _update_spectrum_plot(self, z=None, flux_override=None, _no_draw=False):
            ax = self._ax_spec
            ax.clear()
            self._dark_ax(ax)

            if self.spec is None:
                self._canvas.draw()
                return

            import astropy.units as u
            wave = self.spec.wavelength.to(u.AA).value

            # Raw flux — always available
            if flux_override is not None:
                raw = flux_override
            elif self._raw_flux is not None:
                raw = self._raw_flux
            else:
                raw = self.spec.flux.value
                self._raw_flux = raw.copy()

            # --- raw spectrum ---
            raw_color = self._theme_color('spectrum_raw')
            # Dim raw slightly when a smooth overlay will be drawn on top
            raw_alpha = 0.45 if self._smooth_scale > 1 else 0.9
            ax.plot(wave, raw, color=raw_color, lw=0.6, alpha=raw_alpha,
                    label='raw')

            # --- smoothed overlay (when display smooth > 1) ---
            if self._smooth_scale > 1:
                flux_smooth = convolve(raw, Box1DKernel(self._smooth_scale))
                ax.plot(wave, flux_smooth,
                        color=self._theme_color('spectrum_smooth'),
                        lw=1.0, alpha=0.9,
                        label=f'smooth×{self._smooth_scale}')
                flux_for_ylim = flux_smooth
            else:
                flux_for_ylim = raw

            # --- error envelope (on raw) ---
            if self._raw_err is not None:
                ax.fill_between(wave, raw - self._raw_err, raw + self._raw_err,
                                color=self._theme_color('error_fill'),
                                alpha=0.30, lw=0)

            # --- continuum (if stored on spec) ---
            if self.spec.co_is_set:
                ax.plot(wave, self.spec.co.value,
                        color=self._theme_color('continuum'),
                        lw=0.8, ls='--', alpha=0.7, label='continuum')

            ax.set_xlabel('Wavelength  (Å)', fontsize=8)
            ax.set_ylabel('Flux', fontsize=8)
            ax.tick_params(labelsize=7)

            finite = flux_for_ylim[np.isfinite(flux_for_ylim)]
            if len(finite) > 10:
                ylo = np.percentile(finite, 1)
                yhi = np.percentile(finite, 99)
                pad = 0.15 * (yhi - ylo) if yhi > ylo else 0.1 * abs(yhi)
                ax.set_ylim(ylo - pad, yhi + pad)

            # --- template / PCA model overlay ---
            if z is not None and self._result is not None:
                if (self._last_run_method == 'Template'
                        and self._last_run_template_name not in (None, 'Best of all')):
                    try:
                        from rbcodes.GUIs.zfind.engine import (
                            _load_template, _preprocess)
                        t_wave, t_flux, t_pc = _load_template(
                            self._last_run_template_name)
                        t_obs = t_wave * (1.0 + z)

                        # ivar from spectrum error
                        if self.spec.sig_is_set:
                            sig = self.spec.sig.value
                            ivar_ov = np.where(sig > 0, 1.0 / sig**2, 0.0)
                        else:
                            ivar_ov = np.ones(len(wave))

                        data_norm  = self._last_run_data_norm
                        model_norm = self._last_run_model_norm

                        if data_norm != 'raw' or model_norm != 'raw':
                            # Need data continuum to fit in normalized space
                            # and convert back to flux units.
                            _w, _f, _iv, data_cont, _ = _preprocess(self.spec)
                            data_cont = np.interp(wave, _w, data_cont,
                                                  left=data_cont[0],
                                                  right=data_cont[-1])
                            cont_safe = np.where(
                                np.abs(data_cont) > 1e-10, data_cont, 1.0)

                            if data_norm == 'normalize':
                                F = raw / cont_safe
                            elif data_norm == 'subtract':
                                F = raw - data_cont
                            else:
                                F = raw

                            # template in model_norm space
                            if model_norm == 'normalize':
                                T_work = t_flux / t_pc   # both already array
                            elif model_norm == 'subtract':
                                T_work = t_flux - t_pc
                            else:
                                T_work = t_flux
                            T = np.interp(wave, t_obs, T_work,
                                          left=0.0, right=0.0)
                            in_range = ((wave >= t_obs.min())
                                        & (wave <= t_obs.max()))

                            good  = in_range & np.isfinite(F) & (ivar_ov > 0)
                            denom = np.sum(T[good]**2 * ivar_ov[good])
                            if good.sum() > 50 and denom > 0:
                                A = np.sum(F[good] * T[good] * ivar_ov[good]) / denom
                                if data_norm == 'normalize':
                                    model = A * T * cont_safe
                                elif data_norm == 'subtract':
                                    model = A * T + data_cont
                                else:
                                    model = A * T
                                model[~in_range] = np.nan
                                ax.plot(wave, model,
                                        color=self._theme_color('template'),
                                        lw=1.1, alpha=0.85, ls='-',
                                        label=self._last_run_template_name)
                        else:
                            # raw mode
                            T = np.interp(wave, t_obs, t_flux,
                                          left=0.0, right=0.0)
                            in_range = ((wave >= t_obs.min())
                                        & (wave <= t_obs.max()))
                            good  = in_range & np.isfinite(raw) & (ivar_ov > 0)
                            denom = np.sum(T[good]**2 * ivar_ov[good])
                            if good.sum() > 50 and denom > 0:
                                A = np.sum(raw[good] * T[good] * ivar_ov[good]) / denom
                                model = A * T
                                model[~in_range] = np.nan
                                ax.plot(wave, model,
                                        color=self._theme_color('template'),
                                        lw=1.1, alpha=0.85, ls='-',
                                        label=self._last_run_template_name)
                    except Exception as exc:
                        self._status.setText(f'Template overlay failed: {exc}')

                elif (self._last_run_method == 'PCA'
                        and self._last_run_pca_set is not None):
                    try:
                        from rbcodes.GUIs.zfind.engine import _load_pca
                        t_wave, eigenvecs = _load_pca(self._last_run_pca_set)
                        # Apply same model_norm as the run
                        mn = self._pca_model_norm.currentText()
                        if mn == 'normalize':
                            norms = np.sqrt(
                                np.sum(eigenvecs**2, axis=1, keepdims=True))
                            eigenvecs = eigenvecs / np.where(norms > 0, norms, 1.0)
                        elif mn == 'subtract':
                            eigenvecs = (eigenvecs
                                         - eigenvecs.mean(axis=1, keepdims=True))
                        # Build basis matrix at this z
                        t_obs = t_wave * (1.0 + z)
                        n_c = eigenvecs.shape[0]
                        M = np.column_stack([
                            np.interp(wave, t_obs, eigenvecs[j],
                                      left=0.0, right=0.0)
                            for j in range(n_c)
                        ])
                        coverage = np.any(M != 0.0, axis=1)
                        # ivar for weighting
                        if self.spec.sig_is_set:
                            sig = self.spec.sig.value
                            ivar_d = np.where(sig > 0, 1.0 / sig**2, 0.0)
                        else:
                            ivar_d = np.ones(len(wave))
                        good = coverage & np.isfinite(raw) & (ivar_d > 0)
                        if good.sum() > 50:
                            Mg = M[good]; Fg = raw[good]; IVg = ivar_d[good]
                            WMg = Mg * IVg[:, np.newaxis]
                            c, _, _, _ = np.linalg.lstsq(
                                WMg.T @ Mg, WMg.T @ Fg, rcond=None)
                            model = M @ c
                            model[~coverage] = np.nan
                            ax.plot(wave, model,
                                    color=self._theme_color('template'),
                                    lw=1.1, alpha=0.85, ls='-',
                                    label=f'PCA:{self._last_run_pca_set}')
                    except Exception as exc:
                        self._status.setText(f'PCA overlay failed: {exc}')

            # --- line overlay ---
            lc = self._theme_color('line_marker')
            if z is not None and self._df is not None:
                wmin, wmax = wave.min(), wave.max()
                obs   = self._df['wave'].values * (1.0 + z)
                names = self._df['name'].values
                mask  = (obs > wmin) & (obs < wmax)

                ylo_ax, yhi_ax = ax.get_ylim()
                label_y = yhi_ax - 0.04 * (yhi_ax - ylo_ax)

                for ow, nm in zip(obs[mask], names[mask]):
                    ax.axvline(ow, color=lc, lw=0.9, ls='--', alpha=0.75)
                    short = nm.split()[0] if ' ' in nm else nm
                    ax.text(ow, label_y, short, rotation=90,
                            fontsize=6, color=lc,
                            va='top', ha='right', alpha=0.85)

                title_parts = [f'z = {z:.5f}  ({int(mask.sum())} lines in range)']
                if self._smooth_scale > 1:
                    title_parts.append(f'smooth×{self._smooth_scale}')
                title = 'Spectrum — ' + '   '.join(title_parts)
            else:
                title = 'Spectrum'
                if self._smooth_scale > 1:
                    title += f'   smooth×{self._smooth_scale}'
                if self._result is None:
                    title += '   (run search then click chi2 plot to set z)'

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
                xlim = self._ax_spec.get_xlim()
                ylim = self._ax_spec.get_ylim()
                self._update_spectrum_plot(z=self._selected_z, _no_draw=True)
                self._ax_spec.set_xlim(xlim)
                self._ax_spec.set_ylim(ylim)
                self._canvas.draw()
                return
            if key == 'U':
                self._smooth_scale = max(1, self._smooth_scale - 2)
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
                        curve = self._result.chi2_curves[0]['chi2']
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
                z = float(event.xdata)
                self._selected_z = z

                # If click is within 20 % of the visible z range from a
                # solution, snap the table selection to that solution.
                if self._result and self._result.solutions:
                    zlim   = self._ax.get_xlim()
                    ax_width_px = self._ax.get_window_extent().width or 1
                    thresh = 5.0 * (zlim[1] - zlim[0]) / ax_width_px
                    dists  = [abs(s.z - z)
                              for s in self._result.solutions]
                    nearest_idx = int(np.argmin(dists))
                    if dists[nearest_idx] <= thresh:
                        # Snap selected_z to the solution z
                        self._selected_z = self._result.solutions[nearest_idx].z
                        z = self._selected_z
                        # Update table without re-triggering spectrum update
                        self._em_table.blockSignals(True)
                        self._em_table.selectRow(nearest_idx)
                        self._em_table.blockSignals(False)
                        self._selected_em_idx = nearest_idx

                self._show_z_status(z)
                self._update_chi2_plot(highlight_z=z)
                # Preserve whatever zoom/pan the user has on the spectrum panel
                xlim = self._ax_spec.get_xlim()
                ylim = self._ax_spec.get_ylim()
                self._update_spectrum_plot(z=z, _no_draw=True)
                self._ax_spec.set_xlim(xlim)
                self._ax_spec.set_ylim(ylim)
                self._canvas.draw()
                return

        # ------------------------------------------------------------------
        # Table population
        # ------------------------------------------------------------------

        def _update_table(self):
            if self._result is None:
                return
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

        # ------------------------------------------------------------------
        # Status bar z display
        # ------------------------------------------------------------------

        def _show_z_status(self, z):
            """Update status bar with current z and nearest z_err."""
            z_err = float('nan')
            if self._result and self._result.solutions:
                dists = [abs(s.z - z) for s in self._result.solutions]
                near  = self._result.solutions[int(np.argmin(dists))]
                z_err = near.z_err
            z_err_str = f'{z_err:.5f}' if np.isfinite(z_err) else 'nan'
            self._status.setText(
                f'z = {z:.5f}   z_err = {z_err_str}   '
                f'[click Accept to send to GUI]')

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
                self._show_z_status(z)
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
            if self._result is None or self._selected_z is None:
                return
            z_out = self._selected_z

            # Get z_err from the nearest solution (may differ if user clicked
            # between solutions on the chi² plot).
            z_err_out = float('nan')
            if self._result.solutions:
                dists = [abs(s.z - z_out) for s in self._result.solutions]
                near  = self._result.solutions[int(np.argmin(dists))]
                z_err_out = near.z_err

            payload = [z_out, z_err_out]
            self.accepted_z.emit(payload)
            print(f'rb_zfind accepted z = {z_out:.6f}  z_err = {z_err_out:.6f}')
            self.accept()

except ImportError:
    pass
