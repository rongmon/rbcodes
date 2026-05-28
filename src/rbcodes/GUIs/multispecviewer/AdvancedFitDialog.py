"""
AdvancedFitDialog.py — Advanced multi-Gaussian line fitting dialog for rb_multispec.

Opened with G (shift+g) from the main canvas (cursor panel).

Fits 1–5 Gaussian components with a shared redshift parameter, returning
z ± error from the covariance matrix of scipy curve_fit.

Key interactions inside the dialog canvas
------------------------------------------
  Shift+C   Set z_guess from cursor position (uses Ion 1 rest wavelength)
  d → d     Define two manual continuum anchor points
  x / X     Set left / right x-limit
  [ / ]     Pan left / right
  t / b     Set y max / min
  r         Reset canvas view
"""

import numpy as np
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QComboBox, QPushButton, QCheckBox,
    QRadioButton, QButtonGroup, QSpinBox, QTextEdit,
    QFileDialog, QApplication, QSizePolicy, QWidget,
    QDialogButtonBox, QFrame
)
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

from astropy.convolution import convolve, Box1DKernel
from rbcodes.GUIs.multispecviewer.utils import read_line_options
from rbcodes.IGM import rb_setline as rb_setline


# ─────────────────────────────────────────────────────────────────────────────
# Dark-theme stylesheet shared across dialogs
# ─────────────────────────────────────────────────────────────────────────────
_DARK_STYLE = """
    QDialog, QWidget { background-color: #353535; color: #F2F2F7; }
    QLabel  { color: #F2F2F7; }
    QLineEdit, QTextEdit, QSpinBox {
        background-color: #3C3C3C; color: white;
        border: 1px solid #555; border-radius: 2px; padding: 2px;
        selection-background-color: #2A82DA;
    }
    QComboBox {
        background-color: #3A3A3C; color: #F2F2F7;
        border: 1px solid #636366; border-radius: 2px; padding: 2px;
    }
    QComboBox QAbstractItemView {
        background-color: #3A3A3C; color: #F2F2F7;
        selection-background-color: #0A84FF;
    }
    QComboBox::drop-down { background-color: #48484A; border-left: 1px solid #636366; width: 20px; }
    QPushButton {
        background-color: #474747; color: #F2F2F2;
        border: none; border-radius: 4px; padding: 5px 12px;
    }
    QPushButton:hover  { background-color: #555; }
    QPushButton:pressed { background-color: #2A2A2A; }
    QPushButton:disabled { background-color: #333; color: #666; }
    QCheckBox, QRadioButton { color: #F2F2F7; }
    QFrame[frameShape="4"] { color: #555; }
"""


# ─────────────────────────────────────────────────────────────────────────────
# Multi-Gaussian model
# ─────────────────────────────────────────────────────────────────────────────

class MultiGauss:
    """
    N-component Gaussian model in the rest frame with a shared redshift.

    Parameter vector layout (what curve_fit sees):
      default:           [z, sig_1, …, sig_N, amp_1, …, amp_N]
      tie_sigma=True:    [z, sig,   amp_1, …, amp_N]
      fix_z=<value>:     [sig_1, …, sig_N, amp_1, …, amp_N]
      both:              [sig,   amp_1, …, amp_N]

    sigma units: Å in the rest frame
    """

    def __init__(self, rest_waves, tie_sigma=False, fix_z=None):
        self.rest_waves = np.asarray(rest_waves, dtype=float)
        self.n = len(self.rest_waves)
        self.tie_sigma = tie_sigma
        self.fix_z = fix_z

    @staticmethod
    def _gauss(wave_rest, sig, amp, mu_rest):
        return amp * np.exp(-(wave_rest - mu_rest) ** 2 / (2.0 * sig ** 2))

    def compile_model(self, wave_obs, *params):
        params = list(params)
        idx = 0
        z = params[idx] if self.fix_z is None else self.fix_z
        if self.fix_z is None:
            idx += 1
        sigs = ([params[idx]] * self.n) if self.tie_sigma else list(params[idx:idx + self.n])
        idx += 1 if self.tie_sigma else self.n
        amps = list(params[idx:idx + self.n])
        wave_rest = wave_obs / (1.0 + z)
        total = np.zeros(len(wave_obs))
        for i in range(self.n):
            total += self._gauss(wave_rest, sigs[i], amps[i], self.rest_waves[i])
        return total

    def n_params(self):
        return ((0 if self.fix_z is not None else 1) +
                (1 if self.tie_sigma else self.n) +
                self.n)

    def make_initial_params(self, z_guess, sig_guess, amp_guess):
        p = []
        if self.fix_z is None:
            p.append(z_guess)
        p.append(sig_guess[0] if self.tie_sigma else None)
        if self.tie_sigma:
            pass
        else:
            p = p + list(sig_guess[:self.n])
            p = [x for x in p if x is not None]
        if self.tie_sigma:
            p = ([z_guess] if self.fix_z is None else []) + [sig_guess[0]] + list(amp_guess[:self.n])
        else:
            p = ([z_guess] if self.fix_z is None else []) + list(sig_guess[:self.n]) + list(amp_guess[:self.n])
        return p

    def make_bounds(self, z_guess, amp_guess, v_uncer_kms=1000.0,
                    sig_lo=0.1, sig_hi=200.0):
        beta = v_uncer_kms / 299792.458
        dz = np.sqrt((1 + beta) / (1 - beta)) - 1
        lo, hi = [], []
        if self.fix_z is None:
            lo.append(z_guess - dz); hi.append(z_guess + dz)
        n_sig = 1 if self.tie_sigma else self.n
        lo.extend([sig_lo] * n_sig); hi.extend([sig_hi] * n_sig)
        lo.extend([0.5 * abs(a) for a in amp_guess[:self.n]])
        hi.extend([2.0 * abs(a) + 1e-6 for a in amp_guess[:self.n]])
        return lo, hi

    def param_names(self, ion_names):
        names = []
        if self.fix_z is None:
            names.append('z')
        if self.tie_sigma:
            names.append('σ (Å, rest)')
        else:
            for i, n in enumerate(ion_names):
                names.append(f'σ_{i+1} {n} (Å, rest)')
        for i, n in enumerate(ion_names):
            names.append(f'Amp_{i+1} {n}')
        return names


# ─────────────────────────────────────────────────────────────────────────────
# Fitting constraints dialog
# ─────────────────────────────────────────────────────────────────────────────

class FittingConstraintDialog(QDialog):
    """Edit initial guesses and bounds; toggle tie_sigma / fix_z."""

    def __init__(self, init_guess, bounds_lo, bounds_hi, param_names,
                 tie_sigma=False, fix_z=False, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Fitting Constraints')
        self.setModal(True)
        self.setStyleSheet(_DARK_STYLE)

        layout = QGridLayout()
        for col, hdr in enumerate(['Parameter', 'Guess', 'Lower', 'Upper']):
            lbl = QLabel(f'<b>{hdr}</b>')
            lbl.setStyleSheet('color: #F2F2F7;')
            layout.addWidget(lbl, 0, col)

        self._fields = []
        for row, (name, g, lo, hi) in enumerate(
                zip(param_names, init_guess, bounds_lo, bounds_hi), start=1):
            layout.addWidget(QLabel(name), row, 0)
            g_le  = QLineEdit(f'{g:.5g}');  g_le.setValidator(QDoubleValidator())
            lo_le = QLineEdit(f'{lo:.5g}'); lo_le.setValidator(QDoubleValidator())
            hi_le = QLineEdit(f'{hi:.5g}'); hi_le.setValidator(QDoubleValidator())
            for col, le in enumerate([g_le, lo_le, hi_le], start=1):
                layout.addWidget(le, row, col)
            self._fields.append((g_le, lo_le, hi_le))

        row = len(param_names) + 1
        self.tie_cb = QCheckBox('Tie sigmas — all lines share one σ')
        self.tie_cb.setChecked(tie_sigma)
        layout.addWidget(self.tie_cb, row, 0, 1, 4)

        self.fixz_cb = QCheckBox('Fix z at z_guess — only fit amplitudes and sigmas')
        self.fixz_cb.setChecked(fix_z)
        layout.addWidget(self.fixz_cb, row + 1, 0, 1, 4)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons, row + 2, 0, 1, 4)
        self.setLayout(layout)

    def get_values(self):
        """Return (new_guess, new_lo, new_hi, tie_sigma, fix_z)."""
        new_guess, new_lo, new_hi = [], [], []
        for g_le, lo_le, hi_le in self._fields:
            try:
                new_guess.append(float(g_le.text()))
                new_lo.append(float(lo_le.text()))
                new_hi.append(float(hi_le.text()))
            except ValueError:
                pass
        return new_guess, new_lo, new_hi, self.tie_cb.isChecked(), self.fixz_cb.isChecked()


# ─────────────────────────────────────────────────────────────────────────────
# Main dialog
# ─────────────────────────────────────────────────────────────────────────────

class AdvancedFitDialog(QDialog):
    """
    Non-modal dialog for advanced multi-Gaussian line fitting.

    Parameters
    ----------
    wave, flux, error : ndarray
        Original *unsmoothed* spectrum from the cursor panel.
        Pass error=None if unavailable.
    current_z        : float
    current_linelist : str
    parent_window    : MainWindow
    panel_index      : int
    filename         : str  (for the window title)
    """

    def __init__(self, wave, flux, error, current_z, current_linelist,
                 parent_window, panel_index=0, filename='', xlim=None):
        super().__init__(parent_window)
        self.setWindowTitle(
            f'Advanced Line Fit — Panel {panel_index + 1}:  {filename}')
        self.setMinimumSize(1050, 900)
        self.setModal(False)
        self.setStyleSheet(_DARK_STYLE)

        # ── original arrays (never modified) ──────────────────────────────
        self._wave_orig  = np.asarray(wave,  dtype=float)
        self._flux_orig  = np.asarray(flux,  dtype=float)
        self._error_orig = (np.asarray(error, dtype=float)
                            if error is not None else None)

        self._parent      = parent_window
        self._panel_index = panel_index

        # ── displayed arrays (may be smoothed) ────────────────────────────
        self._wave  = self._wave_orig.copy()
        self._flux  = self._flux_orig.copy()
        self._error = (self._error_orig.copy()
                       if self._error_orig is not None
                       else np.ones_like(self._flux_orig))

        # ── state ─────────────────────────────────────────────────────────
        self._smooth_kernel  = 1
        self._cont_mode      = 'auto'  # 'auto' | 'manual'
        self._n_edge         = 15
        self._cont_anchor    = None    # (x1,y1,x2,y2) for manual mode
        self._d_pending      = None    # first d-keystroke position
        self._cont_artist    = None    # yellow dashed continuum line
        self._fit_artist     = None    # orange fit overlay
        self._guess_artists  = []      # blue dashed lines from C keystroke

        self._xlim           = xlim  # (x1, x2) from G+G keypresses, or None
        self._z_guess        = float(current_z)
        self._linelist_name  = current_linelist
        self._linelist_data  = []      # list of dicts from rb_setline
        self._ion_widgets    = []      # 5 QComboBox widgets
        self._selected_rest  = []      # rest wavelengths of chosen ions
        self._selected_names = []      # ion label strings

        self._fit_result   = None  # (z, z_err)
        self._popt = self._perr = None
        self._fit_cont = None
        self._fit_wave = self._fit_flux_model = None
        self._model      = None
        self._init_guess = None
        self._bounds_lo  = self._bounds_hi = None
        self._tie_sigma  = False
        self._fix_z      = False
        self._param_names = []

        self._initUI()
        self._load_linelist(current_linelist)
        self._plot_spectrum()
        self._update_continuum()

    # ─────────────────────────────────────────────────────────────────
    # UI construction
    # ─────────────────────────────────────────────────────────────────

    def _initUI(self):
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(6)
        main_layout.setContentsMargins(8, 8, 8, 8)

        # 1. Canvas ──────────────────────────────────────────────────
        self._fig = Figure(figsize=(11, 4), dpi=100)
        self._fig.patch.set_facecolor('#353535')
        self._canvas = FigureCanvas(self._fig)
        self._canvas.setFocusPolicy(Qt.ClickFocus)
        self._canvas.mpl_connect('key_press_event', self._on_canvas_key)
        self._canvas.mpl_connect('button_press_event', self._on_canvas_click)
        self._ax = self._fig.add_subplot(111)
        self._ax.set_facecolor('#353535')
        self._fig.tight_layout(pad=0.5)

        toolbar = NavigationToolbar(self._canvas, self)

        smooth_row = QHBoxLayout()
        smooth_row.addWidget(QLabel('Smooth (kernel):'))
        self._smooth_spin = QSpinBox()
        self._smooth_spin.setRange(1, 51)
        self._smooth_spin.setSingleStep(2)
        self._smooth_spin.setValue(1)
        self._smooth_spin.setFixedWidth(60)
        self._smooth_spin.valueChanged.connect(self._on_smooth_changed)
        smooth_row.addWidget(self._smooth_spin)
        smooth_row.addStretch()

        main_layout.addWidget(toolbar)
        main_layout.addWidget(self._canvas, stretch=5)
        main_layout.addLayout(smooth_row)

        main_layout.addWidget(self._hline())

        # 2. Linelist + z_guess + status ─────────────────────────────
        ctrl_row = QHBoxLayout()
        ctrl_row.addWidget(QLabel('Linelist:'))
        self._linelist_combo = QComboBox()
        self._linelist_combo.addItems(read_line_options())
        self._linelist_combo.setCurrentText(self._linelist_name)
        self._linelist_combo.currentTextChanged.connect(self._on_linelist_changed)
        self._linelist_combo.setFixedWidth(130)
        ctrl_row.addWidget(self._linelist_combo)

        ctrl_row.addSpacing(20)
        ctrl_row.addWidget(QLabel('z guess:'))
        self._z_guess_le = QLineEdit(f'{self._z_guess:.6f}')
        self._z_guess_le.setValidator(QDoubleValidator(0, 20, 8))
        self._z_guess_le.setFixedWidth(100)
        ctrl_row.addWidget(self._z_guess_le)
        ctrl_row.addWidget(QLabel('(Shift+C on canvas to set from cursor)'))

        ctrl_row.addStretch()
        self._status_label = QLabel('Select ≥ 1 ion')
        self._status_label.setStyleSheet('color: #FFA500; font-weight: bold;')
        ctrl_row.addWidget(self._status_label)
        main_layout.addLayout(ctrl_row)

        # 3. Ion dropdowns with Lines spinbox ───────────────────────
        ion_header = QHBoxLayout()
        ion_header.addWidget(QLabel('Lines to fit:'))
        self._n_lines_spin = QSpinBox()
        self._n_lines_spin.setRange(1, 5)
        self._n_lines_spin.setValue(2)
        self._n_lines_spin.setFixedWidth(45)
        self._n_lines_spin.valueChanged.connect(self._on_n_lines_changed)
        ion_header.addWidget(self._n_lines_spin)
        ion_header.addStretch()
        main_layout.addLayout(ion_header)

        ion_grid = QGridLayout()
        ion_grid.setSpacing(4)
        self._ion_labels = []
        for col in range(5):
            lbl = QLabel(f'Ion {col+1}')
            ion_grid.addWidget(lbl, 0, col)
            self._ion_labels.append(lbl)
            cb = QComboBox()
            cb.setFixedWidth(160)
            cb.addItem('NONE')
            cb.currentIndexChanged.connect(
                lambda idx, c=col: self._on_ion_changed(idx, c))
            ion_grid.addWidget(cb, 1, col)
            self._ion_widgets.append(cb)
        main_layout.addLayout(ion_grid)
        # Apply initial visibility (default 2 lines)
        self._apply_n_lines_visibility(2)

        main_layout.addWidget(self._hline())

        # 4. Continuum ───────────────────────────────────────────────
        cont_row = QHBoxLayout()
        self._auto_rb   = QRadioButton('Auto  (median of edge pixels)')
        self._manual_rb = QRadioButton('Manual  (press d → d on canvas)')
        self._auto_rb.setChecked(True)
        grp = QButtonGroup(self)
        grp.addButton(self._auto_rb)
        grp.addButton(self._manual_rb)
        self._auto_rb.toggled.connect(self._on_cont_mode_changed)

        cont_row.addWidget(self._auto_rb)
        cont_row.addWidget(QLabel('Edge pixels:'))
        self._edge_spin = QSpinBox()
        self._edge_spin.setRange(3, 100)
        self._edge_spin.setValue(self._n_edge)
        self._edge_spin.setFixedWidth(55)
        self._edge_spin.valueChanged.connect(self._on_edge_changed)
        cont_row.addWidget(self._edge_spin)

        cont_row.addSpacing(20)
        cont_row.addWidget(self._manual_rb)

        cont_row.addStretch()
        self._cont_status = QLabel(f'Continuum: auto — {self._n_edge} px median')
        self._cont_status.setStyleSheet('color: #FFD700;')
        cont_row.addWidget(self._cont_status)
        main_layout.addLayout(cont_row)

        main_layout.addWidget(self._hline())

        # 5. Error weighting + Fit + Advanced ────────────────────────
        fit_row = QHBoxLayout()
        self._weight_cb = QCheckBox('Weight by flux errors')
        self._weight_cb.setChecked(True)
        fit_row.addWidget(self._weight_cb)
        fit_row.addSpacing(30)

        self._fit_btn = QPushButton('Fit')
        self._fit_btn.setFixedWidth(90)
        self._fit_btn.setEnabled(False)
        self._fit_btn.clicked.connect(self._on_fit_clicked)
        fit_row.addWidget(self._fit_btn)

        self._adv_btn = QPushButton('Advanced...')
        self._adv_btn.setFixedWidth(110)
        self._adv_btn.setEnabled(False)
        self._adv_btn.clicked.connect(self._on_advanced_clicked)
        fit_row.addWidget(self._adv_btn)
        fit_row.addStretch()
        main_layout.addLayout(fit_row)

        main_layout.addWidget(self._hline())

        # 6. Results ─────────────────────────────────────────────────
        main_layout.addWidget(QLabel('Results:'))
        self._results_te = QTextEdit()
        self._results_te.setReadOnly(True)
        self._results_te.setFixedHeight(130)
        self._results_te.setFont(QFont('Courier New', 10))
        self._results_te.setPlaceholderText(
            'Fit results will appear here after clicking Fit.')
        main_layout.addWidget(self._results_te, stretch=2)

        results_btn_row = QHBoxLayout()
        copy_btn = QPushButton('Copy to clipboard')
        copy_btn.clicked.connect(self._on_copy_clicked)
        export_btn = QPushButton('Export to file')
        export_btn.clicked.connect(self._on_export_clicked)
        results_btn_row.addWidget(copy_btn)
        results_btn_row.addWidget(export_btn)
        results_btn_row.addStretch()
        main_layout.addLayout(results_btn_row)

        main_layout.addWidget(self._hline())

        # 7. Action buttons ──────────────────────────────────────────
        action_row = QHBoxLayout()
        self._apply_btn = QPushButton('Apply z + linelist to main')
        self._apply_btn.setEnabled(False)
        self._apply_btn.clicked.connect(self._on_apply_clicked)
        action_row.addWidget(self._apply_btn)

        self._add_abs_btn = QPushButton('Add to absorbers')
        self._add_abs_btn.setEnabled(False)
        self._add_abs_btn.clicked.connect(self._on_add_absorber_clicked)
        action_row.addWidget(self._add_abs_btn)

        action_row.addStretch()
        close_btn = QPushButton('Close')
        close_btn.clicked.connect(self.close)
        action_row.addWidget(close_btn)
        main_layout.addLayout(action_row)

    @staticmethod
    def _hline():
        f = QFrame()
        f.setFrameShape(QFrame.HLine)
        f.setStyleSheet('color: #555;')
        return f

    # ─────────────────────────────────────────────────────────────────
    # Lines-to-fit spinbox
    # ─────────────────────────────────────────────────────────────────

    def _apply_n_lines_visibility(self, n):
        """Show the first n ion dropdowns; hide and reset the rest."""
        for col in range(5):
            visible = col < n
            self._ion_labels[col].setVisible(visible)
            self._ion_widgets[col].setVisible(visible)
            if not visible:
                self._ion_widgets[col].blockSignals(True)
                self._ion_widgets[col].setCurrentIndex(0)
                self._ion_widgets[col].blockSignals(False)

    def _on_n_lines_changed(self, n):
        self._apply_n_lines_visibility(n)
        # Reset fit state so stale guesses aren't carried over
        self._init_guess = None
        self._update_ion_selection()

    # ─────────────────────────────────────────────────────────────────
    # Linelist / ion management
    # ─────────────────────────────────────────────────────────────────

    def _load_linelist(self, name):
        self._linelist_name = name
        if name == 'None' or not name:
            self._linelist_data = []
            for cb in self._ion_widgets:
                cb.blockSignals(True)
                cb.clear(); cb.addItem('NONE')
                cb.blockSignals(False)
            return
        try:
            self._linelist_data = rb_setline.read_line_list(name)
        except Exception:
            self._linelist_data = []
        entries = ['NONE'] + [
            f"{d['ion']}  {d['wrest']:.2f}" for d in self._linelist_data]
        for cb in self._ion_widgets:
            cb.blockSignals(True)
            cb.clear()
            cb.addItems(entries)
            cb.blockSignals(False)
        self._update_ion_selection()

    def _on_linelist_changed(self, name):
        self._load_linelist(name)

    def _on_ion_changed(self, idx, col):
        # Auto-populate subsequent *visible* ions when Ion 1 is set
        n = self._n_lines_spin.value()
        if col == 0 and idx > 0:
            n_entries = self._ion_widgets[0].count()
            for k in range(1, n):
                next_idx = idx + k
                self._ion_widgets[k].blockSignals(True)
                self._ion_widgets[k].setCurrentIndex(
                    next_idx if next_idx < n_entries else 0)
                self._ion_widgets[k].blockSignals(False)
        self._update_ion_selection()

    def _update_ion_selection(self):
        rest, names = [], []
        for cb in self._ion_widgets:
            idx = cb.currentIndex()
            if idx > 0 and idx - 1 < len(self._linelist_data):
                d = self._linelist_data[idx - 1]
                rest.append(d['wrest'])
                names.append(d['ion'])
        self._selected_rest  = rest
        self._selected_names = names
        ready = len(rest) >= 1
        self._fit_btn.setEnabled(ready)
        self._status_label.setText(
            f'Ready  ({len(rest)} line{"s" if len(rest) != 1 else ""})'
            if ready else 'Select ≥ 1 ion')
        self._status_label.setStyleSheet(
            'color: #00FF7F; font-weight: bold;' if ready
            else 'color: #FFA500; font-weight: bold;')

    # ─────────────────────────────────────────────────────────────────
    # Spectrum display
    # ─────────────────────────────────────────────────────────────────

    def _plot_spectrum(self):
        self._ax.clear()
        self._ax.step(self._wave, self._flux, where='mid',
                      color='#FFFFFF', lw=0.7, label='flux')
        if self._error_orig is not None:
            self._ax.step(self._wave, self._error, where='mid',
                          color='#FF6B6B', lw=0.5, alpha=0.6, label='error')
        self._ax.set_facecolor('#353535')
        self._ax.tick_params(colors='white')
        self._ax.spines['bottom'].set_color('white')
        self._ax.spines['left'].set_color('white')
        self._ax.set_xlabel('Wavelength (Å)', color='white')
        self._ax.set_ylabel('Flux', color='white')

        # Apply window defined by the two G keypresses
        if self._xlim is not None:
            x1, x2 = self._xlim
            self._ax.set_xlim(x1, x2)
            # Autoscale y to the visible window
            mask = (self._wave >= x1) & (self._wave <= x2) & np.isfinite(self._flux)
            if mask.sum() > 0:
                vis = self._flux[mask]
                pad = 0.15 * (vis.max() - vis.min()) if vis.max() != vis.min() else 0.1
                self._ax.set_ylim(vis.min() - pad, vis.max() + pad)

        self._cont_artist   = None
        self._fit_artist    = None
        self._guess_artists = []
        self._canvas.draw_idle()

    def _on_smooth_changed(self, val):
        # Enforce odd kernel
        if val > 1 and val % 2 == 0:
            val += 1
            self._smooth_spin.blockSignals(True)
            self._smooth_spin.setValue(val)
            self._smooth_spin.blockSignals(False)
        self._smooth_kernel = val
        if val <= 1:
            self._flux  = self._flux_orig.copy()
            self._error = (self._error_orig.copy()
                           if self._error_orig is not None
                           else np.ones_like(self._flux_orig))
        else:
            self._flux  = convolve(self._flux_orig,  Box1DKernel(val))
            if self._error_orig is not None:
                self._error = convolve(self._error_orig, Box1DKernel(val))
        # Preserve current x/y limits
        xlim = self._ax.get_xlim()
        ylim = self._ax.get_ylim()
        self._plot_spectrum()
        self._ax.set_xlim(xlim)
        self._ax.set_ylim(ylim)
        self._update_continuum()

    # ─────────────────────────────────────────────────────────────────
    # Continuum management
    # ─────────────────────────────────────────────────────────────────

    def _compute_auto_continuum(self, wave, flux):
        n = min(self._n_edge, max(3, len(wave) // 5))
        x_l = float(np.median(wave[:n]))
        y_l = float(np.median(flux[:n]))
        x_r = float(np.median(wave[-n:]))
        y_r = float(np.median(flux[-n:]))
        cont = y_l + (y_r - y_l) / (x_r - x_l) * (wave - x_l)
        return cont, x_l, y_l, x_r, y_r

    def _compute_manual_continuum(self, wave):
        x1, y1, x2, y2 = self._cont_anchor
        if x1 > x2:
            x1, x2, y1, y2 = x2, x1, y2, y1
        cont = y1 + (y2 - y1) / (x2 - x1) * (wave - x1)
        return cont

    def _update_continuum(self):
        """Recompute and redraw the continuum on the canvas."""
        xlim = self._ax.get_xlim()
        mask = (self._wave >= xlim[0]) & (self._wave <= xlim[1])
        if mask.sum() < 4:
            return

        w = self._wave[mask]
        f = self._flux[mask]

        if self._cont_mode == 'auto':
            cont, x_l, y_l, x_r, y_r = self._compute_auto_continuum(w, f)
            # Draw across the full view
            cont_full = y_l + (y_r - y_l) / (x_r - x_l) * (self._wave - x_l)
            self._cont_anchor = (x_l, y_l, x_r, y_r)
            self._cont_status.setText(
                f'Continuum: auto — {self._n_edge} px median  '
                f'[{x_l:.1f} Å, {x_r:.1f} Å]')
        else:
            if self._cont_anchor is None:
                return
            cont_full = self._compute_manual_continuum(self._wave)
            x1, y1, x2, y2 = self._cont_anchor
            self._cont_status.setText(
                f'Continuum: manual  [{x1:.1f} Å → {x2:.1f} Å]')

        if self._cont_artist is not None:
            try:
                self._cont_artist.remove()
            except Exception:
                pass
        line, = self._ax.plot(self._wave, cont_full,
                              color='#FFD700', lw=1.0,
                              linestyle='--', alpha=0.8, zorder=4)
        self._cont_artist = line
        self._canvas.draw_idle()

    def _on_cont_mode_changed(self, auto_checked):
        self._cont_mode = 'auto' if auto_checked else 'manual'
        if self._cont_mode == 'auto':
            self._edge_spin.setEnabled(True)
            self._d_pending = None
            self._update_continuum()
        else:
            self._edge_spin.setEnabled(False)
            if self._cont_artist is not None:
                try:
                    self._cont_artist.remove()
                except Exception:
                    pass
                self._cont_artist = None
                self._canvas.draw_idle()
            self._cont_status.setText(
                'Continuum: manual — press d → d on canvas to set anchors')

    def _on_edge_changed(self, val):
        self._n_edge = val
        if self._cont_mode == 'auto':
            self._update_continuum()

    # ─────────────────────────────────────────────────────────────────
    # Canvas interaction
    # ─────────────────────────────────────────────────────────────────

    def _on_canvas_click(self, event):
        if event.button == 1:
            self._canvas.setFocus()

    def _on_canvas_key(self, event):
        x, y = event.xdata, event.ydata
        ax = self._ax

        if event.key == 'C':  # shift+c — set z_guess from cursor
            if not self._selected_rest or x is None:
                return
            z_new = x / self._selected_rest[0] - 1.0
            self._z_guess = z_new
            self._z_guess_le.setText(f'{z_new:.6f}')
            self._draw_guess_lines(z_new)
            return

        if event.key == 'd':  # manual continuum anchor
            if self._cont_mode != 'manual' or x is None:
                return
            if self._d_pending is None:
                self._d_pending = (x, y)
                self._cont_status.setText(
                    f'Continuum: manual — left anchor at {x:.2f} Å, '
                    f'press d again at right edge')
            else:
                x1, y1 = self._d_pending
                self._d_pending = None
                self._cont_anchor = (x1, y1, x, y)
                self._update_continuum()
            return

        if x is None or y is None:
            return

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        if event.key == 'x':
            ax.set_xlim(x, xlim[1]); self._update_continuum()
        elif event.key == 'X':
            ax.set_xlim(xlim[0], x); self._update_continuum()
        elif event.key == 't':
            ax.set_ylim(ylim[0], y)
        elif event.key == 'b':
            ax.set_ylim(y, ylim[1])
        elif event.key == '[':
            dx = xlim[1] - xlim[0]
            ax.set_xlim(xlim[0] - dx, xlim[0]); self._update_continuum()
        elif event.key == ']':
            dx = xlim[1] - xlim[0]
            ax.set_xlim(xlim[1], xlim[1] + dx); self._update_continuum()
        elif event.key == 'r':
            ax.set_xlim(self._wave.min(), self._wave.max())
            finite = self._flux[np.isfinite(self._flux)]
            if len(finite):
                pad = 0.1 * (finite.max() - finite.min())
                ax.set_ylim(finite.min() - pad, finite.max() + pad)
            self._update_continuum()

        self._canvas.draw_idle()

    def _draw_guess_lines(self, z):
        for art in self._guess_artists:
            try: art.remove()
            except Exception: pass
        self._guess_artists = []
        ylim = self._ax.get_ylim()
        for rest, name in zip(self._selected_rest, self._selected_names):
            obs = rest * (1.0 + z)
            vl = self._ax.axvline(obs, color='#4FC3F7', linestyle='--',
                                  lw=1.0, alpha=0.8)
            tx = self._ax.text(obs, ylim[1] * 0.85, name,
                               rotation=90, color='#4FC3F7', fontsize=8,
                               ha='right', va='top')
            self._guess_artists.extend([vl, tx])
        self._canvas.draw_idle()

    # ─────────────────────────────────────────────────────────────────
    # Fitting
    # ─────────────────────────────────────────────────────────────────

    def _on_fit_clicked(self):
        try:
            self._z_guess = float(self._z_guess_le.text())
        except ValueError:
            self._status_label.setText('Invalid z_guess')
            return

        if not self._selected_rest:
            self._status_label.setText('Select ≥ 1 ion')
            return

        xlim = self._ax.get_xlim()
        mask = ((self._wave >= xlim[0]) & (self._wave <= xlim[1]) &
                np.isfinite(self._flux))
        if mask.sum() < 5:
            self._status_label.setText('Too few pixels in view')
            return

        w = self._wave[mask]
        f = self._flux[mask]
        e = self._error[mask]

        # Continuum
        if self._cont_mode == 'auto':
            cont, *_ = self._compute_auto_continuum(w, f)
        else:
            if self._cont_anchor is None:
                self._status_label.setText('Set manual continuum first (d→d)')
                return
            cont = self._compute_manual_continuum(w)

        residual = f - cont
        direction = -1 if residual.sum() < 0 else 1
        ydata = direction * residual  # always positive for fitting

        # Amplitude guess at each line's wavelength
        rest_arr = np.array(self._selected_rest)
        obs_guess = rest_arr * (1.0 + self._z_guess)
        amp_guess = []
        for wg in obs_guess:
            idx = np.argmin(np.abs(w - wg))
            amp_guess.append(float(np.clip(ydata[idx], 1e-6, None)))

        sig_guess = [20.0] * len(self._selected_rest)

        model = MultiGauss(rest_arr, tie_sigma=self._tie_sigma,
                           fix_z=self._z_guess if self._fix_z else None)
        self._model = model

        if self._init_guess is None:
            p0 = model.make_initial_params(self._z_guess, sig_guess, amp_guess)
            lo, hi = model.make_bounds(self._z_guess, amp_guess)
        else:
            p0  = self._init_guess
            lo  = self._bounds_lo
            hi  = self._bounds_hi

        self._param_names = model.param_names(self._selected_names)

        use_weights = self._weight_cb.isChecked()
        sigma_arg = e if use_weights else None

        try:
            popt, pcov = curve_fit(
                model.compile_model, w, ydata,
                p0=p0, bounds=(lo, hi),
                sigma=sigma_arg, absolute_sigma=use_weights,
                maxfev=5000)
        except Exception as exc:
            self._status_label.setText(f'Fit failed: {exc}')
            self._status_label.setStyleSheet('color: #FF4444; font-weight: bold;')
            return

        perr = np.sqrt(np.diag(pcov))
        self._popt = popt
        self._perr = perr
        self._fit_cont = cont
        self._init_guess  = list(popt)
        self._bounds_lo   = lo
        self._bounds_hi   = hi

        # Extract z and z_err
        if self._fix_z:
            z_fit   = self._z_guess
            z_err   = 0.0
        else:
            z_fit   = float(popt[0])
            z_err   = float(perr[0])

        self._fit_result = (z_fit, z_err)

        # Draw fit overlay
        fitted_model = direction * model.compile_model(w, *popt) + cont
        if self._fit_artist is not None:
            try: self._fit_artist.remove()
            except Exception: pass
        line, = self._ax.plot(w, fitted_model, color='#FF6B35',
                              lw=1.8, alpha=0.95, zorder=6)
        self._fit_artist = line
        self._canvas.draw_idle()

        # Enable downstream actions
        self._adv_btn.setEnabled(True)
        self._apply_btn.setEnabled(True)
        self._add_abs_btn.setEnabled(True)
        self._status_label.setText('Fit succeeded')
        self._status_label.setStyleSheet('color: #00FF7F; font-weight: bold;')

        # Format results
        txt = self._format_results(
            z_fit, z_err, popt, perr, w, ydata, cont,
            use_weights, direction)
        self._results_te.setPlainText(txt)

    def _format_results(self, z, z_err, popt, perr, wave, ydata, cont,
                        weighted, direction):
        lines = ['─' * 60,
                 'Advanced Fit Results',
                 '─' * 60]
        wt_tag = 'weighted' if weighted else 'unweighted'
        if self._smooth_kernel > 1:
            wt_tag += f'  ⚠ smoothed (kernel={self._smooth_kernel})'
        lines.append(f'Lines  : {", ".join(self._selected_names)}')
        lines.append(f'z      = {z:.8f}  ±  {z_err:.2e}  ({wt_tag})')
        lines.append('')

        rest_arr = np.array(self._selected_rest)
        idx = 0
        if not self._fix_z:
            idx += 1   # skip z parameter
        n_sig = 1 if self._tie_sigma else len(rest_arr)
        sigs = ([popt[idx]] * len(rest_arr) if self._tie_sigma
                else list(popt[idx:idx + n_sig]))
        idx += n_sig
        amps = list(popt[idx:idx + len(rest_arr)])
        sig_errs = ([perr[idx - n_sig]] * len(rest_arr) if self._tie_sigma
                    else list(perr[idx - n_sig:idx]))
        amp_errs = list(perr[idx:idx + len(rest_arr)])

        c_kms = 2.998e5
        ew_list = []
        for i, (name, rest, sig, amp, sig_e, amp_e) in enumerate(
                zip(self._selected_names, rest_arr, sigs, amps,
                    sig_errs, amp_errs)):
            obs_cen = rest * (1.0 + z)
            sig_obs = sig * (1.0 + z)
            fwhm_ang = 2.3548 * sig_obs
            fwhm_kms = fwhm_ang / obs_cen * c_kms
            # Equivalent width (observed frame, trapezoidal)
            gauss_i = direction * MultiGauss._gauss(
                wave / (1.0 + z), sig, amp, rest)
            cont_at_cen = np.interp(obs_cen, wave, cont)
            ew = float(np.trapz(np.abs(gauss_i / (cont_at_cen + 1e-30)), wave))
            ew_list.append(ew)
            lines.append(
                f'Line {i+1}: {name}  (rest {rest:.2f} Å)')
            lines.append(
                f'  σ = {sig:.2f} ± {sig_e:.2f} Å (rest)   '
                f'FWHM = {fwhm_ang:.2f} Å  ({fwhm_kms:.0f} km/s)')
            lines.append(
                f'  Amp = {direction*amp:+.4f} ± {amp_e:.4f}   '
                f'EW ≈ {ew:.3f} Å (obs)')

        if len(ew_list) > 1:
            lines.append('')
            for j in range(1, len(ew_list)):
                r_amp = amps[j] / (amps[0] + 1e-30)
                r_ew  = ew_list[j] / (ew_list[0] + 1e-30)
                lines.append(
                    f'Ratio Amp[{j+1}]/Amp[1] = {r_amp:.4f}   '
                    f'EW[{j+1}]/EW[1] = {r_ew:.4f}')

        lines.append('─' * 60)
        return '\n'.join(lines)

    # ─────────────────────────────────────────────────────────────────
    # Advanced bounds dialog
    # ─────────────────────────────────────────────────────────────────

    def _on_advanced_clicked(self):
        if self._init_guess is None:
            return
        dlg = FittingConstraintDialog(
            self._init_guess, self._bounds_lo, self._bounds_hi,
            self._param_names,
            tie_sigma=self._tie_sigma, fix_z=self._fix_z,
            parent=self)
        if dlg.exec_():
            new_g, new_lo, new_hi, tie, fixz = dlg.get_values()
            if new_g:
                self._init_guess = new_g
                self._bounds_lo  = new_lo
                self._bounds_hi  = new_hi
            changed = (tie != self._tie_sigma) or (fixz != self._fix_z)
            self._tie_sigma = tie
            self._fix_z     = fixz
            if changed:
                # Reset guesses so new model structure is used
                self._init_guess = None

    # ─────────────────────────────────────────────────────────────────
    # Action buttons
    # ─────────────────────────────────────────────────────────────────

    def _on_apply_clicked(self):
        if self._fit_result is None:
            return
        z_new, _ = self._fit_result
        linelist  = self._linelist_combo.currentText()
        pw = self._parent

        # Save last z for Z-revert in main canvas
        try:
            pw.canvas._last_z = float(
                pw.redshift_widget.redshift_input.text().strip())
        except Exception:
            pass

        pw.redshift_widget.set_redshift(z_new)
        pw.redshift_widget.linelist_combo.setCurrentText(linelist)
        color = pw.redshift_widget.color_combo.currentText()
        pw.canvas.set_redshift_data(z_new, linelist, color)

        if hasattr(pw, 'message_box'):
            pw.message_box.on_sent_message(
                f'[Advanced Fit] Applied  z = {z_new:.6f}  '
                f'linelist = {linelist}', '#FF6B35')

    def _on_add_absorber_clicked(self):
        if self._fit_result is None:
            return
        z_new, _ = self._fit_result
        linelist  = self._linelist_combo.currentText()
        pw = self._parent
        if hasattr(pw, 'absorber_manager'):
            pw.absorber_manager.add_absorber(z_new, linelist)
        if hasattr(pw, 'message_box'):
            pw.message_box.on_sent_message(
                f'[Advanced Fit] Added absorber  z = {z_new:.6f}  '
                f'linelist = {linelist}', '#FF6B35')

    def _on_copy_clicked(self):
        txt = self._results_te.toPlainText()
        if txt:
            QApplication.clipboard().setText(txt)

    def _on_export_clicked(self):
        txt = self._results_te.toPlainText()
        if not txt:
            return
        path, ok = QFileDialog.getSaveFileName(
            self, 'Export fit results', '',
            'Text files (*.txt);;All files (*)')
        if ok and path:
            with open(path, 'w') as fh:
                fh.write(txt)
