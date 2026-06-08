"""
ImageControls — colormap, scale, normalization, min/max bar for the image panel.

Phase 4.
"""
import numpy as np
from astropy.visualization import ZScaleInterval

from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QLabel,
                              QComboBox, QLineEdit, QPushButton, QCheckBox)
from PyQt5.QtCore import Qt


COLORMAPS = [
    'gray', 'gray_r',
    'viridis', 'plasma', 'inferno', 'magma',
    'cividis', 'turbo',
    'hot', 'cubehelix', 'gnuplot2',
    'coolwarm', 'RdBu_r', 'seismic',
    'twilight', 'jet',
]

SCALES = ['Linear', 'Log', 'Sqrt', 'Square']

NORMS = ['ZScale', 'MinMax', '99.5%', '99%', '98%', '97%', '95%', 'Manual']


class ImageControls(QWidget):
    """
    Horizontal control bar below the image canvas.

    Layout:  Scale | Norm | Min [____] Max [____] | Colormap | [Reset]
    """

    def __init__(self, canvas, parent=None):
        super().__init__(parent)
        self._canvas   = canvas
        self._updating = False
        self._data_initial = None   # raw data at load time, for Reset

        self._build_ui()
        self._canvas.clim_changed.connect(self._on_clim_changed_external)

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _build_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 3, 6, 3)
        layout.setSpacing(6)

        # Scale
        layout.addWidget(_label('Scale:'))
        self._scale_box = QComboBox()
        self._scale_box.addItems(SCALES)
        self._scale_box.setFixedWidth(80)
        self._scale_box.currentIndexChanged.connect(self._on_scale_changed)
        layout.addWidget(self._scale_box)

        # Normalization
        layout.addWidget(_label('Norm:'))
        self._norm_box = QComboBox()
        self._norm_box.addItems(NORMS)
        self._norm_box.setFixedWidth(80)
        self._norm_box.currentIndexChanged.connect(self._on_norm_changed)
        layout.addWidget(self._norm_box)

        # Min / Max
        layout.addWidget(_label('Min:'))
        self._min_box = QLineEdit()
        self._min_box.setFixedWidth(85)
        self._min_box.setPlaceholderText('min')
        self._min_box.returnPressed.connect(self._on_minmax_entered)
        layout.addWidget(self._min_box)

        layout.addWidget(_label('Max:'))
        self._max_box = QLineEdit()
        self._max_box.setFixedWidth(85)
        self._max_box.setPlaceholderText('max')
        self._max_box.returnPressed.connect(self._on_minmax_entered)
        layout.addWidget(self._max_box)

        # Colormap (rightmost, less frequently changed)
        layout.addWidget(_label('Colormap:'))
        self._cmap_box = QComboBox()
        self._cmap_box.addItems(COLORMAPS)
        self._cmap_box.setFixedWidth(100)
        self._cmap_box.setMaxVisibleItems(10)
        self._cmap_box.currentTextChanged.connect(self._on_cmap_changed)
        layout.addWidget(self._cmap_box)

        # Invert checkbox
        self._invert_cb = QCheckBox("Invert")
        self._invert_cb.setChecked(False)
        self._invert_cb.toggled.connect(self._on_invert_toggled)
        layout.addWidget(self._invert_cb)

        layout.addStretch()

        # Reset button (far right)
        self._reset_btn = QPushButton("Reset")
        self._reset_btn.setFixedWidth(60)
        self._reset_btn.setToolTip("Restore Linear scale, ZScale norm, gray colormap")
        self._reset_btn.clicked.connect(self._on_reset)
        layout.addWidget(self._reset_btn)

    # ------------------------------------------------------------------
    # Called by MainWindow when a new dataset loads
    # ------------------------------------------------------------------

    def reset(self, data2d):
        """Full reset: restore Linear/ZScale/gray defaults (new cube load only)."""
        self._data_initial = data2d
        self._apply_defaults(data2d)

    def refresh(self, data2d):
        """
        Swap in new image data while keeping all current display settings
        (colormap, scale type, normalization).  Use this for mode switches
        (Whitelight → Narrowband, channel changes, etc.) so the user's
        chosen settings are preserved.
        """
        self._data_initial = data2d   # update Reset target to current frame
        scaled = _apply_scale(data2d, self._scale_box.currentIndex())
        # update_data sets ZScale; _apply_norm immediately corrects to user's norm
        self._canvas.update_data(scaled)
        self._apply_norm(scaled)
        # Reapply colormap (update_data doesn't touch it, but be explicit)
        self._canvas.update_cmap(
            _effective_cmap(self._cmap_box.currentText(), self._invert_cb.isChecked()))

    @property
    def current_cmap(self):
        """Effective colormap string currently displayed (includes _r suffix if inverted)."""
        return _effective_cmap(self._cmap_box.currentText(), self._invert_cb.isChecked())

    def set_cmap(self, cmap_name):
        """
        Programmatically set the colormap (e.g. from moment map auto-selection).
        Updates the UI controls and the canvas.
        """
        inverted = cmap_name.endswith('_r')
        base     = cmap_name[:-2] if inverted else cmap_name
        self._updating = True
        if base in COLORMAPS:
            self._cmap_box.setCurrentText(base)
        self._invert_cb.setChecked(inverted)
        self._updating = False
        self._canvas.update_cmap(cmap_name)

    def get_display_state(self):
        """
        Return a dict capturing the current display settings for later restore.
        Keys: scale, norm, cmap, vmin, vmax.
        """
        return {
            'scale': self._scale_box.currentText(),
            'norm':  self._norm_box.currentText(),
            'cmap':  self.current_cmap,
            'vmin':  self._min_box.text(),
            'vmax':  self._max_box.text(),
        }

    def set_display_state(self, state):
        """
        Restore display settings from a dict produced by get_display_state().
        Updates controls and canvas without triggering a full reset.
        """
        self._updating = True
        scale = state.get('scale', 'Linear')
        norm  = state.get('norm',  'ZScale')
        cmap  = state.get('cmap',  'gray')
        vmin  = state.get('vmin',  '')
        vmax  = state.get('vmax',  '')

        if scale in SCALES:
            self._scale_box.setCurrentText(scale)
        if norm in NORMS:
            self._norm_box.setCurrentText(norm)
        self._min_box.setText(vmin)
        self._max_box.setText(vmax)

        inverted = cmap.endswith('_r')
        base     = cmap[:-2] if inverted else cmap
        if base in COLORMAPS:
            self._cmap_box.setCurrentText(base)
        self._invert_cb.setChecked(inverted)
        self._updating = False

        # Apply to canvas
        self._canvas.update_cmap(cmap)
        if norm == 'Manual' and vmin and vmax:
            try:
                self._canvas.update_clim(float(vmin), float(vmax))
            except ValueError:
                pass

    # ------------------------------------------------------------------
    # Slots — controls → canvas
    # ------------------------------------------------------------------

    def _on_cmap_changed(self, name):
        if self._updating:
            return
        self._canvas.update_cmap(_effective_cmap(name, self._invert_cb.isChecked()))

    def _on_invert_toggled(self, checked):
        if self._updating:
            return
        self._canvas.update_cmap(
            _effective_cmap(self._cmap_box.currentText(), checked)
        )

    def _on_scale_changed(self, idx):
        if self._updating or self._data_initial is None:
            return
        scaled = _apply_scale(self._data_initial, idx)
        self._canvas.update_data(scaled)
        self._apply_norm(scaled)

    def _on_norm_changed(self, idx):
        if self._updating or self._data_initial is None:
            return
        # Always apply scale to the original raw data, not the already-scaled display data
        scaled = _apply_scale(self._data_initial,
                              self._scale_box.currentIndex())
        self._apply_norm(scaled, norm_idx=idx)

    def _on_minmax_entered(self):
        try:
            vmin = float(self._min_box.text())
            vmax = float(self._max_box.text())
        except ValueError:
            return
        if vmin > vmax:
            vmin, vmax = vmax, vmin
            self._set_minmax_text(vmin, vmax)
        self._updating = True
        self._norm_box.setCurrentText('Manual')
        self._updating = False
        self._canvas.update_clim(vmin, vmax)

    def _on_clim_changed_external(self, vmin, vmax):
        """Right-click drag on canvas updated clim — sync boxes."""
        self._updating = True
        self._norm_box.setCurrentText('Manual')
        self._set_minmax_text(vmin, vmax)
        self._updating = False

    def _on_reset(self):
        """Restore Linear / ZScale / gray and redisplay original data."""
        if self._data_initial is None:
            return
        self._canvas.update_data(self._data_initial)
        self._apply_defaults(self._data_initial)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _apply_defaults(self, data2d):
        self._updating = True
        self._scale_box.setCurrentIndex(0)       # Linear
        self._norm_box.setCurrentIndex(0)        # ZScale
        self._cmap_box.setCurrentText('gray')
        self._invert_cb.setChecked(False)
        vmin, vmax = _zscale_limits(data2d)
        self._set_minmax_text(vmin, vmax)
        self._updating = False
        self._canvas.update_clim(vmin, vmax)
        self._canvas.update_cmap('gray')

    def _apply_norm(self, scaled, norm_idx=None):
        if norm_idx is None:
            norm_idx = self._norm_box.currentIndex()
        norm_name = NORMS[norm_idx]

        # Always work on a plain float64 array — avoids MaskedArray warnings
        data = np.asarray(scaled, dtype=np.float64)
        finite = data[np.isfinite(data)]
        if finite.size == 0:
            return

        if norm_name == 'ZScale':
            vmin, vmax = _zscale_limits(data)
        elif norm_name == 'MinMax':
            vmin, vmax = float(finite.min()), float(finite.max())
        elif norm_name == 'Manual':
            try:
                vmin = float(self._min_box.text())
                vmax = float(self._max_box.text())
            except ValueError:
                vmin, vmax = _zscale_limits(data)
        else:
            pct = _NORM_PERCENTILES[norm_name]
            vmin, vmax = np.percentile(finite, pct)

        self._canvas.update_clim(vmin, vmax)
        self._set_minmax_text(vmin, vmax)

    def _set_minmax_text(self, vmin, vmax):
        self._min_box.setText(f'{vmin:.4g}')
        self._max_box.setText(f'{vmax:.4g}')


# ---------------------------------------------------------------------------
# Scale helpers
# ---------------------------------------------------------------------------

def _apply_scale(data, scale_idx):
    d = data.copy().astype(np.float64)
    if scale_idx == 0:
        return d
    elif scale_idx == 1:     # Log
        d[d <= 0] = 1e-99
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.log10(d)
    elif scale_idx == 2:     # Sqrt
        out = np.empty_like(d)
        pos = d >= 0
        out[pos]  =  np.sqrt(d[pos])
        out[~pos] = -np.sqrt(-d[~pos])
        return out
    elif scale_idx == 3:     # Square
        return np.sign(d) * d ** 2
    return d


def _zscale_limits(data):
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return 0.0, 1.0
    try:
        vmin, vmax = ZScaleInterval().get_limits(finite)
    except Exception:
        vmin, vmax = np.nanpercentile(finite, [1, 99])
    return float(vmin), float(vmax)


_NORM_PERCENTILES = {
    '99.5%': [0.25, 99.75],
    '99%':   [0.5,  99.5],
    '98%':   [1.0,  99.0],
    '97%':   [1.5,  98.5],
    '95%':   [2.5,  97.5],
}


def _effective_cmap(name, invert):
    """Return the matplotlib colormap name, inverting by toggling the _r suffix."""
    if not invert:
        return name
    if name.endswith('_r'):
        return name[:-2]
    return name + '_r'


def _label(text):
    lbl = QLabel(text)
    lbl.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
    return lbl
