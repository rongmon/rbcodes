"""
SpectrumCanvas — 1D spectrum display panel.

Phase 6: collapsed spectrum.
Phase 7: real-time spaxel hover.
Phase 8: SpanSelector — on-band (green), continuum windows (red/orange).
Phase 9: locked extracted spectra with colour cycle.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

from PyQt5.QtCore import pyqtSignal, Qt

# Colour cycle for locked/extracted spectra (Catppuccin Mocha palette)
_EXTRACT_COLORS = [
    '#f9e2af',  # yellow
    '#94e2d5',  # teal
    '#cba6f7',  # mauve
    '#fab387',  # peach
    '#74c7ec',  # sky
    '#a6e3a1',  # green
]


class SpectrumCanvas(FigureCanvasQTAgg):
    """
    Matplotlib canvas showing a 1-D spectrum.

    Signals
    -------
    onband_changed(float, float)              — on-band wmin, wmax
    contband_changed(float, float, float, float) — c1min,c1max, c2min,c2max
                                                   (c2 = nan,nan when unused)
    spaxel_locked(int, int)                   — click-locked spaxel [Phase 9]
    """

    onband_changed   = pyqtSignal(float, float)
    contband_changed = pyqtSignal(float, float, float, float)
    window_cleared   = pyqtSignal(str)    # 'onband' | 'cont1' | 'cont2'
    spaxel_locked    = pyqtSignal(int, int)
    cursor_info      = pyqtSignal(str)

    def __init__(self, parent=None, dpi=100):
        self._fig = Figure(dpi=dpi, facecolor='#1e1e2e')
        self._fig.patch.set_facecolor('#1e1e2e')
        super().__init__(self._fig)
        self.setParent(parent)

        self._ax   = self._fig.add_subplot(111)
        self._style_axes()

        self._line        = None   # collapsed spectrum
        self._wave        = None
        self._spaxel_line = None   # hover line (Phase 7)

        # SpanSelector widgets (created in show_spectrum)
        self._onband_sel  = None
        self._cont1_sel   = None
        self._cont2_sel   = None

        # Current span extents (nan = not set)
        self._onband_ext  = (np.nan, np.nan)
        self._cont1_ext   = (np.nan, np.nan)
        self._cont2_ext   = (np.nan, np.nan)

        # Active mode — controls which selectors respond
        self._mode = 'Whitelight'

        # Which selector the user is currently drawing into (Cont-sub mode)
        # None | 'onband' | 'cont1' | 'cont2'
        self._drawing_sel = None

        # Legend visibility state
        self._legend_visible = True

        # In-axes status label (created in _build_selectors)
        self._mode_text = None

        # Key-press handler for c band drawing
        self.setFocusPolicy(Qt.StrongFocus)
        self.mpl_connect('key_press_event',    self._on_key_press)
        self.mpl_connect('motion_notify_event', self._on_mouse_move)

        # Locked / extracted spectra  [Phase 9]
        self._locked_lines    = []   # matplotlib Line2D objects
        self._locked_spectra  = []   # list of (label, rb_spectrum or None)
        self._color_idx       = 0

        self._fig.tight_layout(pad=0.8)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def show_spectrum(self, wave, flux, label='Collapsed'):
        """Plot a new spectrum and (re)create SpanSelectors."""
        self._wave = wave
        self._onband_ext = (np.nan, np.nan)
        self._cont1_ext  = (np.nan, np.nan)
        self._cont2_ext  = (np.nan, np.nan)

        self._ax.cla()
        self._style_axes()

        finite = flux[np.isfinite(flux)]
        ymin = np.nanpercentile(finite, 1)  if finite.size else 0
        ymax = np.nanpercentile(finite, 99) if finite.size else 1
        pad  = (ymax - ymin) * 0.1 or 0.1

        self._line, = self._ax.plot(wave, flux, color='#6c7086',
                                    lw=0.8, alpha=0.6, label=label)

        # Hover line — invisible until first mouse move
        self._spaxel_line, = self._ax.plot(wave, flux, color='#f38ba8',
                                           lw=0.9, alpha=0.85,
                                           label='Spaxel', visible=False)

        self._ax.set_xlim(wave[0], wave[-1])
        self._ax.set_ylim(ymin - pad, ymax + pad)
        self._ax.set_xlabel('Wavelength (Å)', fontsize=8)
        self._ax.set_ylabel('Flux', fontsize=8)
        self._rebuild_legend()

        self._build_selectors()
        self.set_mode(self._mode)   # re-apply current mode's selector state

        self._fig.tight_layout(pad=0.8)
        self.draw_idle()

    def update_spectrum(self, flux):
        """Fast y-data swap — no replot."""
        if self._line is None:
            return
        self._line.set_ydata(flux)
        self.draw_idle()

    def on_spaxel_hovered(self, x, y, spectrum):
        """Overlay spaxel spectrum in real time."""
        if self._spaxel_line is None or self._wave is None:
            return
        self._spaxel_line.set_ydata(spectrum)
        self._spaxel_line.set_visible(True)
        self._spaxel_line.set_label(f'Spaxel ({x},{y})')
        self._rebuild_legend()
        self.draw_idle()

    def set_mode(self, mode):
        """
        Activate the right SpanSelectors for *mode*.

        Narrowband  — on-band selector auto-active (only one, no conflict)
        Cont-sub    — on-band auto-active (mouse drag = emission); c key cycles to cont selectors
        Others      — all selectors inactive
        """
        self._mode = mode
        if self._onband_sel is None:
            return

        in_narrowband = mode == 'Narrowband'
        in_contsub    = mode == 'Cont-sub'

        # Neutral drawing state on mode switch — on-band is active by selector
        # logic but no key has been explicitly pressed yet (_drawing_sel = None).
        # Pressing c enters explicit on-band mode; pressing c again clears it.
        self._drawing_sel = None

        self._onband_sel.set_active(in_narrowband or in_contsub)
        self._cont1_sel.set_active(False)
        self._cont2_sel.set_active(False)

        # Show/hide existing shaded spans
        self._onband_sel.set_visible(in_narrowband or in_contsub)
        self._cont1_sel.set_visible(in_contsub)
        self._cont2_sel.set_visible(in_contsub)

        self._update_drawing_status()
        self.draw_idle()

    def enterEvent(self, event):
        """Grab keyboard focus when mouse enters the spectrum panel."""
        self.setFocus()
        super().enterEvent(event)

    def lock_spectrum(self, wave, flux, err, label, rb_spec=None):
        """
        Add a permanently visible extracted spectrum to the axes.

        Parameters
        ----------
        wave    : ndarray
        flux    : ndarray
        err     : ndarray or None  (not plotted yet, stored in rb_spec)
        label   : str              display label
        rb_spec : rb_spectrum or None  — stored for save/export
        """
        if self._wave is None:
            return

        color = _EXTRACT_COLORS[self._color_idx % len(_EXTRACT_COLORS)]
        self._color_idx += 1

        line, = self._ax.plot(wave, flux, color=color, lw=0.9,
                              alpha=0.9, label=label, zorder=5)
        self._locked_lines.append(line)
        self._locked_spectra.append((label, rb_spec))

        self._rebuild_legend()
        self.draw_idle()

    def remove_locked_at(self, idx):
        """Remove a single extracted spectrum by index."""
        if idx < 0 or idx >= len(self._locked_lines):
            return
        try:
            self._locked_lines[idx].remove()
        except Exception:
            pass
        del self._locked_lines[idx]
        del self._locked_spectra[idx]

        self._rebuild_legend()
        self.draw_idle()

    def clear_locked(self):
        """Remove all locked/extracted spectra from the axes."""
        for line in self._locked_lines:
            try:
                line.remove()
            except Exception:
                pass
        self._locked_lines.clear()
        self._locked_spectra.clear()
        self._color_idx = 0

        self._rebuild_legend()
        self.draw_idle()

    @property
    def locked_spectra(self):
        """List of (label, rb_spectrum) for all locked extractions."""
        return list(self._locked_spectra)

    def apply_band_ranges(self, onband, cont1, cont2):
        """
        Programmatically set band extents and update the SpanSelector visuals.

        Parameters
        ----------
        onband : (float, float)  — on-band wmin/wmax
        cont1  : (float, float)  — continuum 1 wmin/wmax  (nan, nan = clear)
        cont2  : (float, float)  — continuum 2 wmin/wmax  (nan, nan = clear)
        """
        self._onband_ext = tuple(onband)
        self._cont1_ext  = tuple(cont1)
        self._cont2_ext  = tuple(cont2)

        if self._onband_sel is not None and not np.isnan(onband[0]):
            try:
                self._onband_sel.extents = (float(onband[0]), float(onband[1]))
            except Exception:
                pass

        if self._cont1_sel is not None and not np.isnan(cont1[0]):
            try:
                self._cont1_sel.extents = (float(cont1[0]), float(cont1[1]))
            except Exception:
                pass

        if self._cont2_sel is not None and not np.isnan(cont2[0]):
            try:
                self._cont2_sel.extents = (float(cont2[0]), float(cont2[1]))
            except Exception:
                pass

        self.draw_idle()

    def set_yscale(self, scale):
        """Set Y axis scale: 'linear' or 'log'."""
        try:
            self._ax.set_yscale(scale)
            self.draw_idle()
        except Exception:
            pass

    def autoscale_y(self):
        """
        Rescale the Y axis to fit all visible data within the current X limits.
        Accounts for collapsed line, hover line, and any locked spectra.
        """
        xmin, xmax = self._ax.get_xlim()
        yvals = []

        def _collect(line):
            if line is None or not line.get_visible():
                return
            xd = line.get_xdata()
            yd = line.get_ydata()
            if xd is None or yd is None or len(xd) == 0:
                return
            xd = np.asarray(xd); yd = np.asarray(yd)
            mask = (xd >= xmin) & (xd <= xmax) & np.isfinite(yd)
            if mask.any():
                yvals.extend(yd[mask].tolist())

        _collect(self._line)
        for line in self._locked_lines:
            _collect(line)

        if not yvals:
            return
        lo = np.percentile(yvals, 1)
        hi = np.percentile(yvals, 99)
        pad = (hi - lo) * 0.08 or 0.1
        self._ax.set_ylim(lo - pad, hi + pad)
        self.draw_idle()

    def set_wave_limits(self, wmin, wmax):
        """Zoom the spectrum x-axis to [wmin, wmax]."""
        self._ax.set_xlim(wmin, wmax)
        self.draw_idle()

    def reset_wave_limits(self):
        """Restore full wavelength range."""
        if self._wave is not None:
            self._ax.set_xlim(self._wave[0], self._wave[-1])
            self.draw_idle()

    def clear(self):
        self._ax.cla()
        self._style_axes()
        self._line        = None
        self._wave        = None
        self._spaxel_line = None
        self._onband_sel  = None
        self._cont1_sel   = None
        self._cont2_sel   = None
        self._locked_lines.clear()
        self._locked_spectra.clear()
        self._color_idx = 0
        self.draw_idle()

    # ------------------------------------------------------------------
    # SpanSelector construction
    # ------------------------------------------------------------------

    def _build_selectors(self):
        """(Re)create all three SpanSelectors on the current axes."""
        common = dict(
            direction='horizontal',
            useblit=True,
            interactive=True,
            drag_from_anywhere=False,   # only edge handles resize; no accidental whole-span moves
        )

        self._onband_sel = SpanSelector(
            self._ax, self._on_onband_select,
            props=dict(facecolor='#a6e3a1', alpha=0.25, edgecolor='#a6e3a1'),
            **common,
        )

        self._cont1_sel = SpanSelector(
            self._ax, self._on_cont1_select,
            props=dict(facecolor='#f38ba8', alpha=0.20, edgecolor='#f38ba8'),
            **common,
        )

        self._cont2_sel = SpanSelector(
            self._ax, self._on_cont2_select,
            props=dict(facecolor='#fab387', alpha=0.20, edgecolor='#fab387'),
            **common,
        )

        # Start all inactive — set_mode() will enable the right ones
        for sel in (self._onband_sel, self._cont1_sel, self._cont2_sel):
            sel.set_active(False)
            sel.set_visible(False)

        # Status annotation — shows which band the user is currently drawing
        self._mode_text = self._ax.text(
            0.01, 0.97, '',
            transform=self._ax.transAxes,
            fontsize=7, color='#cdd6f4', va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#313244', alpha=0.85),
            visible=False, zorder=10,
        )

    # ------------------------------------------------------------------
    # SpanSelector callbacks
    # ------------------------------------------------------------------

    def _on_onband_select(self, wmin, wmax):
        if wmax <= wmin:
            return
        self._onband_ext = (wmin, wmax)
        self.onband_changed.emit(wmin, wmax)

    def _on_cont1_select(self, wmin, wmax):
        if wmax <= wmin:
            return
        self._cont1_ext = (wmin, wmax)
        self._emit_contband()

    def _on_cont2_select(self, wmin, wmax):
        if wmax <= wmin:
            return
        self._cont2_ext = (wmin, wmax)
        self._emit_contband()

    def _emit_contband(self):
        c1min, c1max = self._cont1_ext
        c2min, c2max = self._cont2_ext
        self.contband_changed.emit(
            float(c1min), float(c1max),
            float(c2min), float(c2max),
        )

    # ------------------------------------------------------------------
    # Key-driven band drawing (Cont-sub mode)
    # ------------------------------------------------------------------

    def _on_key_press(self, event):
        """
        Navigation keys (all modes):
          x / X   — cursor x → new xmin / xmax
          t / b   — cursor y → new ymax (top) / ymin (bottom)
          r       — reset full wavelength range
          a       — autoscale Y to visible data
          [ / ]   — pan left / right 10%
          - / =   — zoom out / in X by 20%
          l       — toggle legend

        Cont-sub band keys (Cont-sub mode only):
          c       — on-band window  (press twice to clear)
          1       — cont window 1   (press twice to clear)
          2       — cont window 2   (press twice to clear)
          Esc     — back to neutral (keep all spans)
        """
        key = event.key

        # ---- Navigation — active in all modes ----
        if key == 'r':
            self.reset_wave_limits()
            return
        if key == 'a':
            self.autoscale_y()
            return
        if key == 'l':
            self._toggle_legend()
            return
        if key in ('[', ']'):
            self._pan_x(-0.1 if key == '[' else 0.1)
            return
        if key in ('-', '='):
            self._zoom_x(1.2 if key == '-' else 1.0 / 1.2)
            return
        if key in ('x', 'X') and event.xdata is not None:
            xmin, xmax = self._ax.get_xlim()
            if key == 'x':
                self._ax.set_xlim(event.xdata, xmax)
            else:
                self._ax.set_xlim(xmin, event.xdata)
            self.draw_idle()
            return
        if key in ('t', 'b') and event.ydata is not None:
            ymin, ymax = self._ax.get_ylim()
            if key == 't':
                self._ax.set_ylim(ymin, event.ydata)
            else:
                self._ax.set_ylim(event.ydata, ymax)
            self.draw_idle()
            return

        # ---- Cont-sub band keys — Cont-sub mode only ----
        if self._mode != 'Cont-sub':
            return

        if key == 'c':
            if self._drawing_sel == 'onband':
                self._clear_window('onband')
            else:
                self._set_drawing_selector('onband')
        elif key == '1':
            if self._drawing_sel == 'cont1':
                self._clear_window('cont1')
            else:
                self._set_drawing_selector('cont1')
        elif key == '2':
            if self._drawing_sel == 'cont2':
                self._clear_window('cont2')
            else:
                self._set_drawing_selector('cont2')
        elif key == 'escape':
            self._set_drawing_selector(None)

    def _clear_window(self, name):
        """
        Clear one band window, rebuild selectors to remove the span visually,
        restore remaining spans, and emit window_cleared so main.py updates image.
        """
        if name == 'onband':
            self._onband_ext = (np.nan, np.nan)
        elif name == 'cont1':
            self._cont1_ext = (np.nan, np.nan)
        elif name == 'cont2':
            self._cont2_ext = (np.nan, np.nan)

        # Snapshot what remains after clearing
        onband = self._onband_ext
        cont1  = self._cont1_ext
        cont2  = self._cont2_ext

        # Remove mode_text before rebuilding to avoid duplication
        if self._mode_text is not None:
            try:
                self._mode_text.remove()
            except Exception:
                pass
            self._mode_text = None

        # Disconnect old selectors so their event handlers don't linger
        for sel in (self._onband_sel, self._cont1_sel, self._cont2_sel):
            if sel is not None:
                try:
                    sel.disconnect_events()
                except Exception:
                    pass

        # Rebuild fresh selectors (no visual spans)
        self._build_selectors()

        # Restore non-cleared spans visually
        self.apply_band_ranges(onband, cont1, cont2)

        # Re-apply mode state (activates right selectors, resets _drawing_sel)
        self.set_mode(self._mode)

        # Notify main.py to update the image
        self.window_cleared.emit(name)
        self.draw_idle()

    def _set_drawing_selector(self, name):
        """Switch which selector is active. On-band is the default; cont selectors replace it temporarily."""
        self._drawing_sel = name
        if self._onband_sel is None:
            return
        # On-band is active whenever we are NOT drawing a cont selector
        self._onband_sel.set_active(name not in ('cont1', 'cont2'))
        self._cont1_sel.set_active(name == 'cont1')
        self._cont2_sel.set_active(name == 'cont2')
        self._update_drawing_status()
        self.draw_idle()

    def _update_drawing_status(self):
        """Show a mode label whenever a key has been explicitly pressed in Cont-sub."""
        if self._mode_text is None:
            return
        if self._mode != 'Cont-sub' or self._drawing_sel is None:
            self._mode_text.set_visible(False)
            return
        labels = {
            'onband': 'On-band  [c×2=clear | 1=Cont1 | 2=Cont2 | Esc]',
            'cont1':  'Cont 1   [1×2=clear | c=On-band | 2=Cont2 | Esc]',
            'cont2':  'Cont 2   [2×2=clear | c=On-band | 1=Cont1 | Esc]',
        }
        self._mode_text.set_text(labels.get(self._drawing_sel, ''))
        self._mode_text.set_visible(True)

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _rebuild_legend(self):
        """Rebuild legend and apply current visibility state."""
        leg = self._ax.legend(fontsize=7, loc='upper right',
                              facecolor='#313244', labelcolor='#cdd6f4',
                              edgecolor='#45475a')
        if leg is not None:
            leg.set_visible(self._legend_visible)

    def _toggle_legend(self):
        self._legend_visible = not self._legend_visible
        leg = self._ax.get_legend()
        if leg is not None:
            leg.set_visible(self._legend_visible)
            self.draw_idle()

    def _pan_x(self, fraction):
        """Pan X axis by *fraction* of current range (negative = left)."""
        xmin, xmax = self._ax.get_xlim()
        delta = (xmax - xmin) * fraction
        self._ax.set_xlim(xmin + delta, xmax + delta)
        self.draw_idle()

    def _zoom_x(self, factor):
        """Zoom X axis by *factor* around its centre (>1 = zoom out)."""
        xmin, xmax = self._ax.get_xlim()
        mid = (xmin + xmax) / 2.0
        half = (xmax - xmin) / 2.0 * factor
        self._ax.set_xlim(mid - half, mid + half)
        self.draw_idle()

    def _on_mouse_move(self, event):
        if event.inaxes != self._ax or event.xdata is None or event.ydata is None:
            self.cursor_info.emit('')
            return
        self.cursor_info.emit(f"λ={event.xdata:.2f} Å    flux={event.ydata:.4g}")

    def _style_axes(self):
        c = '#cdd6f4'
        self._ax.set_facecolor('#1e1e2e')
        self._ax.tick_params(colors=c, labelsize=7, which='both', labelcolor=c)
        self._ax.xaxis.label.set_color(c)
        self._ax.yaxis.label.set_color(c)
        for spine in self._ax.spines.values():
            spine.set_edgecolor('#45475a')
