"""
ImageCanvas — 2D matplotlib image display panel.

Phase 3: whitelight / plain 2D image display with colorbar.
Phase 4: raw data stored, update_data(), clim_changed signal, right-click drag.
Phase 9: left-click → spaxel_locked, left-drag → aperture_drawn.
"""
import numpy as np
import matplotlib as mpl
mpl.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from astropy.visualization import ZScaleInterval

from PyQt5.QtCore import pyqtSignal

# Pixel drag threshold: fewer than this → treated as a click, not a drag
_DRAG_THRESHOLD_PX = 5

# Colors for aperture / spaxel markers (Catppuccin palette)
_AP_COLOR = '#a6e3a1'   # green dashed box
_SP_COLOR = '#f9e2af'   # yellow cross


class ImageCanvas(FigureCanvasQTAgg):
    """
    Matplotlib canvas that displays a single 2-D image.

    Signals
    -------
    spaxel_hovered(int, int, object) — (x, y, spectrum) under cursor      [Phase 7]
    spaxel_locked(int, int, object)  — (x, y, spectrum) on left-click     [Phase 9]
    aperture_drawn(object)           — boolean mask on left-drag           [Phase 9]
    cursor_info(str)                 — formatted status-bar string
    clim_changed(float, float)       — vmin/vmax after right-drag          [Phase 4]
    """

    spaxel_hovered  = pyqtSignal(int, int, object)  # x, y, flux array
    spaxel_locked   = pyqtSignal(int, int, object)  # x, y, flux array
    aperture_drawn  = pyqtSignal(object)             # boolean mask
    region_selected = pyqtSignal(int, int, int, int) # x1, y1, x2, y2
    aperture_picked = pyqtSignal(object)             # artist that was clicked
    cursor_info          = pyqtSignal(str)
    clim_changed         = pyqtSignal(float, float)
    crop_requested          = pyqtSignal()   # right-click context menu → Crop to selection
    sky_region_requested    = pyqtSignal()   # right-click context menu → Set as sky region
    clear_sky_requested     = pyqtSignal()   # right-click context menu → Clear sky region
    clear_overlays_requested = pyqtSignal()  # right-click context menu → Clear region overlays

    def __init__(self, parent=None, dpi=100):
        self._fig = Figure(dpi=dpi, facecolor='#1e1e2e')
        self._fig.patch.set_facecolor('#1e1e2e')
        super().__init__(self._fig)
        self.setParent(parent)

        self._ax = self._fig.add_subplot(111)
        self._style_axes()

        self._im         = None   # current AxesImage
        self._cbar       = None   # current Colorbar
        self._cube       = None   # active IFUCube (set by MainWindow)
        self._wcs        = None   # 2D WCS for sky coordinate lookup
        self.coord_fmt   = 'deg'  # 'deg' = decimal degrees, 'sex' = sexagesimal
        self._data_shape = None   # (ny, nx)
        self._data_raw   = None   # raw (unscaled) float64 array

        # Right-click drag state (contrast/bias)
        self._drag_start       = None   # (x_pixel, y_pixel) at button-press
        self._drag_vmin        = None
        self._drag_vmax        = None
        self._right_drag_moved = False  # True once mouse moves during right-drag

        # Left-click / aperture drag state
        self._left_press       = None   # (canvas_x, canvas_y, data_x, data_y)
        self._ap_rect          = None   # Rectangle patch while dragging
        self._extraction_marks = []     # permanent markers (lines/patches)
        self._last_rect        = None   # (x1, y1, x2, y2) last drag bounds
        self._sky_patch        = None   # cyan rectangle for sky region
        self._preview_artists  = []     # dashed preview circle (cursor-following)
        self.preview_mode      = None   # None or 'circle'
        self.preview_radius    = 3
        self.preview_bg_inner  = None
        self.preview_bg_outer  = None

        # Reference to NavigationToolbar2QT — set by MainWindow after creation.
        # When toolbar is in zoom/pan mode we skip left-click handling.
        self._nav_toolbar      = None

        self._fig.tight_layout(pad=0.5)
        self.mpl_connect('motion_notify_event',  self._on_mouse_move)
        self.mpl_connect('button_press_event',   self._on_button_press)
        self.mpl_connect('button_release_event', self._on_button_release)
        self.mpl_connect('pick_event',           self._on_pick)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def show_image(self, data2d, header=None, cmap='gray', norm=None):
        """
        Display *data2d* (ny × nx).

        Parameters
        ----------
        data2d : np.ndarray, shape (ny, nx)
        header : astropy fits Header or None
        cmap   : str
        norm   : matplotlib Normalize or None  — ZScale used when None
        """
        data2d = np.asarray(data2d, dtype=np.float64)
        self._data_raw   = data2d
        self._data_shape = data2d.shape

        if norm is None:
            norm = _zscale_norm(data2d)

        self._wcs = None
        self._ax.cla()
        self._style_axes()

        if header is not None:
            try:
                from astropy.wcs import WCS
                wcs2d = WCS(header, naxis=2)
                self._wcs = wcs2d
                self._fig.delaxes(self._ax)
                self._ax = self._fig.add_subplot(111, projection=wcs2d)
                self._style_axes()
            except Exception:
                pass

        self._im = self._ax.imshow(
            data2d,
            origin='lower',
            interpolation='nearest',
            cmap=cmap,
            norm=norm,
            aspect='auto',
        )

        if self._cbar is not None:
            try:
                self._cbar.remove()
            except Exception:
                pass

        self._cbar = self._fig.colorbar(self._im, ax=self._ax,
                                        pad=0.02, fraction=0.046)
        self._cbar.ax.yaxis.set_tick_params(color='#cdd6f4', labelcolor='#cdd6f4')
        self._cbar.ax.set_facecolor('#1e1e2e')
        self._cbar.outline.set_edgecolor('#45475a')

        self._apply_label_colors()
        self._fig.tight_layout(pad=0.5)
        self._ax.autoscale(False)   # prevent overlay artists from expanding the view
        self.draw_idle()

    def show_slice(self, i_wave):
        """Display a single wavelength slice from the stored cube (Channel mode)."""
        if self._cube is None or self._im is None:
            return
        data2d = self._cube.flux[i_wave].astype(np.float64)
        self._data_raw = data2d
        self._im.set_data(data2d)
        norm = _zscale_norm(data2d)
        self._im.set_clim(norm.vmin, norm.vmax)
        self.draw_idle()

    def update_data(self, data2d):
        """Replace displayed array without rebuilding axes (scale/method changes).

        Uses ZScale for the new clim so the window position and size are
        unaffected — only pixel values change.
        """
        if self._im is None:
            return
        data2d = np.asarray(data2d, dtype=np.float64)
        self._data_raw = data2d
        self._im.set_data(data2d)
        norm = _zscale_norm(data2d)
        self._im.set_clim(norm.vmin, norm.vmax)
        self.draw_idle()

    def update_clim(self, vmin, vmax):
        """Change display limits without replotting."""
        if self._im is None:
            return
        self._im.set_clim(vmin, vmax)
        self.draw_idle()

    def update_cmap(self, cmap):
        """Change colormap without replotting."""
        if self._im is None:
            return
        self._im.set_cmap(cmap)
        self.draw_idle()

    def draw_sky_region(self, x1, y1, x2, y2):
        """Draw a persistent cyan rectangle marking the sky/background region."""
        self.clear_sky_region()
        self._sky_patch = Rectangle(
            (x1 - 0.5, y1 - 0.5), x2 - x1, y2 - y1,
            linewidth=1.5, edgecolor='#94e2d5', facecolor='#94e2d5',
            alpha=0.15, linestyle='-', zorder=10,
        )
        self._ax.add_patch(self._sky_patch)
        # Label
        self._sky_label = self._ax.text(
            x1, y2, 'sky', color='#94e2d5', fontsize=7, va='bottom', zorder=12)
        self.draw_idle()

    def clear_sky_region(self):
        """Remove the sky region overlay."""
        for attr in ('_sky_patch', '_sky_label'):
            artist = getattr(self, attr, None)
            if artist is not None:
                try:
                    artist.remove()
                except Exception:
                    pass
                setattr(self, attr, None)
        self.draw_idle()

    def clear_extraction_marks(self):
        """Remove all spaxel/aperture markers drawn on the image."""
        for artist in self._extraction_marks:
            try:
                artist.remove()
            except Exception:
                pass
        self._extraction_marks.clear()
        self.draw_idle()

    # ------------------------------------------------------------------
    # Mouse: hover info
    # ------------------------------------------------------------------

    def _on_mouse_move(self, event):
        # Right-click drag → contrast/bias
        if self._drag_start is not None:
            self._handle_drag(event)
            return

        # Left-click drag → draw aperture rectangle
        if self._left_press is not None and event.inaxes:
            px0, py0, dx0, dy0 = self._left_press
            dist = ((event.x - px0) ** 2 + (event.y - py0) ** 2) ** 0.5
            if dist >= _DRAG_THRESHOLD_PX and event.xdata is not None:
                w = event.xdata - dx0
                h = event.ydata - dy0
                if self._ap_rect is None:
                    self._ap_rect = Rectangle(
                        (dx0, dy0), w, h,
                        linewidth=1.2, edgecolor=_AP_COLOR, facecolor='none',
                        linestyle='--', transform=self._ax.transData,
                        zorder=10,
                    )
                    self._ax.add_patch(self._ap_rect)
                else:
                    self._ap_rect.set_xy((dx0, dy0))
                    self._ap_rect.set_width(w)
                    self._ap_rect.set_height(h)
                self.draw_idle()
            return

        if event.inaxes is None or self._im is None:
            self.cursor_info.emit("")
            return

        x, y = event.xdata, event.ydata
        if x is None or y is None:
            self.cursor_info.emit("")
            return

        xi, yi = int(round(x)), int(round(y))
        ny, nx = self._data_shape if self._data_shape else (0, 0)

        parts = [f"x={xi:d}  y={yi:d}"]

        if 0 <= xi < nx and 0 <= yi < ny:
            try:
                val = self._im.get_array()[yi, xi]
                if not np.ma.is_masked(val) and np.isfinite(val):
                    parts.append(f"val={val:.4g}")
            except Exception:
                pass

            # Emit spaxel spectrum for real-time display (Phase 7)
            if self._cube is not None and hasattr(self._cube, 'flux'):
                try:
                    spectrum = self._cube.flux[:, yi, xi]
                    self.spaxel_hovered.emit(xi, yi, spectrum)
                except Exception:
                    pass

        if self._wcs is not None:
            try:
                sky = self._wcs.pixel_to_world(x, y)
                if self.coord_fmt == 'sex':
                    ra_str  = sky.ra.to_string(unit='hour', sep=':', precision=2, pad=True)
                    dec_str = sky.dec.to_string(sep=':', precision=1, alwayssign=True)
                    parts.append(f"RA={ra_str}  Dec={dec_str}")
                else:
                    parts.append(f"RA={sky.ra.deg:.5f}°  Dec={sky.dec.deg:.5f}°")
            except Exception:
                pass

        self.cursor_info.emit("    ".join(parts))

        # Live circular aperture preview — follows cursor
        if self.preview_mode == 'circle' and 0 <= xi < nx and 0 <= yi < ny:
            self.draw_preview_circle(xi, yi, self.preview_radius,
                                     self.preview_bg_inner, self.preview_bg_outer)

    # ------------------------------------------------------------------
    # Mouse: right-click drag → contrast / bias
    # ------------------------------------------------------------------

    def _on_button_press(self, event):
        if event.button == 3 and self._im is not None:  # right click
            self._drag_start = (event.x, event.y)
            vmin, vmax = self._im.get_clim()
            self._drag_vmin = vmin
            self._drag_vmax = vmax
            self._right_drag_moved = False   # reset drag flag

        elif event.button == 1 and self._cube is not None and event.inaxes:
            # Left-click: only handle when nav toolbar is idle
            if self._nav_toolbar is not None and self._nav_toolbar.mode != '':
                return
            self._left_press = (event.x, event.y, event.xdata, event.ydata)

    def _on_button_release(self, event):
        if event.button == 3:
            was_drag = getattr(self, '_right_drag_moved', False)
            self._drag_start = None
            self._right_drag_moved = False
            if not was_drag:
                self._show_context_menu(event)
            return

        elif event.button == 1 and self._left_press is not None:
            px0, py0, dx0, dy0 = self._left_press
            self._left_press = None

            # Remove in-progress rectangle
            if self._ap_rect is not None:
                try:
                    self._ap_rect.remove()
                except Exception:
                    pass
                self._ap_rect = None
                self.draw_idle()

            if event.xdata is None or event.ydata is None:
                return

            dist = ((event.x - px0) ** 2 + (event.y - py0) ** 2) ** 0.5

            if dist < _DRAG_THRESHOLD_PX:
                # --- Single spaxel click ---
                self._handle_spaxel_click(event.xdata, event.ydata)
            else:
                # --- Aperture drag ---
                self._handle_aperture_drag(dx0, dy0, event.xdata, event.ydata)

    def _handle_drag(self, event):
        if self._im is None or event.x is None:
            return

        self._right_drag_moved = True   # mark that a real drag occurred
        dx = event.x  - self._drag_start[0]   # horizontal → bias (shift)
        dy = event.y  - self._drag_start[1]   # vertical   → contrast (stretch)

        vmin0, vmax0 = self._drag_vmin, self._drag_vmax
        span = vmax0 - vmin0 if (vmax0 - vmin0) != 0 else 1.0

        # dy > 0 (drag down) → increase contrast (narrow range)
        # dx > 0 (drag right) → shift brighter
        contrast_factor = 1.0 - dy / 300.0
        contrast_factor = max(0.05, contrast_factor)

        new_span   = span * contrast_factor
        mid        = (vmin0 + vmax0) / 2.0 + dx / 300.0 * span
        new_vmin   = mid - new_span / 2.0
        new_vmax   = mid + new_span / 2.0

        self._im.set_clim(new_vmin, new_vmax)
        self.draw_idle()
        self.clim_changed.emit(float(new_vmin), float(new_vmax))

    def _show_context_menu(self, event):
        """Show right-click context menu at current cursor position."""
        from PyQt5.QtWidgets import QMenu
        from PyQt5.QtGui import QCursor

        has_rect = self._last_rect is not None
        menu = QMenu(self)

        crop_action = menu.addAction("Crop to selection")
        crop_action.setEnabled(has_rect)
        crop_action.setToolTip("Crop the cube to the drawn rectangle  (Ctrl+K)")

        sky_action = menu.addAction("Set as sky region")
        sky_action.setEnabled(has_rect)
        sky_action.setToolTip("Use the drawn rectangle as the sky/background region")

        has_sky = self._sky_patch is not None
        menu.addSeparator()
        clear_sky_action = menu.addAction("Clear sky region")
        clear_sky_action.setEnabled(has_sky)
        clear_sky_action.setToolTip("Remove the current sky region")

        has_overlays = len(self._extraction_marks) > 0
        clear_overlays_action = menu.addAction("Clear region overlays")
        clear_overlays_action.setEnabled(has_overlays)
        clear_overlays_action.setToolTip("Remove all drawn region overlays from this image")

        chosen = menu.exec_(QCursor.pos())
        if chosen == crop_action:
            self.crop_requested.emit()
        elif chosen == sky_action:
            self.sky_region_requested.emit()
        elif chosen == clear_sky_action:
            self.clear_sky_requested.emit()
        elif chosen == clear_overlays_action:
            self.clear_overlays_requested.emit()

    # ------------------------------------------------------------------
    # Left-click / aperture handlers
    # ------------------------------------------------------------------

    def _handle_spaxel_click(self, xdata, ydata):
        """Emit spaxel_locked on left-click — marker drawn by MainWindow via draw_aperture_marker."""
        if self._data_shape is None or self._cube is None:
            return
        ny, nx = self._data_shape
        xi, yi = int(round(xdata)), int(round(ydata))
        if not (0 <= xi < nx and 0 <= yi < ny):
            return

        try:
            spectrum = self._cube.flux[:, yi, xi]
        except Exception:
            return

        self.spaxel_locked.emit(xi, yi, spectrum)

    def draw_aperture_marker(self, cx, cy, shape, size,
                             bg_inner=None, bg_outer=None, color=None):
        """
        Draw a persistent aperture marker on the image after extraction.

        Parameters
        ----------
        cx, cy   : int   — centre pixel
        shape    : str   — 'Circle' or 'Square'
        size     : int   — radius (circle) or half-width (square)
        bg_inner : int or None  — inner radius of background annulus
        bg_outer : int or None  — outer radius of background annulus
        color    : str or None  — marker color; defaults to _SP_COLOR (yellow)
        """
        from matplotlib.patches import Circle, Rectangle as MplRect

        src_color = color if color is not None else _SP_COLOR
        bg_color  = '#74c7ec'    # sky blue — background annulus

        if shape == 'Circle':
            src = Circle((cx, cy), size, fill=False,
                         edgecolor=src_color, lw=2.0, zorder=11, picker=5)
        else:
            s = size + 0.5
            src = MplRect((cx - s, cy - s), 2 * s, 2 * s,
                          fill=False, edgecolor=src_color, lw=2.0, zorder=11, picker=5)
        self._ax.add_patch(src)
        self._extraction_marks.append(src)

        # Centre cross
        mk, = self._ax.plot(cx, cy, '+', color=src_color,
                            ms=10, mew=1.8, zorder=12, picker=5)
        self._extraction_marks.append(mk)

        # Background annulus
        if bg_inner is not None and bg_outer is not None:
            if shape == 'Circle':
                in_patch  = Circle((cx, cy), bg_inner, fill=False,
                                   edgecolor=bg_color, lw=1.4, linestyle='--', zorder=11)
                out_patch = Circle((cx, cy), bg_outer, fill=False,
                                   edgecolor=bg_color, lw=1.4, linestyle='--', zorder=11)
            else:
                si = bg_inner + 0.5
                so = bg_outer + 0.5
                in_patch  = MplRect((cx - si, cy - si), 2*si, 2*si,
                                    fill=False, edgecolor=bg_color, lw=1.4,
                                    linestyle='--', zorder=11)
                out_patch = MplRect((cx - so, cy - so), 2*so, 2*so,
                                    fill=False, edgecolor=bg_color, lw=1.4,
                                    linestyle='--', zorder=11)
            self._ax.add_patch(in_patch)
            self._ax.add_patch(out_patch)
            self._extraction_marks.extend([in_patch, out_patch])

        self.draw_idle()

    def draw_preview_circle(self, cx, cy, radius, bg_inner=None, bg_outer=None):
        """
        Draw a dashed preview circle at (cx, cy) with given radius.

        Replaces any existing preview without an intermediate draw_idle so
        the cursor-following update is a single redraw per mouse move.
        Not added to _extraction_marks — cleared by clear_preview_circle().
        """
        from matplotlib.patches import Circle
        # Remove old artists without redrawing yet
        for a in self._preview_artists:
            try:
                a.remove()
            except Exception:
                pass
        self._preview_artists.clear()

        clr = '#cdd6f4'   # lavender-white
        src = Circle((cx, cy), radius, fill=False,
                     edgecolor=clr, lw=1.4, linestyle='--', zorder=13)
        self._ax.add_patch(src)
        mk, = self._ax.plot(cx, cy, '+', color=clr, ms=8, mew=1.2, zorder=14)
        self._preview_artists.extend([src, mk])
        if bg_inner is not None and bg_outer is not None:
            for r in (bg_inner, bg_outer):
                ring = Circle((cx, cy), r, fill=False,
                              edgecolor='#89dceb', lw=1, linestyle=':', zorder=13)
                self._ax.add_patch(ring)
                self._preview_artists.append(ring)
        self.draw_idle()

    def clear_preview_circle(self):
        """Remove the dashed preview circle from the image."""
        if not self._preview_artists:
            return
        for a in self._preview_artists:
            try:
                a.remove()
            except Exception:
                pass
        self._preview_artists.clear()
        self.draw_idle()

    def draw_region_shape(self, region, wcs, color):
        """
        Draw the actual ds9 region shape on the image canvas.

        Uses ``spatial_mask.draw_region_overlay`` which returns a list of
        matplotlib artists.  Every returned artist is registered in
        ``_extraction_marks`` so it participates in pick events and highlight.

        Parameters
        ----------
        region : dict  — from parse_ds9_regions()
        wcs    : WCS or None
        color  : str

        Returns
        -------
        list of artists
        """
        from rbcodes.GUIs.ifuviewer.processing.spatial_mask import draw_region_overlay

        artists = draw_region_overlay(region, self._ax, wcs,
                                      color=color, draw_label=True)
        for a in artists:
            if a not in self._extraction_marks:
                self._extraction_marks.append(a)
        return artists

    def _on_pick(self, event):
        """Emit aperture_picked when user clicks on a marker artist."""
        if event.artist in self._extraction_marks:
            self.aperture_picked.emit(event.artist)

    def _handle_aperture_drag(self, x0, y0, x1, y1):
        """Lock an aperture rectangle on left-drag."""
        if self._data_shape is None or self._cube is None:
            return
        ny, nx = self._data_shape

        xi0, xi1 = sorted([int(round(x0)), int(round(x1))])
        yi0, yi1 = sorted([int(round(y0)), int(round(y1))])
        xi0 = max(0, xi0); xi1 = min(nx - 1, xi1)
        yi0 = max(0, yi0); yi1 = min(ny - 1, yi1)

        if xi1 < xi0 or yi1 < yi0:
            return

        mask = np.zeros((ny, nx), dtype=bool)
        mask[yi0:yi1 + 1, xi0:xi1 + 1] = True

        # Store bounds for crop / sky region use (exclusive upper bound)
        self._last_rect = (xi0, yi0, xi1 + 1, yi1 + 1)

        # Draw a permanent dashed rectangle
        rect = Rectangle(
            (xi0 - 0.5, yi0 - 0.5),
            xi1 - xi0 + 1, yi1 - yi0 + 1,
            linewidth=2.0, edgecolor=_AP_COLOR, facecolor='none',
            linestyle='--', zorder=11, picker=5,
        )
        self._ax.add_patch(rect)
        self._extraction_marks.append(rect)
        self.draw_idle()

        self.aperture_drawn.emit(mask)
        self.region_selected.emit(xi0, yi0, xi1 + 1, yi1 + 1)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _style_axes(self):
        self._ax.set_facecolor('#1e1e2e')
        self._ax.tick_params(colors='#cdd6f4', labelsize=8, which='both',
                             labelcolor='#cdd6f4')
        for spine in self._ax.spines.values():
            spine.set_edgecolor('#45475a')
        self._apply_label_colors()

    def _apply_label_colors(self):
        c = '#cdd6f4'
        self._ax.xaxis.label.set_color(c)
        self._ax.yaxis.label.set_color(c)
        self._ax.xaxis.label.set_fontsize(8)
        self._ax.yaxis.label.set_fontsize(8)
        self._ax.tick_params(axis='both', colors=c, labelcolor=c, labelsize=8)
        if hasattr(self._ax, 'coords'):
            import astropy
            _new_astropy = int(astropy.__version__.split('.')[0]) >= 6
            for coord in self._ax.coords:
                if _new_astropy:
                    try:
                        coord.set_axislabel(coord.get_axislabel(), color=c)
                    except Exception:
                        pass
                    try:
                        coord.set_ticklabel(color=c)
                    except Exception:
                        pass
                    try:
                        coord.set_ticks(color=c)
                    except Exception:
                        pass
                else:
                    try:
                        coord.axislabels.set_color(c)
                    except Exception:
                        pass
                    try:
                        coord.ticklabels.set_color(c)
                    except Exception:
                        pass
                    try:
                        coord.ticks.set_color(c)
                    except Exception:
                        pass


# ---------------------------------------------------------------------------
# Normalization helper
# ---------------------------------------------------------------------------

def _zscale_norm(data2d):
    finite = data2d[np.isfinite(data2d)]
    if finite.size == 0:
        return mpl.colors.Normalize(vmin=0, vmax=1)
    try:
        vmin, vmax = ZScaleInterval().get_limits(finite)
    except Exception:
        vmin, vmax = np.nanpercentile(finite, [1, 99])
    return mpl.colors.Normalize(vmin=vmin, vmax=vmax)
