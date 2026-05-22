# IFU Viewer — Build Plan

## What It Does (High Level)

1. **Multi-dataset sidebar** — load any number of FITS cubes or 2D images, switch active dataset with one click
2. **Internal IFU viewer** — qfitsview-style: 2D image panel + channel slider + real-time spaxel spectrum on cursor move, fully self-contained, no ds9 required
3. **Image modes** — Channel (single slice), Whitelight, Narrowband, Continuum-subtracted — driven by SpanSelector on spectrum panel with live shaded overlays
4. **Image display controls** — colormap, scale (linear/log/sqrt/square), normalization (ZScale, percentile, manual min/max), right-click drag for contrast/bias
5. **Spectral extraction** — click spaxel, draw aperture, or fetch ds9 regions → extract 1D → display → save or send to other rbcodes GUIs
6. **ds9 bridge (optional)** — push any image to any ds9 frame, WCS lock/match, fetch regions, filter annotation shapes, extract selected or all
7. **Modular IO** — base class + instrument subclasses, easily extensible to new IFUs
8. **rb_spectrum integration** — all extracted spectra become rb_spectrum objects → plugs into specgui, multispecviewer, abstools

---

## Design Principles

- **Modular**: each panel/widget is its own file and class. No file depends on another except through clean interfaces (signals/slots or explicit method calls). Changing one panel must not break others.
- **Reuse rbcodes**: use `rb_spectrum` for all spectral I/O, `spec_advanced2d.ShowAdvanced` logic for image controls, `rb_setline` for line lists, existing GUI dark theme. Do not reinvent what already exists.
- **Simple logic first**: implement the simplest working version of each element before adding options. No speculative features.
- **Test each step**: every phase ends with a runnable, testable state. **Do not proceed to the next phase until the current one is verified working by the user.** If a phase test fails, fix it before moving on — never skip ahead.
- **Not rigid**: use signals/slots for inter-panel communication. Panels do not call each other directly. Adding a new panel = connect its signals, nothing else changes.
- **ds9 is optional**: all core functionality works without ds9. Wrap all pyds9 calls in try/except, disable ds9 controls gracefully if not available.
- **No kcwitools dependency**: copy and generalize relevant functions into `rbcodes/GUIs/ifuviewer/processing/`. Instrument-specific logic lives only in IO subclasses.

---

## File Structure

```
src/rbcodes/GUIs/ifuviewer/
    __init__.py
    main.py                  # MainWindow — top-level QMainWindow, assembles panels
    sidebar.py               # DatasetSidebar — list of loaded files, active selection
    image_panel.py           # ImageCanvas — 2D matplotlib image display
    spectrum_panel.py        # SpectrumCanvas — 1D spectrum + SpanSelectors
    image_controls.py        # ImageControls — scale/norm/colormap/min-max bar
    channel_slider.py        # ChannelSlider — QSlider + wavelength label
    io/
        __init__.py
        base.py              # IFUCube base class — defines the interface
        auto_cube.py         # GenericCube — auto-detect instrument from header (default loader)
        kcwi.py              # KCWICube subclass
        muse.py              # MUSECube subclass
        image2d.py           # FITSImage — wraps a plain 2D FITS image
    processing/
        __init__.py
        cube_collapse.py     # build_whitelight, build_narrowband, build_continuum_sub
        aperture_extract.py  # extract_aperture, extract_variance_weighted
        spatial_mask.py      # parse_ds9_regions, region_to_mask, filter_annotation_shapes
        moment_maps.py       # moment0, moment1, moment2, velocity_array, moment_map
    help.py                  # show_help_dialog() — tabbed help dialog (rb_multispec style)
    HELP.md                  # full user documentation page
    ds9/
        __init__.py
        bridge.py            # DS9Bridge — wraps pyds9, all methods fail gracefully
    TODO.md                  # this file
```

---

## Inter-panel Communication (Signals)

Panels communicate only via Qt signals. No panel imports another panel.

```
DatasetSidebar
    dataset_changed(IFUCube)        → MainWindow → ImageCanvas, SpectrumCanvas, ChannelSlider

ChannelSlider
    channel_changed(int, float)     → ImageCanvas (update slice display)

SpectrumCanvas
    onband_changed(float, float)    → ImageCanvas (recompute narrowband)
    contband_changed(float, float, float, float)  → ImageCanvas (recompute cont-sub)
    spaxel_locked(int, int)         → MainWindow (enable extraction button)

ImageCanvas
    spaxel_hovered(int, int)        → SpectrumCanvas (update spectrum display)
    aperture_drawn(mask)            → MainWindow (enable extraction)

MainWindow
    connects all signals
    owns extraction logic (calls processing/)
    owns ds9 bridge
```

---

## Build Phases

### Phase 1 — IO Module
**Goal**: load a FITS cube, inspect it in a Python shell. No GUI yet.

Files: `io/base.py`, `io/auto_cube.py`, `io/kcwi.py`, `io/image2d.py`

**`IFUCube` base class** (`io/base.py`):
```python
class IFUCube:
    flux    # np.ndarray (n_wave, ny, nx)
    var     # np.ndarray same shape, or None
    wave    # np.ndarray (n_wave,) in Angstroms
    header  # astropy fits Header
    wcs     # astropy WCS, 2D spatial only
    name    # filename stem

    def load(self, path): raise NotImplementedError
    def spatial_header(self): ...   # returns 2D WCS header for ds9/regions
```

**`GenericCube`** (`io/auto_cube.py`):
- Try extension 0 for flux (NAXIS=3), fall back to extension 1
- Build wave from CRVAL3, CDELT3 (or CD3_3), CRPIX3, NAXIS3
- Variance: look for extension named VAR, VARIANCE, STAT, ERR — or None
- Strip 3rd WCS axis from header to get 2D spatial WCS
- Auto-detect instrument from INSTRUME keyword, return subclass

**`FITSImage`** (`io/image2d.py`):
- Wraps a plain 2D FITS image (HST, ground-based, anything)
- Same interface as IFUCube but flux is 2D, wave/var are None

**Test**:
```python
from rbcodes.GUIs.ifuviewer.io.auto_cube import load_fits
cube = load_fits('my_cube.fits')
assert cube.flux.ndim == 3
assert len(cube.wave) == cube.flux.shape[0]
print(cube.wave[0], cube.wave[-1])   # wavelength range
```

---

### Phase 2 — MainWindow Skeleton + Sidebar
**Goal**: window opens, user can load files, they appear in sidebar. No display yet.

Files: `main.py`, `sidebar.py`

**`DatasetSidebar`** (`sidebar.py`):
- `QListWidget` showing loaded datasets
- Each item: filename, type (cube/image), wavelength range if cube
- `[+ Add]` button → `QFileDialog` → `load_fits()` → add to list
- `[- Remove]` button → remove selected
- Click item → emit `dataset_changed(cube)` signal
- Active item shown with filled marker

**`MainWindow`** (`main.py`):
- `QMainWindow` with dark theme (copy palette from zgui/main.py)
- Left: `DatasetSidebar` (fixed width ~200px)
- Center/right: placeholder `QLabel("Load a cube to begin")`
- Menubar: File > Open, File > Quit
- Status bar: shows filename + shape + wave range of active dataset

**Test**: load 2 cubes and 1 image, verify all appear in sidebar, clicking switches active, removing works.

---

### Phase 3 — Image Panel (Whitelight only)
**Goal**: load a cube → whitelight image appears in image panel.

Files: `image_panel.py`, `processing/cube_collapse.py`

**`build_whitelight(flux, wave, wmin=None, wmax=None, method='mean')`** (`processing/cube_collapse.py`):
```python
# simplest version first
mask = np.ones(len(wave), dtype=bool)
if wmin: mask &= wave >= wmin
if wmax: mask &= wave <= wmax
return np.nanmean(flux[mask, :, :], axis=0)   # or nansum, nanmedian
```

**`ImageCanvas`** (`image_panel.py`):
- `FigureCanvasQTAgg` with single `imshow` axes
- `show_image(data2d, header=None)` — displays image, uses WCSAxes if header provided
- Default: whitelight on dataset load
- Colorbar on the side

**Connect**: `DatasetSidebar.dataset_changed` → `MainWindow` → compute whitelight → `ImageCanvas.show_image`

**Test**: load cube → whitelight appears. Load 2D image → image appears. Switch datasets → image updates.

---

### Phase 4 — Image Controls
**Goal**: colormap, scale, normalization, min/max boxes all work on the image panel.

Files: `image_controls.py`

**`ImageControls`** (`image_controls.py`):
- Reuse `ShowAdvanced` logic from `spec_advanced2d.py` directly
- Add colormap `QComboBox`: gray, gray_r, viridis, plasma, inferno, hot, RdBu_r, cubehelix
- Add ZScale option to normalization dropdown using `astropy.visualization.ZScaleInterval`
- Min/Max `QLineEdit` boxes — auto-populate, manual entry with Enter key
- Right-click drag on `ImageCanvas` for contrast/bias:
  ```python
  # track drag delta → adjust imshow vmin/vmax → update min/max boxes
  ```
- All changes call `imshow_obj.set_clim()` + `imshow_obj.set_cmap()` + `canvas.draw_idle()` — no replot

**Test**: load cube → change colormap → image updates. Change normalization → image updates. Type min/max → image updates. Right-click drag → image updates, boxes update.

---

### Phase 5 — Channel Slider
**Goal**: slider scrolls through cube wavelength axis, image updates.

Files: `channel_slider.py`

**`ChannelSlider`** (`channel_slider.py`):
- `QSlider(Qt.Horizontal)` range 0 to N_wave-1
- Label showing: `Channel: 234/1024   λ = 6543.2 Å`
- On value change: emit `channel_changed(i_wave, wavelength)`
- `[◄]` `[►]` buttons for step-by-step

**Image mode button group** (add to `image_controls.py` or `main.py`):
```
[Channel] [Whitelight] [Narrowband] [Cont-sub]
```
- Channel mode: slider active, image shows `flux[i_wave, :, :]`
- Other modes: slider disabled

**Connect**: `ChannelSlider.channel_changed` → `ImageCanvas.show_slice(i_wave)`

**Test**: drag slider → image updates to that wavelength slice. Step buttons work. Mode buttons enable/disable slider.

---

### Phase 6 — Spectrum Panel (Collapsed Spectrum)
**Goal**: collapsed 1D spectrum of active cube displayed below image panel.

Files: `spectrum_panel.py`

**`SpectrumCanvas`** (`spectrum_panel.py`):
- `FigureCanvasQTAgg` with single spectrum axes
- `show_spectrum(wave, flux, label='')` — plots spectrum, stores as `self._line`
- On dataset load: show collapsed spectrum (nanmean over all spaxels)
- `update_spectrum(flux)` — calls `self._line.set_ydata(flux)` + `draw_idle()` — fast, no replot
- x-axis: wavelength in Å, y-axis: flux

**Test**: load cube → collapsed spectrum appears. Switch datasets → spectrum updates. Axis labels correct.

---

### Phase 7 — Real-time Spaxel Spectrum
**Goal**: move cursor over image → spectrum panel updates in real time.

**`ImageCanvas`** additions:
```python
self.canvas.mpl_connect('motion_notify_event', self._on_mouse_move)

def _on_mouse_move(self, event):
    if event.inaxes and self._cube is not None:
        x, y = int(round(event.xdata)), int(round(event.ydata))
        if 0 <= x < nx and 0 <= y < ny:
            spectrum = self._cube.flux[:, y, x]
            self.spaxel_hovered.emit(x, y, spectrum)
```

**`SpectrumCanvas`** additions:
- Slot `on_spaxel_hovered(x, y, spectrum)` → `update_spectrum(spectrum)`
- Status label: `Spaxel: (32, 45)`

**Test**: hover over image → spectrum updates smoothly. Edge spaxels don't crash. Works with variance cube absent.

---

### Phase 8 — SpanSelector + Image Mode Update
**Goal**: drag on spectrum → narrowband image updates in image panel.

**`SpectrumCanvas`** additions:
- `matplotlib.widgets.SpanSelector` for on-band (green, alpha=0.3)
- Second `SpanSelector` for continuum (red, alpha=0.2) — only visible in Cont-sub mode
- On drag end: emit `onband_changed(wmin, wmax)`
- Show shaded spans persistently (they stay after drag ends)

**`processing/cube_collapse.py`** additions:
```python
def build_narrowband(flux, wave, wmin, wmax, method='mean'):
    mask = (wave >= wmin) & (wave <= wmax)
    if method == 'mean':   return np.nanmean(flux[mask], axis=0)
    if method == 'sum':    return np.nansum(flux[mask], axis=0)
    if method == 'median': return np.nanmedian(flux[mask], axis=0)

def build_continuum_sub(flux, wave, wmin, wmax, c1min, c1max, c2min, c2max, method='mean'):
    nb   = build_narrowband(flux, wave, wmin, wmax, method)
    cont = build_narrowband(flux, wave, c1min, c1max, method)
    if c2min is not None:
        cont2 = build_narrowband(flux, wave, c2min, c2max, method)
        cont  = (cont + cont2) / 2.0
    return nb - cont
```

**Connect**: `SpectrumCanvas.onband_changed` → `MainWindow` → `build_narrowband` → `ImageCanvas.show_image`

**Test**: drag span → image updates. Change method dropdown → image updates. Two cont windows work. Shaded spans visible and correct color.

---

### Phase 8 — Status: COMPLETE ✓
Band range dialog (Ctrl+B) implemented:
- `BandRangeDialog` in `main.py` — on-band (green) + continuum 1 (red) + continuum 2 (orange, optional), pre-fills from current SpanSelector extents
- `SpectrumCanvas.apply_band_ranges()` — pushes extents into SpanSelectors visually
- Auto-switches image mode (Narrowband / Cont-sub) and recomputes image on OK

---

### Phase 9 — Status: COMPLETE ✓
Spectral extraction implemented with QFitsView-inspired UI:

**`processing/aperture_extract.py`**
- `make_circular_mask(ny, nx, cx, cy, radius)` — pixels within radius of centre
- `make_annulus_mask(ny, nx, cx, cy, inner_r, outer_r)` — ring background mask
- `extract_with_method(flux, mask, method)` — sum / mean / median dispatch
- `extract_aperture(flux, var, mask)` — summed extraction + propagated error
- `extract_variance_weighted(flux, var, mask)` — 1/σ² weighted extraction per channel
- `subtract_background(spec, flux, bg_mask, method)` — per-pixel bg subtraction

**`utils/rb_optimal_extract.py`** (new utility, no kcwitools/linetools dependency)
- `extract_optimal_weighted(flux, var, mask, method)` — static spatial-profile optimal extraction
  - `method='data'`: weights = collapsed whitelight image within aperture (fast, no fitting)
  - `method='gaussian'`: weights = 2D Gaussian fit to whitelight (smoother, astropy.modeling)
  - Works with any boolean mask shape (circular, rectangular, freehand)
  - Returns `(fl, sig)` — sig=None when var=None; flux renormalized to boxcar-sum scale
  - Gaussian fit falls back to data mode on failure

**`aperture_controls.py`** — `ApertureControls` widget (always visible with cube)
- Three QFitsView-style modes: **Single pixel** | **Circular** | **Circular-Annular**
- **Weighting** dropdown (replaces old checkbox): None | Var-weighted | Optimal (Data) | Optimal (Gaussian)
  - Method (sum/mean/median) hidden when any weighted mode is active — irrelevant then
  - Only Var-weighted disabled when no variance cube; Optimal modes always selectable
- Fields appear/disappear per mode (Radius, BG inner/outer, Method, Weighting)

**`image_panel.py`** additions
- `spaxel_locked(int, int, object)` signal — left-click (< 5px drag threshold)
- `aperture_drawn(object)` signal — left-drag rectangle; dashed green box drawn on image
- `draw_aperture_marker(cx, cy, shape, size, bg_inner, bg_outer)` — persistent yellow circle + cross (+ cyan dashed annulus if bg enabled)
- `clear_extraction_marks()` — removes all markers
- Nav toolbar mode respected (no conflict with zoom/pan)
- Fixed `FutureWarning`: masked array pixel values skipped in coord label rather than formatted

**`spectrum_panel.py`** additions
- `lock_spectrum(wave, flux, err, label, rb_spec)` — colour-cycled overlay line (6-colour Catppuccin cycle)
- `remove_locked_at(idx)` — delete one by index
- `clear_locked()` — remove all
- `locked_spectra` property — list of `(label, rb_spectrum)`
- `set_yscale(scale)` — 'linear' or 'log'
- `autoscale_y()` — fits Y to visible data within current X window (1st–99th percentile)

**`sidebar.py`** additions
- Right-click context menu on any cube item → "Use as variance for '[target]'"
- Shape-mismatch check before assigning; status bar confirmation
- Multiple target cubes → submenu

**`main.py`** additions
- `_on_spaxel_locked` / `_on_aperture_drawn` — dispatch on `extraction_weighting` string
  - extraction labels: `[var]`, `[opt]`, `[opt-g]` for weighted modes
  - warning only for Var-weighted when var=None; Optimal modes run silently without variance
- `_on_use_as_variance` — assigns var_cube.flux → target_cube.var with shape check
- Per-extraction marker tracking (`_extraction_marker_groups`) so Delete removes only the right artists
- **Extract strip** (always visible when cube loaded, disabled when empty):
  - Combo box of all extractions
  - `[Save…]` → `rb_spectrum.write()` via `QFileDialog`
  - `[→ rb_multispec]` → `MultispecSendDialog`
  - `[Delete]` + **Del / Backspace** key → removes selected extraction
  - `[Clear all]` → removes all
- **Spectrum Range & Scale dialog** (`Ctrl+R` / `Ctrl+W`)
- **Image controls**: Reset button moved to far right (addStretch before it)

**`processing/cube_collapse.py`** fixes
- `np.errstate(all='ignore')` around nanmean/nansum/nanmedian — suppresses `Mean of empty slice` on all-NaN channels

**Bug fixes**
- `multispec.py` `LineSelectionDialog`: stored as `self._line_dialog` instead of local variable — fixes non-modal dialog being garbage-collected after nth selection (window disappearing without crash)
- `multispec.py`: `spec.sig.value if hasattr(spec, 'sig') else None` → proper None check (5 occurrences)
- Layout stability: extract strip always-visible (disabled when empty) — no height jitter

---

### Phase 10 — Status: COMPLETE ✓

**`processing/moment_maps.py`** (new file)
- `moment0(flux, wave, wmin, wmax)` — integrated flux (flux·Å)
- `moment1(flux, wave, wmin, wmax, lambda_rest)` — flux-weighted velocity centroid (km/s); NaN where M0 ≤ 0
- `moment2(flux, wave, wmin, wmax, lambda_rest)` — flux-weighted velocity dispersion (km/s); NaN where M0 ≤ 0
- `velocity_array(wave, lambda_rest)` — wavelength → velocity in km/s
- `moment_map(flux, wave, wmin, wmax, order, lambda_rest)` — dispatch wrapper

**`main.py`** additions
- `Analysis > Moment Map…` menu (Ctrl+M) → `MomentMapDialog`
- `MomentMapDialog`: moment order radio buttons (M0/M1/M2), λ window spinboxes, rest wavelength (shown only for M1/M2), `[Use current on-band]` button pre-fills from SpanSelector
- Result replaces current image panel; colormap auto-selected (gray / RdBu_r / plasma)
- Status bar shows moment order, window, rest wavelength, units; persists until next image update
- Clicking any mode button (Whitelight etc.) naturally restores normal display
- `File > Save Current Frame…` (Ctrl+S) → saves `ImageCanvas._data_raw` as FITS with 2D spatial WCS header and `IFUVMODE` keyword tagging what was displayed

### Phase 10 — Moment Maps (original spec)
**Goal**: compute zeroth, first, second moment maps (and higher) from a cube, display them as images, and let the user define the line center + velocity window interactively.

**What a moment map is**:
```
M0(x,y) = Σ F(λ,x,y) · Δλ                              # integrated flux  [flux × Å]
M1(x,y) = Σ F(λ,x,y) · v(λ) / M0                       # velocity centroid [km/s]
M2(x,y) = sqrt( Σ F(λ,x,y) · (v(λ)-M1)² / M0 )        # velocity dispersion [km/s]
```
where `v(λ) = c · (λ - λ_rest) / λ_rest` and the sum runs over the on-band window.

**`processing/moment_maps.py`** (new file):
```python
import numpy as np

C_KMS = 2.998e5  # km/s

def velocity_array(wave, lambda_rest):
    """Convert wavelength array to velocity in km/s relative to lambda_rest."""
    return C_KMS * (wave - lambda_rest) / lambda_rest

def moment0(flux, wave, wmin, wmax):
    """Zeroth moment — integrated flux (flux · Å)."""
    mask = (wave >= wmin) & (wave <= wmax)
    dw   = np.gradient(wave[mask])
    return np.nansum(flux[mask] * dw[:, None, None], axis=0)

def moment1(flux, wave, wmin, wmax, lambda_rest):
    """First moment — flux-weighted velocity centroid (km/s)."""
    mask = (wave >= wmin) & (wave <= wmax)
    vel  = velocity_array(wave[mask], lambda_rest)
    dw   = np.gradient(wave[mask])
    m0   = np.nansum(flux[mask] * dw[:, None, None], axis=0)
    m1   = np.nansum(flux[mask] * vel[:, None, None] * dw[:, None, None], axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(m0 != 0, m1 / m0, np.nan)

def moment2(flux, wave, wmin, wmax, lambda_rest):
    """Second moment — velocity dispersion (km/s)."""
    mask = (wave >= wmin) & (wave <= wmax)
    vel  = velocity_array(wave[mask], lambda_rest)
    dw   = np.gradient(wave[mask])
    m0   = np.nansum(flux[mask] * dw[:, None, None], axis=0)
    m1   = moment1(flux, wave, wmin, wmax, lambda_rest)
    dv2  = (vel[:, None, None] - m1[None]) ** 2
    m2   = np.nansum(flux[mask] * dv2 * dw[:, None, None], axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(m0 > 0, np.sqrt(m2 / m0), np.nan)

def moment_map(flux, wave, wmin, wmax, order, lambda_rest=None):
    """Dispatch to the right function. order in {0, 1, 2}."""
    if order == 0: return moment0(flux, wave, wmin, wmax)
    if order == 1: return moment1(flux, wave, wmin, wmax, lambda_rest)
    if order == 2: return moment2(flux, wave, wmin, wmax, lambda_rest)
    raise ValueError(f"Unsupported moment order {order}")
```

**UI additions**:
- New menu: `Analysis > Compute Moment Map…`  (or a button in the toolbar)
- `MomentDialog` popup:
  - **Moment order**: radio buttons  0 / 1 / 2
  - **Wavelength window**: λ min / λ max (pre-fill from current on-band span)
  - **Rest wavelength** (required for M1/M2): `QDoubleSpinBox` with Å suffix
  - **[Use current on-band]** button — fills λ from SpanSelector
  - OK → compute and display in image panel; image panel tab OR separate result panel TBD

**Display**:
- M0: use existing `ImageCanvas.show_image()` — auto ZScale
- M1 (velocity map): RdBu_r colormap (red=receding, blue=approaching); symmetric clim around 0
- M2 (dispersion map): plasma or inferno colormap; positive-only

**Key design notes**:
- All math is 3–5 numpy lines per moment order — keep in `processing/moment_maps.py`, no logic in GUI files
- `lambda_rest` must be in Å (same units as wave). If user enters in nm, multiply by 10.
- Mask out spaxels where M0 ≤ 0 (or below a user-settable S/N threshold) before computing M1/M2 to avoid noise-dominated velocity pixels
- Optional S/N mask: compute continuum noise in a line-free window, only show M1/M2 where M0 > N_sigma × noise × sqrt(N_channels) × Δλ

**Test**: run on a KCWI Hβ cube at z~0 with λ_rest = 4861.3 Å. M0 should show emission morphology, M1 should show rotation, M2 broadening. Compare with QFitsView result.

---

### Phase 11 — Help System
**Goal**: add a `?` toolbar button and `Help > Help` menu entry that opens a tabbed help dialog, matching the rb_multispec style. Include a dedicated `HELP.md` document for the viewer.

**Pattern**: follow `rbcodes/GUIs/multispecviewer/utils.py → show_help_dialog()` exactly.

**Files to create/modify**:
- `help.py` (new) — `show_help_dialog(parent)` function
- `HELP.md` (new) — full documentation page for the viewer
- `main.py` — add `?` toolbar button + `Help > Help` menu item (shortcut `h`)

**`help.py`**:
- Import pattern identical to multispec utils: lazy PyQt5 imports inside the function so `help.py` has no hard Qt dep at import time
- Tabbed `QDialog` (min 750×520):
  - **Overview** tab — `QTextEdit` (monospace), plain-text description of the viewer: panels, workflow, modes
  - **Keyboard Shortcuts** tab — `QTableWidget` (Key / Action), two columns, alternating rows, key in bold Courier
  - **Extraction** tab — shortcuts + mode descriptions (Single pixel / Circular / Circular-Annular)
  - **Image Modes** tab — Channel / Whitelight / Narrowband / Cont-sub with SpanSelector workflow
  - **Dialogs & Menus** tab — Band Range (Ctrl+B), Spectrum Range (Ctrl+R), Moment Maps, ds9 Bridge
- Bottom row of buttons:
  - `[View README on GitHub]` → opens `docs/GUIs/ifuviewer/HELP.md` on GitHub via `webbrowser.open`
  - `[Open Local README]` → opens `HELP.md` from `importlib.resources` path (if it exists)
  - `[Close]`
- Dark theme: inherit parent palette + same stylesheet as multispec dialog

**Keyboard shortcuts table** (at minimum):
| Key | Action |
|---|---|
| h | Show this help |
| Del / Backspace | Delete selected extraction |
| Ctrl+B | Band range dialog |
| Ctrl+R | Spectrum range & scale dialog |
| Ctrl+O | Open FITS file |
| Ctrl+Q | Quit |

**`main.py` additions**:
```python
from rbcodes.GUIs.ifuviewer.help import show_help_dialog

# In _build_toolbar() — after nav toolbar:
help_btn = QPushButton("?")
help_btn.setToolTip("Show help (h)")
help_btn.setFixedWidth(28)
help_btn.setStyleSheet("QPushButton { font-weight: bold; font-size: 13px; }")
help_btn.clicked.connect(lambda: show_help_dialog(self))
self._mpl_toolbar.addWidget(help_btn)

# In menubar:
help_menu = menubar.addMenu("Help")
help_act = QAction("Help…", self)
help_act.setShortcut("h")
help_act.triggered.connect(lambda: show_help_dialog(self))
help_menu.addAction(help_act)
```

**`HELP.md`** contents (outline):
- Title + one-line description
- Interface overview (sidebar, image panel, spectrum panel, controls)
- Loading files (Open, drag-and-drop if supported)
- Image modes with step-by-step workflow for Narrowband and Cont-sub
- Spectral extraction (all three modes with tips)
- Extractions list — save, send to rb_multispec, delete, clear
- Spectrum Range & Scale dialog
- Band Range dialog
- Moment Maps (Phase 10)
- ds9 Bridge (Phase 11) — optional, what to install
- Full keyboard shortcut reference table
- Troubleshooting (common issues: variance cube not found, pyds9 not available)

**Test**: click `?` button → dialog opens with correct tabs. Click GitHub link → browser opens. Shortcut `h` → same dialog. Dark theme matches rest of GUI.

---

### Phase 12 — ds9 Bridge (Optional)
**Goal**: push images to ds9, fetch regions, extract from regions.

Files: `ds9/bridge.py`

**`DS9Bridge`** (`ds9/bridge.py`):
```python
class DS9Bridge:
    def __init__(self):
        self._ds9 = None   # pyds9.DS9() or None

    def connect(self):
        try:
            import pyds9
            self._ds9 = pyds9.DS9()
            return True
        except:
            return False

    @property
    def available(self):
        return self._ds9 is not None

    def send_image(self, image2d, header, frame=1):
        if not self.available: return
        from astropy.io import fits
        import tempfile, os
        hdu = fits.PrimaryHDU(image2d, header=header)
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            hdu.writeto(f.name, overwrite=True)
            self._ds9.set(f"frame {frame}")
            self._ds9.set(f"fits {f.name}")
        os.unlink(f.name)

    def match_wcs(self):
        if not self.available: return
        self._ds9.set("frame match wcs")
        self._ds9.set("frame lock wcs")

    def get_regions(self, selected_only=False):
        if not self.available: return ""
        cmd = "regions selected sky" if selected_only else "regions sky"
        return self._ds9.get(cmd)
```

**`processing/spatial_mask.py`**:
```python
EXTRACTION_SHAPES = {'circle', 'box', 'ellipse', 'polygon', 'annulus'}
ANNOTATION_SHAPES = {'text', 'compass', 'ruler', 'projection', 'vector', 'point'}

def parse_ds9_regions(reg_text):
    # parse .reg format text → list of dicts with shape, coords, size
    # return only EXTRACTION_SHAPES, skip excluded regions (lines starting with -)

def region_to_mask(region, wcs, shape):
    # sky coords → pixel coords via wcs → rasterize to boolean mask
    # use skimage.draw for circle, rectangle, polygon
```

**Sidebar additions**:
- `[→ ds9 Frame N]` button per dataset
- `[Connect ds9]` button in toolbar
- `[Match WCS]` / `[Lock frames]` buttons

**Extraction panel additions**:
- `[Extract selected region]` → `bridge.get_regions(selected_only=True)`
- `[Extract all regions]` → `bridge.get_regions(selected_only=False)`
- Both filter through `parse_ds9_regions` automatically

**Test**: connect to running ds9 → image appears in correct frame. Match WCS works. Draw region in ds9, fetch selected, spectrum appears. Text/compass regions filtered out.

---

## Reuse Checklist

| Need | Reuse from rbcodes |
|---|---|
| Dark theme palette | `GUIs/zgui/main.py` — copy `set_style()` |
| Image scale/norm controls | `GUIs/zgui/spec_advanced2d.py` — `ShowAdvanced` logic |
| Spectrum container + I/O | `utils/rb_spectrum.py` — `rb_spectrum` |
| Line identification | `IGM/rb_setline.py` — `rb_setline` |
| Spectral rebinning | `IGM/rb_specbin.py` — `rb_specbin` |
| SNR estimation | `utils/compute_SNR_1d.py` — `estimate_snr` |
| Send to specgui | `GUIs/launch_specgui.py` |
| Send to multispecviewer | `GUIs/multispecviewer/` |

---

## External Dependencies

| Package | Purpose | Required? |
|---|---|---|
| `astropy` | FITS I/O, WCS, ZScaleInterval | Yes |
| `numpy` | All array operations | Yes |
| `matplotlib` | Image + spectrum display | Yes |
| `PyQt5` | GUI framework | Yes |
| `skimage.draw` | Rasterize region shapes to masks | Yes (in scipy stack) |
| `pyds9` | ds9 bridge | No — optional |

No kcwitools. No linetools. No other external IFU packages.

---

## ds9 + pyds9 Installation & Verification

ds9 is optional but must be installed and working **before** Phase 10. Verify early.

### Step 1 — Install SAOImageDS9

Download from https://ds9.si.edu/site/Download.html for your platform.

On macOS (recommended via Homebrew):
```bash
brew install --cask saoimage-ds9
```

Verify ds9 is on your PATH:
```bash
ds9 &   # should open the ds9 window
```

### Step 2 — Install pyds9

pyds9 requires XPA, which must also be on your PATH.

```bash
# Install XPA first (needed by pyds9 to communicate with ds9)
brew install xpa        # macOS

# Then install pyds9
pip install pyds9
```

### Step 3 — Verify the full chain works

Run this in a Python shell **while ds9 is already open**:

```python
import pyds9
d = pyds9.DS9()           # connect to running ds9
print(d.get('version'))   # should print ds9 version string

# Send a test image
import numpy as np
from astropy.io import fits
import tempfile, os

data = np.random.random((100, 100))
hdu = fits.PrimaryHDU(data)
with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
    hdu.writeto(f.name, overwrite=True)
    d.set(f'fits {f.name}')
    os.unlink(f.name)

# Verify image appeared in ds9 — you should see a noisy image
print('ds9 + pyds9: OK')
```

**If this fails**, common issues:
- `XPA_METHOD` not set — try `export XPA_METHOD=local` before launching ds9 and Python
- ds9 not on PATH — add `/Applications/SAOImageDS9.app/Contents/MacOS` to PATH
- pyds9 version mismatch — try `pip install pyds9==1.0.3`

**Do not proceed to Phase 12 until this verification script runs without errors.**

---

## Testing Strategy

Each phase has a standalone test. Run from the repo root:

```bash
# Phase 1
python -c "from rbcodes.GUIs.ifuviewer.io.auto_cube import load_fits; c = load_fits('test.fits'); print(c.flux.shape, c.wave[0], c.wave[-1])"

# Phase 2
python -m rbcodes.GUIs.ifuviewer.main   # window opens, load files manually

# Phases 3+ — integrate into main and test each panel interactively
```

Write a small `tests/test_gui/test_ifuviewer_io.py` after Phase 1 covering:
- Load KCWI cube
- Load MUSE cube
- Load plain 2D FITS image
- Load file with missing VAR extension (var=None, no crash)
- Wave array length matches flux axis 0

---

## What NOT to Build (yet)

- RGB false-color (nice to have, not core)
- Optimal extraction (add after variance-weighted works)
- Astrometry correction / WCS tweaking
- Cube arithmetic (add/subtract cubes)
- Any feature not listed above

Add these only when the core phases are complete and tested.
