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

### Phase 11 — Status: COMPLETE ✓

**ds9 Bridge + Region I/O + Per-Dataset State Persistence**

#### What was implemented

**`ds9/bridge.py`**
- `DS9Bridge.get_regions()` rewritten: tries XPA `get` commands first (`regions wcs fk5 degrees`, `regions fk5`, `regions wcs`) then falls back to `regions save {tmpfile}` (no coord args — those are invalid in XPA `set`)
- All methods fail gracefully; GUI works fully without pyds9

**`processing/spatial_mask.py`**
- `parse_ds9_regions()` — uses `astropy-regions` package first, fallback parser for physical/pixel coords
- Fixed: `if result:` instead of `if result is not None:` — empty list from regions package now falls through to fallback
- `_reg_shape_name` — ordered list of tuples, `circleannulus` checked before `circle`
- `_draw_fallback_marker()` — rewritten to draw actual `Circle`/`Ellipse`/`Rectangle`/`Polygon` matplotlib patches from region args (not just crosses)
- `draw_region_overlay()` — dispatches to `astropy-regions` artists or fallback patches

**`image_panel.py`**
- `draw_region_shape(region, wcs, color)` — draws imported ds9 region as real shape on image

**`main.py` — Region import**
- `_extract_from_region_text()` — `FITSImage` path: overlay only (no spectrum extraction); cube path: uses `extract_optimal_weighted(flux, var, mask, method='data')` directly (same pattern as rest of GUI)
- `_on_load_region_file()` / `_on_ds9_import_regions()` — both work for cubes and 2D images
- `marker_info` dicts persisted per extraction: `{'type': 'spaxel'|'circle'|'rect'|'region', ...params}`

**`main.py` — Per-dataset state persistence**

Every time you switch datasets, the current state is saved and restored when you return:

| Field saved | What it covers |
|---|---|
| `image_data` | Whatever 2D array is on screen (whitelight / narrowband / cont-sub / moment map / SNR map) |
| `image_norm` | Colorscale (ZScale / percentile / manual) |
| `slider_mode` | Which channel-slider button was active (Whitelight / Narrowband / Cont-sub / Channel) |
| `onband_ext`, `cont1_ext`, `cont2_ext` | Green/red shaded band spans on the spectrum |
| `spec_mode` | Spectrum SpanSelector mode |
| `spec_ylim` | Y-axis limits |
| `locked` | All extracted spectra (wave, flux, color, label, rb_spec) |
| `color_idx` | Position in the color cycle |
| `sky_mask`, `sky_rect` | Sky region boolean mask + pixel bounds |
| `moment_params` | All MomentMapDialog last-used values (per dataset) |

On restore:
- Image redrawn with saved norm (no recompute)
- Channel slider button restored (setChecked, does not trigger image rebuild)
- Spectrum bands redrawn via `apply_band_ranges()`
- All locked spectra re-plotted with exact saved colors
- All image markers redrawn (`spaxel` cross / `circle` aperture / `rect` / `region` shape)
- Sky region cyan box redrawn
- Y-axis autoscaled to fit restored spectra
- MomentMapDialog class-level `_last_*` attributes swapped to match this dataset

**`spectrum_panel.py`**
- `lock_spectrum(..., color=None)` — `color` param allows restoring exact saved colors
- `_rebuild_legend()` — guards with `get_legend_handles_labels()` before calling `legend()` (eliminates "No artists with labels" warning)

**`processing/cube_collapse.py`**
- Wrapped `nanmean/nansum/nanmedian` with `warnings.catch_warnings()` + `warnings.simplefilter('ignore', RuntimeWarning)` — eliminates "Mean of empty slice" warning

#### Known limitation
Moment maps: the image pixel data restores correctly (the moment map IS shown). However the channel slider button says "Whitelight" — if you click any mode button it will rebuild from the cube. The moment map parameters (wavelength range, continuum, SNR) are saved per-dataset in `moment_params` so reopening the dialog shows the correct values.

---

### Phase 11 + 12 — Testing Checklist — Status: ALL PASS ✓

Run the GUI:
```bash
python -m rbcodes.GUIs.ifuviewer.main
```

#### A — Basic launch ✓
#### B — Image modes ✓
#### C — Spectral extraction ✓
#### D — Multiple extractions ✓
#### E — Per-dataset state: basic ✓
#### F — Per-dataset state: image mode ✓
#### G — Per-dataset state: moment maps ✓
#### H — Sky region ✓
#### I — ds9 region import ✓
#### J — Region file loading ✓
#### K — Cross-instrument region (HST → IFU) ✓
#### L — Warning/error sanity ✓

**Bugs fixed during testing:**
- `subtract_linear_continuum`: replaced `np.errstate` with `warnings.catch_warnings()` to suppress `RuntimeWarning: Mean of empty slice`
- `rb_spectrum._rb_parse_multi_extension_fits`: rewrote to use name-based extension lookup (`FLUX`, `WAVELENGTH`, `ERROR`/`ERR`) instead of positional — fixes misread when no error array present
- `io_manager.py`: added `warnings.warn()` to terminal when 5% flux error is assumed (no error array in file)
- `_restore_dataset_state`: passes `cube.spatial_header()` to `show_image` so RA/Dec WCS is restored on dataset switch
- `DS9Bridge.get_regions`: sets `regions system wcs / sky fk5 / skyformat degrees` before file-save fallback — fixes cross-instrument region import (HST → IFU)
- `_extract_from_region_text` (2D image path): region overlays saved in `_image_overlays` and persisted in `_dataset_states` — overlays survive dataset switch
- `_on_clear_extractions`: added `_from_switch` flag; manual "Clear all" now saves cleared state so overlays don't reappear on next switch
- Extract strip for 2D images: shows only "Clear all" button (hidden when empty, enabled after overlays drawn); spectrum buttons hidden
- `image_panel.py`: added `clear_overlays_requested` signal → right-click "Clear region overlays" menu item

---

### Phase 11 — ds9 Bridge + Region I/O + Batch Extraction (original spec)
**Goal**: send images to ds9, import/export region files, batch-extract spectra from regions, highlight apertures interactively, add RA/Dec provenance to extracted spectra.

All ds9 functionality is optional — the GUI must work fully without pyds9 installed.

---

#### 11.1 — `rb_spectrum` RA/Dec provenance (`utils/rb_spectrum.py`)

Add RA/Dec and extraction metadata to `rb_spectrum.meta` and persist as explicit FITS header keywords.

**`rb_write_fits` additions** — after the existing `METADATA` keyword:
```python
if 'ra' in self.meta:
    hdu_list[0].header['RA']      = (self.meta['ra'],  '[deg] J2000 right ascension')
    hdu_list[0].header['DEC']     = (self.meta['dec'], '[deg] J2000 declination')
    hdu_list[0].header['EQUINOX'] = 2000.0
    hdu_list[0].header['RADESYS'] = 'ICRS'
if 'instrume' in self.meta:
    hdu_list[0].header['INSTRUME'] = self.meta['instrume']
if 'object' in self.meta:
    hdu_list[0].header['OBJECT'] = self.meta['object']
```

**Provenance fields written into `meta` by IFU extractor** (`main.py`) after every extraction:
```python
meta = {}
if cube.wcs is not None:
    ra, dec = cube.wcs.all_pix2world([[cx, cy]], 0)[0]
    meta['ra']  = float(ra)
    meta['dec'] = float(dec)
meta['extr_x']   = cx          # pixel center x
meta['extr_y']   = cy          # pixel center y
meta['extr_rad'] = radius      # arcsec for circular, 'rect WxH' for rect, 1 for single
meta['object']   = cube.header.get('OBJECT', '')
meta['instrume'] = cube.header.get('INSTRUME', '')
rb_spec.meta.update(meta)
```

Apply this in `_on_spaxel_locked`, `_on_aperture_drawn`, and the new batch extractor.

---

#### 11.2 — `processing/spatial_mask.py` (new file)

Region parsing and rasterization. Works without ds9 — used by both live ds9 import and disk file loading.

```python
EXTRACTION_SHAPES = {'circle', 'box', 'ellipse', 'polygon', 'annulus'}
ANNOTATION_SHAPES = {'text', 'compass', 'ruler', 'projection', 'vector', 'point'}

def parse_ds9_regions(reg_text):
    """
    Parse ds9 region format text → list of region dicts.
    Each dict: {'shape': str, 'coords': [...], 'system': 'fk5'|'image', 'name': str}
    Skips annotation shapes and excluded regions (lines starting with -).
    Supports: circle, box, ellipse, polygon, annulus.
    """

def region_to_mask(region, wcs, ny, nx):
    """
    Convert one region dict → boolean mask (ny, nx).
    Sky coords → pixel via wcs.all_world2pix(); rasterize with skimage.draw.
    circle   → draw.disk(center, radius)
    box      → draw.rectangle(top_left, extent)
    polygon  → draw.polygon(rows, cols)
    ellipse  → draw.ellipse(center, r_radius, c_radius, rotation)
    annulus  → outer disk & ~inner disk
    Clips to image bounds; returns all-False mask (not an error) if region outside FOV.
    """

def region_center_sky(region, wcs):
    """
    Return (ra_deg, dec_deg) of region center.
    Used for IAU filename generation and meta['ra'/'dec'].
    """

def iau_name(ra_deg, dec_deg, prefix='spec_'):
    """
    Build IAU-format filename: prefix + J101056.2+213042.fits
    Uses astropy.coordinates.SkyCoord.to_string('hmsdms').
    """
```

---

#### 11.3 — `ds9/bridge.py` (new file)

All methods fail gracefully if pyds9 is unavailable or ds9 is not running.

```python
class DS9Bridge:
    def __init__(self):
        self._ds9 = None

    # --- connection ---
    def connect(self):
        """
        Try to connect to a running ds9 instance.
        Returns (True, '') on success.
        Returns (False, reason) where reason is:
          'no_pyds9'  — pyds9 not installed
          'no_ds9'    — pyds9 installed but no ds9 running
        """
        try:
            import pyds9
        except ImportError:
            return False, 'no_pyds9'
        try:
            self._ds9 = pyds9.DS9()
            return True, ''
        except Exception:
            return False, 'no_ds9'

    def disconnect(self):
        self._ds9 = None

    @property
    def available(self):
        return self._ds9 is not None

    # --- send image ---
    def send_image(self, image2d, header, frame=1):
        """Write 2D image to temp FITS, send to ds9 frame N, delete temp file."""
        if not self.available: return
        from astropy.io import fits
        import tempfile, os
        hdu = fits.PrimaryHDU(image2d, header=header)
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            hdu.writeto(f.name, overwrite=True)
            self._ds9.set(f'frame {frame}')
            self._ds9.set(f'fits {f.name}')
        os.unlink(f.name)

    # --- WCS ---
    def match_wcs(self):
        if not self.available: return
        self._ds9.set('frame match wcs')
        self._ds9.set('frame lock wcs')

    # --- regions ---
    def get_regions(self, selected_only=False):
        """Return region text string from ds9 (sky coordinates)."""
        if not self.available: return ''
        cmd = 'regions selected sky' if selected_only else 'regions sky'
        try:
            return self._ds9.get(cmd)
        except Exception:
            return ''
```

---

#### 11.4 — ds9 Frame Queue panel (`main.py`)

A small collapsible panel docked **below the sidebar**, hidden until ds9 is connected.

```
┌──────────────────────────────┐
│  ds9 Frames            [▲]  │
│  ────────────────────────── │
│  ☑  1  Whitelight           │
│  ☑  2  NB 6560–6580 Å      │
│  ☑  3  M0 6540–6600 Å      │
│                       [× ]  │  ← delete selected row
│  [+ Add current]            │
│  [Send all]  [Match WCS]    │
└──────────────────────────────┘
```

- **[+ Add current]**: adds current `ImageCanvas._data_raw` + label (auto-named from image mode: "Whitelight", "NB 6560–6580", "M0 …") to queue; frame number auto-incremented; frame spinbox inline-editable
- **[Send all]**: iterates checked rows, calls `bridge.send_image(data, header, frame)` for each
- **[Match WCS]**: calls `bridge.match_wcs()` after sending
- **[× ]**: removes selected row from queue
- Queue is in-memory only (not persisted); cleared on dataset switch

---

#### 11.5 — ds9 toolbar buttons (`main.py`)

Add to main toolbar (right side, separated by spacer):

```
[Connect ds9]  [→ Send image]  [← Import regions]
```

- **[Connect ds9]**: calls `bridge.connect()`
  - If `no_pyds9` → `QMessageBox`: "pyds9 is not installed.\nInstall with: pip install pyds9"
  - If `no_ds9` → `QMessageBox`: "No running ds9 instance found.\nStart ds9 first, then connect."
  - On success → button text becomes "Disconnect"; status bar shows "ds9: connected ●"
- **[→ Send image]**: shortcut to add current image to frame queue (same as [+ Add current])
- **[← Import regions]**: opens `ImportRegionsDialog` (see 11.6)

All three buttons disabled when no cube loaded. Send/Import also disabled when not connected.

---

#### 11.6 — Import regions from live ds9 (`main.py`)

`ImportRegionsDialog`:
```
  Source: ● Selected region(s) only
          ○ All regions
  [OK]  [Cancel]
```
On OK: `bridge.get_regions(selected_only)` → `parse_ds9_regions()` → for each → `region_to_mask()` → `extract_aperture()` → add to extract strip exactly as a manual extraction.

---

#### 11.7 — File menu additions (`main.py`)

```
File
  ├── Open FITS…               Ctrl+O
  ├── Load Region File…                  ← new
  ├── Save Apertures as Region File…     ← new
  ├── Save Current Frame…      Ctrl+S
  └── Quit                     Ctrl+Q
```

**Load Region File…** (no ds9 required):
- `QFileDialog` → `.reg` file → read text → `parse_ds9_regions()` → ask user: "Extract immediately?" or "Load into batch dialog?"
- If extract immediately: same path as Import from ds9 above

**Save Apertures as Region File…** (no ds9 required):
- Iterates `_extraction_marker_groups`, converts each stored aperture (pixel center + shape + size) back to ds9 region syntax
- Single pixel → `circle(x, y, 1")`, rect → `box(...)`, circular → `circle(...)`, annular → `annulus(...)`
- Writes `.reg` file with `# Region file format: DS9 version 4.1` header and `global color=green` defaults

---

#### 11.8 — Batch extraction from regions (`main.py`)

`Analysis > Batch Extract from Regions…` (also accessible after loading a region file).

**`BatchExtractThread(QThread)`**:
```python
progress = pyqtSignal(int, int, str)   # (current, total, region_name)
finished = pyqtSignal(list)            # list of (label, rb_spectrum)
error    = pyqtSignal(str)             # non-fatal per-region errors accumulated
cancelled = pyqtSignal()
```
- Loop: for each region → `region_to_mask()` → `extract_aperture()` → inject RA/Dec meta → store
- Checks `self._cancel` flag between regions; emits `cancelled()` if set
- Per-region errors: logged and collected, do not abort the loop
- Emits `finished(results)` when done; caller loads into strip and/or saves

**`BatchExtractDialog`**:
```
┌──────────────────────────────────────────┐
│  Batch Extract from Regions              │
│                                          │
│  Source:  ● Region file  [Browse…]       │
│           ○ Live ds9 (selected only)     │
│           ○ Live ds9 (all regions)       │
│                                          │
│  Method:    [Sum ▾]                      │
│  Weighting: [None ▾]                     │
│                                          │
│  Output:                                 │
│  ☑  Load into extract strip             │
│     (warn if > 10 regions)              │
│  ☑  Save to folder  [Browse…]           │
│     Filename:  ● IAU coords (J…)        │
│                ○ Region index (001…)    │
│                ○ Region name (.reg)     │
│     Prefix: [spec_]                     │
│  ☑  Open in rb_multispec when done      │
│                                          │
│  ☑  Run in background                   │
│                                          │
│         [Extract N regions]  [Cancel]    │
└──────────────────────────────────────────┘
```

- **N regions** button label updates after file selected: "Extract 47 regions"
- If "Load into strip" checked and N > 10: show warning "Loading 47 spectra into the strip may be slow. Continue?"
- If background: dialog closes, status bar shows progress: `Extracting regions… 12/47 [■■■□□□] [Cancel]`
- Cancel button in status bar sets `thread._cancel = True`
- On `finished`: load into strip (if checked), save files (if checked), open multispec (if checked)
- On `error` signal: collect all per-region errors, show summary dialog at end

**IAU filename helper** (in `spatial_mask.py`):
```python
from astropy.coordinates import SkyCoord
import astropy.units as u

def iau_name(ra_deg, dec_deg, prefix='spec_'):
    c = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
    s = c.to_string('hmsdms', sep='', precision=1).replace(' ', '')
    return f'{prefix}{s}.fits'
```

---

#### 11.9 — Aperture highlighting (`image_panel.py`, `main.py`)

**Color coding (always on)**:
- Each extraction's aperture marker on the image is drawn in the **same color** as its locked spectrum line in the spectrum panel (Catppuccin 6-color cycle, already used for spectrum lines)
- Color index tracked in `_extraction_marker_groups` alongside the artists

**Highlight on combo selection** (`main.py`):
- `_extract_combo.currentIndexChanged` → `_on_extract_selected(idx)`
- Selected extraction's artists: `linewidth *= 2`, `alpha = 1.0`
- All other extractions: `linewidth` back to normal, `alpha = 0.5`
- Calls `image_canvas.canvas.draw_idle()`

**Click aperture → select in combo** (`image_panel.py`):
- Make aperture artists pickable: `artist.set_picker(5)` (5px tolerance)
- Connect `pick_event` → `_on_artist_picked(event)`
- `_on_artist_picked`: look up which extraction group owns `event.artist` → emit new signal `aperture_picked(int)` with extraction index
- `main.py` connects `aperture_picked` → sets `_extract_combo.setCurrentIndex(idx)` → triggers highlight

---

#### 11.10 — Tests

- Connect with no pyds9 → message box shown, no crash
- Connect with pyds9 but no ds9 running → correct message
- Connect with ds9 running → status bar updates, buttons enable
- Add 3 images to frame queue → Send all → 3 ds9 frames populated
- Match WCS → frames locked
- Draw circle in ds9 → Import selected → spectrum appears in extract strip
- Load `.reg` file with 5 circles → Batch extract → 5 spectra with RA/Dec in meta
- Save spectra → open one in astropy: `fits.getheader(f)['RA']` returns correct value
- Save apertures as region file → reload in ds9 → regions appear correctly
- Background batch of 20 regions → progress bar updates → cancel midway → stops cleanly
- Click aperture on image → combo selects correct extraction
- Select extraction in combo → correct aperture highlights on image

---

### Phase 12 — Status: COMPLETE ✓

All deferred items implemented and tested:

- **Save Subcube…** (`File > Save Subcube…`): writes 3D FITS with updated header (flux + optional VAR extension)
- **Live circular aperture preview**: dashed lavender circle follows cursor in Circular/Circular-Annular mode; hover also shows live extracted spectrum; `(x, y, radius, mode)` cache prevents lag
- **FITS Header viewer** (`View > Show FITS Header…`): extension dropdown, monospace text area, search bar with highlight-all / Next / Prev / match count, Copy button; uses system fixed-width font
- **Right-click context menu** on image canvas: "Crop to selection", "Set as sky region", "Clear sky region", "Clear region overlays" — all enabled/disabled contextually
- **Sexagesimal RA/Dec** toggle (`View > Coordinates: Sexagesimal`)

---

### Phase 13 — Help System (LAST — implement after all features stable)
**⚠ Do not implement until all other phases are complete and tested on real data.**
**Before writing, audit every feature listed below is actually in the code.**

**Goal**: add a `?` toolbar button and `Help > Help` menu entry that opens a tabbed
help dialog matching the rb_multispec style. Include a dedicated `HELP.md` document.

**Pattern**: follow `rbcodes/GUIs/multispecviewer/utils.py → show_help_dialog()` exactly.

**Files to create/modify**:
- `help.py` (new) — `show_help_dialog(parent)` function
- `HELP.md` (new) — full documentation page for the viewer
- `main.py` — add `?` toolbar button + `Help > Help` menu item (shortcut `F1`)

**`help.py`** — tabbed `QDialog` (min 800×560):
- **Overview** tab — interface description: sidebar, image panel, spectrum panel, controls
- **Image Modes** tab — Channel / Whitelight / Narrowband / Cont-sub workflow
- **Extraction** tab — Single pixel / Rectangle / Circular / Circular-Annular modes, weighting options
- **Moment Maps** tab — M0/M1/M2, continuum subtraction, SNR mask, sky region workflow
- **Crop & Analysis** tab — crop by drag, crop by coordinates (pixel/RA-Dec), sky region, save subcube
- **Keyboard Shortcuts** tab — complete reference table (see below)
- **ds9 Bridge** tab — optional, install instructions, workflow (omit if Phase 11 not implemented)
- Bottom: `[Close]`

**Complete keyboard shortcuts to document**:

*Spectrum panel (cursor must be in panel):*
| Key | Action |
|-----|--------|
| `x` | Cursor x → new xmin |
| `X` | Cursor x → new xmax |
| `t` | Cursor y → new ymax |
| `b` | Cursor y → new ymin |
| `r` | Reset to full wavelength range |
| `a` | Autoscale Y to visible data |
| `[` | Pan left 10% |
| `]` | Pan right 10% |
| `-` | Zoom out X 20% |
| `=` | Zoom in X 20% |
| `l` | Toggle legend |
| `c` | On-band window (Cont-sub mode only; press twice to clear) |
| `1` | Cont window 1 (Cont-sub mode only; press twice to clear) |
| `2` | Cont window 2 (Cont-sub mode only; press twice to clear) |
| `Esc` | Return to neutral (Cont-sub mode only) |

*Application-level:*
| Key | Action |
|-----|--------|
| `Ctrl+O` | Open FITS file(s) |
| `Ctrl+S` | Save current frame as FITS |
| `Ctrl+M` | Moment Map dialog |
| `Ctrl+K` | Crop to drawn rectangle |
| `Ctrl+B` | Band Range dialog |
| `Ctrl+R` / `Ctrl+W` | Spectrum Range & Scale dialog |
| `Del` / `Backspace` | Delete selected extraction |
| `Ctrl+Q` | Quit |
| `F1` | Help |

**`HELP.md`** contents:
- Title + one-line description
- Interface overview with panel diagram (ASCII)
- Loading files — Open dialog, multiple files, sidebar switching
- Image controls — colormap, scale, normalization, right-click drag contrast
- Image modes — step-by-step for each mode including Cont-sub band drawing
- Spectral extraction — all four modes, weighting options (None/Var-weighted/Optimal)
- Extraction list — save, send to rb_multispec, delete, clear
- Aperture controls — Rectangle mode drag, Circular mode Extract button
- Moment maps — full workflow: set on-band, optional continuum subtraction, SNR mask
- Crop — drag method and coordinate dialog (pixel and RA/Dec)
- Sky region — how to set, what it's used for, clearing
- Spectrum navigation keystrokes — full table
- Variance cube — how to assign, what features require it
- ds9 Bridge — optional section
- Troubleshooting — variance not found, WCS warnings, pyds9 not available

**Test**: click `?` → dialog opens. All tabs present. Shortcut table complete. Dark theme matches GUI.

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
