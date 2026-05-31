# rb_align — Astrometry Alignment Tool for `rbcodes`

## Overview

`rb_align` is a modular astrometry alignment class for the `rbcodes` package.
It aligns 2D images and/or IFU cubes (e.g. KCWI, MUSE, NIRSpec IFU, JWST/NIRISS) to a reference
frame using a cascade of source detection strategies, with full WCS fitting via
`astropy`. Supports image↔image, IFU↔IFU, and image↔IFU alignment.
Designed to work standalone, in batch pipelines, and embedded in GUIs.

---

## Placement in rbcodes

Following the existing `src/rbcodes` layout, `rb_align` lives at the same
level as `GUIs/`:

```
src/rbcodes/
├── GUIs/
│   ├── ifuviewer/              # existing IFU viewer GUI
│   │   ├── io/                 # ← REUSE: FITS read/write, WCS loading
│   │   │   ├── auto_cube.py    #   load_fits(), _make_2d_wcs()
│   │   │   └── image2d.py      #   load_image()
│   │   └── processing/         # ← REUSE: whitelight/narrowband image creation
│   │       └── cube_collapse.py#   build_whitelight(), build_narrowband()
├── rb_align/              # ← NEW
│   ├── __init__.py
│   ├── core.py
│   ├── io.py              # thin wrapper — imports from GUIs.ifuviewer.io
│   ├── preprocess.py      # thin wrapper — imports from GUIs.ifuviewer.processing
│   ├── sources.py
│   ├── align.py
│   └── qa.py
└── ...
```

## Reuse of Existing rbcodes Modules

**Do not reimplement** — import directly from existing code:

```python
# rb_align/io.py — reuse ifuviewer IO loaders
from rbcodes.GUIs.ifuviewer.io.auto_cube import load_fits   # IFU cubes (auto-detects KCWI/MUSE/generic)
from rbcodes.GUIs.ifuviewer.io.image2d import load_image    # 2D reference images (HST, ground-based)

# rb_align/preprocess.py — reuse ifuviewer collapse functions
from rbcodes.GUIs.ifuviewer.processing.cube_collapse import build_whitelight, build_narrowband
```

`load_fits(path)` auto-detects the instrument from the FITS header (`INSTRUME` keyword) and
returns the appropriate `IFUCube` subclass with `.flux`, `.wave`, `.wcs`, `.header` already
populated. Do not re-read the FITS file manually — populate `Frame` fields from the returned object.

`preprocess.py` wraps `build_whitelight` and `build_narrowband` from ifuviewer to collapse
IFU cubes to 2D images for alignment. No new implementation needed. If a user injects their
own function via `preprocess(func=...)`, it bypasses these wrappers entirely.

---

## Module Structure

```
src/rbcodes/rb_align/
├── __init__.py
├── core.py          # wcs_align class, pipeline orchestration, GUI hooks
├── io.py            # thin wrapper around GUIs.ifuviewer.io
├── preprocess.py    # thin wrapper around GUIs.ifuviewer.processing
├── sources.py       # source detection and centroiding (photutils optional)
├── align.py         # WCS fitting via fit_wcs_from_points, fallbacks
└── qa.py            # residuals, per-source vectors, QA figures
```

---

## Dependencies

### Required
- `numpy`
- `astropy`
- `scipy`
- `matplotlib`

### Optional
- `photutils` — enables DAOStarFinder, segmentation, centroid_2dg
- `astroquery` — enables Gaia DR3 catalog matching

If `photutils` is not installed, centroiding falls back to
`scipy.ndimage.center_of_mass` seamlessly. All core functionality
works without optional dependencies.

---

## Core Classes and Data Structures

### `Frame` dataclass (`core.py`)
Internal container for each image/cube:
```python
@dataclass
class Frame:
    data: np.ndarray            # 2D image or 3D IFU cube (= IFUCube.flux after loading)
    wcs: astropy.wcs.WCS        # always the 2D spatial WCS — never the 3D cube WCS
    header: fits.Header         # original FITS header (3D for IFU cubes)
    input_type: str             # 'image' | 'ifu'
    wave: np.ndarray | None     # wavelength array (= IFUCube.wave); None for 2D images
    image2d: np.ndarray | None  # 2D collapsed image, filled by preprocess(); None until then
```

When constructing from `load_fits()` / `load_image()` results:
- `frame.data   = cube.flux`    (3D for IFU, 2D for image)
- `frame.wave   = cube.wave`    (None for 2D images)
- `frame.wcs    = cube.wcs`     already the 2D spatial WCS — `load_fits()` calls `_make_2d_wcs` internally
- `frame.header = cube.header`  (original, unmodified; patched in-place by `write_output()`)

### `wcs_align` class (`core.py`)
Main class. All intermediate results stored as attributes for inspection.
Reference is always a single `Frame`. Targets is always a list of `Frame`
objects, even if only one.

---

## API Reference

### Constructors (`io.py`)

```python
# From files — input_type per item: 'image' | 'ifu'
# from_file calls load_fits() for IFU and load_image() for 2D images internally.
# load_fits() auto-detects instrument (KCWI/MUSE/generic) from the FITS header.
c = wcs_align.from_file(
      reference='hst_ref.fits',
      targets=['kcwi_01.fits', 'kcwi_02.fits'],
      input_type='ifu'       # applies to all targets; reference always auto-detected
)

# From data arrays — reference and targets as (array, wcs, header) tuples
c = wcs_align.from_data(
      reference=(arr_ref, wcs_ref, header_ref),
      targets=[(arr1, wcs1, hdr1), (arr2, wcs2, hdr2)],
      input_type='ifu'
)
```

### Preprocessing (`preprocess.py`)

```python
# Internal white light image (default for IFU)
c.preprocess()

# Internal narrow band image — wl_range=[wmin, wmax] unpacked as build_narrowband(flux, wave, wmin, wmax)
c.preprocess(method='narrowband', wl_range=[6550, 6570])

# Inject external function (e.g. from rbcodes or your own tools)
# func must accept a cube (3D array) and return a 2D array
c.preprocess(func=my_narrowband_func)
```

- No-op for 2D images — preprocess detects input type automatically
- In batch mode, the method set here applies to **all frames** automatically
- User sets it once; every IFU cube in the target list is processed the same way
- 2D WCS is extracted from the 3D cube header for alignment work (spatial axes only)

### Source Detection (`sources.py`)

```python
# Strategies: 'interactive' | 'gaia' | 'dao' | 'knots' | 'cross_corr' | 'batch' | 'auto'
c.find_sources(strategy='interactive')

# stretch= controls display normalization for both panels independently.
# Options: 'zscale' (default) | '99.5%' | '99%' | '98%' | '97%' | '95%' | 'minmax' | (vmin, vmax)
c.find_sources(strategy='interactive', stretch='zscale')
c.find_sources(strategy='interactive', stretch='99%')
c.find_sources(strategy='interactive', stretch=(100, 5000))   # manual vmin, vmax

# box= sets centroid refinement half-width in pixels (default 15)
c.find_sources(strategy='interactive', box=15)

# Save catalog after interactive run for batch reuse
c.find_sources(strategy='interactive', save_catalog='sources.fits')

# Reuse saved catalog for batch run — no user input needed
c.find_sources(strategy='batch', catalog='sources.fits')

# Auto cascade (see order below)
c.find_sources(strategy='auto')
```

**Source strategies explained:**

| Strategy | Requires | Best for | Notes |
|---|---|---|---|
| `interactive` | matplotlib | any field | click → centroid refinement, always works |
| `gaia` | astroquery | wide fields with stars | absolute astrometry, may miss small IFU fields |
| `dao` | photutils | point source fields | blind DAOStarFinder detection |
| `knots` | photutils | arcs, extended sources | segmentation → bright knot centroids |
| `cross_corr` | scipy only | same-instrument IFU↔IFU | translation only, no rotation/scale |
| `batch` | saved catalog | 10+ exposures same field | reproject + recentroid per frame |
| `auto` | varies | unknown fields | cascade in priority order below |

**Auto cascade order (reliability first):**
1. `gaia` — absolute astrometry (if astroquery available and sources found)
2. `dao` — point source detection (if photutils available)
3. `knots` — extended/arc fields (if photutils available)
4. `cross_corr` — only if same instrument or no sources found; translation only
5. `interactive` — always available as final fallback

> **Cross-instrument warning:** Do NOT use `cross_corr` between HST and MUSE/KCWI
> without PSF-matching first. PSF mismatch (HST ~0.05" vs MUSE ~0.6-1") will bias
> the correlation peak. Use `dao`, `knots`, or `interactive` for cross-instrument work.

### Alignment (`align.py`)

```python
# Single target
c.align()

# Batch — handles partial field overlap per frame automatically
# on_fail: 'skip' | 'raise' | 'degrade'
c.align_many(on_fail='skip', min_sources=2)
```

**WCS fitting logic:**
- 3+ matched pairs → full affine (shift + rotation + scale) via `fit_wcs_from_points`
- 2 pairs → shift + rotation only
- 1 pair → delta RA/Dec shift only (CRVAL1/CRVAL2)
- 0 pairs → fail per `on_fail` policy

**Partial overlap handling in `align_many`:**
For each frame, catalog sky positions are projected to pixels. Only sources
that fall within image bounds are used. Frames with fewer than `min_sources`
are flagged and handled per `on_fail`.

### QA (`qa.py`)

```python
c.qa()                       # print per-frame residuals summary
c.qa(plot=True)              # show overlay figure
c.qa(save='qa_fig.png')      # save figure without displaying
c.qa_summary()               # compact table of all frames (useful after align_many)
```

**Stored attributes after `align()`:**
```python
c.rms_residuals      # per-frame RMS in arcsec
c.source_residuals   # per-source residual vectors
c.fit_type           # 'full' | 'shift+rot' | 'shift' per frame
c.n_sources_used     # number of matched sources per frame
c.flagged_frames     # list of frames that failed or degraded
```

### Output (`io.py`)

```python
# Write corrected files — default suffix _wcsfix, originals untouched
c.write_output()
# → kcwi_01_wcsfix.fits, kcwi_02_wcsfix.fits, ...

# Custom suffix
c.write_output(suffix='_wcsfix')

# Patch header of one existing file in place (no copy)
c.update_header(target_idx=2)
```

Output filename convention: `{original_name}_wcsfix.fits`
Original files are never overwritten by default.

---

## Typical Workflows

### Single IFU aligned to HST reference (interactive)
```python
from rbcodes.rb_align import wcs_align

c = wcs_align.from_file(reference='hst_ref.fits', targets=['kcwi_obs.fits'],
                      input_type='ifu')
c.preprocess(method='narrowband', wl_range=[4860, 4870])
c.find_sources(strategy='interactive')
c.align()
c.qa(plot=True)
c.write_output()
# → kcwi_obs_wcsfix.fits
```

### Batch IFU survey — 10 exposures, reuse interactive catalog
```python
from rbcodes.rb_align import wcs_align

targets = [f'kcwi_{i:02d}.fits' for i in range(10)]
c = wcs_align.from_file(reference='hst_ref.fits', targets=targets, input_type='ifu')
c.preprocess()
c.find_sources(strategy='interactive', save_catalog='field_sources.fits')
c.align_many(on_fail='skip', min_sources=2)
c.qa_summary()
c.write_output()
# → kcwi_00_wcsfix.fits ... kcwi_09_wcsfix.fits
```

### KCWI blue + red datacubes (same target, cross-correlation)
```python
# Same instrument, same pointing, same PSF — cross_corr is ideal
c = wcs_align.from_file(
      reference='kcwi_blue.fits',
      targets=['kcwi_red.fits'],
      input_type='ifu')
c.preprocess()   # white light image works for both
c.find_sources(strategy='cross_corr')   # same PSF, same morphology
c.align()
c.write_output()
# → kcwi_red_wcsfix.fits
```
Fall back to `strategy='interactive'` if one cube is dominated by a strong
emission line that makes white light morphology too different from the other.

### IFU↔IFU alignment (same instrument, cross-correlation)
```python
c = wcs_align.from_file(reference='kcwi_ref.fits', targets=['kcwi_obs.fits'],
                      input_type='ifu')
c.preprocess()
c.find_sources(strategy='cross_corr')   # safe: same instrument
c.align()
c.write_output()
```

### Lensed arc — no point sources
```python
c = wcs_align.from_file(reference='hst_arc.fits', targets=['muse_arc.fits'],
                      input_type='ifu')
c.preprocess()
c.find_sources(strategy='knots')   # segmentation → bright knot centroids
c.align()
c.qa(plot=True)
c.write_output()
```

### External narrow band function injection
```python
from rbcodes.GUIs.ifuviewer.processing.cube_collapse import build_narrowband  # or your own function

c = wcs_align.from_file(reference='ref.fits', targets=['cube.fits'], input_type='ifu')
c.preprocess(func=my_narrowband_func)     # class calls your function, expects 2D array back
c.find_sources(strategy='auto')
c.align()
c.write_output()
```

---

## GUI Integration

`rb_align` separates display logic from alignment logic entirely.
The interactive matplotlib canvas is replaceable by any external GUI.

### Exposed hooks (`core.py`)

```python
# Returns 2D numpy array ready for display (summed over wavelength if IFU)
arr = c.get_display_image(idx=0)

# Add a matched pair by pixel coordinates on reference and target.
# ra, dec are computed from ref_x, ref_y via the reference WCS.
c.add_pair(ref_x=245.3, ref_y=312.7, tgt_x=189.1, tgt_y=204.5)

# Assign external callback — called AFTER a pair is stored.
# Signature: callback(ref_x, ref_y, tgt_x, tgt_y, ra, dec) -> None
# Use it to update GUI overlays on both panels.
c.on_pair_callback = my_gui_pair_handler

# Remove the last stored pair (undo)
c.remove_last_pair()

# Clear all pairs and restart
c.clear_pairs()

# Inspect current pair list — list of (ref_x, ref_y, ra, dec, tgt_x, tgt_y) tuples
print(c.pairs)
```

### PyQt/PySide embedding
```python
# The external GUI manages the two-panel state machine itself,
# calling add_pair() once both clicks are confirmed.
aligner.add_pair(ref_x, ref_y, tgt_x, tgt_y)
```

### Jupyter notebook
```python
c.find_sources(strategy='interactive')   # renders inline, click as normal
```

The class never assumes how clicks arrive — it only requires confirmed pairs
via `add_pair(ref_x, ref_y, tgt_x, tgt_y)`. Any GUI that can manage the
two-panel state machine and call a Python function can drive `rb_align`.

---

## Implementation Notes for Claude Code

### `sources.py` — interactive session design

The interactive strategy opens a two-panel matplotlib figure (reference left, target right).
Sources are selected as **matched pairs** — one click on reference followed by one click on
target. The state machine enforces this strictly: a double-click on the same panel is ignored.
A pair is never written to the catalog until both sides are confirmed.

#### Layout

```python
fig, (ax_ref, ax_tgt) = plt.subplots(1, 2, figsize=(12, 6))
fig.suptitle(...)   # updated live to show current state and key bindings
```

Each panel has its own display normalization computed independently via `_compute_norm`.
Zoom and pan are **independent** — do not link axes. HST and MUSE/KCWI have very
different pixel scales (e.g. 0.05"/px vs 0.2"/px), so linked zoom would make one
panel unusable.

#### State machine

```
IDLE    — only left-clicks on ax_ref are accepted
PENDING — only left-clicks on ax_tgt are accepted; right-click anywhere cancels

IDLE → [left-click ax_ref]
         centroid around click within box px
         convert (x_ref, y_ref) → (ra, dec) via reference WCS
         project (ra, dec) → (x_pred, y_pred) via current target WCS
         draw solid red + on ax_ref at centroided position
         draw dashed cyan circle on ax_tgt at predicted position
         → PENDING

PENDING → [left-click ax_tgt]
         centroid around click within box px on target image
         store pair: (x_ref, y_ref, ra, dec, x_tgt, y_tgt)
         replace dashed cyan circle with solid red + on ax_tgt
         add source number label on both panels
         call on_pair_callback if set
         → IDLE

PENDING → [spacebar]
         auto-accept: recentroid around predicted position (x_pred, y_pred)
         store pair as above
         → IDLE

PENDING → [right-click anywhere]
         discard pending reference click — remove red + from ax_ref
         remove dashed cyan circle from ax_tgt
         → IDLE

any state → [u]
         remove last stored pair — remove markers from both panels
         renumber remaining sources

any state → [Enter or close window]
         finalise — return stored pairs
```

The dashed cyan circle on ax_tgt is a **prediction guide only**, not a required click
target. When the initial WCS is very wrong the circle will be far from the star; the user
simply clicks on the actual star position instead.

#### Title bar as state indicator

```python
IDLE:    'Left-click reference to select source  |  u: undo  |  Enter: done'
PENDING: 'Source N pending — click target to confirm  |  Space: accept prediction  |  Right-click: cancel'
```

#### Zoom

Scroll-wheel zoom on each panel independently:
```python
def _on_scroll(event, ax):
    if event.inaxes is not ax:
        return
    factor = 0.85 if event.button == 'up' else 1.15
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    cx, cy = event.xdata, event.ydata
    ax.set_xlim([cx + (x - cx) * factor for x in xlim])
    ax.set_ylim([cy + (y - cy) * factor for y in ylim])
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('scroll_event', lambda e: _on_scroll(e, ax_ref))
fig.canvas.mpl_connect('scroll_event', lambda e: _on_scroll(e, ax_tgt))
```

Toolbar zoom/pan remains available as backup. Before processing any click, check:
```python
if fig.canvas.toolbar.mode != '':   # '' = pointer mode; ignore clicks during zoom/pan
    return
```

#### Display normalization

```python
from astropy.visualization import ZScaleInterval
import matplotlib.colors as mcolors

def _compute_norm(data2d, stretch='zscale'):
    finite = data2d[np.isfinite(data2d)]
    if stretch == 'zscale':
        try:
            vmin, vmax = ZScaleInterval().get_limits(finite)
        except Exception:
            vmin, vmax = np.nanpercentile(finite, [1, 99])
    elif stretch == 'minmax':
        vmin, vmax = float(finite.min()), float(finite.max())
    elif isinstance(stretch, (tuple, list)):
        vmin, vmax = float(stretch[0]), float(stretch[1])   # manual (vmin, vmax)
    else:
        pct = float(str(stretch).rstrip('%'))
        vmin, vmax = np.nanpercentile(finite, [100 - pct, pct])
    return mcolors.Normalize(vmin=vmin, vmax=vmax)
```

Supported `stretch=` values (match ifuviewer `NORMS` convention):
`'zscale'` (default) | `'99.5%'` | `'99%'` | `'98%'` | `'97%'` | `'95%'` | `'minmax'` | `(vmin, vmax)`

Each panel calls `_compute_norm(frame.image2d, stretch)` independently — reference and
target almost always have different flux scales and should never share a norm.

#### Catalog written by interactive session

Each confirmed pair appends one row:
```
x_ref, y_ref, ra, dec, x_target, y_target
```
`ra`, `dec` come from the reference WCS (sky truth). `x_target`, `y_target` come from
centroiding on the target frame. `residual_arcsec` is added after `align()`.

For `align_many`: the interactive session runs once against the **first target only**.
The saved catalog (sky positions `ra`, `dec`) is then used for all remaining targets
via the `batch` strategy — project sky → pixel, recentroid, fit. The user clicks once;
the batch loop handles the rest.

### `sources.py` — centroiding fallback
```python
try:
    from photutils.centroids import centroid_com, centroid_2dg
    HAS_PHOTUTILS = True
except ImportError:
    from scipy.ndimage import center_of_mass
    def centroid_com(cutout):
        return center_of_mass(cutout)   # identical behavior
    HAS_PHOTUTILS = False

# Use centroid_2dg if available (more robust), else centroid_com
def refine_centroid(cutout):
    if HAS_PHOTUTILS:
        return centroid_2dg(cutout)
    return centroid_com(cutout)
```

### `preprocess.py` — IFU 3D→2D WCS handling (important)

IFU cubes have a 3D WCS (RA, Dec, wavelength). All alignment works on the
spatial 2D WCS only. The wavelength axis must never be touched.

**Do not reimplement 3D→2D WCS extraction.** `_make_2d_wcs` already exists in
`rbcodes.GUIs.ifuviewer.io.auto_cube` and is used by all `IFUCube` loaders. It uses
`WCS(header, naxis=2)` which avoids axis-count mismatch warnings from astropy.
Reuse it directly:

```python
# load_fits() already calls _make_2d_wcs internally — IFUCube.wcs is already 2D.
# When building a Frame from a loaded cube, just assign directly:
frame.wcs = cube.wcs   # already the 2D spatial WCS; no need to call _make_2d_wcs again
```

To patch the corrected WCS back into the 3D cube header for output, update only spatial keywords:

```python
def _update_3d_header(header_3d, new_2d_wcs):
    """Patch only spatial WCS keywords back into 3D cube header.
    Wavelength axis keywords (CRVAL3, CDELT3, CTYPE3 etc.) are untouched."""
    new_hdr = new_2d_wcs.to_header()
    spatial_keys = ['CRVAL1','CRVAL2','CRPIX1','CRPIX2',
                    'CD1_1','CD1_2','CD2_1','CD2_2',
                    'CDELT1','CDELT2','CTYPE1','CTYPE2','CUNIT1','CUNIT2']
    for key in spatial_keys:
        if key in new_hdr:
            header_3d[key] = new_hdr[key]
    return header_3d
```

- `Frame.wcs` is **always** the 2D spatial WCS — set at load time via `_make_2d_wcs`, never the raw 3D WCS
- Final corrected WCS is always written back to the **original 3D cube header**
- This applies to both single and batch modes

### `align.py` — WCS fitting
```python
from astropy.wcs.utils import fit_wcs_from_points

new_wcs = fit_wcs_from_points(
    xy=pixel_coords,         # (2, N) array of pixel positions
    world_coords=sky_coords  # SkyCoord object of matched sky positions
)
# Fall back to CRVAL shift when only 1 pair:
# header['CRVAL1'] += delta_ra
# header['CRVAL2'] += delta_dec
```

### `align_many` — partial overlap per frame
```python
for frame in self.targets:
    # frame.wcs is always the 2D spatial WCS (set by _make_2d_wcs at load time).
    # Never call all_world2pix on a 3D WCS — it will fail or return wrong results.
    px, py = frame.wcs.all_world2pix(catalog_ra, catalog_dec, 0)
    # keep only in-bounds sources
    mask = (px > 0) & (px < nx) & (py > 0) & (py < ny)
    if mask.sum() < min_sources:
        handle_fail(frame, on_fail)
        continue
    # recentroid and fit
```

### `io.py` — output filename convention
```python
from pathlib import Path
def _output_path(input_path, suffix='_wcsfix'):
    p = Path(input_path)
    return p.parent / f"{p.stem}{suffix}{p.suffix}"
```

### Catalog format (FITS table)
```
columns: x_ref, y_ref, ra, dec, x_target, y_target, residual_arcsec
```
Save with `astropy.table.Table` as FITS for reuse in batch runs.

---

## Testing Suggestions

- Unit test each module independently with synthetic FITS data
- Test `from_data` and `from_file` constructors separately
- Test `align_many` with mock set where 2 of 5 targets have partial overlap
- Test photutils-absent fallback by mocking the import at test time
- Test GUI hooks with headless backend: `matplotlib.use('Agg')`
- Test cross-correlation on shifted synthetic IFU pair with known offset

---

## Future Extensions (not in v1)
- PSF matching before cross-correlation (needed for cross-instrument use)
- Proper motion correction for Gaia sources
- Distortion correction (SIP coefficients)
