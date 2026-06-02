# rb_align — Astrometry Alignment Tool

[Back to Main Page](../../main_readme.md) | [Back to GUIs](../GUIs_readme.md)

A modular astrometry alignment class for aligning 2D images and IFU datacubes
(KCWI, MUSE, NIRSpec IFU, JWST/NIRISS) to a reference frame with full WCS
fitting via `astropy`. Works standalone, in batch pipelines, and embedded in
GUIs.

---

## Quick start

```python
import rbcodes.rb_align as rb_align

# show all workflow examples
rb_align.help()
```

```python
from rbcodes.rb_align import wcs_align

c = wcs_align.from_file(reference='hst_ref.fits', targets=['cube.fits'],
                         input_type='ifu')
c.preprocess()
c.find_sources(strategy='interactive')
c.align()
c.qa(plot=True)
c.write_output()            # → cube_wcsfix.fits
```

---

## Overview

`rb_align` supports:

- **Image↔image** — two 2D FITS images, any instrument
- **IFU↔image** — IFU datacube aligned to a 2D reference (typical case: KCWI/MUSE → HST)
- **IFU↔IFU** — two IFU cubes aligned to each other (same or different instruments)
- **Batch mode** — one interactive session drives alignment of 10+ exposures automatically

### Supported instruments (auto-detected from FITS header)

| Instrument | Auto-detected |
|------------|--------------|
| KCWI | Yes (`INSTRUME = 'KCWI'`) |
| MUSE | Yes (`INSTRUME = 'MUSE'`) |
| Any generic IFU cube | Yes (any 3-D FITS array) |
| 2-D FITS image (HST, ground-based) | Yes |

---

## Installation / dependencies

`rb_align` is part of `rbcodes`. No separate install step is needed.

### Required
`numpy`, `astropy`, `scipy`, `matplotlib`

### Optional

| Package | Enables |
|---------|---------|
| `photutils` | DAOStarFinder, segmentation (`dao`, `knots` strategies), centroid_2dg |
| `astroquery` | Gaia DR3 catalog matching (`gaia` strategy) |

All core functionality (interactive selection, WCS fitting, QA, output) works
without optional dependencies. Missing packages are detected at runtime and
the appropriate strategy or centroiding fallback is used automatically.

---

## Constructors

Load from files — instrument is auto-detected from the FITS header:

```python
c = wcs_align.from_file(
    reference='hst_f814w.fits',          # 2D reference image
    targets=['kcwi_obs.fits'],           # IFU cube(s) to align
    input_type='ifu',                    # 'ifu' or 'image' for all targets
)
```

Load from data arrays already in memory:

```python
c = wcs_align.from_data(
    reference=(arr_ref, wcs_ref, header_ref),
    targets=[(arr1, wcs1, hdr1), (arr2, wcs2, hdr2)],
    input_type='ifu',
)
```

---

## Step-by-step workflow

### 1 — Preprocess (IFU cubes only)

Collapse the 3D cube to a 2D image for alignment. No-op for 2D images.

```python
c.preprocess()                                     # white-light (default)
c.preprocess(method='narrowband', wl_range=[6550, 6570])   # narrow band
c.preprocess(func=my_collapse_func)                # inject your own function
```

`my_collapse_func` must accept a 3D array and return a 2D array.

### 2 — Find sources

```python
c.find_sources(strategy='interactive')             # two-panel click interface
c.find_sources(strategy='interactive',
               stretch='zscale', box=0.1,          # arcsec
               save_catalog='sources.fits')        # save catalog for reuse
```

**Source detection strategies:**

| Strategy | Requires | Best for |
|----------|----------|---------|
| `interactive` | matplotlib | any field — always works |
| `gaia` | astroquery | wide fields with stars, absolute astrometry |
| `dao` | photutils | point source fields |
| `knots` | photutils | arcs and extended sources |
| `cross_corr` | scipy | same-instrument IFU↔IFU only |
| `batch` | saved catalog | 10+ exposures, same field |
| `auto` | varies | unknown fields — cascade in priority order |

**`auto` cascade order:** gaia → dao → knots → cross_corr → interactive

> **Cross-instrument warning:** do NOT use `cross_corr` between HST and
> MUSE/KCWI without PSF-matching. PSF mismatch (HST ~0.05" vs MUSE ~0.6–1")
> will bias the correlation peak. Use `interactive`, `dao`, or `knots` instead.

**`box`** — centroid refinement half-width in **arcsec** (default `0.1`). Converted to pixels per-frame via the WCS pixel scale (minimum 3 px). Typical values: 0.05–0.5" for HST, 0.3–1.0" for KCWI/MUSE.

**Display stretch options** (`stretch=` parameter):

`'zscale'` (default) | `'99.5%'` | `'99%'` | `'98%'` | `'97%'` | `'95%'` |
`'minmax'` | `(vmin, vmax)` tuple

### 3 — Align

```python
c.align()                                          # single target
c.align_many(on_fail='skip', min_sources=2)        # batch — list of targets
```

**WCS fitting logic (automatic):**

| Matched pairs | Fit type |
|--------------|---------|
| ≥ 3 | Full affine (shift + rotation + scale) |
| 2 | Shift + rotation |
| 1 | RA/Dec shift only |
| 0 | Fail per `on_fail` policy |

`on_fail` options: `'skip'` (log and continue), `'raise'` (exception),
`'degrade'` (fall back to simpler fit).

### 4 — Quality assurance

```python
c.qa()                  # print per-frame residuals table
c.qa(plot=True)         # residuals table + overlay figure
c.qa(save='qa.png')     # save figure without displaying
c.qa_summary()          # compact table across all frames (batch)
```

Attributes available after `align()`:

```python
c.rms_residuals      # per-frame RMS in arcsec
c.source_residuals   # per-source residual vectors
c.fit_type           # 'full' | 'shift+rot' | 'shift' per frame
c.n_sources_used     # number of matched sources per frame
c.flagged_frames     # frames that failed or degraded
```

### 5 — Write output

```python
c.write_output()                    # writes {name}_wcsfix.fits
c.write_output(suffix='_wcsfix')    # explicit suffix
c.update_header(target_idx=0)       # patch header in-place (no copy)
```

Original files are never overwritten.

---

## Interactive alignment window

`find_sources(strategy='interactive')` opens a two-panel matplotlib figure
(reference left, target right) with independent zoom and display stretch.

**Interaction:**

| Action | Effect |
|--------|--------|
| Left-click on reference | Pick source → centroid refined, cyan prediction circle appears on target |
| Left-click on target | Confirm pair → recentroid on target |
| Double-click a numbered marker | Toggle edit mode — cyan circle on target; left-click to re-place |
| `Space` | Auto-accept prediction (useful when WCS is already close) |
| Right-click | IDLE: delete nearest pair — PENDING: cancel |
| `u` | Delete nearest pair (hover cursor near it) |
| `Enter` / close window | Finalise and return |
| Scroll wheel | Zoom in/out on each panel independently |

Aim for **4–6 sources spread across the field** for a stable affine fit.

> **Batch strategy:** loads RA/Dec from a saved FITS catalog, reprojects onto the current WCS, recentroids, then opens the interactive window for inspection and editing before returning.

---

## Worked examples

### Single IFU aligned to HST (interactive)

```python
from rbcodes.rb_align import wcs_align

c = wcs_align.from_file(reference='hst_f814w.fits',
                         targets=['kcwi_obs.fits'],
                         input_type='ifu')
c.preprocess(method='narrowband', wl_range=[4860, 4870])
c.find_sources(strategy='interactive', stretch='zscale', box=0.1,  # arcsec
               save_catalog='sources.fits')
c.align()
c.qa(plot=True)
c.write_output()
# → kcwi_obs_wcsfix.fits
```

### Batch survey — 10 exposures, one interactive run

```python
from rbcodes.rb_align import wcs_align

targets = [f'kcwi_{i:02d}.fits' for i in range(10)]
c = wcs_align.from_file(reference='hst_f814w.fits',
                         targets=targets, input_type='ifu')
c.preprocess()
c.find_sources(strategy='interactive', save_catalog='field_sources.fits')
c.align_many(on_fail='skip', min_sources=2)
c.qa_summary()
c.write_output()
# → kcwi_00_wcsfix.fits ... kcwi_09_wcsfix.fits
```

### Reuse saved catalog (skip interactive step)

```python
c = wcs_align.from_file(reference='hst_f814w.fits',
                         targets=['kcwi_obs.fits'], input_type='ifu')
c.preprocess()
c.find_sources(strategy='batch', catalog='sources.fits')
c.align()
c.write_output()
```

### KCWI blue + red cubes (same instrument, cross-correlation)

```python
c = wcs_align.from_file(reference='kcwi_blue.fits',
                         targets=['kcwi_red.fits'], input_type='ifu')
c.preprocess()
c.find_sources(strategy='cross_corr')
c.align()
c.write_output()
# → kcwi_red_wcsfix.fits
```

### Lensed arc / extended source (no point sources)

```python
c = wcs_align.from_file(reference='hst_arc.fits',
                         targets=['muse_arc.fits'], input_type='ifu')
c.preprocess()
c.find_sources(strategy='knots')
c.align()
c.qa(plot=True)
c.write_output()
```

### Inject a custom collapse function

```python
def my_ha_image(cube):
    """Sum flux between Hα ± 5 Å — returns 2D array."""
    return cube[490:510].sum(axis=0)

c = wcs_align.from_file(reference='hst_ref.fits',
                         targets=['muse_cube.fits'], input_type='ifu')
c.preprocess(func=my_ha_image)
c.find_sources(strategy='interactive')
c.align()
c.write_output()
```

---

## Example scripts

Full runnable scripts are in `rbcodes/rb_align/examples/`:

| Script | Description |
|--------|-------------|
| `example_generic.py` | Generic placeholder-path script covering all five workflows |
| `example_J1004_kcwi_hst.py` | Real-data example: KCWI blue aligned to KCWI red for J1004+4112 |

---

## Module layout

```
src/rbcodes/rb_align/
├── __init__.py     # wcs_align, Frame, help()
├── core.py         # wcs_align class, pipeline orchestration, GUI hooks
├── io.py           # load/write wrappers (reuses ifuviewer IO)
├── preprocess.py   # white-light / narrowband collapse (reuses ifuviewer)
├── sources.py      # source detection and interactive matplotlib session
├── align.py        # WCS fitting via fit_wcs_from_points, fallbacks
├── qa.py           # residuals, overlay figures, QA summary
└── examples/
    ├── example_generic.py
    └── example_J1004_kcwi_hst.py
```
