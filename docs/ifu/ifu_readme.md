# rbcodes.ifu — IFU cube tools
[Back to Main Page](../main_readme.md)

Reading IFU cubes and making images, moment maps, spectra and region masks from them.
These are the functions behind [rb_ifuview](../GUIs/ifuview/rb_ifuview.md) and `rb_align`,
usable directly in scripts and notebooks. (Moved here from `GUIs/ifuviewer/io` and
`GUIs/ifuviewer/processing` in 3.1.0.)

```python
from rbcodes.ifu import load_fits
cube = load_fits('cube.fits')              # optional: load_fits('cube.fits', var='var.fits')
cube.flux    # (n_wave, ny, nx)
cube.var     # same shape, or None
cube.wave    # (n_wave,) Angstrom
cube.wcs     # 2D spatial WCS
```

## Modules

| Module | Contents |
|--------|----------|
| `io` | `load_fits(path, var=None)` — reads the `INSTRUME` keyword and returns a `MUSECube` (flux in `DATA`, variance in `STAT`), `KCWICube` or `GenericCube`; a 2D file gives a `FITSImage`. Variance is taken from a `VAR`/`VARIANCE`/`STAT`/`ERR` extension or a separate file. `load_image(path)` reads a 2D image. New instruments: subclass `GenericCube` and add it to `_INSTRUMENT_MAP` in `io/auto_cube.py`. |
| `cube_collapse` | `build_whitelight(flux, wave, wmin, wmax, method)`, `build_narrowband`, `build_continuum_sub` (on-band minus one or two continuum windows) |
| `moment_maps` | `moment0/1/2(flux, wave, wmin, wmax[, lambda_rest])` (flux, velocity in km/s, dispersion in km/s), `moment_map(..., order)`, `compute_snr_map`, `sky_stats`, `subtract_linear_continuum`, `velocity_array` |
| `aperture_extract` | `make_circular_mask`, `make_square_mask`, `make_annulus_mask`; `extract_aperture` (sum, error from variance), `extract_variance_weighted` (inverse-variance weighted mean), `extract_with_method`, `subtract_background` |
| `spatial_mask` | `parse_ds9_regions(text)`, `region_to_mask(region, wcs, ny, nx)` (pixel-centre rule; image or sky coordinates), `draw_region_overlay`, `regions_to_ds9_text`, `iau_name` |

## Example: velocity map inside a ds9 region

```python
from rbcodes.ifu import load_fits
from rbcodes.ifu.moment_maps import moment1, subtract_linear_continuum
from rbcodes.ifu.spatial_mask import parse_ds9_regions, region_to_mask

cube = load_fits('cube.fits')
ny, nx = cube.flux.shape[1:]
line = subtract_linear_continuum(cube.flux, cube.wave, 6630, 6640, 6670, 6680)
v = moment1(line, cube.wave, 6645, 6665, lambda_rest=6655.0)   # km/s, observed-frame centre

reg = parse_ds9_regions(open('nebula.reg').read())[0]
v[~region_to_mask(reg, cube.wcs, ny, nx)] = float('nan')
```

Sky regions are converted with the cube WCS; give them in the cube's frame (usually `icrs`),
since `fk5` and `icrs` differ by ~0.02".
