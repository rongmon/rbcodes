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
| `noise` | `noise_cube(flux, var=None, line_free=None, sky_mask=None)` — sqrt(variance), or, without a variance cube, a sky spectrum (robust spread over your sky region per channel) times a per-spaxel scale (robust spread over line-free channels). `sky_mask_from_region(reg, wcs, shape)` makes the sky mask from a ds9 region file or text |
| `linefit` | Emission-line fitting → velocity, dispersion, flux, doublet-ratio maps with errors; quality cuts (see below) |
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

## Line fitting (`rbcodes.ifu.linefit`)

```python
from rbcodes.ifu import load_fits
from rbcodes.ifu.linefit import EmissionLine, fit_cube, quality_mask, LineFitResult

cube = load_fits('muse_cube.fits')
line = EmissionLine('OII', z=0.785273, air=True, lsf='muse')
res  = fit_cube(cube, line, exclude=[(6680, 6685)])      # noise from the variance extension
keep, code = quality_mask(res, snr_min=5, verr_max=30, chi2_max=3)
res.save('OII_fit.fits')                                  # all maps + WCS; LineFitResult.load(...)
v = res['v']                                              # km/s relative to z
```

**Model.** All components share one velocity v and one intrinsic dispersion s:
f(λ) = C + F Σ w_k G(λ; λ_k, σ_k), with λ_k = λ_rest,k (1+z)(1+v/c) and
σ_k² = (λ_k s/c)² + σ_LSF(λ_k)². G has unit area, so F is the total line flux. The instrumental
width is part of the model, so `sig` is the intrinsic dispersion.

**Lines.** Presets `OII` (doublet, free ratio R = F3729/F3726), `OIII`, `OIII_doublet` (ratio fixed
to 1/2.98), `Hbeta`, `Halpha`, `NII` (fixed 1/3.05), `MgII` (free R = F2803/F2796). Rest
wavelengths come from the rbcodes line lists (vacuum); `air=True` converts them to air (MUSE
cubes are in air — using vacuum wavelengths on an air cube shifts velocities by ~85 km/s).
Custom lines: `EmissionLine('CIII', z=2, rest=(1906.68, 1908.73), ratio_bounds=(0.5, 2))`.

**Instrumental width (`lsf`).** `'muse'` (Bacon et al. 2017), `lsf_resolution(R)` for a constant
resolving power (e.g. KCWI), a number (constant σ in Å), any function of observed wavelength,
or None.

**Fitting options** (`fit_cube`): `window` (default ±max(900, vmax+300) km/s around the
line), `exclude` ranges (sky residuals), `vmax`, `sig_bounds`, `ratio_bounds`, `detect` /
`detect_snr` (which spaxels to fit; the detection S/N counts continuum too, so pass `detect`
for cubes that are not continuum-subtracted), `overrides={(iy, ix): {...}}` for single spaxels.
Errors are from the fit covariance with the given noise, multiplied by sqrt(max(1, χ²_ν))
(`inflate=False` to switch off). Fitting is done in units of the median noise, so any flux
units work.

**Maps** (`res[key]`): `v`, `v_err`, `sig`, `sig_err`, `F`, `F_err`, `R`, `R_err`, `C`,
`C_err`, `chi2`, `status` (NaN = not fitted, 0 ok, 1 bad data, 2 fit failed), `det_snr`.

**Quality cuts.** `quality_mask` applies S/N, velocity-error, χ²_ν and velocity-bound cuts, then
removes isolated pixels and small islands; `code` gives the first reason each spaxel was
rejected (`REJECT`), with `force_include` / `force_exclude` lists. `spaxel_report(res, code, iy, ix)`
prints a one-line diagnosis and `plot_spaxel(axes, line, wave, flux, noise, iy, ix, res)` plots
data, fit, components and residuals.

**No variance extension?** Define a blank-sky region (no sources, away from the cube edges)
and estimate the noise from it:
```python
from rbcodes.ifu.noise import noise_cube, sky_mask_from_region
from rbcodes.ifu.linefit import fit_domain
sky = sky_mask_from_region('sky.reg', cube.wcs, cube.flux.shape[1:])   # or any bool (ny, nx) array
line_free = ~fit_domain(cube.wave, line.default_window(1000))
noise, _ = noise_cube(cube.flux, line_free=line_free, sky_mask=sky)
res = fit_cube(cube, line, noise=noise)
```
noise(λ, y, x) = σ_sky(λ) × s(y, x): σ_sky is the robust scatter (1.4826 MAD) across the sky
spaxels in each channel (follows sky-line residuals); s is the robust scatter of flux/σ_sky over
the line-free channels of each spaxel (follows exposure-map variations).

## Notes

Sky regions are converted with the cube WCS; give them in the cube's frame (usually `icrs`),
since `fk5` and `icrs` differ by ~0.02".
