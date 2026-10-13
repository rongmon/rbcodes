# rbcodes.ifu — IFU cube tools
[Back to Main Page](../main_readme.md)

Reading IFU cubes and making images, moment maps, spectra and region masks from them.
These are the functions behind [rb_ifuview](../GUIs/ifuview/rb_ifuview.md) and `rb_align`,
usable directly in scripts and notebooks. (Moved here from `GUIs/ifuviewer/io` and
`GUIs/ifuviewer/processing` in 3.1.0.)

```python
from rbcodes.ifu import load_fits
cube = load_fits('cube.fits')              # optional: load_fits('cube.fits', var='var.fits')
cube = load_fits('big_cube.fits', wave_range=(7520, 7600), air_to_vacuum=True)   # read only a slab
cube.flux    # (n_wave, ny, nx)
cube.var     # same shape, or None
cube.wave    # (n_wave,) Angstrom
cube.wcs     # 2D spatial WCS
```

## Modules

| Module | Contents |
|--------|----------|
| `io` | `load_fits(path, var=None, wave_range=None, air_to_vacuum=False)` — `wave_range` reads only the channels in that range (default: all); `air_to_vacuum` converts an air-wavelength cube (CTYPE3 `AWAV`) to vacuum (`cube.airvac` records the convention).  reads the `INSTRUME` keyword and returns a `MUSECube` (flux in `DATA`, variance in `STAT`), `KCWICube` or `GenericCube`; a 2D file gives a `FITSImage`. Variance is taken from a `VAR`/`VARIANCE`/`STAT`/`ERR` extension or a separate file. `load_image(path)` reads a 2D image. New instruments: subclass `GenericCube` and add it to `_INSTRUMENT_MAP` in `io/auto_cube.py`. |
| `cube_collapse` | `build_whitelight(flux, wave, wmin, wmax, method)`, `build_narrowband`, `build_continuum_sub` (on-band minus one or two continuum windows) |
| `moment_maps` | `moment0/1/2(flux, wave, wmin, wmax[, lambda_rest])` (flux, velocity in km/s, dispersion in km/s), `moment_map(..., order)`, `compute_snr_map`, `sky_stats`, `subtract_linear_continuum`, `velocity_array` |
| `aperture_extract` | `make_circular_mask`, `make_square_mask`, `make_annulus_mask`; `extract_aperture` (sum, error from variance), `extract_variance_weighted` (inverse-variance weighted mean), `extract_with_method`, `subtract_background` |
| `noise` | `NoiseModel.from_cube(cube, linefree=..., sky_box= / sky_mask= / sky_region=, spaxel_scale=True)` — sqrt(variance), or the hybrid model σ_sky(λ)·s(y,x) from a user-defined sky region; `.correlation(...)` estimates adjacent-spaxel / adjacent-channel noise correlation |
| `lines` | `LINES` catalogue (air and vacuum rest wavelengths), `LineModel` (single lines, multiplets with fixed or bounded free ratios, several lines tied in v and σ), LSFs (`lsf_muse`, `lsf_resolution(R)`, `lsf_constant`) |
| `measure` | `GaussianFitter`, `MomentMeasurer`, `MCErrors`, `QualityMask`, `SpaxelDecisions`, `FitResult` (maps; FITS save/load), `spaxel_report`, `plot_spaxel` |
| `vmap` | `VelocityMap`: v, σ_v, width, flux, mask, rejection codes; image-plane pixel positions and optional source-plane β (kpc) + image id; kpc/pix, PSF, line, z, provenance; `from_fit`, `from_moments`, `select`, `subtract_bulk`, `points`, `save`/`load`, `plot` |
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

## Velocity measurements (`lines`, `noise`, `measure`)

```python
from rbcodes.ifu import load_fits
from rbcodes.ifu.lines import LineModel, lsf_muse
from rbcodes.ifu.noise import NoiseModel
from rbcodes.ifu.measure import GaussianFitter, QualityMask, SpaxelDecisions, FitResult

cube = load_fits('cube.fits')
lm   = LineModel('OII', z=0.785273, windows=dict(flux=(6640, 6680), fit=(6640, 6680)),
                 air=False, lsf=lsf_muse)                      # air/vacuum must match the cube
nm   = NoiseModel.from_cube(cube, linefree=lm, sky_box=(0, 5, 0, 8))
dec  = SpaxelDecisions('OII_decisions.json')
fit  = GaussianFitter(lm, vmax=600, sig_bounds=(1, 500)).fit_cube(cube, nm, detect_snr=3, decisions=dec)
qm   = QualityMask(snr_min=5, verr_max=30, chi2_max=3).apply(fit, decisions=dec)
fit.save('OII_fit.fits')                                     # FitResult.load(...)
```

**Nothing instrument- or physics-specific is fixed.** Air or vacuum (`air`), the LSF, the
windows, doublet ratio bounds (`ratio_bounds`), fixed ratios (`ratios`) and the line catalogue
(`catalogue`) are all arguments; `LINES` only provides defaults.

**Model.** All lines share one velocity v and one intrinsic dispersion s:
f(λ) = C_g + F_g Σ w_k G(λ; λ_k, σ_k) in the fit window of line g, λ_k = λ_rest,k (1+z)(1+v/c),
σ_k² = (λ_k s/c)² + σ_LSF(λ_k)². G has unit area, so F is the line flux; with an LSF, `sig` is
the intrinsic dispersion. Doublets: `ratio='free'` (R = F of the 2nd catalogue entry / 1st,
within bounds), `'fixed'`, or `'single'`. Several lines: `LineModel(['OII', 'Hbeta'], z,
windows={'OII': dict(fit=...), 'Hbeta': dict(fit=...)})` — independent fluxes and continua,
each line in its own fit window with its own LSF (`lsf={'OII': ..., 'Hbeta': ...}`).

**Continuum:** `LineModel(..., continuum='linear')` fits C + S (λ − λ_c) for each line (λ_c = centre
of its fit window) — for cubes that are not continuum-subtracted. It starts from a straight-line
fit to the fit-window channels more than `vmax` km/s from the lines and is unbounded; maps `S`,
`S_err`. The default `'constant'` keeps C within ±`cont_bound` × noise of `cont_p0`.

**Windows** are observed-frame Å per line: `fit` (required), `flux` (detection S/N; default
`fit`), `exclude` (ranges dropped from the fit), `linefree` (noise estimate; default: more
than 5 Å outside `fit`).

**GaussianFitter.** `curve_fit` with `absolute_sigma=True`; errors × sqrt(max(1, χ²_ν)) unless
`inflate=False`. Bounds: ±`vmax`, `sig_bounds`, ratio bounds, continuum ±`cont_bound` × median
noise around `cont_p0` (defaults 5 and 0 assume continuum-subtracted cubes). Fits are done in
units of the median noise, so any flux units work. Per-spaxel `overrides={(iy, ix): dict(fit=,
exclude=, vmax=, sig_bounds=, ratio_bounds=, p0_v=, p0_sig=)}`; `refit(...)` fits one spaxel.

**MomentMeasurer.** Fixed velocity window moments (v, width including the instrument, flux) with
analytic errors; the components of tied lines are stacked in velocity with inverse-variance
weights (gains S/N for lines of similar strength). Choose `vwin` close to the line extent.

**MCErrors.** Injects noise (white, or correlated with Gaussian smoothing of `corr_spatial`
spaxels / `corr_spectral` channels) into the best-fit models, refits, and returns the ratio of
the scatter in v to the reported σ_v. The seed is a required argument.

**QualityMask** codes (`REJECT`): 0 kept, 1 not attempted / bad data, 2 fit failed, 3 low S/N,
4 large σ_v, 5 high χ²_ν, 6 at velocity bound, 7 isolated pixel, 8 small island, 9 forced
include, 10 forced exclude. No hole filling.

**SpaxelDecisions** (JSON file): `accept(fit, iy, ix)` keeps a spaxel; `accept(fit, iy, ix,
fitter.refit(...))` writes a refit into the maps and stores its override; `reject(line, iy, ix)`;
`forget(line, iy, ix)`. Re-applied by `fit_cube(..., decisions=)` and `QualityMask.apply(...,
decisions=)`; they change the maps without refitting the cube.

**No variance extension?** Give a blank-sky region (`sky_box`, `sky_mask`, or a ds9
`sky_region`) and line-free channels: σ(λ, y, x) = σ_sky(λ) × s(y, x), where σ_sky is the
robust scatter across the sky spaxels per channel and s the robust scatter of f/σ_sky over the
line-free channels of each spaxel (`spaxel_scale=False` drops s).

## VelocityMap (`vmap`)

Data and mocks, lensed or unlensed, all become a `VelocityMap`; everything downstream
(`rbcodes.turbulence`: bulk flows, structure functions, summaries, inference) only uses it.

```python
from rbcodes.ifu.vmap import VelocityMap
vm = VelocityMap.from_fit(fit, qm, kpc_per_pix=1.49, psf_fwhm_kpc=5.2)
vm.save('OII_vmap.fits')                      # VelocityMap.load(...)
inner = vm.select(region='inner.reg')         # ds9 region / bool mask; select(image=1) for lensed images
x, y, v, v_err = vm.points()                  # kept spaxels; positions='source' uses beta (kpc)
res = vm.subtract_bulk(model)                 # model.fit(x, y, v, v_err), model.predict(x, y)
```

Rejected spaxels keep their values but are outside `mask`; nothing is filled. Each operation
appends to `provenance`. Source-plane positions are set with `set_source_positions(beta_x,
beta_y, image)` (kpc), e.g. from a lens model.

## Notes

Sky regions are converted with the cube WCS; give them in the cube's frame (usually `icrs`),
since `fk5` and `icrs` differ by ~0.02".
