# rbcodes.turbulence — bulk flows and velocity structure functions
[Back to Main Page](../main_readme.md)

Works on `rbcodes.ifu.vmap.VelocityMap` objects (data or mocks; image or source plane).
Velocity maps are made with `rbcodes.ifu` ([IFU tools](../ifu/ifu_readme.md)).

```python
from rbcodes.ifu.vmap import VelocityMap
from rbcodes.turbulence.bulk import Plane, compare_bulk
from rbcodes.turbulence.vsf import StructureFunction

vm = VelocityMap.load('OII_vmap.fits')
results, best = compare_bulk(vm)                       # none, plane, quad, arctan, twolobe; BIC
res = vm.subtract_bulk(Plane())                        # or results['plane']['resid']
sf  = StructureFunction(res, n_bins=20, estimator='debiased')
fit = sf.fit_powerlaw(rmin=5, rmax=10)                 # alpha, alpha_err, A (r0 = 10)
sf.rc, sf.S, sf.err, sf.S_raw, sf.r_geo
```

## bulk

| Model | Parameters |
|-------|------------|
| `Constant` | mean |
| `Plane` | c0 + c1 x + c2 y (`.gradient` gives |∇v| and its position angle) |
| `Quadratic` | 6-term polynomial |
| `ArctanDisk` | x0, y0, pa, q, V_s, r_t, v0: v0 + (2/π) V_s arctan(R/r_t) cos φ |
| `TwoLobe` | xc, yc, θ, v0, Δv, w: v0 + Δv tanh(s/w) |

Polynomials are least-squares fits about the mean position (`weighted=True` uses the
velocity errors). The disk and two-lobe models use a robust bounded fit (soft L1,
`f_scale`=20 km/s); their length bounds and starting values (`rt_bounds_pix`, `rt0_pix`,
`w_bounds_pix`, `w0_pix`) are in pixels and scaled by the pixel spacing of the positions, so
results do not depend on the units. Every model has `fit(x, y, v, v_err)`, `predict(x, y)`,
`describe()`, so `VelocityMap.subtract_bulk` accepts any object with these methods.
`compare_bulk` returns rms, BIC and ΔBIC per model and the residual maps.

## vsf

`StructureFunction(vm, n_bins=14, rmin=None, min_pairs=100, block=None, positions='image',
estimator='debiased', other=None, cross_images=False)`

* **Pairs** between the kept spaxels; separations from the image plane (pixels × kpc/pix) or
  the source plane (`positions='source'`, β from a lens model). Pairs between different
  lensed images are left out unless `cross_images=True`.
* **Bins:** `n_bins` logarithmic bins from `rmin` (default one pixel) to the largest
  separation; bins with fewer than `min_pairs` pairs are dropped (`sf.ok`).
* **Estimators:** `raw` (Δv²), `debiased` (Δv² − σ_i² − σ_j²), `ivw` (debiased, weights
  1/(σ_i² + σ_j²)²), `cross` ((vA_i − vA_j)(vB_i − vB_j) for two maps with independent noise,
  `other=`). `sf.s2(estimator)` evaluates any of them on the same pairs.
* **Errors:** delete-one-block jackknife (`block` pixels, default ≈ 16 blocks). It only sees
  variance from structure smaller than the region; when the field has power on the scale of
  the region it underestimates the scatter between independent regions — use mocks for that.
* **Directional:** `sf.directional(theta, tol=π/6)` → S2 and covariance from pairs within
  `tol` of the axis at angle `theta` and within `tol` of the perpendicular.
* **Fits:** `fit_powerlaw(rmin, rmax, r0=10, noise_floor=False, psf_fwhm=None)` fits
  A (r/r0)^α in linear space with the (shrunk) jackknife covariance; `noise_floor=True` adds a
  constant N (use on the raw S2; N and α trade off when the range is short); `psf_fwhm` fits the
  power law smoothed by a Gaussian PSF (uniform-brightness approximation, averaged over the
  pair separations in each bin). `fit_range()` gives [2 px, 0.5 × largest separation].
* **Geometry:** `major`, `minor` (4√eigenvalues of the position covariance), `r_geo` =
  0.3 × minor; `r_turn(rmax, frac=0.9)`; `increments(sep)` for increment distributions.
* **Regions / images:** `StructureFunction(vm.select(region=...))`, `vm.select(image=1)`.

## forward

```python
from rbcodes.turbulence.forward import TurbulenceField, Observation
obs  = Observation.from_vmap(vm, psf_fwhm_kpc=5.2)                     # mask, flux, errors, scale of the data
mock = obs(TurbulenceField(sigma=10., alpha=2/3, L_kpc=20.), np.random.default_rng(1))
mock = obs.simulate((10., 2/3, np.log10(20.), 0.), rng)               # theta = (sigma, alpha, log10 L, gamma)
```

* **Field** (`grf`): P(k) ∝ (k² + k₀²)^−(α+D)/2, k₀ = 1/L, made in a box `pad`× larger and cut
  at a random offset (large modes present, not periodic). S2 ∝ r^α for r ≪ L.
* **Modes:** `'2D'` — σ is the rms of the observed (PSF-weighted) field in the mask; `'3D'` — a
  slab of depth `depth_kpc`, σ is the rms of the 3-D field, the map is the PSF-weighted
  line-of-sight mean and the width map σ_obs² = W(σ_los² + v_c²) − W(v_c)² + σ_th².
* **Intermittency** (toy): g → g exp(γ g₂ − γ²/2), g₂ a smooth field of scale `L_env`.
* **Observation:** flux-weighted PSF average (F v ⊗ P)/(F ⊗ P), mean removed over the mask,
  optional bulk model, Gaussian noise from the velocity errors; returns a `VelocityMap`.
  Seeds: pass a `numpy.random.Generator`.
* **Lensed systems:** `Observation.from_vmap(vm_traced, positions='source')` (or `beta=(bx, by)`)
  makes the field on a source-plane grid (`src_kpp`, default the median spacing of distinct
  source positions) covering the spaxels' β, samples it at each β (bilinear), then applies the
  PSF, mask and noise in the image plane. Spaxels of different images at the same β get the same
  value; spaxels without β get zero weight. Mocks carry β and image ids, so
  `StructureFunction(mock, positions='source')` and `match_images(mock, ...)` work on them.

```python
traced = lens.raytrace(vm, z=1.70, images={1: 'arc1.reg', 2: 'arc2.reg'})
obs    = Observation.from_vmap(traced, positions='source')
mock   = obs.simulate((10., 2/3, 1.3, 0.), rng)
sf     = StructureFunction(mock, positions='source')
```

## demo

`demo_velocity_field()` → (flux, velocity, turbulent part, galaxy mask) with `DEMO_TRUTH`
(α = 2/3, L = 10 kpc, σ = 20 km/s, gradient (−7.0, −5.8) km/s/pix).
`demo_cube(line, z)` → (cube, truth): a continuum-subtracted cube with that field (MUSE LSF,
σ_int = 60 km/s, R = 1.35 for free-ratio doublets, exact variance), ready for
`NoiseModel` / `GaussianFitter`.

## summaries

```python
from rbcodes.turbulence.summaries import Summaries, S2, Scattering, Kurtosis, LineWidth
sf    = StructureFunction(vm_data, n_bins=20)
summ  = Summaries(sf, blocks=['S2', 'WST', 'KURT', 'SIG'], removal='plane', transform=True)
s_obs = summ(vm_data); s_mock = summ(mock)
```

Every map is measured on the data's spaxels, pairs, bins and errors (`sf`) after fitting and
removing the same bulk model. Blocks: `S2` (debiased S2 in bins `kfit`), `Scattering` /
`'WST'` (reduced wavelet scattering: log S1(j), s2(j1,j2) = S2/S1; J = 4 scales, L = 4
orientations), `Kurtosis` / `'KURT'` (excess kurtosis of increments at 2, 4, 8 px),
`LineWidth` / `'SIG'` (16/50/84th percentiles of width², needs a width map such as 3-D mocks).
`transform=True`: asinh(S2/s2_ref), asinh(kurtosis), log(width²). `sf.s2(u=sf.values(mock))`
gives S2 of any map on the data's pairs.

## inference

```python
from rbcodes.turbulence.inference import (latin_hypercube, Emulator, validate_emulator,
                                          adaptive_metropolis, gelman_rubin, fisher,
                                          moped_vectors, moped_loglike, synthetic_loglike)
design = latin_hypercube(400, bounds, rng)                         # theta = (log10 sigma, alpha, log10 L, gamma)
sims   = np.array([[summ(obs.simulate((10**t[0], *t[1:]), rng)) for _ in range(nrep)] for t in design])
emu_var, report = validate_emulator(design, sims, bounds, rng)     # held-out std(z), emulator error
emu    = Emulator(bounds).fit(design, sims).set_emulator_error(emu_var)
chain, lp, acc = adaptive_metropolis(logpost, x0, 20000, rng)
F, err = fisher(emu, theta_map)                                    # per summary block with idx=
```

* `synthetic_loglike`: Gaussian synthetic likelihood (Wood 2010) with the Hartlap factor.
* `Emulator`: polynomial mean (degree 3) and log-variance (degree 2) in θ scaled to [−1, 1], pooled
  correlation matrix; `set_emulator_error` adds the held-out emulator variance.
* `adaptive_metropolis` (Haario et al. 2001), `gelman_rubin` (R̂), `fisher` (J^T Σ^−1 J),
  `moped_vectors` / `moped_loglike` (score compression; lossless for a linear Gaussian model),
  `credible_interval` for gridded posteriors. Random generators are explicit arguments.

## results

```python
from rbcodes.turbulence.results import Results, pm, asym, num, savefig
R = Results('paper/results.json', 'paper/results.tex')
R.update({'alpha-plane': pm(1.08, 0.46)}, section='nb03')         # a rerun replaces its section
```
In LaTeX: `\input{results.tex}`, `\newcommand{\R}[1]{\csname R@#1\endcsname}`, then `\R{alpha-plane}`.
