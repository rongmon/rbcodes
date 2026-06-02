# rb_spec — Absorption Line Analysis Pipeline

[Back to Main Page](../../main_readme.md)

```python
from rbcodes.GUIs.rb_spec import rb_spec as r

spec = r.from_file('spectrum.fits', filetype='linetools')
spec.shift_spec(zabs=0.348)
spec.slice_spec(1548.2, -1500, 1500, use_vel=True)
spec.fit_continuum(mask=[-200, 200], Legendre=3)
spec.compute_EW(1548.2, vmin=-150, vmax=150, plot=True)
spec.save_slice('results.json', file_format='json')
```

---

## Loading Spectra

```python
# From file
spec = r.from_file('spectrum.fits', filetype='linetools')

# From numpy arrays
spec = r.from_data(wave, flux, error)
```

**Supported file types:** `'linetools'`, `'fits'`, `'ascii'`, `'p'` (pickle), `'xfits'`

```python
spec.shift_spec(zabs=0.511)
spec.slice_spec(2796.3, -1500, 1500, use_vel=True)  # velocity limits in km/s
```

---

## Continuum Fitting

### Interactive spline

```python
spec.fit_continuum(Interactive=True)
```

| Key / Click | Action |
|-------------|--------|
| Left-click | Add continuum point (median of ±2.5 px) |
| Right-click | Delete nearest point |
| `b` | Add point at exact cursor position |
| `Enter` | Fit spline through current points |
| `n` | Show normalized spectrum |
| `w` | Save continuum and exit |
| `q` | Quit without saving |

<img src="images/interactive_continuum.png" width="700" alt="Interactive Continuum Fitting">

### Interactive masking GUI

```python
spec.fit_continuum_interactive(order=3, use_weights=False, domain=[-1000, 1000])
```

| Key / Click | Action |
|-------------|--------|
| Left-click pairs | Define a mask region |
| Right-click pairs | Remove masks within selected range |
| `f` | Fit continuum |
| `a` | Accept and exit |
| `r` | Reset all masks |
| `z` | Undo last mask |
| `m` | Manual mask entry by velocity values |
| `c` / `Esc` | Cancel |

<img src="images/interactive_masking_gui.png" width="700" alt="Interactive Masking GUI">

### Polynomial fitting

```python
# Fixed order
spec.fit_continuum(mask=[-200, 200], domain=[-1500, 1500], Legendre=3)

# BIC-optimized order
spec.fit_continuum(mask=[-200, 200], domain=[-1500, 1500],
                   Legendre=True, optimize_cont=True,
                   min_order=1, max_order=7)
```

<img src="images/polynomial_fit_comparison.png" width="700" alt="Polynomial Fit Comparison">

### RANSAC fitting

```python
spec.fit_continuum_ransac(window=149)          # window-based
spec.fit_polynomial_ransac(degree=4, residual_threshold=0.1)
```

Inspect any continuum fit:

```python
spec.plot_continuum_fit(verbose=True)
```

<img src="images/continuum_fitting_example.png" width="700" alt="Continuum Fitting Example">

---

## Equivalent Width and Column Density

```python
spec.compute_EW(transition, vmin=-80, vmax=55, plot=True, SNR=True, _binsize=3)

print(f"W_λ   = {spec.W:.3f} ± {spec.W_e:.3f} Å")
print(f"log N = {spec.logN:.2f} ± {spec.logN_e:.2f}")
print(f"v_cen = {spec.vel_centroid:.1f} km/s")
print(f"SNR   = {spec.SNR:.1f}")
```

Saturation:

```python
result = spec.compute_EW(transition, vmin=-100, vmax=100)
if result['line_saturation']:
    print(f"Saturated — {result['saturation_fraction']:.1%} of pixels affected")
```

<img src="images/ew_measurement.png" width="700" alt="Equivalent Width Measurement">

---

## Visualization

```python
spec.plot_continuum_fit(verbose=True)
spec.plot_spec()
spec.plot_slice()
spec.plot_doublet(1548.2, 1550.8, vmin=-600, vmax=600)
```

<img src="images/spectrum_visualization.png" width="700" alt="Spectrum Visualization">

---

## Saving and Loading

```python
# Save
spec.save_slice('results.json', file_format='json')  # preferred
spec.save_slice('results.p')                          # pickle

# Load
from rbcodes.GUIs.rb_spec import load_rb_spec_object
spec = load_rb_spec_object('results.json')
```

---

## Output Attributes

| Attribute | Description |
|-----------|-------------|
| `wave_slice` | Sliced wavelength array |
| `flux_slice` | Sliced flux array |
| `error_slice` | Sliced error array |
| `velo` | Velocity array (km/s) |
| `cont` | Fitted continuum |
| `fnorm` | Normalized flux |
| `enorm` | Normalized error |
| `cont_mask` | Velocity masks used for continuum fitting |
| `continuum_fit_params` | Method, order, and fit metadata |
| `W` | Rest-frame equivalent width (Å) |
| `W_e` | EW uncertainty (Å) |
| `logN` | log AOD column density |
| `logN_e` | log column density uncertainty |
| `vel_centroid` | EW-weighted velocity centroid (km/s) |
| `vel_disp` | Velocity dispersion (km/s) |
| `SNR` | Signal-to-noise ratio |

---

## Examples

### CGM absorption — multiple systems

```python
for zabs in [0.511, 1.026, 1.564]:
    spec.shift_spec(zabs)
    spec.slice_spec(2796.3, -1500, 1500, use_vel=True)
    spec.fit_continuum_interactive(order=3)
    spec.compute_EW(2796.3, vmin=-100, vmax=100)
    spec.save_slice(f'mgii_z{zabs:.3f}.json', file_format='json')
```

### Lyman series

```python
for lam in [1215.67, 1025.72, 972.54]:   # Lyα, Lyβ, Lyγ
    spec.shift_spec(2.354)
    spec.slice_spec(lam, -2000, 2000, use_vel=True)
    spec.fit_continuum_interactive(order=4, use_weights=True)
    spec.compute_EW(lam, vmin=-300, vmax=300)
```

### Batch processing

```python
results = []
for file in ['qso1.fits', 'qso2.fits', 'qso3.fits']:
    spec = r.from_file(file, filetype='linetools')
    spec.shift_spec(zabs=0.5)
    spec.slice_spec(1548.2, -1000, 1000, use_vel=True)
    spec.fit_continuum(domain=[-1000, 1000], Legendre=3)
    spec.compute_EW(1548.2, vmin=-200, vmax=200)
    results.append({'file': file, 'EW': spec.W, 'logN': spec.logN})
```
