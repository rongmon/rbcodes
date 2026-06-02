# interactive_continuum_fit

[Back to Main Page](../main_readme.md)

Interactive GUI for fitting a continuum to spectral data. Choose between
polynomial (Legendre) or spline methods with interactive masking.

```python
from rbcodes.GUIs.interactive_continuum_fit import launch_interactive_continuum_fit
from rbcodes.utils import rb_spectrum as r

sp = r.rb_spectrum.from_file('spectrum.fits')

result = launch_interactive_continuum_fit(
    wave=sp.wavelength.value,
    flux=sp.flux.value,
    error=sp.sig.value
)

continuum = result['continuum']
```

Also accessible from within `launch_specgui` via the continuum fitting controls.

---

## Workflow

1. Select method: **Polynomial** or **Spline** (tab or `p` / `s` keys)
2. Add masks over absorption/emission features to exclude from the fit
3. Set order or enable Auto Optimize (BIC), then click **Fit Continuum**
4. Check the normalized spectrum in the bottom panel
5. Click **Accept & Return** when satisfied

---

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `f` | Fit continuum |
| `a` | Accept and return |
| `r` | Reset masks |
| `R` | Reset everything (masks, spline points, fit) |
| `z` | Undo last mask |
| `m` | Manual mask entry |
| `p` | Switch to polynomial mode |
| `s` | Switch to spline mode |
| `b` | Add spline point at exact cursor position |
| `c` | Clear all spline points |
| `h` | Help dialog |
| `+` / `-` | Zoom in / out |
| `0` | Reset zoom |
| `Esc` | Cancel / close without saving |

---

## Masks

| Interaction | Action |
|-------------|--------|
| Left-click × 2 | Define mask region (start → end) |
| Right-click × 2 | Remove masks overlapping this range |
| Auto-Mask button | Detect features automatically (sigma threshold) |
| Manual Entry button | Enter mask boundaries numerically |

Overlapping masks are merged automatically.

---

## Polynomial options

| Option | Description |
|--------|-------------|
| Order | Polynomial degree — start low (2–3), increase only if needed |
| Use Weights | Weight fit by error spectrum |
| Auto Optimize | BIC-based order selection over a specified range |

## Spline options

| Option | Description |
|--------|-------------|
| Degree | Polynomial piece order (3 = cubic, standard) |
| Smoothing | 0 = exact interpolation; higher = smoother curve |
| Median Window | Half-width used when placing points by left-click |

Place anchor points on true continuum regions, not on spectral features.
Requires at least 3 points.
