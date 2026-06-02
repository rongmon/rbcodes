# rb_cont

[Back to Main Page](../main_readme.md)

> **This tool is superseded by [`interactive_continuum_fit`](interactive_continuum_fit.md)**, which provides the same spline fitting plus polynomial fitting, interactive masking, and BIC optimization. Use that instead.

If you still need the legacy spline-only fitter:

```bash
python /path/to/rbcodes/GUIs/rb_cont.py spectrum.fits
```

```python
from rbcodes.GUIs.rb_fit_interactive_continuum import rb_fit_interactive_continuum

fitter = rb_fit_interactive_continuum(wave, flux, error)
cont = fitter.cont   # available after pressing 'w' to save
```

| Key | Action |
|-----|--------|
| Left-click | Add continuum point (median ±2.5 px) |
| Right-click | Delete nearest point |
| `b` | Add point at exact cursor position |
| `Enter` | Fit spline |
| `n` | Show normalized spectrum |
| `w` | Save and exit |
| `q` | Quit without saving |

Output is written to `{name}_norm.fits`.
