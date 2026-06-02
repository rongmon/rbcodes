# rb_multispec

Multi-spectrum viewer for visualizing and analyzing multiple spectra simultaneously —
line identification, redshift overlay, Gaussian fitting, and velocity-space analysis.

```bash
rb_multispec                            # launch empty
rb_multispec file.fits                  # load one file
rb_multispec file1.fits file2.fits      # load multiple files
rb_multispec *.fits                     # wildcard
rb_multispec -c file.fits               # classic mode
rb_multispec --examples                 # detailed usage examples
```

<img src="images/main_interface.png" width="700" alt="Main Interface">

---

## Navigation

| Key | Action |
|-----|--------|
| `x` / `X` | Set left / right x-limit to cursor |
| `t` / `b` | Set y-max / y-min to cursor (per panel) |
| `Y` | Enter y-limits manually |
| `[` / `]` | Pan left / right one viewport width |
| `o` | Zoom out x-axis |
| `r` | Reset view |
| `R` | Clear all line IDs and fit overlays |
| `a` | Autoscale y of all panels to visible x range |
| `A` | Autoscale y of current panel only |
| `L` | Toggle line labels on/off |
| `S` / `U` | Increase / decrease smoothing |
| `h` / `H` / `?` | Open tabbed help dialog |
| `q` | Quit |

---

## Quick line identification

| Key | Lines |
|-----|-------|
| `C` | CIV doublet (1548, 1550 Å) |
| `M` | MgII doublet (2796, 2803 Å) |
| `F` | FeII multiplet (2600, 2586, 2382 Å) |
| `6` | OVI doublet (1031, 1037 Å) |
| `4` | SiIV doublet (1393, 1402 Å) |
| `8` | NeVIII doublet (778, 770 Å) |
| `2` | Lyβ/Lyα |
| `1` | Lyα/Lyβ |

Right-click anywhere to open a full line list at the clicked wavelength.

<img src="images/quick_id.png" width="700" alt="Quick Line Identification">
<img src="images/right_click_id.png" width="700" alt="Right-Click Line Identification">

---

## Quick line fitting

Two clicks define the window and continuum level (flux at each click sets the
continuum on that side — tilted continua are fully supported).

| Keys | Action |
|------|--------|
| `g` → `g` | Gaussian fit. Draws orange overlay; updates z if a line list is active. |
| `c` → `c` | Centre-of-mass centroid. Draws cyan marker; updates z. |
| `Z` | Undo last z (toggle back with a second `Z`) |
| `R` | Clear fit overlay and cancel any pending fit |

<img src="images/quick_fit_gaussian.png" width="700" alt="Quick Gaussian fit">
<img src="images/quick_fit_com.png" width="700" alt="Quick CoM centroid">

---

## Advanced fitting dialog (`G` → `G`)

Two `G` keypresses define the wavelength window. The dialog opens zoomed to
that region.

| Control | Description |
|---------|-------------|
| **Lines to fit** spinbox | Number of Gaussian components (1–5) |
| **Ion 1…N** dropdowns | Transitions to fit; Ion 1 auto-populates remaining slots |
| **z guess** field | Starting redshift; or use **Shift+C** on canvas to set from cursor |
| **Auto continuum** | Median of first/last N edge pixels |
| **Manual continuum** | `d` → `d` on canvas to set two anchor points |
| **Weight by flux errors** | Inverse-variance weighting |
| **Fit** button | Run the fit |
| **Advanced…** | Edit bounds, tie sigmas, fix z |

### Canvas keys inside dialog

| Key | Action |
|-----|--------|
| `Shift+C` | Set z_guess from cursor wavelength |
| `d` → `d` | Set manual continuum anchor points |
| `x` / `X` | Set x-limits |
| `[` / `]` | Pan |
| `t` / `b` | Set y-limits |
| `r` | Reset view |

### Result buttons

| Button | Action |
|--------|--------|
| **Copy to clipboard** | Copy results text |
| **Export to file** | Save results to `.txt` |
| **Apply z + linelist to main** | Push fitted z back to main GUI |
| **Add to absorbers** | Add to absorber manager |

<img src="images/advanced_fit_dialog.png" width="700" alt="Advanced Fit Dialog">

---

## Redshift and absorber manager

Enter redshift, select line list and color, click **Submit**. Click **Catalog**
to add to the absorber manager. Use checkboxes to toggle system visibility.

<img src="images/redshift_controls.png" width="700" alt="Redshift Controls">

---

## vStack — velocity-space analysis

| Key | Action |
|-----|--------|
| `v` | Launch vStack (default ±1000 km/s) |
| `V` | Launch vStack with custom velocity limits |
| `>` / `<` | Next / previous page of lines |
| `w` | Cycle flag: Non-Detection → Detection → Blended → Low-Confidence |
| `Y` | Set y-limits for current panel |
| `S` | Save flags and return to main display |

<img src="images/vstack.png" width="700" alt="vStack Interface">

---

## Action buttons

| Button | Action |
|--------|--------|
| Clear All Lines | Remove all line IDs |
| Load | Load saved line IDs and absorber systems |
| Save | Save current IDs and absorbers (JSON recommended) |
| Show | Toggle visibility of all lines |
| List | Show table of all identified lines |

<img src="images/line_list_dialog.png" width="700" alt="Line List Dialog">

---

## Reconcile line catalogs from multiple sessions

```python
from rbcodes.GUIs.multispecviewer.utils import reconcile_linelists

reconciled, absorbers = reconcile_linelists(
    ['session1.json', 'session2.json'],
    velocity_threshold=20,           # km/s
    output_file='combined.json',
    create_absorber_df=True
)
```

Lines with the same transition name within `velocity_threshold` km/s of the
anchor entry are merged (mean z, mean wavelength).

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| No error spectrum | Assumes 5% of flux |
| Lines not visible after applying z | Check line list covers your wavelength range |
| vStack not launching | Need a redshift and line list selected first |
| Absorbers missing after loading JSON | Check checkboxes in absorber manager |
| Quick fit "Only N finite pixels" | Fit window too narrow or in a masked region — widen it |
| Advanced Fit "Fit failed" | Relax bounds in **Advanced…** or adjust z_guess with Shift+C |
