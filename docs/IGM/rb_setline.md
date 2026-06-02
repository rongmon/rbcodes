# rb_setline

[Back to Main Page](../main_readme.md)

Look up atomic transition data by wavelength or species name.

```python
from rbcodes.IGM.rb_setline import rb_setline

# Closest line to a wavelength
result = rb_setline(2796.3, 'closest')
print(result['name'], result['wave'], result['fval'])

# Exact wavelength match
result = rb_setline(1215.67, 'Exact')

# By species name
result = rb_setline(0, 'Name', target_name='HI 1215')

# Different line list
result = rb_setline(1550.0, 'closest', linelist='LBG')
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `lambda_rest` | float | Rest wavelength to match (Å) — ignored when `method='Name'` |
| `method` | str | `'closest'`, `'Exact'`, or `'Name'` |
| `linelist` | str | Line list to search (default: `'atom'`) |
| `target_name` | str | Species name — required when `method='Name'` |

## Available line lists

| List | Contents |
|------|---------|
| `atom` | Full atomic list with oscillator strengths and damping constants |
| `LLS` | Lines common in Lyman Limit Systems |
| `DLA` | Lines common in Damped Lyman-alpha systems |
| `LBG` | Lines common in Lyman Break Galaxies |
| `Gal` | Galaxy vacuum wavelength lines |
| `Gal_Em` | Galaxy emission lines |
| `Gal_Abs` | Galaxy absorption lines |
| `Gal_long` | Extended galaxy emission and absorption catalog |
| `AGN` | Active Galactic Nuclei lines |
| `Eiger_Strong` | Strong lines in the Eiger catalog |
| `HI_recomb` | Hydrogen recombination lines |
| `HI_recomb_light` | Subset of hydrogen recombination lines |

## Output keys

| Key | Description |
|-----|-------------|
| `wave` | Rest wavelength(s) (Å) — scalar for `'closest'`, array otherwise |
| `fval` | Oscillator strength(s) |
| `name` | Species name(s) |
| `gamma` | Radiation damping constant — `atom` list only |
