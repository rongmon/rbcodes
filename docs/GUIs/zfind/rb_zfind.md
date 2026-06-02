# rb_zfind — Semi-Automated Redshift Finder

[Back to Main Page](../../main_readme.md) | [Back to GUIs](../GUIs_readme.md)

Chi2-based redshift search tool. Runs three search methods against a spectrum, ranks solutions, and feeds the best redshift back into `rb_zgui` or `rb_multispec`.

---

## Quick start

```bash
rb_zfind                                       # open empty dialog
rb_zfind spectrum.fits                         # load file and open
rb_zfind spectrum.fits -l zfind_galaxy         # with default linelist
```

```python
from rbcodes.GUIs.zfind.main import launch_zfind
launch_zfind('spectrum.fits', linelist='zfind_em')
```

---

## Search methods

| Method | Description | Best for |
|--------|-------------|---------|
| **PCA** | Chi2 scan against DESI Redrock-style PCA eigenvectors | Galaxies, QSOs with good S/N |
| **Template** | Chi2 scan against 1D FITS templates (MARZ set) | Galaxies, passive/early-type |
| **Picket Fence** | Matched-filter chi2 using a curated linelist | Emission-line galaxies, quick scans |

---

## Dialog layout

**Top panel** — chi2 vs redshift. Click anywhere to jump to that z and update the spectrum overlay.

**Bottom panel** — spectrum with linelist overlaid at the selected z.

**Solutions table** — top-N ranked solutions with z, z_err, chi2. Click a row to select; **Accept z** sends it to the parent GUI.

---

## Parameters

| Parameter | Description |
|-----------|-------------|
| **Method** | PCA / Template / Picket Fence |
| **z min / z max** | Redshift search range; **↺** resets to template defaults |
| **n steps** | Redshift grid points |
| **Overlay lines** | Linelist drawn on spectrum at selected z |
| **λ min / λ max** | Mask pixels outside this observed-frame range (blank = full) |
| **Resolution** | Instrument LSF — R, FWHM(Å/km/s/pix); 0 = no convolution |
| **Smooth (pix)** | Boxcar pre-smoothing before chi2 scan; 1 = off |

**Picket Fence only:** Linelist, Win (pix), Data norm (`subtract`/`normalize`/`raw`), Mode (`Direct` or `Detect+Match`).

---

## Curated linelists

| Name | Contents |
|------|---------|
| `zfind_em` | Galaxy nebular emission: Lyα, [OII], Hβ, [OIII], Hα, [NII], [SII] |
| `zfind_galaxy` | Emission + stellar absorption (one-stop galaxy preset) |
| `zfind_stellar` | Stellar absorption lines |
| `zfind_igm` | IGM/CGM absorbers: OVI, HI, SiIV, CIV, MgII |
| `zfind_qso` | QSO/AGN broad emission lines |

All `rb_setline` linelists (`Gal_Em`, `LLS`, `DLA`, `AGN`, etc.) are also available.

---

## Spectrum panel keyboard shortcuts

| Key | Action |
|-----|--------|
| `x` / `X` | Set left / right x limit to cursor wavelength |
| `t` / `b` | Set top / bottom y limit |
| `r` | Reset axes to full range |
| `[` / `]` | Pan left / right |
| `S` / `U` | Increase / decrease display smoothing |
| `a` | Autoscale y to visible x range |

---

## Integration with rb_zgui

Click **Find z** in the `rb_zgui` toolbar. After running a search, select a solution and click **Accept z** — the redshift is passed directly to `rb_zgui`'s redshift widget.

---

## Python API

```python
from rbcodes.GUIs.zfind.io import to_rb_spectrum
from rbcodes.GUIs.zfind.engine import line_search
from rbcodes.GUIs.zfind.linelists import get_curated_df

spec   = to_rb_spectrum('spectrum.fits')
df     = get_curated_df('zfind_em')
result = line_search(spec, df, z_min=0.0, z_max=2.0, mode='emission')
print(result.solutions[0].z, result.solutions[0].z_err)
```
