# zfind Template Files

Template FITS files used by `template_search()` (Phase 7) and `pca_search()` (Phase 8).

---

## `marz/` — 1D spectral templates

**Status: to be populated before Phase 7**

MARZ (Samreay et al.) stores templates as JavaScript arrays — no FITS files exist in the
repository. Templates must be extracted from the JavaScript source.

**To extract MARZ templates:**

```python
# extract_marz_templates.py — run once
# Clone MARZ: git clone https://github.com/Samreay/Marz.git /tmp/marz
# Templates are embedded in /tmp/marz/js/templates.js as float32 arrays.
# Parse the JS, extract wave + flux arrays, save as FITS.

# Template names to extract:
#   'Early Type Absorption Galaxy'  → passive.fits
#   'Late Type Emission Galaxy'     → elg.fits
#   'Intermediate Type Galaxy'      → intermediate.fits
#   'High Redshift Star Forming Galaxy' → lbg.fits
#   'Quasar'                        → qso.fits
```

**Alternatively:** use SDSS eigenspectra (see `sdss/` section below) — actual FITS files,
broader wavelength coverage.

**Template FITS format (for all sources):**
```
HDU 0: PRIMARY  — minimal header
HDU 1: WAVE     — rest-frame wavelength array (Angstrom, vacuum)
HDU 2: FLUX     — template flux array (arbitrary units, normalised to median=1)
```

---

## `sdss/` — SDSS eigenspectra

**Status: to be populated before Phase 7**

SDSS spectroscopic eigenspectra cover rest-frame ~1000–10000 Å.
Available from SDSS data release SAS server.

Files needed:
- `galaxy_eigen.fits` — galaxy PCA eigenspectra (~4 components)
- `qso_eigen.fits`    — QSO PCA eigenspectra (~4 components)

**Format:** multi-extension FITS, each extension is one eigenvector (same wave grid).

---

## `pca/` — Redrock PCA eigenvectors

**Status: to be populated before Phase 8**

Extracted from the desihub/redrock-templates HDF5 files.
Run `extract_redrock_templates.py` once to generate these.

```bash
# Download redrock templates:
pip install redrock-templates   # or clone desihub/redrock-templates
# Then run:
python extract_redrock_templates.py
```

**FITS format (one file per object class):**
```
HDU 0: PRIMARY  — header: n_components, wavelength range, source info
HDU 1: WAVE     — rest-frame wavelength grid (Angstrom, vacuum)
HDU 2: COMP0    — first eigenvector
HDU 3: COMP1    — second eigenvector
...
```

Files: `galaxy_pca.fits`, `qso_pca.fits`
