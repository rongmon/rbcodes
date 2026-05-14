# rb_zfind — Semi-Automated Redshift & Absorber Finder
## Implementation Plan

**Status:** Design complete, not yet implemented  
**Branch:** `zfind-feature` (create from master before starting)  
**Rule:** test every phase standalone before moving to the next  

---

## What This Is

A self-contained GUI tool + engine for:
1. **Emission mode** — find galaxy/QSO redshift from 1D spectrum (integrates with rb_zgui)
2. **Absorption mode** — find intervening absorber systems along a QSO sightline (integrates with rb_multispec)

Both modes share identical IO and engine architecture. Only the adapters and dialog callers differ.

---

## Package Location

```
rbcodes/GUIs/zfind/
├── __init__.py
├── main.py          ← CLI entry point: rb_zfind command
├── engine.py        ← ALL computation. No Qt, no matplotlib, no GUI imports ever.
├── dialog.py        ← PyQt5 popup. mode='emission' or mode='absorption' parameter.
├── io.py            ← Input normalization + ZFindResult / AbsorberResult dataclasses
├── adapters.py      ← Thin translators: ZFindResult→zgui, AbsorberResult→multispec
└── templates/
    ├── README.md    ← sources, licenses, how to regenerate PCA files
    ├── marz/        ← MARZ 1D template FITS (MIT, from Samreay/Marz on GitHub)
    │   ├── ELG.fits          (emission-line galaxy)
    │   ├── Passive.fits      (absorption/elliptical galaxy)
    │   ├── QSO.fits          (broad-line QSO)
    │   ├── LRG.fits          (luminous red galaxy)
    │   └── LBG.fits          (Lyman-break galaxy)
    ├── sdss/        ← SDSS eigenspectra FITS (wider rest-frame coverage than MARZ)
    │   ├── galaxy_eigen.fits
    │   └── qso_eigen.fits
    └── pca/         ← redrock eigenvectors extracted to FITS (no redrock dep needed)
        ├── galaxy_pca.fits   (HDU0: header, HDU1: WAVE, HDU2+: eigenvectors)
        ├── qso_pca.fits
        └── extract_redrock_templates.py  ← run once manually to generate these
```

**Template download (do before coding):**
```bash
git clone --depth 1 --filter=blob:none --sparse https://github.com/Samreay/Marz.git /tmp/marz
cd /tmp/marz
git sparse-checkout set web/data/templates
cp web/data/templates/*.fits rbcodes/GUIs/zfind/templates/marz/
```
Start with: ELG, Passive, QSO, LRG, LBG. Stars optional for v1.

---

## IO Design (`io.py`)

### Input: all three forms accepted, all normalized to `rb_spectrum`

```python
# All equivalent inside the engine:
to_rb_spectrum('spectrum.fits')           # filepath string
to_rb_spectrum(spec)                      # rb_spectrum object
to_rb_spectrum(wave, flux)               # raw arrays, minimal
to_rb_spectrum(wave, flux, error)        # with error
to_rb_spectrum(wave, flux, error, cont)  # with continuum
to_rb_spectrum((wave, flux, error))      # tuple form
```

`to_rb_spectrum()` normalizes everything. Engine always receives `rb_spectrum`.

### Air/vacuum: automatic
- `spec.meta.get('airvac', 'vac')` — default vacuum
- If `'air'`: call `spec.air2vac()` (already implemented in rb_spectrum)
- Dialog shows yellow warning label if header has no wavelength frame info

### Continuum handling (engine internal, NEVER displayed in dialog)
Priority order:
1. `spec.co_is_set` → use `spec.co.value` directly (user provided, don't override)
2. Checkbox "Fit continuum" checked → `fit_optimal_polynomial(wave, flux, error, use_weights=False, silent=True, plot=False)` from `rbcodes.IGM.rb_iter_contfit` — BIC selects polynomial order automatically
3. Checkbox unchecked (skip) → `continuum = np.zeros_like(flux)` for emission mode (not needed — looking for peaks above noise, not above continuum)

### Output dataclasses

```python
# Emission mode
@dataclass
class ZSolution:
    z            : float
    z_err        : float        # from chi2 curvature: sqrt(1/|d²chi2/dz²|)
    chi2_dof     : float
    method       : str          # 'LineSearch', 'Template:Passive', 'PCA:Galaxy'
    template_type: str          # 'Galaxy', 'QSO', 'Star'
    n_features   : int          # lines or pixels used in fit

@dataclass
class ZFindResult:
    z_array     : np.ndarray
    chi2_curves : list          # [{'label': 'LineSearch:ISM', 'chi2': ndarray}, ...]
    solutions   : list          # list of ZSolution, sorted by chi2_dof ascending
    input_spec  : rb_spectrum   # stored for re-plotting, no recomputation needed
    warnings    : list          # e.g. ['No error array — using MAD-STD IVAR']

# Absorption mode  
@dataclass
class AbsorberCandidate:
    z            : float
    significance : float        # sigma above noise floor
    n_lines      : int
    is_doublet   : bool         # doublet match = far more reliable than single line
    linelist_name: str          # 'CIV', 'MgII', 'SiIV'
    lines_matched: list         # e.g. ['CIV 1548', 'CIV 1550']

@dataclass
class AbsorberResult:
    z_array           : np.ndarray
    significance_curve: np.ndarray   # significance vs z (not chi2)
    candidates        : list         # AbsorberCandidate list, sorted by significance desc
    input_spec        : rb_spectrum
    warnings          : list
```

Both objects are self-contained — carry everything needed to re-plot or pass downstream.

---

## Engine Design (`engine.py`)

**Rule: never import PyQt5, matplotlib, or any GUI library here.**

### Shared preprocessing (runs once before z scan)

```python
def _preprocess(spec, fit_continuum=True):
    # 1. Air→vac if needed
    if spec.meta.get('airvac', 'vac') == 'air':
        spec = spec.copy(); spec.air2vac()
    wave = spec.wavelength.to(u.AA).value
    flux = spec.flux.value
    
    # 2. IVAR
    if spec.sig_is_set:
        ivar = 1.0 / spec.sig.value**2
        ivar[~np.isfinite(ivar)] = 0.0
    else:
        sigma = 1.4826 * median_abs_deviation(flux, nan_policy='omit')
        ivar = np.full_like(flux, 1.0 / sigma**2)
        # flag this in warnings
    
    # 3. Continuum
    if spec.co_is_set:
        continuum = spec.co.value
    elif fit_continuum:
        from rbcodes.IGM.rb_iter_contfit import fit_optimal_polynomial
        result = fit_optimal_polynomial(wave, flux, 
                     error=spec.sig.value if spec.sig_is_set else None,
                     use_weights=False, silent=True, plot=False)
        continuum = result['continuum']
    else:
        continuum = np.zeros_like(flux)
    
    return wave, flux, ivar, continuum
```

### `line_search(spec, linelist_df, z_min, z_max, n_steps, mode, fwhm_ang, fit_continuum)`

```
mode = 'emission' or 'absorption'
```

**Algorithm:**
1. Preprocess → wave, flux, ivar, continuum
2. If `fwhm_ang > 0`: build Gaussian LSF kernel (for emission line broadening)
3. z grid: `z_array = np.linspace(z_min, z_max, n_steps)`
4. For each z:
   - `obs_wave = linelist_df['wave'].values * (1 + z)`
   - `in_range = obs_wave[(obs_wave > wave.min()) & (obs_wave < wave.max())]`
   - if `len(in_range) == 0`: `chi2[i] = np.nan`; continue
   - For each in-range line: extract window ±W pixels around `obs_wave`
     - **emission:** fit Gaussian amplitude A analytically on `flux - continuum`; accumulate chi2
     - **absorption:** measure `continuum - flux` decrement; accumulate significance
   - Normalize by `len(in_range)` so sparse coverage doesn't penalize
5. Top-N minima: find local minima with min z-separation 0.01; extract z_err from curvature
6. Return `ZFindResult` (emission) or `AbsorberResult` (absorption)

**Absorption mode extra notes:**
- Use doublet linelists preferentially (CIV 1548/1550, MgII 2796/2803, SiIV 1393/1402)
- Doublet match: significance is product of both line detections — much harder to fake
- Restrict z_max < z_QSO (user-provided; no absorbers beyond the background source)
- Output is `AbsorberResult` with candidate list, NOT a single z
- This output goes to rb_multispec adapter only — NOT connected to zgui dialog

### `template_search(spec, template_name, template_set, z_min, z_max, n_steps, fwhm_ang, fit_continuum)`

1. Load template FITS from `templates/marz/` or `templates/sdss/`
2. If `fwhm_ang > 0`: convolve template with Gaussian
3. At each z: `template_wave_shifted = template_wave * (1+z)`, interpolate onto observed grid
4. Fit amplitude: `A = Σ(flux·template·ivar) / Σ(template²·ivar)` — overlap region only
5. `chi2 = Σ((flux - A·template)² · ivar)` — overlap region only
6. Top-N minima + z_err from curvature → `ZFindResult`

### `pca_search(spec, template_set, z_min, z_max, n_steps, fwhm_ang, fit_continuum)`

1. Load PCA FITS: `pca_wave` (HDU1) + eigenvectors `E` (HDU2+)
2. At each z: shift + interpolate each component onto observed grid → matrix `M` (n_obs × n_comp)
3. Solve: `A = (MᵀΛM)⁻¹ Mᵀ Λ flux` where `Λ = diag(ivar)` — `np.linalg.solve`, fast
4. `chi2 = Σ((flux - M@A)² · ivar)`
5. Top-N minima + z_err → `ZFindResult`

---

## Adapters (`adapters.py`)

**zgui adapter (emission mode → zgui signals):**
```python
def zfind_to_zgui_z(result: ZFindResult, idx=0) -> list:
    s = result.solutions[idx]
    return [s.z, s.z_err]    # matches widget_z._on_estZ_changed() exactly
```

**rb_multispec adapter (absorption mode → AbsorberManager format):**
```python
def absorbers_to_multispec(result: AbsorberResult, accepted_indices: list) -> list:
    # AbsorberManager expects list of {'zabs': float, 'name': str, 'label': str}
    return [{'zabs': c.z, 'name': c.linelist_name, 'label': f'z={c.z:.4f}'}
            for i, c in enumerate(result.candidates) if i in accepted_indices]

def zfind_to_multispec_z(result: ZFindResult, idx=0) -> float:
    return result.solutions[idx].z    # for single-z update in RedshiftInputWidget
```

---

## Dialog (`dialog.py`)

Single `ZFindDialog` class. `mode` parameter controls which UI is shown.

**Instantiation:**
```python
# From zgui (emission mode):
ZFindDialog(spec, mode='emission', default_linelist='ISM').exec_()

# From rb_multispec (absorption mode, future):
ZFindDialog(spec, mode='absorption', z_qso=2.5, default_linelist='CIV').exec_()
```

**Type:** `QDialog.exec_()` — modal. Main window stays visible behind it.

**Layout:**
```
┌──────────────────────────────────────────────────────────────┐
│  Redshift Estimator                                    [X]   │
│                                                              │
│  Mode: ☑ Line Search  ☐ Template Match  ☐ PCA             │
│  Linelist: [ISM ▼]  Emission ● / Absorption ○  (mode lock) │
│  z range: [0.0] to [6.0]   FWHM: [0 Å]  R: [0]           │
│  Wav frame: [Vacuum ▼]  ☑ Fit continuum                    │
│  [Run ▶]  ░░░░░░░ progress                                  │
│                                                              │
│  [yellow label if wavelength frame ambiguous]               │
│                                                              │
│  ┌─ Spectrum + line overlay at selected z ───────────────┐  │
│  │  NavigationToolbar2QT (pan/zoom/home)                 │  │
│  │  Smooth: [1 ▲▼]  (display only, Gaussian kernel)     │  │
│  └────────────────────────────────────────────────────────┘  │
│                                                              │
│  ┌─ Chi-square vs z  (click to select z) ───────────────┐  │
│  │  Vertical dashed line = selected z                   │  │
│  │  Multiple colored curves if multiple modes ran       │  │
│  │  Minima auto-labeled: "z=1.23 Galaxy"                │  │
│  └────────────────────────────────────────────────────────┘  │
│                                                              │
│  EMISSION:  z=1.2345  χ²=1.12  Galaxy  ← click row         │
│             z=0.6123  χ²=1.45  QSO     ← click row         │
│             z=2.8901  χ²=1.89  Galaxy  ← click row         │
│                                                              │
│  ABSORPTION (future):  table of candidate systems           │
│  with checkboxes, significance, doublet flag                │
│                                                              │
│  Selected: z [________]  z_err [________]  χ²/dof [______] │
│                                                              │
│  [Accept → Send to caller]                    [Cancel]      │
└──────────────────────────────────────────────────────────────┘
```

**Signals:**
```python
accepted_z         = pyqtSignal(list)   # [z, z_err] → zgui widget_z._on_estZ_changed()
accepted_absorbers = pyqtSignal(list)   # [{zabs,name,...}] → rb_multispec AbsorberManager
```

**QThread workers:** one thread per selected engine mode; curves appear on chi-square plot as each completes. Progress shown via spinner or progress bar.

**Chi-square plot click:** `mpl_connect('button_press_event')` on chi-square axes → read x-coordinate as z → update z field + redraw line overlay on spectrum panel.

---

## CLI (`main.py`)

```python
def launch_zfind(spec_or_file=None, linelist=None, mode='emission'):
    app = QApplication.instance() or QApplication(sys.argv[:1])
    spec = io.to_rb_spectrum(spec_or_file) if spec_or_file else None
    ZFindDialog(spec=spec, mode=mode, default_linelist=linelist).exec_()

def main():
    parser = argparse.ArgumentParser(description='rb_zfind: Redshift and absorber finder')
    parser.add_argument('fitsfile', nargs='?', help='FITS file to analyze')
    parser.add_argument('-l', '--linelist', default=None, help='Default linelist')
    parser.add_argument('-m', '--mode', default='emission',
                        choices=['emission', 'absorption'], help='Search mode')
    args = parser.parse_args()
    launch_zfind(args.fitsfile, linelist=args.linelist, mode=args.mode)
```

**pyproject.toml entry point:** `rb_zfind = rbcodes.GUIs.zfind.main:main`

**From Python/notebook:**
```python
from rbcodes.GUIs.zfind.main import launch_zfind
from rbcodes.utils.rb_spectrum import rb_spectrum

launch_zfind('spectrum.fits')                           # from file
launch_zfind(rb_spectrum.from_file('spec.fits'))        # from rb_spectrum
launch_zfind((wave, flux, error))                       # from tuple
launch_zfind('spec.fits', mode='absorption')            # absorption mode
```

---

## Changes to Existing Files

### `rbcodes/GUIs/zgui/menu_toolbars.py`
- Add one button: `"Find z"`
- Enabled only when spectrum is loaded (`fitsobj.flux is not None`)
- On click:
  ```python
  from rbcodes.GUIs.zfind.dialog import ZFindDialog
  from rbcodes.GUIs.zfind.io import to_rb_spectrum
  # prefer filepath (gets continuum too); fallback to fitsobj arrays
  spec = to_rb_spectrum(current_filepath) if current_filepath else \
         to_rb_spectrum(self.fitsobj.wave, self.fitsobj.flux, self.fitsobj.error)
  dlg = ZFindDialog(spec, mode='emission',
                    default_linelist=main_window.widget_z.l_lln.currentText())
  dlg.accepted_z.connect(main_window.widget_z._on_estZ_changed)
  dlg.exec_()
  ```

### `rbcodes/GUIs/zgui/main.py`
- No changes needed — connection is made inline when dialog is created (see above)

### `rbcodes/GUIs/multispecviewer/multispec.py` (FUTURE — separate session)
- Add "Find Absorbers" button
- Same pattern but `mode='absorption'`, connect `accepted_absorbers` to `AbsorberManager`

### `pyproject.toml`
- Add: `rb_zfind = rbcodes.GUIs.zfind.main:main`

---

## Phase Plan — Test Every Step Before Proceeding

### Phase 0 — Setup (before any code)
- [ ] Create branch: `git checkout -b zfind-feature`
- [ ] Download MARZ templates into `templates/marz/` (see download command above)
- [ ] Inspect template FITS structure (read headers, check wavelength/flux extensions)
- [ ] Create `rbcodes/GUIs/zfind/__init__.py` and empty stubs for each file
- [ ] Add `rb_zfind` entry point to `pyproject.toml`

### Phase 1 — IO Layer (`io.py`)
- [ ] Implement `to_rb_spectrum()` — all three input forms
- [ ] Implement `ZFindResult` and `ZSolution` dataclasses
- [ ] Implement `AbsorberResult` and `AbsorberCandidate` dataclasses
- [ ] **Test:** round-trip each input form (filepath, rb_spectrum, arrays)
- [ ] **Test:** `ZFindResult` and `AbsorberResult` create/access correctly

### Phase 2 — Engine: `line_search` emission mode (`engine.py`)
- [ ] Implement `_preprocess()` — IVAR fallback, continuum via `rb_iter_contfit`, air→vac
- [ ] Implement `line_search()` emission mode
- [ ] **Test standalone (no GUI):** synthetic spectrum with 3 known emission lines at z=1.5
  - Verify chi2 minimum is at z=1.5
  - Verify z_err is reasonable
  - Verify degenerate second minimum appears for sparse line case
- [ ] **Test:** MAD-STD IVAR fallback (pass flux-only spec)
- [ ] **Test:** narrow wavelength coverage (only 2 lines in range)
- [ ] **Test:** `fit_continuum=False` skip option

### Phase 3 — Engine: `line_search` absorption mode (`engine.py`)
- [ ] Implement absorption mode branch in `line_search()`
- [ ] **Test standalone:** synthetic QSO spectrum + 2 known absorbers at different z
  - Verify chi2/significance curve shows two minima at correct positions
  - Verify doublet matching gives higher significance than single line
- [ ] **Test:** `rb_iter_contfit` continuum used correctly for absorption depth measurement

### Phase 4 — Adapters (`adapters.py`)
- [ ] Implement `zfind_to_zgui_z()`
- [ ] Implement `absorbers_to_multispec()`
- [ ] Implement `zfind_to_multispec_z()`
- [ ] **Test:** adapter output matches expected zgui signal format `[z, z_err]`
- [ ] **Test:** adapter output matches rb_multispec AbsorberManager dict format

### Phase 5 — Dialog emission mode (`dialog.py`)
- [ ] Build `ZFindDialog` with emission mode UI
- [ ] Spectrum panel: matplotlib canvas + NavigationToolbar + Smooth spinbox
- [ ] Chi-square panel: clickable, vertical dashed line, minima labels
- [ ] Controls: mode checkboxes (Line Search only for now), linelist dropdown,
      z range inputs, FWHM/R inputs, continuum checkbox, Run button, wav frame dropdown
- [ ] QThread worker wiring for `line_search`
- [ ] **Test:** open dialog with a test spectrum, click Run, verify chi2 plot appears
- [ ] **Test:** click chi2 plot → z field updates + line overlay redraws
- [ ] **Test:** click solution row → same update
- [ ] **Test:** Accept button → `accepted_z` signal emitted with correct `[z, z_err]`

### Phase 6 — zgui integration
- [ ] Add "Find z" button to `menu_toolbars.py`
- [ ] Wire `accepted_z` signal → `widget_z._on_estZ_changed()`
- [ ] **Test end-to-end:** open zgui → load spectrum → click "Find z" → run line search →
      accept z → verify estZ field populated and line overlay appears in main canvas
- [ ] **Test:** CLI `rb_zfind spectrum.fits` opens dialog standalone
- [ ] **Test:** `launch_zfind((wave, flux, error))` from Python

### Phase 7 — Template Search (`engine.py` + `dialog.py`)
- [ ] Implement `template_search()` in engine
- [ ] **Test standalone:** known SDSS galaxy spectrum vs Passive/ELG templates
- [ ] Add Template Match checkbox + template set dropdown to dialog
- [ ] **Test:** both Line Search + Template run simultaneously, two curves on chi2 plot

### Phase 8 — PCA Search (`engine.py` + `dialog.py`)
- [ ] Extract redrock eigenvectors → FITS files in `templates/pca/`
      (run `extract_redrock_templates.py` manually once)
- [ ] Implement `pca_search()` in engine
- [ ] **Test standalone:** compare PCA result vs template_search on same spectrum
- [ ] Add PCA checkbox + type dropdown to dialog
- [ ] **Test:** all three modes running simultaneously

### Phase 9 — rb_multispec absorption integration (SEPARATE SESSION)
- [ ] Add "Find Absorbers" button to `multispec.py`
- [ ] Wire `accepted_absorbers` signal → `AbsorberManager`
- [ ] **Test:** QSO spectrum → absorption mode → confirm multiple absorbers detected,
      accepted systems populate rb_multispec absorber overlays

---

## Notes for Future Claude Sessions

- **Read this file first**, then `git log --oneline -20` to see what phases are done
- Engine (`engine.py`) must NEVER import PyQt5, matplotlib, or GUI libraries
- All engine functions accept `rb_spectrum` — never raw numpy arrays directly
- `rb_spectrum` is in `rbcodes/utils/rb_spectrum.py`
  - `.wavelength.to(u.AA).value` — wavelength
  - `.flux.value` — flux
  - `.sig` — error (None if absent)
  - `.co` — continuum (None if absent)
  - `.sig_is_set`, `.co_is_set` — boolean properties
  - `.air2vac()` — converts in place, already implemented
  - `.meta['airvac']` — 'vac' or 'air'
- `rb_iter_contfit` is in `rbcodes/IGM/rb_iter_contfit.py`
  - Use `fit_optimal_polynomial(wave, flux, error, use_weights=False, silent=True, plot=False)`
  - Returns dict with `result['continuum']`
- MARZ templates in `templates/marz/` — simple 1D FITS, rest-frame wavelength + flux
- `ZFindResult.chi2_curves` is a list of dicts `{'label': str, 'chi2': ndarray}`
  so multiple engine modes can be overlaid on the same chi2 plot
- Absorption mode output (`AbsorberResult`) is NOT connected to zgui dialog
  It goes to rb_multispec adapter only (Phase 9)
- Only one existing file changes for zgui integration: `menu_toolbars.py` (one button)
