# rb_zfind — Semi-Automated Redshift & Absorber Finder
## Status & TODO

**Phases 0–8 complete.** Tool is functional and integrated with rb_zgui.
Branch: `zfind-feature`.

---

## File Structure (current + planned)

```
rbcodes/GUIs/zfind/
├── __init__.py
├── main.py            ← CLI: rb_zfind command + launch_zfind()
├── engine.py          ← thin wrappers: preprocess → delegate to core → package result
├── dialog.py          ← ZFindDialog PyQt5 popup
├── io.py              ← ZFindResult, ZSolution, AbsorberResult, AbsorberCandidate
├── adapters.py        ← zfind_to_zgui_z, absorbers_to_multispec
├── linelists.py       ← 5 curated presets with weights + type: zfind_em, zfind_stellar,
│                          zfind_igm, zfind_galaxy, zfind_qso  ✓ DONE
├── picket_fence.py    ← [NEW] PicketFenceZ — weighted line SNR scan, plain arrays
├── template_fit.py    ← [NEW] TemplateFitZ — resolution-matched template chi2 scan, plain arrays
└── templates/
    ├── download_marz_templates.py
    ├── download_desi_pca.py
    ├── marz/   ← EarlyType, Intermediate, LateTypeEmission, Composite, QSO, HighZSFG
    └── pca/    ← rrtemplate-GALAXY-None-v2.6.fits, QSO-LOZ-v1.1, QSO-HIZ-v1.1
```

### Three search methods (engine.py + dialog):
1. **Picket Fence** (was: Line Search) — weighted SNR scan against a linelist, emission or absorption
2. **Template Fit** (was: Template) — resolution-matched chi2 scan against MARZ 1D templates
3. **PCA Fit** (was: PCA) — weighted least-squares against DESI redrock PCA eigenvectors

---

## Root Cause Diagnosis: Why Methods Work Poorly

### Template fit — resolution mismatch (primary bug)
The template is interpolated pixel-by-pixel onto the data grid and compared directly.
This is wrong when the template and data have different spectral resolutions:

- MARZ templates are at ~2 Å/px native resolution, no noise.
- SDSS spectra: R~2000 (FWHM ~3 Å at 6000 Å), similar pixel scale but with real noise.
- Pixel-by-pixel chi2 is dominated by noise residuals, not spectral features.
- Fix: **degrade both to the same (coarser) resolution** before comparing.
  Smooth both to a common FWHM (e.g., 10 Å for SDSS) using a Gaussian kernel.
  Fewer effective pixels, but each pixel now carries real spectral information.

### Picket fence — no weights, fragile sign check  [linelists.py FIXED]
- All lines contribute equally regardless of strength — doublet ratios and f-values ignored.
- Sign check (`if S > 0`) discards valid detections when noise flips the sign at weak lines.
- Fix: per-line weights + type column now in linelist DataFrame. PicketFenceZ accumulates
  weighted SNR normalised by sum(weights), applies sign check per-line from type column.

### Both methods — rb_spectrum coupling
- All engine functions require an `rb_spectrum` object to enter via `_preprocess()`.
- Cumbersome for standalone testing or REPL experimentation.
- Fix: extract core math into `picket_fence.py` and `template_fit.py` that work on plain
  numpy arrays `(wave, flux, ivar)`. Runnable standalone with synthetic data.

---

## Refactor Plan: Two Standalone Core Modules

### picket_fence.py — `PicketFenceZ`

Pure numpy/scipy. No GUI, no rb_spectrum.

```python
class PicketFenceZ:
    """
    Weighted picket-fence redshift scan. Pure numpy — no GUI, no rb_spectrum.

    Preprocessing (at construction):
      - Resolution specified as any of: R=, fwhm_ang=, fwhm_pix=, fwhm_kms=
        → converted internally to fwhm_pix using pixel scale from wave array.
      - Gaussian smooth of flux at smooth_fwhm_pix (any size); None = skip.
        Rebinning (rb_specbin) is the caller's responsibility before passing arrays.
      - Error handling via use_error='auto'|True|False:
          True  → use provided ivar directly
          False → ivar from MAD-STD of smoothed flux (uniform)
          auto  → use ivar if it looks valid (non-uniform, no zeros/negatives),
                  else fall back to MAD-STD with a warning

    Scan (Mode A — direct, recommended):
      For each z:
        obs_wave_i = rest_wave_i * (1+z)
        window = ± window_fwhm * fwhm_pix pixels around obs_wave_i
        S_i = sum(flux_smoothed * ivar_eff) over window  [signed]
        sign check per line from type column: emission keeps S>0, absorption keeps S<0
        score(z) = sum_i[ w_i * |S_i| / sqrt(iv_sum_i) ] / sum(w_i)
      Returns score_array; minima = candidate redshifts.

    Scan (Mode B — detect then match, optional):
      - find_peaks(flux_smoothed, prominence=N*sigma_noise) via scipy
        (negate flux for absorption lines)
      - sigma_noise = 1.4826 * mad(flux_smoothed)
      - user sets prominence_sigma (default 3.0)
      - match detected peaks against linelist to propose z candidates
      - cross-check: how many other lines align at each candidate z?
      Returns ranked z list + detected peak positions for dialog overlay.
    """
    def __init__(self, wave, flux, ivar, line_df,
                 smooth_fwhm_pix=None,    # Gaussian FWHM in pixels; None=skip
                 window_fwhm=1.5,         # window half-width in units of fwhm_pix
                 use_error='auto',        # True | False | 'auto'
                 prominence_sigma=3.0,    # Mode B only: peak detection threshold
                 **resolution):           # R=, fwhm_ang=, fwhm_pix=, fwhm_kms=
        ...

    def run(self, z_array, mode='direct') -> np.ndarray:
        """Return score_array (mode='direct') or ranked z list (mode='detect')."""
        ...
```

**Linelists — DONE.** Five presets in `linelists.py`, each with columns
['wave', 'name', 'weight', 'type']:

  zfind_em      10 lines  pure emission     — Lyα→[NII], weights from fval×4.31 or empirical
  zfind_stellar  8 lines  pure absorption   — MgII, CaII H&K, G-band, MgI b, NaI D
  zfind_igm     11 lines  pure absorption   — OVI, HI, SiII, CII, SiIV, CIV, MgII; all fval
  zfind_galaxy  14 lines  mixed em+abs      — [OII]→[NII] emission + CaII/NaI absorption
  zfind_qso      8 lines  pure emission     — Lyα, CIV, CIII], MgII blend, Hβ, [OIII]

Weight sources: absorption = fval from atom_full.dat; HI recombination emission = fval
normalised to Hα=3.0 (scale ×4.31); forbidden emission = empirical 1–3 scale.
To add a line: add one dict to the relevant list in linelists.py, no other file needed.
Default weight=1.0 if column absent (backwards compatible).

### template_fit.py — `TemplateFitZ`

Pure numpy/scipy. No GUI, no rb_spectrum.

```python
class TemplateFitZ:
    """
    Resolution-matched chi2 template scan.

    At construction: smooths data and template to a common resolution
    (Gaussian convolution to the coarser of data_fwhm_ang / template_fwhm_ang).
    This is done ONCE before the z loop — not inside the loop.

    At each trial z:
      1. Redshift template to observed frame: t_wave_obs = t_wave * (1+z)
      2. Interpolate smoothed template onto smoothed data wavelength grid
      3. Optimal amplitude:  A = sum(F * T * ivar) / sum(T^2 * ivar)
      4. Reduced chi2:       chi2 = sum((F - A*T)^2 * ivar) / N_good

    Returns chi2_array. Minimum = best z.

    get_model(z) — returns (wave_obs, scaled_flux) for the best-fit template
    at redshift z, for plotting the template overlay on the data.
    """
    def __init__(self, wave, flux, ivar,
                 t_wave, t_flux,
                 data_fwhm_ang=3.0,       # FWHM of data (Å); SDSS ~ 3 Å
                 template_fwhm_ang=5.0,   # native FWHM of template (Å)
                 smooth_pixels=1):        # additional boxcar on data; 1 = off
        ...

    def run(self, z_array) -> np.ndarray:
        """Return chi2_array. Minima = best z."""
        ...

    def get_model(self, z) -> tuple[np.ndarray, np.ndarray]:
        """
        Return (wave_obs, model_flux) — the scaled template at redshift z,
        interpolated onto the data wavelength grid. Used for dialog overlay.
        wave_obs is self.wave; model_flux is A * T_interp.
        """
        ...
```

**Resolution matching** — computed once at construction:
```
sigma_extra = sqrt(max(data_fwhm, template_fwhm)^2 - min(data_fwhm, template_fwhm)^2) / 2.355
```
Convolve the finer one to match the coarser. Both go into the z-scan at the same resolution.

### Standalone test blocks

Both files have a `if __name__ == '__main__':` block:

- `picket_fence.py`: build a synthetic spectrum (boxcar emission lines at known z=0.5),
  run `PicketFenceZ.run()`, assert chi2 minimum is within dz < 0.01 of truth, plot.
- `template_fit.py`: load a real MARZ template, build a noisy redshifted version of it
  as fake data, run `TemplateFitZ.run()`, assert minimum at truth z, plot chi2 + overlay.

### engine.py after refactor

`line_search()`, `template_search()`, `pca_search()` become thin wrappers:
1. `_preprocess(spec)` → wave, flux, ivar, continuum
2. Apply data_norm (subtract/normalize/raw)
3. Instantiate `PicketFenceZ` or `TemplateFitZ` or `PCAFitZ`
4. Call `.run(z_array)`
5. Package into `ZFindResult` / `AbsorberResult`

All chi2 math lives in the standalone modules, not in `engine.py`.

---

## Template Overlay in Dialog

When method = Template Fit and user clicks a z on the chi2 panel:
1. Dialog calls `TemplateFitZ.get_model(z)` → `(wave_obs, model_flux)`
2. Overplot `model_flux` on the spectrum panel as a distinct colour (e.g., orange, dashed)
3. Label it: `"{template_name} z={z:.4f}"`
4. Remove the previous overlay when user clicks a different z or runs a new search

**ZFindResult change needed**: store the `TemplateFitZ` instance (or enough to call
`get_model`) on the result object so the dialog can call it after the search completes.
Add `result.searcher` attribute (set by engine.py wrapper). For line search, `searcher=None`
(no overlay needed — line markers already serve that role).

---

## TODO — Prioritised

### HIGH: Core refactor + template overlay (do first)

- [ ] **Create `picket_fence.py`** with `PicketFenceZ` class.
      Weighted SNR accumulation, `run(z_array)` method, standalone test block.

- [ ] **Create `template_fit.py`** with `TemplateFitZ` class.
      Resolution matching at construction, optimal-amplitude chi2 scan, `get_model(z)` method,
      standalone test block using a synthetic noisy spectrum built from a MARZ template.

- [x] **Add weights + type to `linelists.py`** — DONE. Five presets with ['wave','name',
      'weight','type'] columns. Old zfind_abs split into zfind_stellar + zfind_igm.
      Dialog dropdown needs updating to reflect new preset names.

- [ ] **Refactor `engine.py`** — `line_search` → uses `PicketFenceZ`,
      `template_search` → uses `TemplateFitZ`, `pca_search` → uses `PCAFitZ`.
      Store `searcher` on the returned `ZFindResult`.

- [ ] **Template overlay in dialog** — when template search result is displayed and user
      clicks a z on the chi2 panel, call `result.searcher.get_model(z)` and overplot
      on the spectrum panel. Clear overlay on new search or new z click.

- [ ] **Test on real spectra** — SDSS star-forming galaxy, passive galaxy, QSO.
      Confirm template chi2 minimum lands at correct z after resolution matching.
      Compare to old behaviour without resolution matching.

### MEDIUM

- [ ] **Fix "Best of all" chi2 comparison** — normalise each template's chi2 curve by
      `(chi2 - median) / mad` before cross-template comparison. Pick winner on
      normalised scale, not raw chi2.

- [ ] **QThread for engine** — run search in background thread, progress bar in dialog.
      Implement after core refactor is stable.

- [ ] **Update dialog dropdown** — replace `zfind_abs` with `zfind_stellar`, `zfind_igm`,
      `zfind_galaxy`, `zfind_qso` in the linelist combo box.

- [ ] **Overlay combo auto-sync** — when user changes linelist combo in Picket Fence mode,
      auto-update the overlay combo to match.

- [ ] **Delta-chi2 plot mode** — show `(chi2 - median) / mad` instead of raw chi2.
      Makes local minima more visible when baseline is large.

- [ ] **Save/load results** — save ZFindResult to JSON sidecar; reload without rerunning.

### LOW / FUTURE

- [ ] **Phase 9: rb_multispec absorption integration** — "Find Absorbers" button in
      multispec.py, wire accepted_absorbers → AbsorberManager. Separate session.

- [ ] **Absorption mode validation** — test PicketFenceZ absorption on a real QSO
      sightline with known MgII/CIV absorbers. Tune significance threshold and weights.

- [ ] **SDSS eigenspectra** — add as additional template set.

- [ ] **Redshift uncertainty** — z_err from chi2 curvature unreliable for broad minima.
      Consider parabolic fit around minimum or bootstrap.

- [ ] **Star templates** — DESI star PCA templates (A–M, WD) for contamination ID.

---

## Notes for Future Claude Sessions

- Read `engine.py`, `picket_fence.py`, `template_fit.py`, `dialog.py` — they are ground truth.
- `git log --oneline -20` to see recent changes.
- `picket_fence.py` and `template_fit.py` must NEVER import PyQt5, matplotlib, or GUI libs.
- `engine.py` must NEVER import PyQt5 or matplotlib.
- `_preprocess()` in engine.py returns `(wave, flux, ivar, continuum, warn_list)` — wave in Å vacuum.
- MARZ templates live in `templates/marz/`, DESI PCA in `templates/pca/`.
- Resolution matching: MARZ templates ~2–5 Å/px native; always degrade to coarser resolution.
  `TemplateFitZ` does this at construction — do NOT do it inside the z-scan loop.
- `result.searcher` carries the `TemplateFitZ` instance; dialog uses it to call `get_model(z)`.
- Dialog row 2: overlay combo | λmin | λmax | data_norm combo | model_norm combo.
- `model_norm` is greyed out for Picket Fence — only meaningful for Template/PCA.
- `multi_template_search` and `multi_pca_search` exist for "Best of all" runs.
- zgui integration: one "Find z" button in `menu_toolbars.py`; signal wired inline.

### linelists.py — key design decisions (finalised)
- Five presets: zfind_em, zfind_stellar, zfind_igm, zfind_galaxy, zfind_qso
- Columns: wave, name, weight, type ('emission'|'absorption')
- Absorption weights = fval from atom_full.dat (authoritative, physically motivated)
- HI recombination emission weights = fval normalised to Hα=3.0 (scale factor 4.31)
- Forbidden emission weights = empirical 1–3 scale (no fval exists for forbidden lines)
- zfind_galaxy: no Lyα (resonant scattering complex), no MgII (ambiguous em/abs in ISM)
- zfind_qso: MgII listed as blend centre 2799 Å (doublet unresolved in broad-line QSOs)
- atom_full.dat and euv.lst are ABSORPTION line references only — not for emission lines

### PicketFenceZ — key design decisions (finalised)
- Resolution input: any of R=, fwhm_ang=, fwhm_pix=, fwhm_kms= → converted to fwhm_pix
- Smoothing: Gaussian kernel of user-specified FWHM (any size); rebinning via rb_specbin
  is caller's responsibility before passing arrays in
- Error handling: use_error='auto'|True|False; auto validates ivar before using it,
  falls back to MAD-STD (sigma = 1.4826 * mad(flux_smoothed)) if suspect
- Sign check is per-line from type column — not a global mode flag
- Mode A (direct scan): accumulates weighted SNR at predicted line positions
- Mode B (detect+match): scipy find_peaks with prominence=N*sigma_noise, then match peaks
  to linelist; prominence_sigma is user-facing parameter
