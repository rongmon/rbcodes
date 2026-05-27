# rb_zfind — Semi-Automated Redshift & Absorber Finder
## Status & TODO

**Phases 0–8 complete.** Tool is functional and integrated with rb_zgui.
Branch: `master` (merged).

---

## What Exists Now

```
rbcodes/GUIs/zfind/
├── __init__.py
├── main.py          ← CLI: rb_zfind command + launch_zfind()
├── engine.py        ← line_search, template_search, pca_search + multi_* variants
├── dialog.py        ← ZFindDialog PyQt5 popup
├── io.py            ← ZFindResult, ZSolution, AbsorberResult, AbsorberCandidate
├── adapters.py      ← zfind_to_zgui_z, absorbers_to_multispec
├── linelists.py     ← curated presets: zfind_em, zfind_abs
└── templates/
    ├── download_marz_templates.py   ← run once to fetch MARZ templates
    ├── download_desi_pca.py         ← run once to fetch DESI PCA templates
    ├── marz/   ← EarlyType, Intermediate, LateTypeEmission, Composite, QSO, HighZSFG
    └── pca/    ← rrtemplate-GALAXY-None-v2.6.fits, QSO-LOZ-v1.1, QSO-HIZ-v1.1
```

### Three search methods (all in engine.py, all in dialog):
1. **Line Search** — matched-filter chi2 scan against a linelist (emission or absorption mode)
2. **Template** — chi2 scan against MARZ 1D templates (6 types)
3. **PCA** — weighted least-squares against DESI redrock PCA eigenvectors (galaxy, qso_loz, qso_hiz)

### Engine parameters (all three methods):
- `z_min`, `z_max`, `n_steps` — redshift grid
- `wave_min`, `wave_max` — observed-frame Å mask; pixels outside → ivar=0
- `data_norm` — how to preprocess observed spectrum: `normalize` (flux/cont), `subtract` (flux−cont), `raw`
- `model_norm` — how to preprocess template/PCA: `normalize`, `subtract`, `raw` (N/A for line search)
- `smooth_pixels` — boxcar pre-smoothing of flux_work before scan; 1=off
- `fit_continuum` — run polynomial continuum fitter if no continuum in spec

### Dialog UI (Parameter panel rows):
- Row 0: Linelist | z_min | z_max | n_steps
- Row 1: FWHM | Win (pix) | Smooth (pix) | Method | Template/PCA combo | Run
- Row 2: Overlay lines | λ min | λ max | Data norm | Model norm

---

## Known Limitations / Open Issues

### Performance / Accuracy
- **All three methods give mediocre results** — they work but don't reliably identify the
  correct redshift without user guidance. More testing and refinement needed.
- **QSO redshifts** are particularly problematic across all methods. Broad lines confuse
  both template and line search. High-z QSOs need better template coverage.
- **Line search** is sensitive to `window_pixels` and `smooth_pixels` — defaults often
  need manual tuning per spectrum type. No good universal default found yet.
- **Template search** works decently when user selects the right template type manually.
  "Best of all" often picks wrong template due to chi2 scale differences between templates.
- **PCA** runs correctly (DESI redrock eigenvectors, weighted least-squares) but galaxy
  template at 10 eigenvectors is slow (~2.5s at 5000 steps). Results not yet validated
  on real spectra.

### Line Search Specific
- **No line weights** — all lines contribute equally regardless of expected strength.
  Doublet ratios (5007:4960=3:1, MgII 2796:2803=2:1) and f-value weighting would
  meaningfully improve emission and absorption ranking respectively. Deferred.
- **Absorption mode** untested on real QSO spectra. Significance threshold (>3.0 for
  auto-accept) not validated. `zfind_abs` linelist (OVI, HI, SiIV, CIV, MgII) not
  tuned. Needs testing on a QSO sightline with known absorbers.
- **zfind_agn** curated preset discussed but not implemented. Broad-line AGN needs
  different approach (wider windows, Hα+Hβ+MgII+CIV broad lines).

### Template / PCA Specific
- **Continuum mismatch** is the main failure mode. When `data_norm='normalize'` and
  the continuum fit is poor (emission lines biasing the polynomial), normalized features
  are distorted. `data_norm='raw'` + `fit_continuum=False` sometimes works better.
- **Template "Best of all"** picks lowest absolute chi2/dof across all templates.
  This is unreliable because different templates have different chi2 scales at their
  best-fit z. Need to normalize chi2 curves before cross-template comparison, or use
  a delta-chi2 relative to baseline.
- **MARZ QSO template** covers 900–5748 Å rest-frame — fine for low-z, but for
  z>3 the template becomes too short to cover the optical. HIZ template not tried.
- **PCA galaxy template** has been downsampled 10× (0.1→1 Å/px) to reduce runtime.
  May affect accuracy for narrow-line galaxies.

### UX
- **No progress indicator** while search runs — UI freezes on long runs (5000 steps,
  galaxy PCA). Should run engine in a QThread with progress signal.
- **Wave range fields** (λ min/λ max) are new — not yet tested interactively.
- **Overlay combo** default does not auto-sync when user changes linelist combo —
  the two are independent by design but this is sometimes confusing.

---

## TODO — Prioritised

### High priority (next session)
- [ ] **Interactive testing on diverse spectra** — SDSS galaxies (star-forming, passive,
      AGN), z>1 galaxy, QSO at z~1 and z~3. For each, note which method+settings works.
      Build up a best-practices guide.
- [ ] **Fix "Best of all" chi2 comparison** — normalise each template's chi2 curve by
      its own baseline/median before picking the winner. Currently picks the template
      with lowest absolute chi2 regardless of fit quality.
- [ ] **Line weights for line_search** — add optional `weight` column to linelist DataFrames.
      zfind_em: use doublet ratios (Hα=3, [OIII]5007=3, [OIII]4959=1, [OII]=2, etc.).
      zfind_abs: use f-values (CIV stronger than OVI, MgII doublet ratio 2:1).
      Engine: `delta_chi2 *= weight` in accumulation loop.
- [ ] **Absorption mode validation** — test on a real QSO sightline with known MgII/CIV
      absorbers. Adjust significance threshold and zfind_abs linelist as needed.

### Medium priority
- [ ] **QThread for engine** — run search in background thread so UI stays responsive.
      Emit progress signal after each z step or each 10% chunk. Show progress bar.
- [ ] **zfind_agn curated preset** — Hα, Hβ, MgII 2798, CIV 1549, CIII] 1909, Lyα.
      Wider `window_pixels` default for broad-line matching.
- [ ] **Overlay combo auto-sync** — when user changes linelist combo in Line Search mode,
      auto-update overlay combo to match (user can still override manually).
- [ ] **Delta-chi2 plot mode** — option to show (chi2 - baseline) / std instead of raw chi2.
      Makes local minima more visible when baseline chi2 is large.
- [ ] **Save/load results** — save ZFindResult to a JSON or FITS sidecar file so user
      can reload a previous search without rerunning.

### Low priority / future
- [ ] **Phase 9: rb_multispec absorption integration** — add "Find Absorbers" button to
      multispec.py, wire accepted_absorbers → AbsorberManager. Separate session.
- [ ] **SDSS eigenspectra** — add SDSS galaxy/QSO eigenspectra as additional template set
      (broader rest-frame coverage than MARZ for some types).
- [ ] **Redshift uncertainty** — current z_err from chi2 curvature is unreliable for
      broad/flat chi2 minima. Consider bootstrap or MCMC-like error estimation.
- [ ] **Star templates** — DESI star PCA templates available (A, B, F, G, K, M, WD).
      Useful for contamination identification. Low priority for current science case.

---

## Notes for Future Claude Sessions

- Read engine.py and dialog.py directly — they are the ground truth. This TODO may lag.
- `git log --oneline -20` to see recent changes.
- Engine (`engine.py`) must NEVER import PyQt5, matplotlib, or GUI libraries.
- `_preprocess()` returns `(wave, flux, ivar, continuum, warn_list)` — wave in Å vacuum.
- All engine functions: `wave_min`/`wave_max` zero-out ivar before the z scan.
- `data_norm` applies to observed spectrum; `model_norm` applies to template/PCA eigenvectors.
- MARZ templates live in `templates/marz/`, DESI PCA in `templates/pca/`.
- Dialog row 2 has: overlay combo | λmin | λmax | data_norm combo | model_norm combo.
- `model_norm` is greyed out (disabled) for Line Search — only meaningful for Template/PCA.
- PCA galaxy template is downsampled 10× inside `_load_pca()` to keep runtime reasonable.
- `multi_template_search` and `multi_pca_search` exist for "Best of all" runs.
- zgui integration: one "Find z" button in `menu_toolbars.py`; signal wired inline.
