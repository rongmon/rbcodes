"""
picket_fence.py — Weighted picket-fence redshift scanner.

Standalone module: pure numpy/scipy, no GUI, no rb_spectrum.
Accepts plain arrays (wave, flux, ivar) from any instrument.

Algorithm (Mode A — direct scan)
---------------------------------
For each trial redshift z:
  1. Shift line rest wavelengths into observed frame: λ_obs = λ_rest * (1+z)
  2. For each in-range line, extract a pixel window centred on λ_obs
  3. Matched-filter signal: S = sum(flux_smooth * ivar_eff) over the window
  4. Per-line sign check from the 'type' column:
       emission   → keep only S > 0  (feature above continuum)
       absorption → keep only S < 0  (feature below continuum)
  5. Weighted accumulation:
       score(z) = -sum_i[ w_i * |S_i| / sqrt(iv_sum_i) ] / sum(w_used)
     Negative so that minima = best redshift (consistent with chi2 convention).

Algorithm (Mode B — detect then match)
----------------------------------------
  1. Find emission/absorption peaks in smoothed flux via scipy.signal.find_peaks
     with prominence = prominence_sigma * sigma_noise
     (sigma_noise = 1.4826 * MAD of smoothed flux)
  2. For each detected peak, generate z candidates by pairing peak wavelength
     with every same-type line in the linelist
  3. Score each candidate z by summing weights of all linelist lines that have
     a detected peak within 2-pixel tolerance at that z
  4. Return ranked list of candidates

Usage
-----
    from rbcodes.GUIs.zfind.picket_fence import PicketFenceZ
    import numpy as np

    pf = PicketFenceZ(wave, flux, ivar, line_df, fwhm_ang=3.0, smooth_fwhm_pix=3.0)
    z_array = np.linspace(0.0, 1.5, 5000)
    score = pf.run(z_array)          # minima = best z

    peaks = pf.detect_peaks()        # {'emission': array, 'absorption': array}
    candidates = pf.match_peaks()    # Mode B ranked list

Run standalone for a self-test on synthetic data:
    python picket_fence.py

Continuum handling
------------------
PicketFenceZ expects flux that is already continuum-subtracted or
continuum-normalised before being passed in.  It does NOT fit or remove
the continuum internally.

When called through engine.py (the normal path):
  - engine._preprocess() fits a polynomial continuum via fit_optimal_polynomial
    (BIC selects order automatically, iterative sigma-clipping rejects lines)
  - engine.line_search() applies data_norm='subtract' (default) before
    instantiating PicketFenceZ, so the class always receives flux_work = flux - cont.

When used standalone (CLI / notebook):
  - Pass continuum-subtracted flux: flux_sub = flux - continuum
  - Or normalised flux: flux_norm = flux / continuum - 1.0
  - A polynomial fit via fit_optimal_polynomial works well for smooth continua.
  - Raw flux (no continuum removal) will degrade performance: a non-zero
    continuum baseline makes S = sum(flux*ivar) positive everywhere, so the
    emission sign check passes spuriously and every z accumulates signal.
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import median_abs_deviation


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _to_fwhm_pix(wave, **resolution):
    """
    Convert any resolution specification to FWHM in pixels.

    Accepted keyword arguments (exactly one should be provided):
        R          — resolving power (dimensionless): FWHM = λ_centre / R
        fwhm_ang   — FWHM in Angstroms
        fwhm_pix   — FWHM directly in pixels (passthrough)
        fwhm_kms   — FWHM in km/s

    Returns float or None if no resolution keyword was given.
    """
    if not resolution:
        return None

    dwave   = float(np.median(np.diff(wave)))   # Å per pixel
    lam_c   = float(np.median(wave))            # central wavelength
    C_KMS   = 299792.458                        # speed of light km/s

    if 'fwhm_pix' in resolution:
        return float(resolution['fwhm_pix'])
    if 'fwhm_ang' in resolution:
        return float(resolution['fwhm_ang']) / dwave
    if 'R' in resolution:
        return (lam_c / float(resolution['R'])) / dwave
    if 'fwhm_kms' in resolution:
        return (float(resolution['fwhm_kms']) / C_KMS * lam_c) / dwave

    raise ValueError(
        f"Unrecognised resolution keyword(s): {list(resolution.keys())}. "
        "Use one of: R, fwhm_ang, fwhm_pix, fwhm_kms."
    )


def _validate_ivar(ivar):
    """
    Return True if ivar looks like a genuine error spectrum.

    Fails (returns False) for:
      - All-zero or all-negative arrays
      - All-constant arrays (placeholder like all-ones): std/median < 1 %
      - Fewer than 10 finite positive values
    """
    good = ivar[np.isfinite(ivar) & (ivar > 0)]
    if len(good) < 10:
        return False
    if np.std(good) / (np.median(good) + 1e-30) < 0.01:
        return False
    return True


def _mad_ivar(flux_smooth):
    """Uniform ivar from MAD-STD of smoothed flux (single sigma for all pixels)."""
    sigma = 1.4826 * float(
        median_abs_deviation(flux_smooth, nan_policy='omit')
    )
    if sigma == 0.0:
        finite = flux_smooth[np.isfinite(flux_smooth)]
        sigma = float(np.std(finite)) if len(finite) > 0 else 1.0
    sigma = max(sigma, 1e-30)
    return np.full_like(flux_smooth, 1.0 / sigma ** 2)


# ---------------------------------------------------------------------------
# PicketFenceZ
# ---------------------------------------------------------------------------

class PicketFenceZ:
    """
    Weighted picket-fence redshift scanner.

    Parameters
    ----------
    wave : array-like
        Observed wavelengths in Å, monotonically increasing.
    flux : array-like
        Flux array (any units, continuum-subtracted or normalised preferred).
    ivar : array-like
        Inverse-variance array. Use zeros for bad/masked pixels.
    line_df : pandas.DataFrame
        Line list with columns:
          wave   (float) — rest wavelength in Å          [required]
          name   (str)   — label                          [required]
          weight (float) — relative line strength         [optional, default 1.0]
          type   (str)   — 'emission' or 'absorption'     [optional, default 'emission']
    smooth_fwhm_pix : float or None
        FWHM of Gaussian pre-smoothing kernel in pixels.
        None or 0 → no smoothing.  Any positive value is accepted.
        For noisy spectra, set to ~1–2× the resolution element.
        Rebinning (rb_specbin) should be done by the caller before passing arrays.
    window_fwhm : float
        Half-width of the pixel window around each line centre,
        in units of fwhm_pix.  Ignored if no resolution is specified
        (falls back to window_pixels).
    window_pixels : int
        Fallback half-window in pixels when no resolution is given (default 5).
    use_error : True | False | 'auto'
        True  → use provided ivar as-is
        False → always use MAD-STD uniform ivar (ignores input ivar)
        'auto'→ use ivar if it passes validation, else fall back to MAD-STD
    prominence_sigma : float
        Peak detection threshold for Mode B, in units of sigma_noise.
        sigma_noise = 1.4826 * MAD(flux_smooth).  Default 3.0.
    **resolution : keyword
        Exactly one of: R=, fwhm_ang=, fwhm_pix=, fwhm_kms=
        Determines the window half-width when window_fwhm is used.
        Optional — if omitted, window_pixels is used instead.

    Attributes (available after construction)
    ------------------------------------------
    flux_smooth : ndarray   — Gaussian-smoothed flux (or raw if no smoothing)
    ivar        : ndarray   — ivar actually used (validated or MAD-STD)
    fwhm_pix    : float|None — resolution in pixels (None if not specified)
    warnings    : list[str]  — any issues detected during construction
    """

    def __init__(self, wave, flux, ivar, line_df,
                 smooth_fwhm_pix=None,
                 window_fwhm=1.5,
                 window_pixels=5,
                 use_error='auto',
                 prominence_sigma=3.0,
                 **resolution):

        self.wave = np.asarray(wave, dtype=float)
        flux_arr  = np.asarray(flux, dtype=float)
        ivar_arr  = np.asarray(ivar, dtype=float)

        self.warnings = []

        # --- line_df: ensure weight and type columns exist ---
        self._lines = line_df.copy()
        if 'weight' not in self._lines.columns:
            self._lines['weight'] = 1.0
        if 'type' not in self._lines.columns:
            self._lines['type'] = 'emission'
        # cache as arrays for speed
        self._rest  = self._lines['wave'].values.astype(float)
        self._names = self._lines['name'].values
        self._wts   = self._lines['weight'].values.astype(float)
        self._types = self._lines['type'].values

        # --- resolution → fwhm_pix ---
        self.fwhm_pix = _to_fwhm_pix(self.wave, **resolution)

        # --- window half-width in pixels ---
        if self.fwhm_pix is not None and self.fwhm_pix > 0:
            self._wp = max(2, int(round(window_fwhm * self.fwhm_pix)))
        else:
            self._wp = max(2, int(window_pixels))

        # --- pixel scale (precomputed) ---
        self._dwave = float(np.median(np.diff(self.wave)))

        # --- Gaussian smoothing ---
        flux_clean = np.nan_to_num(flux_arr, nan=0.0)
        if smooth_fwhm_pix is not None and smooth_fwhm_pix > 0:
            sigma_pix = float(smooth_fwhm_pix) / 2.355
            self.flux_smooth = gaussian_filter1d(flux_clean, sigma=sigma_pix)
        else:
            self.flux_smooth = flux_clean.copy()

        # --- ivar selection ---
        # zero out ivar for non-finite flux pixels before validation
        ivar_clean = ivar_arr.copy()
        ivar_clean[~np.isfinite(flux_arr)] = 0.0

        if use_error is True:
            self.ivar = ivar_clean
        elif use_error is False:
            self.ivar = _mad_ivar(self.flux_smooth)
            self.warnings.append(
                "use_error=False: using MAD-STD uniform ivar."
            )
        else:   # 'auto'
            if _validate_ivar(ivar_clean):
                self.ivar = ivar_clean
            else:
                self.ivar = _mad_ivar(self.flux_smooth)
                self.warnings.append(
                    "ivar failed validation (may be placeholder) — "
                    "using MAD-STD uniform ivar."
                )

        # --- sigma_noise for Mode B ---
        self._sigma_noise = 1.4826 * float(
            median_abs_deviation(self.flux_smooth, nan_policy='omit')
        )
        if self._sigma_noise == 0.0:
            finite = self.flux_smooth[np.isfinite(self.flux_smooth)]
            self._sigma_noise = float(np.std(finite)) if len(finite) > 0 else 1.0
        self._prominence = prominence_sigma * self._sigma_noise

    # -----------------------------------------------------------------------
    # Mode A: direct scan
    # -----------------------------------------------------------------------

    def run(self, z_array):
        """
        Direct picket-fence scan over a redshift grid (Mode A).

        For each z in z_array, accumulates weighted matched-filter SNR at
        predicted line positions.  Score is negative; minima = best redshift.

        Parameters
        ----------
        z_array : array-like
            Trial redshifts, e.g. np.linspace(0.0, 1.5, 5000).

        Returns
        -------
        score_array : ndarray, shape (len(z_array),)
            NaN where no lines fall in range.  Minima indicate best z.
        """
        z_array   = np.asarray(z_array, dtype=float)
        wave      = self.wave
        flux_s    = self.flux_smooth
        ivar      = self.ivar
        rest      = self._rest
        weights   = self._wts
        types     = self._types
        wp        = self._wp
        dwave     = self._dwave

        wmin = wave[0]  + wp * dwave
        wmax = wave[-1] - wp * dwave

        score_array = np.full(len(z_array), np.nan)

        for i, z in enumerate(z_array):
            obs_waves = rest * (1.0 + z)
            in_range  = (obs_waves > wmin) & (obs_waves < wmax)

            if not np.any(in_range):
                continue

            accum   = 0.0
            w_total = 0.0   # all in-range lines (normalisation denominator)
            w_used  = 0.0   # lines that passed sign check (must be > 0 to record)

            for obs_w, wt, ltype in zip(obs_waves[in_range],
                                         weights[in_range],
                                         types[in_range]):
                w_total += wt   # count every in-range line regardless of sign

                # nearest pixel via searchsorted
                pix = int(np.searchsorted(wave, obs_w))
                pix = min(max(pix, wp), len(wave) - wp - 1)

                fw = flux_s[pix - wp : pix + wp + 1]
                iv = ivar  [pix - wp : pix + wp + 1]

                iv_pos = iv[iv > 0]
                if len(iv_pos) < 2:
                    continue
                iv_sum = float(np.sum(iv_pos))

                S = float(np.sum(fw * iv))

                # per-line sign check
                if ltype == 'emission'   and S <= 0.0:
                    continue
                if ltype == 'absorption' and S >= 0.0:
                    continue

                accum  -= wt * abs(S) / np.sqrt(iv_sum)
                w_used += wt

            # Normalise by sqrt(w_total_in_range).
            # Numerator: strong lines dominate via weights (weighted total SNR).
            # Denominator: sqrt(w_total) gives a soft penalty for missing lines —
            #   at the right z all lines add constructively, score ∝ sqrt(w_total);
            #   an accidental single-line hit scores w_i / sqrt(w_total) < 1.
            # Dividing by w_total (weighted average) over-penalises; no division
            # (raw sum) ignores missing lines. sqrt(w_total) balances both.
            if w_used > 0.0:
                score_array[i] = accum / np.sqrt(w_total)

        return score_array

    # -----------------------------------------------------------------------
    # Peak detection (shared by Mode B and dialog visualization)
    # -----------------------------------------------------------------------

    def detect_peaks(self):
        """
        Detect emission and absorption peaks in the smoothed flux.

        Uses scipy.signal.find_peaks with prominence = prominence_sigma * sigma_noise.

        Returns
        -------
        dict with keys:
            'emission'   : ndarray of peak wavelengths (Å)
            'absorption' : ndarray of trough wavelengths (Å)
        """
        from scipy.signal import find_peaks

        em_idx,  _ = find_peaks( self.flux_smooth, prominence=self._prominence)
        abs_idx, _ = find_peaks(-self.flux_smooth, prominence=self._prominence)

        return {
            'emission':   self.wave[em_idx]  if len(em_idx)  > 0 else np.array([]),
            'absorption': self.wave[abs_idx] if len(abs_idx) > 0 else np.array([]),
        }

    # -----------------------------------------------------------------------
    # Mode B: detect then match
    # -----------------------------------------------------------------------

    def match_peaks(self, z_tol=0.005):
        """
        Detect peaks and match them to the linelist to generate z candidates (Mode B).

        For each detected peak, pairs it with every same-type line in the linelist
        to generate z candidates.  Each candidate is scored by the total weight of
        all linelist lines that have a detected peak within 2 pixels at that z.
        Nearby candidates within z_tol are deduplicated (highest score kept).

        Parameters
        ----------
        z_tol : float
            Minimum separation between returned candidates (default 0.005).

        Returns
        -------
        list of dicts, sorted by score descending:
            {'z', 'n_matches', 'score', 'lines': [{'name', 'wave_obs'}, ...]}
        """
        peaks = self.detect_peaks()
        em_peaks  = peaks['emission']
        abs_peaks = peaks['absorption']

        wave  = self.wave
        tol   = 2.0 * self._dwave   # 2-pixel matching tolerance in Å

        # generate z candidates from peak-line pairings
        z_pool = []
        for pk_arr, ltype in [(em_peaks, 'emission'), (abs_peaks, 'absorption')]:
            mask = self._types == ltype
            if not np.any(mask) or len(pk_arr) == 0:
                continue
            for obs_w in pk_arr:
                for rw in self._rest[mask]:
                    zc = float(obs_w) / float(rw) - 1.0
                    if zc >= 0.0:
                        z_pool.append(zc)

        if not z_pool:
            return []

        # score each unique z candidate
        results = []
        seen = set()
        for zc in z_pool:
            key = round(zc, 3)
            if key in seen:
                continue
            seen.add(key)

            matched = []
            score   = 0.0
            for rw, nm, wt, lt in zip(self._rest, self._names,
                                       self._wts,  self._types):
                obs_w = rw * (1.0 + zc)
                if obs_w < wave[0] or obs_w > wave[-1]:
                    continue
                pk_arr = em_peaks if lt == 'emission' else abs_peaks
                if len(pk_arr) == 0:
                    continue
                if float(np.min(np.abs(pk_arr - obs_w))) <= tol:
                    matched.append({'name': nm, 'wave_obs': obs_w})
                    score += wt

            if matched:
                results.append({
                    'z':         zc,
                    'n_matches': len(matched),
                    'score':     score,
                    'lines':     matched,
                })

        # sort by score descending
        results.sort(key=lambda r: r['score'], reverse=True)

        # deduplicate: keep highest-score candidate within z_tol
        deduped = []
        for r in results:
            if not any(abs(r['z'] - d['z']) < z_tol for d in deduped):
                deduped.append(r)

        return deduped


# ---------------------------------------------------------------------------
# Standalone self-test
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    import pandas as pd

    print("picket_fence.py — self-test on synthetic spectra")
    print("=" * 55)

    RNG      = np.random.default_rng(42)
    FWHM_ANG = 4.0    # instrument resolution FWHM in Å

    # ----------------------------------------------------------------
    # Helper: inject Gaussian lines into a flux array
    # ----------------------------------------------------------------
    def _inject(wave, flux, line_df, z, amp_scale, sign=+1):
        """Add Gaussian lines (sign=+1 emission, sign=-1 absorption)."""
        sigma = FWHM_ANG / 2.355
        for _, row in line_df.iterrows():
            obs_w = row['wave'] * (1.0 + z)
            if wave[0] < obs_w < wave[-1]:
                amp = sign * row['weight'] * amp_scale
                flux += amp * np.exp(-0.5 * ((wave - obs_w) / sigma) ** 2)
        return flux

    def _top3(z_array, score):
        working = score.copy()
        print("  Top 3 minima:")
        for rank in range(3):
            if not np.any(np.isfinite(working)):
                break
            idx = int(np.nanargmin(working))
            print(f"    #{rank+1}  z = {z_array[idx]:.5f}  score = {score[idx]:.4f}")
            working[np.abs(z_array - z_array[idx]) < 0.03] = np.nan

    # ================================================================
    # TEST 1 — Emission lines  (star-forming galaxy, zfind_em subset)
    # ================================================================
    print("\n── TEST 1: Emission lines ──────────────────────────────")

    EM_LINES = pd.DataFrame([
        {'wave': 3727.09, 'name': '[OII] 3727',  'weight': 2.0, 'type': 'emission'},
        {'wave': 3729.87, 'name': '[OII] 3729',  'weight': 2.0, 'type': 'emission'},
        {'wave': 4862.68, 'name': 'Hbeta',        'weight': 0.5, 'type': 'emission'},
        {'wave': 4960.30, 'name': '[OIII] 4960',  'weight': 1.0, 'type': 'emission'},
        {'wave': 5008.24, 'name': '[OIII] 5007',  'weight': 3.0, 'type': 'emission'},
        {'wave': 6564.61, 'name': 'Halpha',        'weight': 3.0, 'type': 'emission'},
        {'wave': 6585.28, 'name': '[NII] 6583',   'weight': 1.0, 'type': 'emission'},
    ])

    Z_EM      = 0.35
    wave_em   = np.linspace(3800.0, 10000.0, 6200)   # ~1 Å/pixel
    flux_em   = np.ones_like(wave_em)
    flux_em   = _inject(wave_em, flux_em, EM_LINES, Z_EM, amp_scale=2.0, sign=+1)
    noise_em  = 0.15
    flux_em  += RNG.normal(0.0, noise_em, size=len(wave_em))
    ivar_em   = np.full_like(flux_em, 1.0 / noise_em ** 2)
    flux_em  -= 1.0    # subtract flat continuum

    print(f"  True z = {Z_EM}  |  {len(wave_em)} px  "
          f"  {wave_em[0]:.0f}–{wave_em[-1]:.0f} Å  |  noise σ = {noise_em}")

    pf_em = PicketFenceZ(
        wave_em, flux_em, ivar_em, EM_LINES,
        smooth_fwhm_pix=2.0, window_fwhm=1.5,
        fwhm_ang=FWHM_ANG, use_error=True,
    )
    print(f"  fwhm_pix = {pf_em.fwhm_pix:.2f} px,  window = ±{pf_em._wp} px")

    z_em    = np.linspace(0.0, 1.0, 10000)
    sc_em   = pf_em.run(z_em)
    z_best_em = float(z_em[np.nanargmin(sc_em)])
    dz_em     = abs(z_best_em - Z_EM)
    print(f"  Mode A best z = {z_best_em:.5f}  Δz = {dz_em:.5f}")
    _top3(z_em, sc_em)
    assert dz_em < 0.01, f"FAIL: Δz = {dz_em:.5f}"
    print(f"  Mode A PASS")

    # Mode B
    cands_em = pf_em.match_peaks()
    if cands_em:
        dz_b = abs(cands_em[0]['z'] - Z_EM)
        lines_str = ', '.join(l['name'] for l in cands_em[0]['lines'])
        print(f"  Mode B best z = {cands_em[0]['z']:.4f}  Δz = {dz_b:.5f}  "
              f"[{lines_str}]")
        status = "PASS" if dz_b < 0.05 else "WARNING (>0.05)"
        print(f"  Mode B {status}")
    else:
        print("  Mode B: no candidates found")

    # ================================================================
    # TEST 2 — Absorption lines  (IGM absorber on QSO sightline)
    # ================================================================
    print("\n── TEST 2: Absorption lines ────────────────────────────")

    ABS_LINES = pd.DataFrame([
        {'wave': 1548.20, 'name': 'CIV 1548',  'weight': 0.191, 'type': 'absorption'},
        {'wave': 1550.78, 'name': 'CIV 1550',  'weight': 0.095, 'type': 'absorption'},
        {'wave': 2796.35, 'name': 'MgII 2796', 'weight': 0.612, 'type': 'absorption'},
        {'wave': 2803.53, 'name': 'MgII 2803', 'weight': 0.305, 'type': 'absorption'},
        {'wave': 1393.76, 'name': 'SiIV 1393', 'weight': 0.514, 'type': 'absorption'},
        {'wave': 1402.77, 'name': 'SiIV 1402', 'weight': 0.255, 'type': 'absorption'},
    ])

    Z_ABS     = 0.72
    wave_abs  = np.linspace(3800.0, 10000.0, 6200)
    # flat QSO continuum = 5.0
    flux_abs  = np.full_like(wave_abs, 5.0)
    flux_abs  = _inject(wave_abs, flux_abs, ABS_LINES, Z_ABS, amp_scale=3.0, sign=-1)
    noise_abs = 0.2
    flux_abs += RNG.normal(0.0, noise_abs, size=len(wave_abs))
    ivar_abs  = np.full_like(flux_abs, 1.0 / noise_abs ** 2)
    flux_abs -= 5.0    # subtract flat QSO continuum

    print(f"  True z = {Z_ABS}  |  {len(wave_abs)} px  "
          f"  {wave_abs[0]:.0f}–{wave_abs[-1]:.0f} Å  |  noise σ = {noise_abs}")

    pf_abs = PicketFenceZ(
        wave_abs, flux_abs, ivar_abs, ABS_LINES,
        smooth_fwhm_pix=2.0, window_fwhm=1.5,
        fwhm_ang=FWHM_ANG, use_error=True,
    )
    print(f"  fwhm_pix = {pf_abs.fwhm_pix:.2f} px,  window = ±{pf_abs._wp} px")

    z_abs    = np.linspace(0.0, 2.0, 10000)
    sc_abs   = pf_abs.run(z_abs)
    z_best_abs = float(z_abs[np.nanargmin(sc_abs)])
    dz_abs     = abs(z_best_abs - Z_ABS)
    print(f"  Mode A best z = {z_best_abs:.5f}  Δz = {dz_abs:.5f}")
    _top3(z_abs, sc_abs)
    assert dz_abs < 0.01, f"FAIL: Δz = {dz_abs:.5f}"
    print(f"  Mode A PASS")

    # Mode B
    cands_abs = pf_abs.match_peaks()
    if cands_abs:
        dz_b = abs(cands_abs[0]['z'] - Z_ABS)
        lines_str = ', '.join(l['name'] for l in cands_abs[0]['lines'])
        print(f"  Mode B best z = {cands_abs[0]['z']:.4f}  Δz = {dz_b:.5f}  "
              f"[{lines_str}]")
        status = "PASS" if dz_b < 0.05 else "WARNING (>0.05)"
        print(f"  Mode B {status}")
    else:
        print("  Mode B: no candidates found")

    # ================================================================
    # TEST 3 — Mixed linelist  (galaxy: emission + stellar absorption)
    # ================================================================
    print("\n── TEST 3: Mixed linelist (galaxy) ─────────────────────")

    MIX_LINES = pd.DataFrame([
        {'wave': 5008.24, 'name': '[OIII] 5007', 'weight': 3.0, 'type': 'emission'},
        {'wave': 6564.61, 'name': 'Halpha',       'weight': 3.0, 'type': 'emission'},
        {'wave': 3727.09, 'name': '[OII] 3727',   'weight': 2.0, 'type': 'emission'},
        {'wave': 3934.78, 'name': 'CaII K',       'weight': 0.635,'type': 'absorption'},
        {'wave': 3969.59, 'name': 'CaII H',       'weight': 0.315,'type': 'absorption'},
        {'wave': 5891.58, 'name': 'NaI D1',       'weight': 0.631,'type': 'absorption'},
    ])

    Z_MIX    = 0.18
    wave_mix = np.linspace(3800.0, 10000.0, 6200)
    flux_mix = np.ones_like(wave_mix)
    # inject emission lines
    em_mask  = MIX_LINES['type'] == 'emission'
    flux_mix = _inject(wave_mix, flux_mix, MIX_LINES[em_mask],  Z_MIX, amp_scale=2.0, sign=+1)
    # inject absorption lines
    ab_mask  = MIX_LINES['type'] == 'absorption'
    flux_mix = _inject(wave_mix, flux_mix, MIX_LINES[ab_mask],  Z_MIX, amp_scale=2.0, sign=-1)
    noise_mix = 0.1
    flux_mix += RNG.normal(0.0, noise_mix, size=len(wave_mix))
    ivar_mix  = np.full_like(flux_mix, 1.0 / noise_mix ** 2)
    flux_mix -= 1.0

    print(f"  True z = {Z_MIX}  |  {len(wave_mix)} px  "
          f"  {wave_mix[0]:.0f}–{wave_mix[-1]:.0f} Å  |  noise σ = {noise_mix}")

    pf_mix = PicketFenceZ(
        wave_mix, flux_mix, ivar_mix, MIX_LINES,
        smooth_fwhm_pix=2.0, window_fwhm=1.5,
        fwhm_ang=FWHM_ANG, use_error=True,
    )
    print(f"  fwhm_pix = {pf_mix.fwhm_pix:.2f} px,  window = ±{pf_mix._wp} px")

    z_mix    = np.linspace(0.0, 0.6, 10000)
    sc_mix   = pf_mix.run(z_mix)
    z_best_mix = float(z_mix[np.nanargmin(sc_mix)])
    dz_mix     = abs(z_best_mix - Z_MIX)
    print(f"  Mode A best z = {z_best_mix:.5f}  Δz = {dz_mix:.5f}")
    _top3(z_mix, sc_mix)
    assert dz_mix < 0.01, f"FAIL: Δz = {dz_mix:.5f}"
    print(f"  Mode A PASS")

    # ================================================================
    # TEST 4 — Raw flux with polynomial continuum fit then subtract
    # ================================================================
    # Demonstrates the correct standalone workflow when continuum is present.
    # fit_optimal_polynomial (BIC order selection + iterative line rejection)
    # is the same fitter used by engine._preprocess().
    # ================================================================
    print("\n── TEST 4: Raw flux + continuum fit ────────────────────")

    Z_CONT   = 0.52
    wave_c   = np.linspace(3800.0, 9000.0, 5200)   # ~1 Å/pixel

    # quadratic continuum — rises and curves significantly across the range
    cont_true = (  3.0
                 + 1.5e-3 * (wave_c - 6000.0)
                 - 8.0e-8 * (wave_c - 6000.0) ** 2 )

    flux_c = cont_true.copy()
    flux_c = _inject(wave_c, flux_c, EM_LINES, Z_CONT, amp_scale=1.5, sign=+1)
    noise_c = 0.12
    flux_c += RNG.normal(0.0, noise_c, size=len(wave_c))
    ivar_c  = np.full_like(flux_c, 1.0 / noise_c ** 2)

    print(f"  True z = {Z_CONT}  |  continuum: quadratic (3–8 flux units across range)")
    print(f"  noise σ = {noise_c}  |  peak line amplitude ≈ {1.5 * 3.0:.1f} × noise")

    # --- fit continuum with fit_optimal_polynomial (same as engine._preprocess) ---
    try:
        from rbcodes.IGM.rb_iter_contfit import fit_optimal_polynomial
        result   = fit_optimal_polynomial(
            wave_c, flux_c, error=None,
            use_weights=False, silent=True, plot=False
        )
        cont_fit  = result['continuum']
        flux_csub = flux_c - cont_fit
        cont_rms  = float(np.std(cont_true - cont_fit))
        print(f"  Continuum fit RMS residual = {cont_rms:.4f} flux units")
    except Exception as e:
        print(f"  fit_optimal_polynomial failed ({e}) — falling back to median")
        cont_fit  = np.full_like(flux_c, np.median(flux_c))
        flux_csub = flux_c - cont_fit

    pf_c = PicketFenceZ(
        wave_c, flux_csub, ivar_c, EM_LINES,
        smooth_fwhm_pix=2.0, window_fwhm=1.5,
        fwhm_ang=FWHM_ANG, use_error=True,
    )

    z_c      = np.linspace(0.2, 0.8, 10000)
    sc_c     = pf_c.run(z_c)
    z_best_c = float(z_c[np.nanargmin(sc_c)])
    dz_c     = abs(z_best_c - Z_CONT)
    print(f"  Mode A best z = {z_best_c:.5f}  Δz = {dz_c:.5f}")
    _top3(z_c, sc_c)
    assert dz_c < 0.01, f"FAIL: Δz = {dz_c:.5f}"
    print(f"  Mode A PASS")

    # ================================================================
    # Optional matplotlib plot — 4 panels
    # ================================================================
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(4, 2, figsize=(14, 13))
        fig.suptitle('PicketFenceZ self-test', fontsize=13)

        def _plot_pair(ax_spec, ax_score, wave, flux_raw, pf_obj,
                       z_arr, sc, z_true, z_found, lines_df, title,
                       flux_raw2=None, label2=None):
            # spectrum
            ax_spec.plot(wave, flux_raw, lw=0.4, color='C0', alpha=0.7,
                         label='flux − cont')
            if flux_raw2 is not None:
                ax_spec.plot(wave, flux_raw2, lw=0.8, color='C4', alpha=0.8,
                             label=label2)
            ax_spec.plot(wave, pf_obj.flux_smooth, lw=1.0, color='C1',
                         label='smoothed')
            for _, row in lines_df.iterrows():
                obs_w = row['wave'] * (1.0 + z_true)
                if wave[0] < obs_w < wave[-1]:
                    col = 'C2' if row['type'] == 'emission' else 'C3'
                    ax_spec.axvline(obs_w, color=col, lw=0.8, alpha=0.7)
            ax_spec.set_ylabel('Flux'); ax_spec.legend(fontsize=7)
            ax_spec.set_title(f'{title}  (em=green, abs=red)')

            # score
            fin = np.isfinite(sc)
            ax_score.plot(z_arr[fin], sc[fin], lw=0.8, color='C0')
            ax_score.axvline(z_true,  color='C2', lw=1.5, ls='--',
                             label=f'true z={z_true}')
            ax_score.axvline(z_found, color='C3', lw=1.5, ls=':',
                             label=f'best z={z_found:.4f}')
            ax_score.set_xlabel('Redshift'); ax_score.set_ylabel('Score')
            ax_score.legend(fontsize=7)

        _plot_pair(axes[0,0], axes[0,1], wave_em, flux_em, pf_em,
                   z_em, sc_em, Z_EM, z_best_em, EM_LINES,
                   'Test 1: Emission')
        _plot_pair(axes[1,0], axes[1,1], wave_abs, flux_abs, pf_abs,
                   z_abs, sc_abs, Z_ABS, z_best_abs, ABS_LINES,
                   'Test 2: Absorption')
        _plot_pair(axes[2,0], axes[2,1], wave_mix, flux_mix, pf_mix,
                   z_mix, sc_mix, Z_MIX, z_best_mix, MIX_LINES,
                   'Test 3: Mixed')
        _plot_pair(axes[3,0], axes[3,1], wave_c, flux_csub, pf_c,
                   z_c, sc_c, Z_CONT, z_best_c, EM_LINES,
                   'Test 4: Raw flux + continuum fit',
                   flux_raw2=flux_c - np.median(flux_c),
                   label2='raw flux (shifted)')

        plt.tight_layout()
        plt.savefig('picket_fence_test.png', dpi=120)
        print("\nPlot saved to picket_fence_test.png")
        plt.show()
    except ImportError:
        print("\nmatplotlib not available — skipping plot")

    print("\n" + "=" * 55)
    print("All 4 tests passed.")
    sys.exit(0)
