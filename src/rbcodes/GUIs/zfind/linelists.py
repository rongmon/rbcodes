"""
linelists.py — Curated line lists for rb_zfind.

Five science-case presets:

  zfind_em      — galaxy nebular emission lines (z~0–7)
  zfind_stellar — stellar absorption lines for passive/SF galaxies
  zfind_igm     — IGM/CGM intervening absorbers on QSO sightlines
  zfind_galaxy  — combined nebular emission + stellar absorption (one-stop galaxy preset)
  zfind_qso     — QSO/AGN broad emission lines

Each line carries four fields:
  wave   (float) — vacuum rest wavelength in Å
  name   (str)   — human-readable label
  weight (float) — relative line strength (see notes below)
  type   (str)   — 'emission' or 'absorption'

Weight convention
-----------------
Absorption lines:
  weight = fval (oscillator strength) taken directly from atom_full.dat.
  Physically motivated — fval determines line depth for a given column density.
  Source column in comments: [atom] = atom_full.dat, [manual] = no fval available.

Emission lines:
  HI recombination lines (Lyα, Hα, Hβ, Hγ): fval from atom_full.dat,
  normalised so Hα = 3.0  (scale factor = 3.0 / 0.6958 = 4.31).
  Forbidden/nebular lines ([OII], [OIII], [NII], [SII], CIII]):
  no fval exists — weights assigned empirically on the same 1–3 scale
  based on typical observed brightness relative to Hα in SF galaxies.
  Source column in comments: [fval×4.31] or [empirical].

To add a line to a preset
-------------------------
Add a dict entry to the relevant list below:
    {'wave': <Å>, 'name': '<label>', 'weight': <float>, 'type': 'emission'|'absorption'}
No other file needs changing.

Public API
----------
get_curated_df(name)  → pandas DataFrame(['wave','name','weight','type'])
                         with df.attrs['name'] set.
CURATED_NAMES         → list of preset name strings
"""

import pandas as pd

# ---------------------------------------------------------------------------
# zfind_em — galaxy nebular emission (pure emission preset)
# ---------------------------------------------------------------------------
# Use for: star-forming galaxies, emission-line galaxies, z~0–7.
# Lyα weight deliberately moderate (1.8) — useful at z<2, IGM-suppressed at z>3.
# ---------------------------------------------------------------------------
_EM_LINES = [
    # wave       name            weight  type        weight source
    {'wave': 1216.00, 'name': 'Lya 1216',    'weight': 1.8,  'type': 'emission'},  # fval×4.31
    {'wave': 1908.73, 'name': 'CIII] 1909',  'weight': 1.0,  'type': 'emission'},  # empirical
    {'wave': 3727.09, 'name': '[OII] 3727',  'weight': 2.0,  'type': 'emission'},  # empirical
    {'wave': 3729.87, 'name': '[OII] 3729',  'weight': 2.0,  'type': 'emission'},  # empirical
    {'wave': 4341.69, 'name': 'Hgamma',      'weight': 0.2,  'type': 'emission'},  # fval×4.31
    {'wave': 4862.68, 'name': 'Hbeta',       'weight': 0.5,  'type': 'emission'},  # fval×4.31
    {'wave': 4960.30, 'name': '[OIII] 4960', 'weight': 1.0,  'type': 'emission'},  # empirical (1/3 of 5007)
    {'wave': 5008.24, 'name': '[OIII] 5007', 'weight': 3.0,  'type': 'emission'},  # empirical
    {'wave': 6564.61, 'name': 'Halpha',      'weight': 3.0,  'type': 'emission'},  # fval×4.31 (anchor)
    {'wave': 6585.28, 'name': '[NII] 6583',  'weight': 1.0,  'type': 'emission'},  # empirical
]

# ---------------------------------------------------------------------------
# zfind_stellar — stellar absorption lines for galaxy redshifts
# ---------------------------------------------------------------------------
# Use for: passive/early-type galaxies, post-starburst, any galaxy with
# strong stellar continuum. Works well even without emission lines.
# G-band and MgI b have no fval (molecular/subordinate) — manual weights.
# ---------------------------------------------------------------------------
_STELLAR_LINES = [
    # wave       name            weight  type          weight source
    {'wave': 2796.35, 'name': 'MgII 2796', 'weight': 0.612, 'type': 'absorption'},  # atom_full.dat
    {'wave': 2803.53, 'name': 'MgII 2803', 'weight': 0.305, 'type': 'absorption'},  # atom_full.dat
    {'wave': 3934.78, 'name': 'CaII K',    'weight': 0.635, 'type': 'absorption'},  # atom_full.dat
    {'wave': 3969.59, 'name': 'CaII H',    'weight': 0.315, 'type': 'absorption'},  # atom_full.dat
    {'wave': 4305.61, 'name': 'G-band',    'weight': 0.300, 'type': 'absorption'},  # manual (CH molecular)
    {'wave': 5175.00, 'name': 'MgI b',     'weight': 0.300, 'type': 'absorption'},  # manual (subordinate)
    {'wave': 5891.58, 'name': 'NaI D1',    'weight': 0.631, 'type': 'absorption'},  # atom_full.dat
    {'wave': 5897.56, 'name': 'NaI D2',    'weight': 0.318, 'type': 'absorption'},  # atom_full.dat
]

# ---------------------------------------------------------------------------
# zfind_igm — IGM/CGM intervening absorbers
# ---------------------------------------------------------------------------
# Use for: finding intervening absorbers on QSO sightlines.
# All fval from atom_full.dat. Lines ordered by rest wavelength.
# SiII 1260 is the strongest single optical IGM line (fval=1.007).
# HI Lyα is almost always saturated for strong absorbers but position is reliable.
# ---------------------------------------------------------------------------
_IGM_LINES = [
    # wave       name            weight  type          weight source
    {'wave': 1031.93, 'name': 'OVI 1031',  'weight': 0.133, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1037.62, 'name': 'OVI 1037',  'weight': 0.066, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1215.67, 'name': 'HI Lya',    'weight': 0.416, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1260.42, 'name': 'SiII 1260', 'weight': 1.007, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1334.53, 'name': 'CII 1334',  'weight': 0.128, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1393.76, 'name': 'SiIV 1393', 'weight': 0.514, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1402.77, 'name': 'SiIV 1402', 'weight': 0.255, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1548.20, 'name': 'CIV 1548',  'weight': 0.191, 'type': 'absorption'},  # atom_full.dat
    {'wave': 1550.78, 'name': 'CIV 1550',  'weight': 0.095, 'type': 'absorption'},  # atom_full.dat
    {'wave': 2796.35, 'name': 'MgII 2796', 'weight': 0.612, 'type': 'absorption'},  # atom_full.dat
    {'wave': 2803.53, 'name': 'MgII 2803', 'weight': 0.305, 'type': 'absorption'},  # atom_full.dat
]

# ---------------------------------------------------------------------------
# zfind_galaxy — combined emission + stellar absorption
# ---------------------------------------------------------------------------
# Use for: general galaxy redshift search when you don't know whether
# emission or absorption lines will dominate. Covers both SF and passive.
# Lyα excluded — complex resonant scattering in galaxy ISM, unreliable.
# MgII excluded — can be emission or absorption in galaxies, ambiguous.
# ---------------------------------------------------------------------------
_GALAXY_LINES = [
    # --- nebular emission ---
    # wave       name            weight  type          weight source
    {'wave': 3727.09, 'name': '[OII] 3727',  'weight': 2.0,  'type': 'emission'},  # empirical
    {'wave': 3729.87, 'name': '[OII] 3729',  'weight': 2.0,  'type': 'emission'},  # empirical
    {'wave': 4341.69, 'name': 'Hgamma',      'weight': 0.2,  'type': 'emission'},  # fval×4.31
    {'wave': 4862.68, 'name': 'Hbeta',       'weight': 0.5,  'type': 'emission'},  # fval×4.31
    {'wave': 4960.30, 'name': '[OIII] 4960', 'weight': 1.0,  'type': 'emission'},  # empirical
    {'wave': 5008.24, 'name': '[OIII] 5007', 'weight': 3.0,  'type': 'emission'},  # empirical
    {'wave': 6564.61, 'name': 'Halpha',      'weight': 3.0,  'type': 'emission'},  # fval×4.31
    {'wave': 6585.28, 'name': '[NII] 6583',  'weight': 1.0,  'type': 'emission'},  # empirical
    # --- stellar absorption ---
    {'wave': 3934.78, 'name': 'CaII K',      'weight': 0.635,'type': 'absorption'}, # atom_full.dat
    {'wave': 3969.59, 'name': 'CaII H',      'weight': 0.315,'type': 'absorption'}, # atom_full.dat
    {'wave': 4305.61, 'name': 'G-band',      'weight': 0.300,'type': 'absorption'}, # manual
    {'wave': 5175.00, 'name': 'MgI b',       'weight': 0.300,'type': 'absorption'}, # manual
    {'wave': 5891.58, 'name': 'NaI D1',      'weight': 0.631,'type': 'absorption'}, # atom_full.dat
    {'wave': 5897.56, 'name': 'NaI D2',      'weight': 0.318,'type': 'absorption'}, # atom_full.dat
]

# ---------------------------------------------------------------------------
# zfind_qso — QSO/AGN broad emission lines
# ---------------------------------------------------------------------------
# Use for: finding redshifts of type-1 QSOs and broad-line AGN.
# All lines are emission from the QSO itself — not intervening absorbers
# (use zfind_igm for those at a different z).
# CIV 1548/1550 kept as two entries — broad profile covers both; picket
# fence window will catch whichever is stronger.
# MgII listed as blend centre 2799 Å — doublet unresolved in broad-line QSOs.
# Hβ and [OIII] useful for lower-z (z<1) type-1 AGN / NLS1.
# ---------------------------------------------------------------------------
_QSO_LINES = [
    # wave       name            weight  type          weight source
    {'wave': 1216.00, 'name': 'Lya 1216',    'weight': 3.0,  'type': 'emission'},  # empirical (brightest QSO line)
    {'wave': 1548.20, 'name': 'CIV 1548',    'weight': 3.0,  'type': 'emission'},  # empirical
    {'wave': 1550.78, 'name': 'CIV 1550',    'weight': 1.5,  'type': 'emission'},  # empirical (~half of 1548)
    {'wave': 1908.73, 'name': 'CIII] 1909',  'weight': 2.0,  'type': 'emission'},  # empirical
    {'wave': 2799.00, 'name': 'MgII 2799',   'weight': 2.0,  'type': 'emission'},  # empirical (blend centre)
    {'wave': 4862.68, 'name': 'Hbeta',       'weight': 0.5,  'type': 'emission'},  # fval×4.31
    {'wave': 4960.30, 'name': '[OIII] 4960', 'weight': 0.5,  'type': 'emission'},  # empirical (weaker in QSOs)
    {'wave': 5008.24, 'name': '[OIII] 5007', 'weight': 1.5,  'type': 'emission'},  # empirical
]

# ---------------------------------------------------------------------------
# Registry — maps preset name → line list
# ---------------------------------------------------------------------------

_REGISTRY = {
    'zfind_em':      _EM_LINES,
    'zfind_stellar': _STELLAR_LINES,
    'zfind_igm':     _IGM_LINES,
    'zfind_galaxy':  _GALAXY_LINES,
    'zfind_qso':     _QSO_LINES,
}

CURATED_NAMES = list(_REGISTRY.keys())


def get_curated_df(name: str) -> pd.DataFrame:
    """
    Return a pandas DataFrame for the named curated linelist.

    Parameters
    ----------
    name : str
        One of: 'zfind_em', 'zfind_stellar', 'zfind_igm',
                'zfind_galaxy', 'zfind_qso'

    Returns
    -------
    pandas.DataFrame
        Columns: ['wave', 'name', 'weight', 'type']
        df.attrs['name'] is set to the preset name.

    Raises
    ------
    KeyError if name is not a known curated preset.
    """
    if name not in _REGISTRY:
        raise KeyError(
            f'{name!r} is not a curated preset. '
            f'Available: {CURATED_NAMES}'
        )
    df = pd.DataFrame(_REGISTRY[name])
    df.attrs['name'] = name
    return df
