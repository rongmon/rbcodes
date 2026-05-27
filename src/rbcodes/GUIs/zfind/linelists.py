"""
linelists.py — Curated line lists for rb_zfind.

Provides two built-in presets that appear at the top of the dialog dropdown:

  zfind_em   — 7 strong emission lines (galaxy / QSO redshift finding)
  zfind_abs  — 9 strong absorption doublets (CGM / IGM absorber finding)

These are defined as plain Python dicts so the module has no file I/O and
no dependency on rb_setline.  All existing rb_setline linelists are still
available in the dialog under the 'Other linelists' section.

Public API
----------
get_curated_df(name)   → pandas DataFrame with columns ['wave', 'name']
                          and df.attrs['name'] set.
CURATED_NAMES          → list of curated preset names
"""

import pandas as pd

# ---------------------------------------------------------------------------
# Curated emission lines  (≈ 7 strongest galaxy/QSO emission features)
# ---------------------------------------------------------------------------
# Selected on the basis of:
#   • High intrinsic luminosity / equivalent width
#   • Commonly detected across 0 < z < 7
#   • Distinctive spacing — minimises accidental chi2 minima
# ---------------------------------------------------------------------------
_EMISSION_LINES = [
    {'wave': 1216.00, 'name': 'Lya 1216'},
    {'wave': 1908.73, 'name': 'CIII] 1909'},
    {'wave': 3728.00, 'name': '[OII] 3728'},    # blend of 3727+3729
    {'wave': 4862.68, 'name': 'Hb 4863'},
    {'wave': 4960.30, 'name': '[OIII] 4960'},
    {'wave': 5008.24, 'name': '[OIII] 5008'},
    {'wave': 6549.85, 'name': '[NII] 6550'},
    {'wave': 6564.61, 'name': 'Ha 6565'},
    {'wave': 6585.28, 'name': '[NII] 6585'},
    {'wave': 6718.29, 'name': '[SII] 6718'},
    {'wave': 6732.67, 'name': '[SII] 6733'},
]

# ---------------------------------------------------------------------------
# Curated absorption lines  (≈ 9 strong CGM/IGM absorber doublets/singlets)
# ---------------------------------------------------------------------------
# Selected on the basis of:
#   • Large oscillator strength (f-value) → deep absorption
#   • Doublet pairs give two corroborating detections
#   • Cover a wide ionisation range: HI → OVI
# ---------------------------------------------------------------------------
_ABSORPTION_LINES = [
    {'wave': 1031.93, 'name': 'OVI 1031'},
    {'wave': 1037.62, 'name': 'OVI 1037'},
    {'wave': 1215.67, 'name': 'HI 1216'},
    {'wave': 1393.76, 'name': 'SiIV 1393'},
    {'wave': 1402.77, 'name': 'SiIV 1402'},
    {'wave': 1548.20, 'name': 'CIV 1548'},
    {'wave': 1550.78, 'name': 'CIV 1550'},
    {'wave': 2796.35, 'name': 'MgII 2796'},
    {'wave': 2803.53, 'name': 'MgII 2803'},
]

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_REGISTRY = {
    'zfind_em':  _EMISSION_LINES,
    'zfind_abs': _ABSORPTION_LINES,
}

CURATED_NAMES = list(_REGISTRY.keys())


def get_curated_df(name: str) -> pd.DataFrame:
    """
    Return a pandas DataFrame for the named curated linelist.

    Parameters
    ----------
    name : str   — 'zfind_em' or 'zfind_abs'

    Returns
    -------
    DataFrame with columns ['wave', 'name'] and df.attrs['name'] set.

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
