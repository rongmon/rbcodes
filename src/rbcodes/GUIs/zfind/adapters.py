"""
adapters.py — Thin translators between zfind output and target GUI formats.

These are the only files that know about both zfind output AND target GUI internals.
The engine and io modules never import anything GUI-specific.
"""

from rbcodes.GUIs.zfind.io import ZFindResult, AbsorberResult


# ---------------------------------------------------------------------------
# zgui adapter  (emission mode → rb_zgui)
# ---------------------------------------------------------------------------

def zfind_to_zgui_z(result: ZFindResult, idx: int = 0) -> list:
    """
    Translate a ZFindResult into the format expected by
    LineListWidget._on_estZ_changed([z, z_err]).

    Parameters
    ----------
    result : ZFindResult
    idx : int  — which solution to use (0 = best)

    Returns
    -------
    [z, z_err]  — matches _on_estZ_changed() signature exactly
    """
    if not result.solutions:
        return [float('nan'), float('nan')]
    idx = min(max(idx, 0), len(result.solutions) - 1)
    s = result.solutions[idx]
    return [s.z, s.z_err]


# ---------------------------------------------------------------------------
# rb_multispec adapter  (emission + absorption → rb_multispec)
# ---------------------------------------------------------------------------

def zfind_to_multispec_z(result: ZFindResult, idx: int = 0) -> float:
    """
    Return single best redshift for rb_multispec RedshiftInputWidget.

    Parameters
    ----------
    result : ZFindResult
    idx : int  — which solution to use

    Returns
    -------
    float — redshift value
    """
    if not result.solutions:
        return float('nan')
    return result.solutions[idx].z


def absorbers_to_multispec(result: AbsorberResult,
                           accepted_indices: list) -> list:
    """
    Translate accepted AbsorberCandidates into the dict format
    expected by rb_multispec AbsorberManager.

    Parameters
    ----------
    result : AbsorberResult
    accepted_indices : list of int — indices into result.candidates to include

    Returns
    -------
    list of {'zabs': float, 'name': str, 'label': str}
    """
    output = []
    for i, candidate in enumerate(result.candidates):
        if i in accepted_indices:
            output.append({
                'zabs' : candidate.z,
                'name' : candidate.linelist_name,
                'label': f'z={candidate.z:.4f} ({candidate.linelist_name})'
            })
    return output
