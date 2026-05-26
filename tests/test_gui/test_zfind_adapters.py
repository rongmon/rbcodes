"""
Phase 4 tests: adapters.py

Verifies the exact output contracts that the target GUIs depend on:

  zfind_to_zgui_z
    - returns [z, z_err] list (length-2, values match solution)
    - selects correct solution via idx
    - empty solutions → [nan, nan]
    - NaN z_err propagated as-is
    - idx out of range raises IndexError

  zfind_to_multispec_z
    - returns a plain float
    - correct value for idx=0 and idx=1
    - empty solutions → nan

  absorbers_to_multispec
    - returns list of dicts with keys 'zabs', 'name', 'label'
    - only accepted_indices included
    - label format: 'z=X.XXXX (linelist_name)'
    - empty accepted_indices → []
    - all indices accepted → full list
    - order preserved (follows candidate order)
    - out-of-range indices silently ignored
"""

import math
import numpy as np
import pytest

from rbcodes.GUIs.zfind.io import (
    ZFindResult, ZSolution, AbsorberResult, AbsorberCandidate,
)
from rbcodes.GUIs.zfind.adapters import (
    zfind_to_zgui_z,
    zfind_to_multispec_z,
    absorbers_to_multispec,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_solution(z, z_err=0.001, chi2=-50.0, method='LineSearch:Test', n=2):
    return ZSolution(
        z=z, z_err=z_err, chi2_dof=chi2,
        method=method, template_type='Unknown', n_features=n,
    )


def _make_zfind_result(solutions):
    """Minimal ZFindResult with given solutions list."""
    z_arr = np.linspace(0.0, 1.0, 100)
    return ZFindResult(
        z_array=z_arr,
        chi2_curves=[{'label': 'Test', 'chi2': np.zeros(100)}],
        solutions=solutions,
        input_spec=None,
        warnings=[],
    )


def _make_absorber_result(candidates):
    """Minimal AbsorberResult with given candidates list."""
    z_arr = np.linspace(0.0, 2.0, 200)
    return AbsorberResult(
        z_array=z_arr,
        significance_curve=np.zeros(200),
        candidates=candidates,
        input_spec=None,
        warnings=[],
    )


def _make_candidate(z, sig=10.0, n_lines=2, linelist='MgII',
                    lines=None):
    return AbsorberCandidate(
        z=z,
        significance=sig,
        n_lines=n_lines,
        is_doublet=(n_lines == 2),
        linelist_name=linelist,
        lines_matched=lines or ['MgII2796', 'MgII2803'],
    )


# ---------------------------------------------------------------------------
# zfind_to_zgui_z
# ---------------------------------------------------------------------------

class TestZfindToZguiZ:

    def test_returns_list(self):
        result = _make_zfind_result([_make_solution(0.5)])
        out = zfind_to_zgui_z(result)
        assert isinstance(out, list)

    def test_list_has_two_elements(self):
        result = _make_zfind_result([_make_solution(0.5)])
        out = zfind_to_zgui_z(result)
        assert len(out) == 2

    def test_z_value_matches_solution(self):
        result = _make_zfind_result([_make_solution(0.5123)])
        z, _ = zfind_to_zgui_z(result)
        assert z == pytest.approx(0.5123)

    def test_z_err_value_matches_solution(self):
        result = _make_zfind_result([_make_solution(0.5, z_err=0.0034)])
        _, z_err = zfind_to_zgui_z(result)
        assert z_err == pytest.approx(0.0034)

    def test_default_idx_is_zero(self):
        sols = [_make_solution(0.5), _make_solution(0.8)]
        result = _make_zfind_result(sols)
        z, _ = zfind_to_zgui_z(result)
        assert z == pytest.approx(0.5)

    def test_idx_one_selects_second_solution(self):
        sols = [_make_solution(0.5), _make_solution(0.8)]
        result = _make_zfind_result(sols)
        z, _ = zfind_to_zgui_z(result, idx=1)
        assert z == pytest.approx(0.8)

    def test_empty_solutions_returns_nan_nan(self):
        result = _make_zfind_result([])
        out = zfind_to_zgui_z(result)
        assert len(out) == 2
        assert math.isnan(out[0])
        assert math.isnan(out[1])

    def test_nan_z_err_propagated(self):
        result = _make_zfind_result([_make_solution(0.5, z_err=float('nan'))])
        _, z_err = zfind_to_zgui_z(result)
        assert math.isnan(z_err)

    def test_idx_out_of_range_raises(self):
        result = _make_zfind_result([_make_solution(0.5)])
        with pytest.raises(IndexError):
            zfind_to_zgui_z(result, idx=5)

    def test_z_is_float_not_array(self):
        result = _make_zfind_result([_make_solution(0.5)])
        z, z_err = zfind_to_zgui_z(result)
        assert isinstance(z, float)
        assert isinstance(z_err, float)

    def test_zero_redshift(self):
        result = _make_zfind_result([_make_solution(0.0)])
        z, _ = zfind_to_zgui_z(result)
        assert z == pytest.approx(0.0)

    def test_high_redshift(self):
        result = _make_zfind_result([_make_solution(5.7)])
        z, _ = zfind_to_zgui_z(result)
        assert z == pytest.approx(5.7)


# ---------------------------------------------------------------------------
# zfind_to_multispec_z
# ---------------------------------------------------------------------------

class TestZfindToMultispecZ:

    def test_returns_float(self):
        result = _make_zfind_result([_make_solution(0.5)])
        out = zfind_to_multispec_z(result)
        assert isinstance(out, float)

    def test_value_matches_best_solution(self):
        result = _make_zfind_result([_make_solution(0.7321)])
        assert zfind_to_multispec_z(result) == pytest.approx(0.7321)

    def test_default_idx_is_zero(self):
        sols = [_make_solution(0.3), _make_solution(0.9)]
        result = _make_zfind_result(sols)
        assert zfind_to_multispec_z(result) == pytest.approx(0.3)

    def test_idx_one_selects_second(self):
        sols = [_make_solution(0.3), _make_solution(0.9)]
        result = _make_zfind_result(sols)
        assert zfind_to_multispec_z(result, idx=1) == pytest.approx(0.9)

    def test_empty_solutions_returns_nan(self):
        result = _make_zfind_result([])
        assert math.isnan(zfind_to_multispec_z(result))

    def test_idx_out_of_range_raises(self):
        result = _make_zfind_result([_make_solution(0.5)])
        with pytest.raises(IndexError):
            zfind_to_multispec_z(result, idx=3)


# ---------------------------------------------------------------------------
# absorbers_to_multispec
# ---------------------------------------------------------------------------

class TestAbsorbersToMultispec:

    def test_returns_list(self):
        result = _make_absorber_result([_make_candidate(0.5)])
        out = absorbers_to_multispec(result, [0])
        assert isinstance(out, list)

    def test_single_accepted_candidate(self):
        result = _make_absorber_result([_make_candidate(0.5)])
        out = absorbers_to_multispec(result, [0])
        assert len(out) == 1

    def test_dict_has_required_keys(self):
        result = _make_absorber_result([_make_candidate(0.5)])
        entry = absorbers_to_multispec(result, [0])[0]
        assert set(entry.keys()) == {'zabs', 'name', 'label'}

    def test_zabs_value_matches_candidate(self):
        result = _make_absorber_result([_make_candidate(0.5123)])
        entry = absorbers_to_multispec(result, [0])[0]
        assert entry['zabs'] == pytest.approx(0.5123)

    def test_name_matches_linelist_name(self):
        result = _make_absorber_result([_make_candidate(0.5, linelist='MgII')])
        entry = absorbers_to_multispec(result, [0])[0]
        assert entry['name'] == 'MgII'

    def test_label_format(self):
        result = _make_absorber_result([_make_candidate(0.5, linelist='MgII')])
        entry = absorbers_to_multispec(result, [0])[0]
        assert entry['label'] == 'z=0.5000 (MgII)'

    def test_label_format_four_decimal_places(self):
        result = _make_absorber_result([_make_candidate(1.23456789)])
        entry = absorbers_to_multispec(result, [0])[0]
        # label must be 4 d.p.
        assert 'z=1.2346' in entry['label']

    def test_empty_accepted_indices_returns_empty_list(self):
        result = _make_absorber_result([_make_candidate(0.5), _make_candidate(0.8)])
        out = absorbers_to_multispec(result, [])
        assert out == []

    def test_all_indices_accepted(self):
        cands = [_make_candidate(z) for z in (0.3, 0.7, 1.1)]
        result = _make_absorber_result(cands)
        out = absorbers_to_multispec(result, [0, 1, 2])
        assert len(out) == 3

    def test_subset_of_indices_accepted(self):
        cands = [_make_candidate(z) for z in (0.3, 0.7, 1.1)]
        result = _make_absorber_result(cands)
        out = absorbers_to_multispec(result, [0, 2])
        assert len(out) == 2
        zabs_vals = [e['zabs'] for e in out]
        assert pytest.approx(0.3) in zabs_vals
        assert pytest.approx(1.1) in zabs_vals
        assert not any(abs(e['zabs'] - 0.7) < 1e-9 for e in out)

    def test_order_preserved(self):
        cands = [_make_candidate(z) for z in (0.3, 0.7, 1.1)]
        result = _make_absorber_result(cands)
        out = absorbers_to_multispec(result, [0, 1, 2])
        assert [e['zabs'] for e in out] == pytest.approx([0.3, 0.7, 1.1])

    def test_out_of_range_index_silently_ignored(self):
        result = _make_absorber_result([_make_candidate(0.5)])
        # index 99 does not exist — should just be skipped
        out = absorbers_to_multispec(result, [0, 99])
        assert len(out) == 1
        assert out[0]['zabs'] == pytest.approx(0.5)

    def test_different_linelist_names_in_label(self):
        cands = [
            _make_candidate(0.5, linelist='MgII'),
            _make_candidate(1.2, linelist='CIV'),
        ]
        result = _make_absorber_result(cands)
        out = absorbers_to_multispec(result, [0, 1])
        assert 'MgII' in out[0]['label']
        assert 'CIV' in out[1]['label']

    def test_empty_candidates_empty_accepted(self):
        result = _make_absorber_result([])
        out = absorbers_to_multispec(result, [])
        assert out == []

    def test_no_extra_keys_in_dict(self):
        result = _make_absorber_result([_make_candidate(0.5)])
        entry = absorbers_to_multispec(result, [0])[0]
        assert len(entry) == 3
