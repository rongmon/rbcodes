"""
Unit tests for rb_setline.py

These tests run against the line lists shipped with the package (no mocking),
so they also check that every packaged list is found and parsed.
"""

import logging
from unittest.mock import patch

import numpy as np
import pytest

from rbcodes.IGM import rb_setline as rs
from rbcodes.IGM.rb_setline import rb_setline, read_line_list

logging.disable(logging.CRITICAL)

ALL_LABELS = ['atom', 'LLS', 'LLS Small', 'DLA', 'LBG', 'Gal', 'Eiger_Strong',
              'Gal_Em', 'Gal_Abs', 'Gal_long', 'AGN', 'HI_recomb',
              'HI_recomb_light', 'HI', 'EUV', 'LLS_EUV']


@pytest.fixture(autouse=True)
def clear_cache():
    rs._LINE_LIST_CACHE.clear()
    yield
    rs._LINE_LIST_CACHE.clear()


class TestRbSetline:

    def test_closest_match(self):
        result = rb_setline(2796.3, 'closest')
        assert result['name'] == 'MgII 2796'
        assert result['wave'] == pytest.approx(2796.352)
        assert result['fval'] == pytest.approx(0.6123)
        assert result['gamma'] == pytest.approx(2.612e8)

        result = rb_setline(1550.7, 'closest')
        assert result['name'] == 'CIV 1550'

    def test_exact_match(self):
        result = rb_setline(1215.6701, 'Exact')
        assert len(result['wave']) == 1
        assert result['name'][0] == 'HI 1215'
        assert result['fval'][0] == pytest.approx(0.4164)
        assert 'gamma' in result

    def test_exact_no_match_returns_empty(self):
        result = rb_setline(1215.0, 'Exact')
        assert len(result['wave']) == 0

    def test_name_match(self):
        result = rb_setline(0, 'Name', target_name='HI 1215')
        assert len(result['wave']) == 1
        assert result['wave'][0] == pytest.approx(1215.6701)

    def test_name_no_match_returns_empty(self):
        result = rb_setline(0, 'Name', target_name='NotAnIon 1234')
        assert len(result['wave']) == 0

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            rb_setline(1215.67, 'invalid_method')

    def test_name_without_target(self):
        with pytest.raises(ValueError):
            rb_setline(1215.67, 'Name')

    def test_different_line_list(self):
        result = rb_setline(1215.67, 'closest', linelist='LLS')
        assert result['wave'] == pytest.approx(1215.6701)
        assert result['name'] == 'HI 1215'
        assert 'gamma' not in result

    def test_line_list_caching(self):
        first = read_line_list('atom')
        assert 'atom' in rs._LINE_LIST_CACHE
        with patch('astropy.io.ascii.read') as mock_read:
            second = read_line_list('atom')
            mock_read.assert_not_called()
        assert second is first


class TestReadLineList:

    @pytest.mark.parametrize('label', ALL_LABELS)
    def test_every_packaged_list_loads(self, label):
        data = read_line_list(label)
        assert len(data) > 0
        for entry in data:
            assert {'wrest', 'ion', 'fval'} <= set(entry)
            assert np.isfinite(entry['wrest']) and entry['wrest'] > 0

    def test_read_atom_line_list(self):
        data = read_line_list('atom')
        assert all('gamma' in d for d in data)
        assert any(d['ion'] == 'HI 1215' for d in data)

    def test_recomb_list_in_angstrom(self):
        # HI_recomb lists are stored in micron and converted to Angstrom
        waves = np.array([d['wrest'] for d in read_line_list('HI_recomb')])
        assert waves.min() > 900.

    def test_invalid_line_list(self):
        with pytest.raises(ValueError):
            read_line_list('invalid_list')

    def test_file_not_found(self, tmp_path):
        with patch.object(rs, 'files', return_value=tmp_path):
            with pytest.raises(FileNotFoundError):
                read_line_list('atom')
