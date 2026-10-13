"""
LoadOptionsDialog — choose a wavelength range (and air -> vacuum) before loading a cube.

Shown automatically for cubes larger than ``LARGE_FILE_BYTES`` and from File -> "Open Wavelength
Range…". Reads only the header to show the cube's wavelength coverage.
"""
import os

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QRadioButton,
                             QDoubleSpinBox, QCheckBox, QDialogButtonBox)

LARGE_FILE_BYTES = 1_000_000_000


def cube_header_info(path):
    """(wmin, wmax, airvac, is_cube) from the headers only."""
    from astropy.io import fits
    from rbcodes.ifu.io.auto_cube import _build_wave
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.header.get('NAXIS', 0) == 3:
                w = _build_wave(hdu.header)
                ctype = str(hdu.header.get('CTYPE3', '')).strip().upper()
                return float(w.min()), float(w.max()), {'AWAV': 'air', 'WAVE': 'vac'}.get(ctype, 'unknown'), True
    return None, None, None, False


def is_large(path):
    try:
        return os.path.getsize(path) > LARGE_FILE_BYTES
    except OSError:
        return False


class LoadOptionsDialog(QDialog):
    """Full cube or a wavelength range; optional air -> vacuum conversion."""

    def __init__(self, path, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Load options')
        wmin, wmax, airvac, _ = cube_header_info(path)
        size = os.path.getsize(path) / 1e9
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(f'{os.path.basename(path)}  ({size:.2f} GB)\n'
                             f'wavelengths {wmin:.1f}–{wmax:.1f} Å ({airvac})'))
        self.full = QRadioButton('Load the full cube')
        self.part = QRadioButton('Load only a wavelength range (Å):')
        self.part.setChecked(is_large(path))
        self.full.setChecked(not is_large(path))
        lay.addWidget(self.full)
        row = QHBoxLayout()
        row.addWidget(self.part)
        self.lo, self.hi = QDoubleSpinBox(), QDoubleSpinBox()
        for sb, val in ((self.lo, wmin), (self.hi, min(wmax, wmin + 200.))):
            sb.setDecimals(2)
            sb.setRange(wmin, wmax)
            sb.setValue(val)
            row.addWidget(sb)
        lay.addLayout(row)
        self.vac = QCheckBox('Convert air → vacuum wavelengths')
        self.vac.setEnabled(airvac == 'air')
        lay.addWidget(self.vac)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def options(self):
        """dict(wave_range=..., air_to_vacuum=...) for ``load_fits``."""
        wr = (self.lo.value(), self.hi.value()) if self.part.isChecked() else None
        return dict(wave_range=wr, air_to_vacuum=self.vac.isChecked())


def ask_load_options(path, parent=None, always=False):
    """
    Options for loading ``path``: {} without asking for small files (unless ``always``);
    None if the user cancelled.
    """
    if not (always or is_large(path)):
        return {}
    try:
        if not cube_header_info(path)[3]:
            return {}
    except Exception:
        return {}
    dlg = LoadOptionsDialog(path, parent)
    return dlg.options() if dlg.exec_() == QDialog.Accepted else None
