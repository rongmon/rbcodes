"""
dialog.py — PyQt5 popup dialog for rb_zfind.

ZFindDialog is the single dialog class for both emission and absorption modes.
The 'mode' parameter controls which UI elements are shown and which adapters are used.

Stub: full implementation in Phase 5 (emission) and Phase 9 (absorption).
"""

# Full implementation: Phase 5 (ZFIND_TODO.md)
# This stub exists so imports don't fail during early phases.

try:
    from PyQt5.QtWidgets import QDialog
    from PyQt5.QtCore import pyqtSignal

    class ZFindDialog(QDialog):
        """
        Self-contained redshift/absorber finder popup.

        Parameters
        ----------
        spec : rb_spectrum or None
        mode : 'emission' or 'absorption'
        default_linelist : str  — linelist name to pre-select
        z_qso : float or None   — for absorption mode, caps z_max

        Signals
        -------
        accepted_z : list          — [z, z_err] for emission mode → zgui
        accepted_absorbers : list  — [{zabs, name, label}, ...] for absorption → multispec
        """

        accepted_z         = pyqtSignal(list)
        accepted_absorbers = pyqtSignal(list)

        def __init__(self, spec=None, mode='emission', default_linelist=None,
                     z_qso=None, parent=None):
            super().__init__(parent)
            self.spec = spec
            self.mode = mode
            self.default_linelist = default_linelist
            self.z_qso = z_qso

            self.setWindowTitle('Redshift Estimator' if mode == 'emission'
                                else 'Absorber Finder')
            self.resize(900, 700)

            # TODO: build full UI in Phase 5
            raise NotImplementedError(
                "ZFindDialog UI not yet implemented. See Phase 5 in ZFIND_TODO.md."
            )

except ImportError:
    # PyQt5 not available — headless/engine-only use is still possible
    pass
