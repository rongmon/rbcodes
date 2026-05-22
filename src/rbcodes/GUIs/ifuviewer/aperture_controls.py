"""
ApertureControls — QFitsView-style extraction aperture settings bar.

Three modes (matching cubespectrum behaviour):
  Single pixel   — extract the clicked pixel only
  Circular       — circular aperture of given radius
  Circular-Annular — circular aperture with annular background subtraction

Phase 9.
"""
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QLabel, QComboBox, QSpinBox)
from PyQt5.QtCore import pyqtSignal, Qt


MODES = ['Single pixel', 'Circular', 'Circular-Annular']

WEIGHTINGS = ['None', 'Var-weighted', 'Optimal (Data)', 'Optimal (Gaussian)']


class ApertureControls(QWidget):
    """
    Horizontal bar for configuring spectral extraction parameters.

    Layout (adaptive — fields appear only when relevant):
      Mode:[Circular▼]  Radius:[3]px  BG inner:[5]px  BG outer:[8]px
      Method:[sum▼]   Weighting:[None▼]

    Method is hidden when any weighted mode is active (it is irrelevant then).
    """

    settings_changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 2, 6, 2)
        layout.setSpacing(6)

        # ---- Mode ----
        layout.addWidget(_lbl("Mode:"))
        self._mode_box = QComboBox()
        self._mode_box.addItems(MODES)
        self._mode_box.setFixedWidth(150)
        self._mode_box.setToolTip(
            "Single pixel: extract one spaxel at cursor\n"
            "Circular: aperture of given radius\n"
            "Circular-Annular: circular aperture minus annular background"
        )
        self._mode_box.currentIndexChanged.connect(self._on_mode_changed)
        layout.addWidget(self._mode_box)

        # ---- Radius ----
        self._radius_label = _lbl("Radius:")
        layout.addWidget(self._radius_label)
        self._radius_spin = QSpinBox()
        self._radius_spin.setRange(1, 200)
        self._radius_spin.setValue(3)
        self._radius_spin.setSuffix(" px")
        self._radius_spin.setFixedWidth(65)
        self._radius_spin.setToolTip("Source aperture radius in pixels")
        self._radius_spin.valueChanged.connect(self.settings_changed)
        layout.addWidget(self._radius_spin)

        # ---- BG inner / outer (Circular-Annular only) ----
        self._bg_inner_label = _lbl("BG inner:")
        layout.addWidget(self._bg_inner_label)
        self._bg_inner_spin = QSpinBox()
        self._bg_inner_spin.setRange(1, 500)
        self._bg_inner_spin.setValue(5)
        self._bg_inner_spin.setSuffix(" px")
        self._bg_inner_spin.setFixedWidth(65)
        self._bg_inner_spin.setToolTip("Inner radius of background annulus in pixels")
        self._bg_inner_spin.valueChanged.connect(self.settings_changed)
        layout.addWidget(self._bg_inner_spin)

        self._bg_outer_label = _lbl("BG outer:")
        layout.addWidget(self._bg_outer_label)
        self._bg_outer_spin = QSpinBox()
        self._bg_outer_spin.setRange(2, 500)
        self._bg_outer_spin.setValue(8)
        self._bg_outer_spin.setSuffix(" px")
        self._bg_outer_spin.setFixedWidth(65)
        self._bg_outer_spin.setToolTip("Outer radius of background annulus in pixels")
        self._bg_outer_spin.valueChanged.connect(self.settings_changed)
        layout.addWidget(self._bg_outer_spin)

        # ---- Method (only shown when Weighting = None) ----
        self._method_label = _lbl("Method:")
        layout.addWidget(self._method_label)
        self._method_box = QComboBox()
        self._method_box.addItems(["sum", "mean", "median"])
        self._method_box.setFixedWidth(72)
        self._method_box.setToolTip(
            "How to collapse spaxels (only used when Weighting = None):\n"
            "  sum    — total flux (scales with aperture area)\n"
            "  mean   — average flux per spaxel\n"
            "  median — robust to bad pixels"
        )
        self._method_box.currentIndexChanged.connect(self.settings_changed)
        layout.addWidget(self._method_box)

        # ---- Extraction weighting ----
        self._weight_label = _lbl("Weighting:")
        layout.addWidget(self._weight_label)
        self._weight_box = QComboBox()
        self._weight_box.addItems(WEIGHTINGS)
        self._weight_box.setFixedWidth(155)
        self._weight_box.setToolTip(
            "None             — simple sum / mean / median\n"
            "Var-weighted     — weight each spaxel by 1/σ² per channel\n"
            "                   (requires variance cube)\n"
            "Optimal (Data)   — profile-weighted: spatial weights from\n"
            "                   the collapsed whitelight image (fast)\n"
            "Optimal (Gauss)  — profile-weighted: 2-D Gaussian fit to\n"
            "                   the whitelight image (smoother weights)\n"
            "\nOptimal modes work without a variance cube\n"
            "(no error spectrum will be returned in that case)."
        )
        self._weight_box.currentIndexChanged.connect(self._on_weight_changed)
        layout.addWidget(self._weight_box)

        layout.addStretch()

        # Start in Single pixel mode
        self._apply_mode_visibility('Single pixel')

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_mode_changed(self):
        self._apply_mode_visibility(self._mode_box.currentText())
        self.settings_changed.emit()

    def _on_weight_changed(self):
        self._update_method_visibility()
        self.settings_changed.emit()

    def _apply_mode_visibility(self, mode):
        """Show/hide widgets depending on the selected aperture mode."""
        is_circular = mode in ('Circular', 'Circular-Annular')
        is_annular  = mode == 'Circular-Annular'

        for w in (self._radius_label, self._radius_spin):
            w.setVisible(is_circular)

        for w in (self._bg_inner_label, self._bg_inner_spin,
                  self._bg_outer_label, self._bg_outer_spin):
            w.setVisible(is_annular)

        for w in (self._weight_label, self._weight_box):
            w.setVisible(is_circular)

        if not is_circular:
            self._weight_box.setCurrentIndex(0)   # reset to None

        # Method visibility depends on both mode and current weighting
        self._update_method_visibility()

    def _update_method_visibility(self):
        """Show Method only in circular modes when Weighting = None."""
        is_circular = self._mode_box.currentText() in ('Circular', 'Circular-Annular')
        show_method = is_circular and (self._weight_box.currentText() == 'None')
        self._method_label.setVisible(show_method)
        self._method_box.setVisible(show_method)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def mode(self):
        """'Single pixel', 'Circular', or 'Circular-Annular'."""
        return self._mode_box.currentText()

    @property
    def radius(self):
        """Source aperture radius in pixels (Circular modes only)."""
        return self._radius_spin.value()

    @property
    def bg_inner(self):
        """Background annulus inner radius (Circular-Annular only)."""
        return self._bg_inner_spin.value()

    @property
    def bg_outer(self):
        """Background annulus outer radius (Circular-Annular only)."""
        return self._bg_outer_spin.value()

    @property
    def method(self):
        """'sum', 'mean', or 'median' (relevant only when weighting is None)."""
        return self._method_box.currentText()

    @property
    def extraction_weighting(self):
        """'None', 'Var-weighted', 'Optimal (Data)', or 'Optimal (Gaussian)'."""
        return self._weight_box.currentText()

    def set_has_variance(self, has_var):
        """
        Enable/disable Var-weighted when no variance cube is loaded.
        Optimal modes remain selectable regardless — they just won't
        return an error spectrum without variance.
        """
        item = self._weight_box.model().item(WEIGHTINGS.index('Var-weighted'))
        if has_var:
            item.setFlags(item.flags() | Qt.ItemIsEnabled)
        else:
            item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
            if self._weight_box.currentText() == 'Var-weighted':
                self._weight_box.setCurrentIndex(0)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _lbl(text):
    lbl = QLabel(text)
    lbl.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
    return lbl
