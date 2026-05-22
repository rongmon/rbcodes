"""
ChannelSlider — wavelength channel navigation + image mode selector.

Phase 5.
"""
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QSlider, QLabel,
                              QPushButton, QButtonGroup, QComboBox)
from PyQt5.QtCore import Qt, pyqtSignal


IMAGE_MODES = ['Channel', 'Whitelight', 'Narrowband', 'Cont-sub']


class ChannelSlider(QWidget):
    """
    Horizontal bar with:
      [Channel] [Whitelight] [Narrowband] [Cont-sub]  ◄  ────slider────  ►  Channel: 42/1024  λ=5432.1 Å

    Signals
    -------
    channel_changed(int, float)   — (channel_index, wavelength_Å)
    mode_changed(str)             — new mode name
    """

    channel_changed = pyqtSignal(int, float)
    mode_changed    = pyqtSignal(str)
    method_changed  = pyqtSignal(str)   # 'mean', 'sum', 'median'

    def __init__(self, parent=None):
        super().__init__(parent)
        self._wave = None    # np.ndarray set by set_cube()
        self._n    = 0

        self._build_ui()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _build_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 2, 6, 2)
        layout.setSpacing(4)

        # Mode buttons
        self._mode_group = QButtonGroup(self)
        self._mode_group.setExclusive(True)
        self._mode_btns = {}
        for mode in IMAGE_MODES:
            btn = QPushButton(mode)
            btn.setCheckable(True)
            btn.setFixedWidth(84)
            self._mode_group.addButton(btn)
            self._mode_btns[mode] = btn
            layout.addWidget(btn)
        self._mode_btns['Whitelight'].setChecked(True)
        self._mode_group.buttonClicked.connect(self._on_mode_clicked)

        # Collapse method (used by Whitelight / Narrowband / Cont-sub)
        layout.addWidget(QLabel("Method:"))
        self._method_box = QComboBox()
        self._method_box.addItems(['mean', 'sum', 'median'])
        self._method_box.setFixedWidth(72)
        self._method_box.currentTextChanged.connect(self._on_method_changed)
        layout.addWidget(self._method_box)

        # Step back
        self._prev_btn = QPushButton("◄")
        self._prev_btn.setFixedWidth(28)
        self._prev_btn.clicked.connect(self._step_back)
        layout.addWidget(self._prev_btn)

        # Slider
        self._slider = QSlider(Qt.Horizontal)
        self._slider.setMinimum(0)
        self._slider.setMaximum(0)
        self._slider.setValue(0)
        self._slider.setTickInterval(1)
        self._slider.valueChanged.connect(self._on_slider_changed)
        layout.addWidget(self._slider, stretch=1)

        # Step forward
        self._next_btn = QPushButton("►")
        self._next_btn.setFixedWidth(28)
        self._next_btn.clicked.connect(self._step_forward)
        layout.addWidget(self._next_btn)

        # Info label
        self._info_label = QLabel("No cube loaded")
        self._info_label.setFixedWidth(240)
        layout.addWidget(self._info_label)

        # Start in Whitelight: slider disabled, method enabled
        self._set_slider_enabled(False)
        self._method_box.setEnabled(True)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_cube(self, cube):
        """Configure slider range from *cube*. Pass None to reset."""
        if cube is None or not hasattr(cube, 'wave') or cube.wave is None:
            self._wave = None
            self._n    = 0
            self._slider.setMaximum(0)
            self._info_label.setText("No cube loaded")
            self._set_slider_enabled(False)
            # Force Whitelight mode (no cube = no channel)
            self._mode_btns['Whitelight'].setChecked(True)
            for mode in ('Channel', 'Narrowband', 'Cont-sub'):
                self._mode_btns[mode].setEnabled(False)
            return

        self._wave = cube.wave
        self._n    = len(cube.wave)

        for mode in IMAGE_MODES:
            self._mode_btns[mode].setEnabled(True)

        self._slider.setMaximum(self._n - 1)
        self._slider.setValue(0)
        self._update_label(0)

        # Enable slider only when in Channel mode
        in_channel = self._mode_btns['Channel'].isChecked()
        self._set_slider_enabled(in_channel)

    @property
    def current_index(self):
        return self._slider.value()

    @property
    def current_mode(self):
        for mode, btn in self._mode_btns.items():
            if btn.isChecked():
                return mode
        return 'Whitelight'

    @property
    def current_method(self):
        return self._method_box.currentText()

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_slider_changed(self, idx):
        if self._wave is None:
            return
        wav = float(self._wave[idx])
        self._update_label(idx)
        self.channel_changed.emit(idx, wav)

    def _on_mode_clicked(self, btn):
        mode = btn.text()
        in_channel = (mode == 'Channel')
        in_band    = mode in ('Narrowband', 'Cont-sub')
        # Slider active in Channel mode (browse slices) AND in band modes
        # (slider position drives on-band center wavelength)
        self._set_slider_enabled(in_channel or in_band)
        self._method_box.setEnabled(not in_channel)

        # Cont-sub: sum is meaningless when on-band and cont windows differ in width
        current = self._method_box.currentText()
        self._method_box.blockSignals(True)
        self._method_box.clear()
        if mode == 'Cont-sub':
            self._method_box.addItems(['mean', 'median'])
            self._method_box.setCurrentText(current if current in ('mean', 'median') else 'mean')
        else:
            self._method_box.addItems(['mean', 'sum', 'median'])
            self._method_box.setCurrentText(current if current in ('mean', 'sum', 'median') else 'mean')
        self._method_box.blockSignals(False)

        self.mode_changed.emit(mode)

    def _on_method_changed(self, method):
        self.method_changed.emit(method)

    def _step_back(self):
        v = self._slider.value()
        if v > 0:
            self._slider.setValue(v - 1)

    def _step_forward(self):
        v = self._slider.value()
        if v < self._slider.maximum():
            self._slider.setValue(v + 1)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _update_label(self, idx):
        if self._wave is None:
            self._info_label.setText("No cube loaded")
            return
        wav = self._wave[idx]
        self._info_label.setText(
            f"Channel: {idx + 1}/{self._n}   λ = {wav:.2f} Å"
        )

    def _set_slider_enabled(self, enabled):
        self._slider.setEnabled(enabled)
        self._prev_btn.setEnabled(enabled)
        self._next_btn.setEnabled(enabled)
