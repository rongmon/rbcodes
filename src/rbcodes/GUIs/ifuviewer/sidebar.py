"""
DatasetSidebar — list of loaded FITS cubes/images, active dataset selection.
"""
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout,
                              QPushButton, QListWidget, QListWidgetItem,
                              QFileDialog, QLabel, QMenu, QAction)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor

from rbcodes.GUIs.ifuviewer.io.auto_cube import load_fits
from rbcodes.GUIs.ifuviewer.io.image2d import FITSImage


class DatasetSidebar(QWidget):
    """
    Left-panel widget: lists loaded datasets, lets the user add/remove them,
    and emits ``dataset_changed`` when the active dataset switches.

    Signals
    -------
    dataset_changed(object)
        Emitted with the newly active IFUCube (or FITSImage) whenever the user
        clicks a different item.  Emits None when the list becomes empty.
    """

    dataset_changed  = pyqtSignal(object)
    use_as_variance  = pyqtSignal(object, object)   # (var_cube, target_cube)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._cubes = []   # parallel list of IFUCube / FITSImage objects

        self._build_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        header = QLabel("Datasets")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)

        self._list = QListWidget()
        self._list.setAlternatingRowColors(True)
        self._list.currentRowChanged.connect(self._on_row_changed)
        self._list.setContextMenuPolicy(Qt.CustomContextMenu)
        self._list.customContextMenuRequested.connect(self._on_context_menu)
        layout.addWidget(self._list)

        btn_row = QHBoxLayout()
        btn_row.setContentsMargins(0, 0, 0, 0)

        self._btn_add = QPushButton("+ Add")
        self._btn_add.setToolTip("Load a FITS cube or 2D image")
        self._btn_add.clicked.connect(self._on_add)

        self._btn_remove = QPushButton("− Remove")
        self._btn_remove.setToolTip("Remove the selected dataset")
        self._btn_remove.clicked.connect(self._on_remove)

        btn_row.addWidget(self._btn_add)
        btn_row.addWidget(self._btn_remove)
        layout.addLayout(btn_row)

        self.setFixedWidth(220)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def add_dataset(self, path, var_path=None):
        """
        Load *path* (and optional *var_path*) and append to the list.
        Selects the new item automatically.
        """
        try:
            cube = load_fits(path, var=var_path)
        except Exception as exc:
            # Caller (or MainWindow) should handle the error; re-raise
            raise RuntimeError(f"Cannot load '{path}': {exc}") from exc

        self._cubes.append(cube)
        item = QListWidgetItem(_item_text(cube))
        item.setToolTip(str(path))
        self._list.addItem(item)
        self._list.setCurrentRow(len(self._cubes) - 1)   # triggers dataset_changed

    @property
    def active_cube(self):
        """Currently selected IFUCube / FITSImage, or None."""
        row = self._list.currentRow()
        if row < 0 or row >= len(self._cubes):
            return None
        return self._cubes[row]

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_add(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Open FITS file(s)",
            "",
            "FITS files (*.fits *.fit *.fits.gz *.fit.gz);;All files (*)",
        )
        for p in paths:
            try:
                self.add_dataset(p)
            except RuntimeError as exc:
                # Surface the error via status bar / message box in MainWindow
                self.dataset_changed.emit(None)   # signal with None = error sentinel
                print(f"[ifuviewer] {exc}")

    def _on_remove(self):
        row = self._list.currentRow()
        if row < 0:
            return
        self._list.takeItem(row)
        self._cubes.pop(row)

        if self._cubes:
            new_row = min(row, len(self._cubes) - 1)
            self._list.setCurrentRow(new_row)
        else:
            self.dataset_changed.emit(None)

    def _on_row_changed(self, row):
        if row < 0 or row >= len(self._cubes):
            return
        self.dataset_changed.emit(self._cubes[row])

    def _on_context_menu(self, pos):
        """Right-click context menu on a list item."""
        item = self._list.itemAt(pos)
        if item is None:
            return
        row = self._list.row(item)
        if row < 0 or row >= len(self._cubes):
            return

        var_cube = self._cubes[row]
        # Only IFU cubes can act as variance (need a flux array)
        if isinstance(var_cube, FITSImage):
            return

        # Build candidate target cubes = all other IFU cubes
        targets = [
            (i, c) for i, c in enumerate(self._cubes)
            if c is not var_cube and not isinstance(c, FITSImage)
        ]
        if not targets:
            return

        menu = QMenu(self)
        if len(targets) == 1:
            _, target = targets[0]
            action = QAction(f"Use as variance for '{target.name}'", self)
            action.triggered.connect(
                lambda checked, v=var_cube, t=target:
                    self.use_as_variance.emit(v, t))
            menu.addAction(action)
        else:
            sub = menu.addMenu("Use as variance for…")
            for _, target in targets:
                action = QAction(target.name, self)
                action.triggered.connect(
                    lambda checked, v=var_cube, t=target:
                        self.use_as_variance.emit(v, t))
                sub.addAction(action)

        menu.exec_(self._list.mapToGlobal(pos))


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _item_text(cube):
    """One-line label for the list widget."""
    is_image = isinstance(cube, FITSImage)
    if is_image:
        kind = "image"
        wave_str = ""
    else:
        kind = "cube"
        w0, w1 = cube.wave_range
        if w0 is not None:
            wave_str = f"  {w0:.0f}–{w1:.0f} Å"
        else:
            wave_str = ""

    return f"{cube.name}  [{kind}]{wave_str}"
