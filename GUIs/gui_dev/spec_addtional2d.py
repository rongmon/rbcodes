from PyQt5.QtWidgets import QWidget, QVBoxLayout
from PyQt5.QtCore import Qt

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np

class ShowStamp(QWidget):
	def __init__(self, stamp):
		super().__init__()
		layout = QVBoxLayout()
		stampobj = Stamp2dCanvas()
		stampobj._imshow(stamp)
		mpl_toolbar = NavigationToolbar(stampobj, self)
		layout.addWidget(mpl_toolbar)
		layout.addWidget(stampobj)
		self.setLayout(layout)
		self.setFixedSize(600,600)

class Stamp2dCanvas(FigureCanvasQTAgg):
	def __init__(self, parent=None, width=6, height=6, dpi=100):
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
		super().__init__(self.fig)
		self.fig.canvas.setFocusPolicy(Qt.ClickFocus)
		self.fig.canvas.setFocus()

	def _imshow(self, stamp):
		ax = self.fig.add_subplot(111)
		ax.cla()
		pos_ax = ax.imshow(stamp, origin='lower',
						vmin=stamp.min(), vmax=stamp.max())
		ax_cb = self.fig.colorbar(pos_ax, ax=ax, location='right')
		ax.set_title('Stamp of Current Object')		
		self.draw()