from PyQt5.QtWidgets import QWidget, QVBoxLayout, QDialog
from PyQt5.QtCore import Qt

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np

class Gaussfit_2d(QDialog):
	def __init__(self,wave, flux1d,error1d):
		super().__init__()
		layout = QVBoxLayout()
		line1d = LineCanvas()
		line1d._plot_hist(flux2d)
		mpl_toolbar = NavigationToolbar(line1d, self)
		layout.addWidget(mpl_toolbar)
		layout.addWidget(line1d)
		self.setLayout(layout)
		self.setFixedSize(600,600)

class LineCanvas(FigureCanvasQTAgg):
	def __init__(self, parent=None, width=5, height=3, dpi=100):
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
		super().__init__(self.fig)
		self.fig.canvas.setFocusPolicy(Qt.ClickFocus)
		self.fig.canvas.setFocus()

	def _plot_line(self, wave,flux1d,error1d):
		self.axline = self.fig.add_subplot(111)
		self.axline.cla()
		self.axline.plot(wave,flux1d,'k')
		self.axline.plot(wave,error1d,'r')


		self.axline.set_xlabel('Wavelength')
		self.axline.set_title('Fit Gaussians')
		
		self.draw()