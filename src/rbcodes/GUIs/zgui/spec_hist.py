from PyQt5.QtWidgets import QWidget, QVBoxLayout, QDialog
from PyQt5.QtCore import Qt

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np

from menu_toolbars import Custom_ToolBar

class FluxHistogram(QDialog):
	def __init__(self, flux2d):
		super().__init__()
		layout = QVBoxLayout()
		hist2d = Histo2dCanvas()
		hist2d._plot_hist(flux2d)
		mpl_toolbar = NavigationToolbar(hist2d, self)
		layout.addWidget(mpl_toolbar)
		layout.addWidget(hist2d)
		self.setLayout(layout)
		self.setFixedSize(600,600)


class Histo2dCanvas(FigureCanvasQTAgg):
	def __init__(self, parent=None, width=5, height=3, dpi=100):
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
		super().__init__(self.fig)
		self.fig.canvas.setFocusPolicy(Qt.ClickFocus)
		self.fig.canvas.setFocus()

	def _plot_hist(self, flux2d):
		self.axhist = self.fig.add_subplot(111)
		self.axhist.cla()
		#print(flux2d.shape)
		#print(len(flux2d.flatten()))
		x = np.linspace(0, len(flux2d), len(flux2d))
		ysum = np.sum(flux2d, axis=1)
		self.axhist.plot(x,ysum)

		self.axhist.set_xlabel('Y-axis')
		self.axhist.set_title('Data Summation along Wavelength Axis')
		
		self.draw()

class PixelHistogram(QDialog):
	def __init__(self, flux2d):
		super().__init__()
		layout = QVBoxLayout()
		hist2d = PixelHisto2dCanvas()
		hist2d._plot_hist(flux2d)
		mpl_toolbar = NavigationToolbar(hist2d, self)
		layout.addWidget(mpl_toolbar)
		layout.addWidget(hist2d)
		self.setLayout(layout)
		self.setFixedSize(600,600)

class PixelHisto2dCanvas(FigureCanvasQTAgg):
	def __init__(self, parent=None, width=5, height=3, dpi=100):
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
		super().__init__(self.fig)
		self.fig.canvas.setFocusPolicy(Qt.ClickFocus)
		self.fig.canvas.setFocus()

	def _plot_hist(self, flux2d):
		psize = len(flux2d.flatten())
		bins = 200 if psize > 200 else 20
		self.axhist = self.fig.add_subplot(111)
		self.axhist.cla()
		self.axhist.hist(flux2d.flatten(), bins=bins, density=True)

		self.axhist.set_xlabel('Pixel Values')
		self.axhist.set_title('Pixel Distribution')
		
		self.draw()