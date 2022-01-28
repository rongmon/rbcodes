from PyQt5 import QtWidgets

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

# not necessary if not displaying images
import matplotlib.image as mpimg

matplotlib.use('Qt5Agg')

class MplCanvas(FigureCanvasQTAgg):
	def __init__(self, parent=None, width=5, height=4, dpi=100):
		fig = Figure(figsize=(width, height), dpi=dpi)
		self.axes = fig.add_subplot(111)
		'''
		img = mpimg.imread('cat.jpg')
		self.axes.imshow(img)
		'''

		#self.wave = None
		#self.flux = None
		#self.error = None
		self.wave=[0,1,2,3,4]
		self.flux = [10, 1, 20, 3, 40]
		self.error = [0.5, 0.1, 0.4, 0.06, 0.6]
		super().__init__(fig)

	def plot(self, wave, flux, error):
		self.axes.cla()
		self.axes.plot(wave, flux, label='Flux')
		self.axes.plot(wave, error, label='Error')
		self.axes.legend(loc='upper right')
		self.axes.set_ylim([-1e-13, 0.01])
		self.axes.set_xlabel('Wavelength')
		self.axes.set_ylabel('Flux')
		self.draw()

