import sys

import numpy as np

from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
from pyqtgraph import PlotWidget, plot


class SpecCanvas(pg.PlotWidget):
	def __init__(self, parent=None):
		super(SpecCanvas, self).__init__(parent, viewBox=CustomViewBox())
		self.wave, self.flux, self.error = [], [], []
		self.clear()
		self.setBackground('w')
		styles = {'color': 'k', 'font-size': '15pt'}
		self.setLabel('left', 'Flux', **styles)
		self.setLabel('bottom', 'Wavelength', **styles)
		self.showGrid(x=True, y=True)
		self.addLegend()

		self.pen_flux = pg.mkPen(color='k')
		self.pen_err = pg.mkPen(color='r')

	def plot_spec(self, wave, flux, error):
		self.plot(wave, flux, name='Flux', pen=self.pen_flux)
		self.plot(wave, error, name='Error', pen=self.pen_err)
		self.setXRange(np.min(wave), np.max(wave), padding=0)
		self.setYRange(np.median(flux)*0.05, np.median(flux)*3, padding=0)

		
class CustomViewBox(pg.ViewBox):
	def __init__(self):
		super().__init__()
		self.setMouseMode(self.RectMode)

		self.firstLeftClick = True
		self.xview, self.yview = [], []

	def mouseClickEvent(self, event):
		if self.firstLeftClick & (event.button() == QtCore.Qt.LeftButton):
			xlim, ylim = self.viewRange()
			self.xview.append(xlim)
			self.yview.append(ylim)
			self.firstLeftClick = False
		if event.button() == QtCore.Qt.RightButton:
			#if len(self.xview) > 1:
			#	self.setRange(xRange=self.xview.pop(), yRange=self.yview.pop())
			#	print('view left: ', len(self.xview))
			#else:
			#	self.setRange(xRange=self.xview[0], yRange=self.yview[0])
			self.setRange(xRange=self.xview[0], yRange=self.yview[0])			


	def mouseDragEvent(self, event):
		if event.button() == QtCore.Qt.RightButton:
			event.ignore()
		else:
			pg.ViewBox.mouseDragEvent(self, event)
			#xzoom, yzoom = self.lastPos()
			#self.xview.append(xzoom)
			#self.yview.append(yzoom)
			#print(len(self.xview))
		#print(event.lastPos(), event.pos())

	def keyPressEvent(self, event):
		print(event.text())

		if event.text() == '-':
			self.scaleHistory(-1)
		elif event.text() in ['+', '=']:
			self.scaleHistory(1)
		elif event.key() == QtCore.Qt.Key.Key_Backspace:
			self.scaleHistory(len(self.axHistory))

	#def getDataItem(self):
	#	self.ObjItemList = pg.GraphicsScene.items(self.scene(), self.ax)
	#	self.dataxy = self.ObjItemList[0].listDataItems()[0].getData()
