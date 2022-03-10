from PyQt5.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QDialog, QComboBox, QLineEdit
from PyQt5.QtCore import Qt, pyqtSignal

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import pandas as pd

from IGM.rb_setline import read_line_list

class Gaussfit_2d(QDialog):
	send_linelist = pyqtSignal(object)
	send_lineindex = pyqtSignal(int)

	def __init__(self,wave, flux1d,error1d):
		super().__init__()
		layout = QVBoxLayout()
		# sublayout for Guess z functions
		lines_layout = QHBoxLayout()
		
		self.line_combo = QComboBox()
		self.line_combo.setFixedWidth(120)
		self.line_combo.addItems(['NONE', 'LBG', 'Gal', 'LLS', 'LLS Small', 'DLA', 'atom'])
		self.line_combo.currentTextChanged.connect(self._linelist_changed)
		
		self.ion_combo = QComboBox()
		self.ion_combo.setFixedWidth(150)
		self.ion_combo.addItem('NONE')
		self.ion_combo.setCurrentIndex(0)
		self.ion_combo.currentIndexChanged.connect(self._index_changed)

		z_guess = QLineEdit()
		z_guess.setPlaceholderText('Guess z')

		lines_layout.addWidget(self.line_combo)
		lines_layout.addWidget(self.ion_combo)
		lines_layout.addWidget(z_guess)

		line1d = LineCanvas()
		line1d._plot_line(wave,flux1d, error1d)
		mpl_toolbar = NavigationToolbar(line1d, self)

		# main layout
		layout.addLayout(lines_layout)
		layout.addWidget(mpl_toolbar)
		layout.addWidget(line1d)
		self.setLayout(layout)
		self.setFixedSize(600,600)

	def _index_changed(self, i): # i is an int
		self.send_lineindex.emit(i)

		#receiver should have

	def _linelist_changed(self, s):
		if s in 'NONE':
			self.send_linelist.emit(s)
			self.ion_combo.clear()
			self.ion_combo.addItem('NONE')
			self.ion_combo.setCurrentIndex(0)
		else:
			llist = self._get_linelist_df(s)
			self.linelist = llist
			self.ion_combo.clear()
			self.ion_combo.addItems(['ALL'] + self.linelist['name'].tolist())
			self.send_linelist.emit(self.linelist)
			self.ion_combo.setCurrentIndex(0)

	def _get_linelist_df(self, linelist_name):
		llist = pd.DataFrame(columns=['wave', 'name'])
		tmp = read_line_list(linelist_name)

		#need a line to append wrest to name if it doesn't have one
		if any(map(str.isdigit, tmp[1]['ion'])):
			# if name column has wrest
			for li in tmp:
				newrow = {'wave': li['wrest'], 'name': li['ion']}
				llist = llist.append(newrow, ignore_index=True)
		else:
			# if name column doesn't have wrest, need to append
			for li in tmp:
				newrow = {'wave': li['wrest'], 'name': li['ion']+' '+str(round(li['wrest']))}
				llist = llist.append(newrow, ignore_index=True)

		return llist

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