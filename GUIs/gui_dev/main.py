import sys
import pandas as pd

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QStatusBar, QVBoxLayout, QGridLayout, QLabel, QComboBox

from PyQt5 import QtCore
from PyQt5 import QtGui
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from menu_toolbars import Custom_ToolBar, Custom_MenuBar
from spec_plot import MplCanvas
from spec_plot_pyqtgraph import SpecCanvas
from linelist_selection import lineListWidget
from utils import FitsObj

use_pyqtgraph = False

class MainWindow(QMainWindow):
	'''This is the main window of this GUI
	This class only assembles different widgets from customed components.
	'''
	def __init__(self):
		super().__init__()

		#----------- External data ----------------------------
		# save a fits copy in the main window
		self.fitsobj = FitsObj([],[],[])
		self.linelist = []
		self.newlinelist = []
		self.z_est = pd.DataFrame(columns=['Name', 'RA', 'DEC', 'Z', 'Z_err', 'Confidence', 'Flag'])

		# --------- Menus and Toolbars -----------------------
		self.setStatusBar(QStatusBar(self))
		toolbar = Custom_ToolBar(self)
		self.addToolBar(toolbar)

		menubar = Custom_MenuBar(self)
		self.setMenuBar(menubar)		

		# --------- Assembling All Widgets ------------------- 
		#placeholder to hold the central widget in main window
		widget = QWidget()
		widget.setMinimumSize(1000, 600)

		widget_z = lineListWidget(menubar)


		l_out = QLabel('GUI Output >> Gaussian profile, EW, etc')
		l_out.setAlignment(QtCore.Qt.AlignCenter)

		layout = QVBoxLayout()
		layout.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
		if use_pyqtgraph:
			self.sc = SpecCanvas()
		else:
			self.sc = MplCanvas(self, width=10, height=8, dpi=100)
			mpl_toolbar = NavigationToolbar(self.sc, self)
			layout.addWidget(mpl_toolbar)
		layout.addWidget(self.sc)
		layout.addWidget(widget_z)
		layout.addWidget(l_out)
		widget.setLayout(layout)

		# MainWindow parameters
		self.setCentralWidget(widget)
		self.setWindowTitle('Spectrum Analysis GUI')
		


		# ----------------- Menubar signal connect ---------------		
		menubar.send_fitsobj.connect(self.on_fitsobj_slot)
		menubar.send_linelist.connect(self.on_linelist_slot)
		menubar.send_z_est.connect(self.on_z_est_slot)
		print(self.fitsobj.wave)

	def on_fitsobj_slot(self, sent_fitsobj):
		self.fitsobj = sent_fitsobj
		self.update()
		print(self.fitsobj.flux)
	def on_linelist_slot(self, sent_linelist):
		self.linelist = sent_linelist
		self.update()
		print(self.linelist)
	def on_newlinelist_slot(self, sent_newlinelist):
		self.newlinelist = sent_newlinelist
		self.update()
	def on_z_est_slot(self, sent_z_est):
		self.z_est = sent_z_est
		self.update()
		print(self.z_est)
		



app = QApplication(sys.argv)
#app.setQuitOnLastWindowClosed(True)
window = MainWindow()
window.show()
app.exec_()
app.quit()