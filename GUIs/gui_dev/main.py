import sys, argparse
import pandas as pd

from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QStatusBar, QDesktopWidget,
							QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox)

from PyQt5 import QtCore
from PyQt5 import QtGui
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from menu_toolbars import Custom_ToolBar, Custom_MenuBar
from spec_plot import MplCanvas
from linelist_selection import LineListWidget
from tableview_pandas import CustomZTable
from message_box import MessageBox
from utils import FitsObj

class MainWindow(QMainWindow):
	'''This is the main window of this GUI
	This class only assembles different widgets from customed components.
	'''
	def __init__(self, xspecio=False):
		super().__init__()
		self.xspecio = xspecio

		#----------- External data ---------------------------------
		# save a fits copy in the main window
		self.fitsobj = FitsObj([],None,None)
		self.linelist = []
		self.newlinelist = []
		self.z_est = [] #pd.DataFrame(columns=['Name', 'RA', 'DEC', 'Z', 'Z_err', 'Confidence', 'Flag'])

		# --------- Menus and Toolbars -----------------------------
		self.setStatusBar(QStatusBar(self))
		toolbar = Custom_ToolBar(self)
		self.addToolBar(toolbar)

		#menubar = Custom_MenuBar(self)
		#self.setMenuBar(menubar)		

		# --------- Define all necessary widgets ------------------- 
		#placeholder to hold the central widget in main window
		widget = QWidget()
		#widget.setMinimumSize(1000, 800)
		#widget.setFixedSize(1600, 900)

		widget_z = LineListWidget()
		#widget_z.setFixedSize(1000,80)
		table_z = CustomZTable()
		mbox = MessageBox()
		sublayout = QHBoxLayout()
		sublayout.addWidget(mbox)
		sublayout.addWidget(table_z)
		sublayout.setContentsMargins(0,0,0,0)

		layout = QVBoxLayout()
		layout.setAlignment(QtCore.Qt.AlignLeft) #QtCore.Qt.AlignTop | 

		self.sc = MplCanvas(width=15, height=9, dpi=100)
		self.sc.setMinimumSize(1000,500)
		sc_layout = QVBoxLayout()
		sc_layout.setContentsMargins(0,0,0,0)
		mpl_toolbar = NavigationToolbar(self.sc, self)
		sc_layout.addWidget(mpl_toolbar)
		sc_layout.addWidget(self.sc)
		
		layout.addLayout(sc_layout)
		layout.addWidget(widget_z)
		#layout.addWidget(table_z)
		#layout.addWidget(mbox)
		layout.addLayout(sublayout)
		layout.setContentsMargins(0,0,0,0)
		widget.setLayout(layout)

		# MainWindow parameters
		self.setCentralWidget(widget)
		self.setWindowTitle('Spectrum Analysis GUI')
		#self.setMinimumSize(1000, 900)
		


		# ----------------- Signal connections ---------------
		# 1. menubar signal exports
		#menubar.send_fitsobj.connect(self.on_fitsobj_slot)
		#menubar.send_linelist.connect(self.on_linelist_slot)
		#menubar.send_z_est.connect(self.on_z_est_slot)
		#menubar.send_filename.connect(widget_z.on_linelist_name_slot)
		#menubar.send_linelist.connect(widget_z.on_linelist_slot)
		#menubar.send_linelist.connect(self.sc.on_linelist_slot)
		#menubar.send_z_est.connect(table_z._on_sent_estZ)

		# 2. toolbar signal exports
		toolbar.send_filename.connect(widget_z._on_sent_filename)
		toolbar.send_fitsobj.connect(widget_z._on_sent_fitsobj)
		toolbar.send_filename.connect(table_z._move_current_filename_top)
		toolbar.send_filename.connect(self.sc._update_lines_for_newfile)
		toolbar.send_fitsobj.connect(self.on_fitsobj_slot)
		toolbar.send_message.connect(lambda s,c='#ff0000': mbox.on_sent_message(s, c))

		# 3. widget_z signal exports
		widget_z.send_linelist.connect(self.sc.on_linelist_slot)
		widget_z.send_lineindex.connect(self.sc.on_lineindex_slot)
		widget_z.send_gauss_num.connect(self.sc._on_sent_gauss_num)
		widget_z.send_z_returnPressed.connect(self.sc._on_estZ_return_pressed)
		widget_z.send_more_linelist.connect(self.sc.on_additional_linelist_slot)
		widget_z.send_more_linelist_z.connect(self.sc.on_additional_linelist_slot_z)
		widget_z.send_linelists2multiG.connect(self.sc._on_sent_linelists2multiG)
		widget_z.send_message.connect(mbox.on_sent_message)
		widget_z.send_data.connect(table_z._on_sent_data)

		# 4. sc (SpecCanvas) signal exports
		self.sc.send_message.connect(mbox.on_sent_message)
		self.sc.send_z_est.connect(widget_z._on_estZ_changed)
		self.sc.send_z_manual.connect(widget_z._on_estZ_manual)
		self.sc.send_scale_limits.connect(toolbar._on_scale_limits_slot)
		self.sc.send_extract1d.connect(toolbar._on_sent_extract1d)
		#self.sc.gauss2d.send_gcenter
		
		# 5. table_z ==> widget_z
		table_z.send_dictdata.connect(widget_z._on_sent_dictdata)

	
		#sizeObj = QDesktopWidget().screenGeometry(-1)
		#print('Screen size:' + str(sizeObj.height())+ 'x' + str(sizeObj.width()))	



	def _location_on_screen(self):
		screen = QDesktopWidget().screenGeometry()
		avail = QDesktopWidget().availableGeometry()
		widget = self.geometry()

		x = avail.width() - widget.width()
		y = 2 * avail.height() - screen.height() - widget.height()
		self.move(5,5)


	def on_fitsobj_slot(self, sent_fitsobj):
		self.fitsobj = sent_fitsobj

	def on_linelist_slot(self, sent_linelist):
		self.linelist = sent_linelist

	def on_newlinelist_slot(self, sent_newlinelist):
		self.newlinelist = sent_newlinelist

	def on_z_est_slot(self, sent_z_est):
		self.z_est = sent_z_est
		self.update()
		print(self.z_est)

	def passing_estZ(self, estZ):
		#print(estZ.text())
		self.sc.estZ = float(estZ.text())
		self.sc._plot_lines(lineindex=self.sc.lineindex, 
							estZ=self.sc.estZ)
		#self.sc.estZ = float(estZ)
		
qss = '''
	.QLabel {font-size: 8pt}
	.QComboBox {font-size: 8pt}
	.QLineEdit {font-size: 8pt}
	.QPushButton {font-size: 8pt}
	.QAction {font-size: 8pt}
'''


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='GUI default IO Initialization',
									add_help=True)
	parser.add_argument('-d', '--default', action='store_true',
		help='Read fits files using default gui_io IO class')
	parser.add_argument('-x', '--xspec', action='store_true',
		help='Read fits files using XSpectrum1D class from linetools')
	args = parser.parse_args()

	app = QApplication(sys.argv[:1])
	if args.xspec:
		print('Enable XSpectrum1D IO')
		window = MainWindow(xspecio=True)
	else:
		print('Enable Default IO')
		window = MainWindow(xspecio=False)


	#app.setQuitOnLastWindowClosed(True)
	window._location_on_screen()
	window.show()

	app.setStyleSheet(qss)
	app.exec_()
	app.quit()




