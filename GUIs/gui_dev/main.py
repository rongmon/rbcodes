import sys, argparse, os
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
		# if using XSpectrum1D io
		self.xspecio = xspecio

		#----------- External data ---------------------------------
		# save a fits copy in the main window
		self.fitsobj = FitsObj([],None,None)
		self.linelist = []
		self.newlinelist = []
		self.z_est = [] #pd.DataFrame(columns=['Name', 'RA', 'DEC', 'Z', 'Z_err', 'Confidence', 'Flag'])

		# --------- Menus and Toolbars -----------------------------
		self.setStatusBar(QStatusBar(self))
		self.toolbar = Custom_ToolBar(self)
		self.addToolBar(self.toolbar)

		#menubar = Custom_MenuBar(self)
		#self.setMenuBar(menubar)		

		# --------- Define all necessary widgets ------------------- 
		#placeholder to hold the central widget in main window
		widget = QWidget()
		#widget.setMinimumSize(1000, 800)
		#widget.setFixedSize(1600, 900)

		# Primary Linelist/Redshift Estimation
		widget_z = LineListWidget()
		#widget_z.setFixedSize(1000,80)

		# GUI Database/Table DataFrame
		table_z = CustomZTable()
		# Message box
		mbox = MessageBox()

		sublayout = QHBoxLayout()
		sublayout.addWidget(mbox)
		sublayout.addWidget(table_z)
		sublayout.setContentsMargins(0,0,0,0)

		layout = QVBoxLayout()
		layout.setAlignment(QtCore.Qt.AlignLeft) #QtCore.Qt.AlignTop | 

		# Main Plotting Canvas
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
		self.toolbar.send_filename.connect(widget_z._on_sent_filename)
		self.toolbar.send_fitsobj.connect(widget_z._on_sent_fitsobj)
		self.toolbar.send_filename.connect(table_z._move_current_filename_top)
		self.toolbar.send_filename.connect(self.sc._update_lines_for_newfile)
		self.toolbar.send_fitsobj.connect(self.on_fitsobj_slot)
		self.toolbar.send_message.connect(lambda s,c='#ff0000': mbox.on_sent_message(s, c))

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
		self.sc.send_scale_limits.connect(self.toolbar._on_scale_limits_slot)
		self.sc.send_extract1d.connect(self.toolbar._on_sent_extract1d)
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

def read_file_from_commandline(filepath, win=None):
	filename = window.toolbar._get_filename(filepath, extension=False)
	win.toolbar.filepaths.append(filepath)
	win.toolbar.filenames.append(filename)		
	win.toolbar.f_combobox.addItem(filename)
	win.toolbar.f_combobox.setCurrentIndex(1)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='GUI default IO Initialization',
									add_help=True)
	# add mutual exclusive arguments to allow users to intialize one IO at a time
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-d', '--default', action='store_true', required=False, default=True,
		help='Read fits files using default gui_io IO class')
	group.add_argument('-x', '--xspec', action='store_true', required=False, default=False,
		help='Read fits files using XSpectrum1D class from linetools')
	parser.add_argument('-f', '--fitsfile', type=str, action='store', required=False, default='',
		help='Feed one FITS file to the internal GUI database')
	args = parser.parse_args()

	app = QApplication(sys.argv[:1])

	# Select IO class
	if args.xspec:
		print('Enable XSpectrum1D IO')
		window = MainWindow(xspecio=True)	
	else:
		print('Enable Default IO')
		window = MainWindow(xspecio=False)

	# Check if filename provided
	if len(args.fitsfile) > 0:
		if args.fitsfile.endswith('fits'):
			filepath = os.path.abspath(args.fitsfile)
			file_exists = os.path.exists(filepath)
			if file_exists:
				print(f'Reading: {filepath}')
				read_file_from_commandline(filepath, win=window)

			else:
				print('Your file does not exists.')
				exit()
		else:
			print('GUI can only read FITS files.')
			exit()

	#app.setQuitOnLastWindowClosed(True)
	window._location_on_screen()
	window.show()

	app.setStyleSheet(qss)
	app.exec_()
	app.quit()	