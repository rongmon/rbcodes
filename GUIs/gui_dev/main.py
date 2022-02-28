import sys
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
	def __init__(self):
		super().__init__()

		#----------- External data ---------------------------------
		# save a fits copy in the main window
		self.fitsobj = FitsObj([],[],[])
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

		layout = QVBoxLayout()
		layout.setAlignment(QtCore.Qt.AlignLeft) #QtCore.Qt.AlignTop | 

		self.sc = MplCanvas(width=15, height=9, dpi=100)
		self.sc.setMinimumSize(1000,500)
		sc_layout = QVBoxLayout()
		mpl_toolbar = NavigationToolbar(self.sc, self)
		sc_layout.addWidget(mpl_toolbar)
		sc_layout.addWidget(self.sc)
		
		layout.addLayout(sc_layout)
		layout.addWidget(widget_z)
		#layout.addWidget(table_z)
		#layout.addWidget(mbox)
		layout.addLayout(sublayout)
		widget.setLayout(layout)

		# MainWindow parameters
		self.setCentralWidget(widget)
		self.setWindowTitle('Spectrum Analysis GUI')
		


		# ----------------- Signal connections ---------------
		# 1. menubr ==> mainWindow		
		#menubar.send_fitsobj.connect(self.on_fitsobj_slot)
		#menubar.send_linelist.connect(self.on_linelist_slot)
		#menubar.send_z_est.connect(self.on_z_est_slot)
		#print(self.fitsobj.wave)
		# 2. menubar ==> widget_z
		#menubar.send_filename.connect(widget_z.on_linelist_name_slot)
		#menubar.send_linelist.connect(widget_z.on_linelist_slot)
		# 3. menubar ==> sc (SpecCanvas)
		#menubar.send_linelist.connect(self.sc.on_linelist_slot)
		# 4. menubar ==> table_z
		#menubar.send_z_est.connect(table_z._on_sent_estZ)
		# 5. widget_z ==> sc (SpecCanvas)
		widget_z.send_linelist.connect(self.sc.on_linelist_slot)
		widget_z.send_lineindex.connect(self.sc.on_lineindex_slot)
		widget_z.send_gauss_num.connect(self.sc._on_sent_gauss_num)
		widget_z.estZ.returnPressed.connect(lambda z=widget_z.estZ: self.passing_estZ(z))
		# 6. sc (SpecCanvas) ==> mbox (MessageBox)
		self.sc.send_message.connect(mbox.on_sent_message)
		# 7. sc (SpecCanvas) ==> widget_z.estZ
		self.sc.send_z_est.connect(widget_z._on_estZ_changed)
		# 8. sc (SpecCanvas) ==> toolbar
		self.sc.send_scale_limits.connect(toolbar._on_scale_limits_slot)
		# 9. toolbar ==> widget_z
		toolbar.send_filename.connect(widget_z._on_sent_filename)
		# 10. widget_z ==> table_z
		widget_z.send_data.connect(table_z._on_sent_data)
		# 11. toolbar ==> table_z
		toolbar.send_filename.connect(table_z._move_current_filename_top)
		# 12. toolbar ==> sc (SpecCanvas)
		toolbar.send_filename.connect(self.sc._update_lines_for_newfile)
		# 13. table_z ==> widget_z
		table_z.send_dictdata.connect(widget_z._on_sent_dictdata)
		# 14. widget_z ==> mbox (MessageBox)
		widget_z.send_message.connect(mbox.on_sent_message)



	def _location_on_screen(self):
		screen = QDesktopWidget().screenGeometry()
		avail = QDesktopWidget().availableGeometry()
		widget = self.geometry()

		x = avail.width() - widget.width()
		y = 2 * avail.height() - screen.height() - widget.height()
		self.move(50,50)


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


app = QApplication(sys.argv)
#app.setQuitOnLastWindowClosed(True)
window = MainWindow()
window._location_on_screen()
window.show()

app.setStyleSheet(qss)
app.exec_()
app.quit()