import sys, argparse, os
import pandas as pd

from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QStatusBar,
							QVBoxLayout, QHBoxLayout)

from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QKeySequence, QDesktopServices
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from rbcodes.GUIs.zgui.menu_toolbars import Custom_ToolBar, Custom_MenuBar
from rbcodes.GUIs.zgui.spec_plot import MplCanvas
from rbcodes.GUIs.zgui.linelist_selection import LineListWidget
from rbcodes.GUIs.zgui.tableview_pandas import CustomZTable
from rbcodes.GUIs.zgui.message_box import MessageBox
from rbcodes.GUIs.zgui.utils import FitsObj

class MainWindow(QMainWindow):
	'''Simplified main window with proper positioning'''
	def __init__(self, xspecio=False, toggle_frames=False):
		super().__init__()
		
		# Store parameters
		self.xspecio = xspecio
		self.toggle_frames = toggle_frames
		
		# External data
		self.fitsobj = FitsObj([],None,None)
		self.linelist = []
		self.newlinelist = []
		self.z_est = []
		
		# Setup UI
		self.setup_ui()
		self.setup_connections()
		
		# Window setup
		self.setWindowTitle('Spectrum Analysis GUI')
		self.resize(1000, 900)
		
	def setup_ui(self):
		"""Setup the user interface"""
		# Status bar and toolbar
		self.setStatusBar(QStatusBar(self))
		self.toolbar = Custom_ToolBar(self)
		self.addToolBar(self.toolbar)
		
		# Central widget
		central_widget = QWidget()
		self.setCentralWidget(central_widget)
		
		# Main layout
		main_layout = QVBoxLayout(central_widget)
		main_layout.setContentsMargins(0,0,0,0)
		
		# Plotting canvas with toolbar
		self.sc = MplCanvas(width=12, height=7, dpi=100)
		self.sc.setMinimumSize(1000,400)
		
		mpl_toolbar = NavigationToolbar(self.sc, self)
		main_layout.addWidget(mpl_toolbar)
		main_layout.addWidget(self.sc)
		
		# Line selection widget
		widget_z = LineListWidget(canvas=self.sc)
		main_layout.addWidget(widget_z)
		
		# Bottom section with table and message box
		bottom_layout = QHBoxLayout()
		bottom_layout.setContentsMargins(0,0,0,0)
		
		table_z = CustomZTable()
		mbox = MessageBox()
		
		bottom_layout.addWidget(mbox)
		bottom_layout.addWidget(table_z)
		main_layout.addLayout(bottom_layout)
		
		# Store references for signal connections
		self.widget_z = widget_z
		self.table_z = table_z
		self.mbox = mbox
		
	def setup_connections(self):
		"""Setup signal connections between widgets"""
		# Toolbar signals
		self.toolbar.send_filename.connect(self.widget_z._on_sent_filename)
		self.toolbar.send_fitsobj.connect(self.widget_z._on_sent_fitsobj)
		self.toolbar.send_filename.connect(self.table_z._move_current_filename_top)
		self.toolbar.send_filename.connect(self.sc._update_lines_for_newfile)
		self.toolbar.send_fitsobj.connect(self.on_fitsobj_slot)
		self.toolbar.send_fitsobj.connect(self.sc._on_sent_fitsobj)
		self.toolbar.send_frame.connect(self.sc._on_frame_selected)
		self.toolbar.send_message.connect(lambda s,c='#ff0000': self.mbox.on_sent_message(s, c))
		
		# Widget_z signals
		self.widget_z.send_linelist.connect(self.sc.on_linelist_slot)
		self.widget_z.send_lineindex.connect(self.sc.on_lineindex_slot)
		self.widget_z.send_gauss_num.connect(self.sc._on_sent_gauss_num)
		self.widget_z.send_z_returnPressed.connect(self.sc._on_estZ_return_pressed)
		self.widget_z.send_more_linelist.connect(self.sc.on_additional_linelist_slot)
		self.widget_z.send_more_linelist_z.connect(self.sc.on_additional_linelist_slot_z)
		self.widget_z.send_linelists2multiG.connect(self.sc._on_sent_linelists2multiG)
		self.widget_z.send_message.connect(self.mbox.on_sent_message)
		self.widget_z.send_data.connect(self.table_z._on_sent_data)
		
		# Canvas signals
		self.sc.send_message.connect(self.mbox.on_sent_message)
		self.sc.send_z_est.connect(self.widget_z._on_estZ_changed)
		self.sc.send_scale_limits.connect(self.toolbar._on_scale_limits_slot)
		self.sc.send_extract1d.connect(self.toolbar._on_sent_extract1d)
		
		# Table signals
		self.table_z.send_dictdata.connect(self.widget_z._on_sent_dictdata)

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
		self.sc.estZ = float(estZ.text())
		self.sc._plot_lines(lineindex=self.sc.lineindex, 
							estZ=self.sc.estZ)

qss = '''
QMainWindow {
    background-color: #f5f5f5;
}

QLabel { 
    font-size: 9pt; 
    color: #2c3e50;
    font-weight: 500;
}

QComboBox { 
    font-size: 9pt; 
    padding: 4px 8px;
    border: 1px solid #bdc3c7;
    border-radius: 4px;
    background: #ffffff;
    color: #2c3e50;
    min-height: 20px;
}

QComboBox:hover { 
    border-color: #3498db;
    background: #ecf0f1;
}

QComboBox:drop-down {
    border: none;
    padding-right: 8px;
}

QLineEdit { 
    font-size: 9pt;
    padding: 5px 8px;
    border: 1px solid #bdc3c7;
    border-radius: 4px;
    background: #ffffff;
    color: #2c3e50;
    selection-background-color: #3498db;
}

QLineEdit:focus { 
    border-color: #3498db;
    background: #ecf0f1;
}

QPushButton { 
    font-size: 9pt;
    padding: 5px 12px;
    border: 1px solid #bdc3c7;
    border-radius: 4px;
    background: #3498db;
    color: white;
    min-height: 20px;
}

QPushButton:hover { 
    background: #2980b9;
    border-color: #2980b9;
}

QPushButton:pressed {
    background: #2472a4;
}

QCheckBox { 
    font-size: 9pt; 
    color: #2c3e50;
    spacing: 8px;
}

QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border: 1px solid #bdc3c7;
    border-radius: 3px;
}

QCheckBox::indicator:checked {
    background: #3498db;
    border-color: #3498db;
}

QStatusBar {
    background: #ecf0f1;
    color: #2c3e50;
}

QToolBar {
    background: #ecf0f1;
    border-bottom: 1px solid #bdc3c7;
    spacing: 6px;
    padding: 3px;
}

QToolBar QToolButton {
    background: transparent;
    border: 1px solid transparent;
    border-radius: 4px;
    padding: 3px;
}

QToolBar QToolButton:hover {
    background: #d5dbdb;
    border-color: #bdc3c7;
}
'''

def read_file_from_commandline(filepath, win):
    filename = win.toolbar._get_filename(filepath, extension=False)
    win.toolbar.filepaths.append(filepath)
    win.toolbar.filenames.append(filename)        
    win.toolbar.f_combobox.addItem(filename)
    win.toolbar.f_combobox.setCurrentIndex(1)

def launch_gui(xspecio=False, toggle_frames=False, fitsfile=''):
    """Launch the zgui GUI application.
    
    Args:
        xspecio (bool): Whether to use XSpectrum1D for IO
        toggle_frames (bool): Enable frame toggling
        fitsfile (str): Optional FITS file to load on startup
    """
    app = QApplication(sys.argv[:1])
    
    # Report IO mode
    if xspecio:
        print('Enable XSpectrum1D IO') 
    else:
        print('Enable Default IO')
        
    # Report toggling mode
    if toggle_frames:
        print('Frame toggling enabled')
    else:
        print('Frame toggling disabled')

    # Create main window
    window = MainWindow(xspecio=xspecio, toggle_frames=toggle_frames)
    window.setAttribute(Qt.WA_DeleteOnClose)

    # Load file if provided
    if fitsfile:
        if fitsfile.endswith('fits'):
            filepath = os.path.abspath(fitsfile)
            if os.path.exists(filepath):
                print(f'Reading: {filepath}')
                read_file_from_commandline(filepath, win=window)
            else:
                print('File does not exist.')
                return
        else:
            print('GUI can only read FITS files.')
            return

    # Show window and start event loop
    window.show()
    app.setStyleSheet(qss)
    app.exec_()
    app.quit()

def main():
    """Entry point for the zgui command line interface"""
    parser = argparse.ArgumentParser(
        description='ZGUI: A GUI tool for analyzing spectroscopic data and estimating redshifts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
---------
rb_zgui                     # Launch with default IO
rb_zgui -tf                 # Launch with frame toggling enabled
rb_zgui -x                  # Launch with XSpectrum1D IO
rb_zgui -f spectrum.fits    # Launch and load specific FITS file

For more information see: https://github.com/rongmon/rbcodes
''')

    # add mutual exclusive arguments to allow users to intialize one IO at a time
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d', '--default', action='store_true', required=False, default=True,
        help='Read FITS files using default gui_io IO class')
    group.add_argument('-x', '--xspec', action='store_true', required=False, default=False,
        help='Read FITS files using XSpectrum1D class from linetools')
    parser.add_argument('-f', '--fitsfile', type=str, action='store', required=False, default='',
        help='Feed one FITS file to the internal GUI database at startup')
    parser.add_argument('-tf', '--toggleframes', action='store_true', required=False, default=False,
        help='Enable toggling between different frames (SCI, EMLINE, etc.) within FITS files')
    args = parser.parse_args()

    launch_gui(
        xspecio=args.xspec,
        toggle_frames=args.toggleframes,
        fitsfile=args.fitsfile
    )

if __name__ == '__main__':
    main()