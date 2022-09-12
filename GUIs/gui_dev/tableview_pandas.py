import sys

import pandas as pd

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, pyqtSignal


class CustomZTable(QtWidgets.QWidget):
	send_dictdata = pyqtSignal(dict)


	def __init__(self):
		super().__init__()
		self.table = QtWidgets.QTableView()
		self.estZ = pd.DataFrame(#[[0,0,0,0,0,0,0,0]],
			columns=['Name', 'RA', 'DEC', 'z', 'z_err', 'Confidence', 'Linelist', 'Flag', 'z_guess'])
		self.filename = ''

		self.model = TableModel(self.estZ)
		self.table.setModel(self.model)
		self.table.setSortingEnabled(True)

		b_load = QtWidgets.QPushButton('Load')
		b_load.clicked.connect(self._load_button_clicked)
		b_save = QtWidgets.QPushButton('Save')
		b_save.clicked.connect(self._save_button_clicked)
		#b_clear = QtWidgets.QPushButton('Clear')
		#b_clear.clicked.connect(self._clear_button_clicked)


		layout_b = QtWidgets.QVBoxLayout()
		layout_b.addWidget(b_load)
		layout_b.addWidget(b_save)
		#layout_b.addWidget(b_clear)


		layout = QtWidgets.QHBoxLayout()
		layout.addWidget(self.table)
		layout.addLayout(layout_b)
		layout.setAlignment(Qt.AlignCenter)
		self.setLayout(layout)
		self.setFixedHeight(200)

	def _on_sent_estZ(self, sent_z_est):
		if self.estZ.shape[0] == 0:
			self.estZ = self.estZ.append(sent_z_est, ignore_index=True)

		self._update_table()
		#print(self.estZ)

	def _on_sent_data(self, sent_data):
		#print(self.estZ.iloc[0].to_dict())
		#print(sent_data['Name'] in self.estZ['Name'])
		if sent_data['Name'] in self.estZ['Name'].values:
			ind = self.estZ[self.estZ['Name'] == sent_data['Name']].index.values[0]
			s = self.estZ.iloc[ind].to_dict()
			#print(s)
			s.update(sent_data)                                                    
			self.estZ.iloc[ind] = s
		else:
			self.estZ = self.estZ.append(sent_data, ignore_index=True)
		self._update_table()

	def _move_current_filename_top(self, sent_filename):
		checklist = self.estZ[self.estZ['Name'] == sent_filename].index.tolist()
		if checklist:
			target_idx = checklist[0]
			new_index = self.estZ.index.tolist()
			# swap index
			new_index[0], new_index[target_idx] = target_idx, 0
			self.estZ = self.estZ.reindex(new_index)
			self.estZ.reset_index(inplace=True, drop=True)
			self._update_table()
			self.send_dictdata.emit(self.estZ.iloc[0].to_dict())
			print(self.estZ.iloc[0].to_dict())
		else:
			self.send_dictdata.emit({})


	#def _clear_button_clicked(self):
	#	self.estZ = pd.DataFrame(#[[0,0,0,0,0,0,0,0]],
	#		columns=['Name', 'RA', 'DEC', 'z', 'z_err', 'Confidence', 'Linelist', 'Flag'])
	#	self._update_table()

	def _load_button_clicked(self):
		#Load estimated redshift working file
		filepath, check = QtWidgets.QFileDialog.getOpenFileName(None,
			'Load estimated redshifts from csv',
			'',
			'CSV Files (*.csv)')
		if check:
			self.estZ = pd.read_csv(filepath, sep=',')
		self._update_table()

	def _save_button_clicked(self):
		filepath, check = QtWidgets.QFileDialog.getSaveFileName(None,
			'Save estimated redshifts to csv',
			'',
			'CSV Files (*.csv)')
		if check:
			self.estZ.to_csv(filepath,
							 sep=',',
							 header=True,
							 index=False)

	def _update_table(self):
		self.model = TableModel(self.estZ)
		self.table.setModel(self.model)


class TableModel(QtCore.QAbstractTableModel):
	def __init__(self, data):
		super().__init__()
		self._data = data

	def data(self, index, role):
		if index.isValid():
			if role == Qt.DisplayRole or role == Qt.EditRole:
				value = self._data.iloc[index.row(), index.column()]
				return str(value)

	def setData(self, index, value, role):
		if role == Qt.EditRole:
			self._data.iloc[index.row(), index.column()] = value
			return True
		return False

	def rowCount(self, index):
		return self._data.shape[0]

	def columnCount(self, index):
		return self._data.shape[1]

	def headerData(self, section, orientation, role):
		if role == Qt.DisplayRole:
			if orientation == Qt.Horizontal:
				return str(self._data.columns[section])

			if orientation == Qt.Vertical:
				return str(self._data.index[section])

	def flags(self, index):
		return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable

	def sort(self, Ncol, order):
		self.layoutAboutToBeChanged.emit()
		self._data = self._data.sort_values(self._data.columns[Ncol],
											ascending=not order)
		self.layoutChanged.emit()
