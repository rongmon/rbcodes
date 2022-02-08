import sys

import pandas as pd

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt


class CustomZTable(QtWidgets.QWidget):
	def __init__(self):
		super().__init__()
		self.table = QtWidgets.QTableView()
		self.estZ = pd.DataFrame(#[[0,0,0,0,0,0,0]],
			columns=['Name', 'RA', 'DEC', 'Z', 'Z_err', 'Confidence', 'Flag'])

		self.model = TableModel(self.estZ)
		self.table.setModel(self.model)
		layout = QtWidgets.QHBoxLayout()
		layout.addWidget(self.table)
		layout.setAlignment(Qt.AlignCenter)
		self.setLayout(layout)

	def on_sent_estZ(self, sent_z_est):
		self.estZ = sent_z_est
		#self.model = TableModel(self.estZ)
		#self.table.setModel(self.model)
		print(self.estZ)


class TableModel(QtCore.QAbstractTableModel):
	def __init__(self, data):
		super().__init__()
		self._data = data

	def data(self, index, role):
		if role == Qt.DisplayRole:
			value = self._data.iloc[index.row(), index.column()]
			return str(value)

	def setData(self, index, value, role):
		row = self._data.index[index.row()]
		col = self._data.columns[index.column()]
		if hasattr(value, 'toPyObject'):
			value = value.toPyObject()
		self._data.set_value(row, col, value)
		return True

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


