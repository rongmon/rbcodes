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
			columns=['Name', 'RA', 'DEC', 'z', 'z_err', 'Confidence', 'Linelist', 'Flag', 'z_guess', 'z_source', 'primary_frame'])
		self.filename = ''
		self.frame_columns = {}  # Track frame-specific columns: {'EMLINEA': ['z_EMLINEA', 'z_err_EMLINEA'], ...}
		self.current_context_row = -1  # Track which row context menu was opened on

		self.model = TableModel(self.estZ)
		self.table.setModel(self.model)
		self.table.setSortingEnabled(True)

		# Enable context menu
		self.table.setContextMenuPolicy(Qt.CustomContextMenu)
		self.table.customContextMenuRequested.connect(self._on_context_menu)

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
		#self.setFixedHeight(200)

	def _ensure_frame_columns(self, frame_name):
		"""Ensure frame-specific columns exist for a given frame.
		Creates z_FRAMENAME and z_err_FRAMENAME columns if they don't exist."""
		z_col = f'z_{frame_name}'
		z_err_col = f'z_err_{frame_name}'

		if z_col not in self.estZ.columns:
			self.estZ[z_col] = None
		if z_err_col not in self.estZ.columns:
			self.estZ[z_err_col] = None

		if frame_name not in self.frame_columns:
			self.frame_columns[frame_name] = [z_col, z_err_col]

	def _on_sent_estZ(self, sent_z_est):
		if self.estZ.shape[0] == 0:
			self.estZ = self.estZ.append(sent_z_est, ignore_index=True)

		self._update_table()
		#print(self.estZ)

	def _on_sent_data(self, sent_data):
		#print(self.estZ.iloc[0].to_dict())
		#print(sent_data['Name'] in self.estZ['Name'])

		# Extract frame_source if present, default to 'DEFAULT'
		frame_source = sent_data.pop('frame_source', 'DEFAULT')

		# Ensure frame-specific columns exist
		self._ensure_frame_columns(frame_source)

		if sent_data['Name'] in self.estZ['Name'].values:
			ind = self.estZ[self.estZ['Name'] == sent_data['Name']].index.values[0]
			s = self.estZ.iloc[ind].to_dict()
			#print(s)
			if float(s['z_guess']) > 0:
				#print(type(sent_data))
				sent_data.pop('z_guess', None)
			s.update(sent_data)

			# Update frame-specific columns
			z_col = f'z_{frame_source}'
			z_err_col = f'z_err_{frame_source}'
			s[z_col] = s['z']
			s[z_err_col] = s['z_err']

			# Set primary_frame to most recent measurement
			s['primary_frame'] = frame_source

			# PRINT OUT measured redshifts and errors in the database
			print(f"Redshift= {s['z']}; Error= {s['z_err']}; Frame= {frame_source}")
			self.estZ.iloc[ind] = s
		else:
			# Ensure all required columns are in sent_data for new row
			if 'z_source' not in sent_data:
				sent_data['z_source'] = 'Manual'
			if 'primary_frame' not in sent_data:
				sent_data['primary_frame'] = frame_source

			# Add frame-specific columns to new row
			z_col = f'z_{frame_source}'
			z_err_col = f'z_err_{frame_source}'
			sent_data[z_col] = sent_data['z']
			sent_data[z_err_col] = sent_data['z_err']

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

	def _on_context_menu(self, position):
		"""Handle right-click context menu on table."""
		index = self.table.indexAt(position)

		# Only show menu if clicked on a valid row
		if not index.isValid():
			return

		row = index.row()
		self.current_context_row = row

		# Get the row data
		row_data = self.estZ.iloc[row]

		# Get available frame columns for this row (columns that have non-NaN z values)
		available_frames = []
		for frame_name in sorted(self.frame_columns.keys()):
			z_col = f'z_{frame_name}'
			# Check if this frame has a measurement (non-NaN)
			if z_col in self.estZ.columns:
				z_val = row_data[z_col]
				# Check if z_val is not NaN and is a number
				try:
					if pd.notna(z_val) and float(z_val) != 0:
						available_frames.append(frame_name)
				except (ValueError, TypeError):
					pass

		# Only show menu if there are multiple frame measurements
		if len(available_frames) < 2:
			return

		# Create context menu
		menu = QtWidgets.QMenu(self.table)
		menu.setTitle("Set Primary Frame")

		current_primary = row_data.get('primary_frame', '')

		for frame_name in available_frames:
			action = menu.addAction(frame_name)
			# Mark current primary with checkmark indicator
			if frame_name == current_primary:
				action.setText(f"✓ {frame_name}")
			action.triggered.connect(lambda checked, f=frame_name: self._set_primary_frame(row, f))

		# Show menu at cursor position
		menu.exec_(self.table.mapToGlobal(position))

	def _set_primary_frame(self, row, frame_name):
		"""Set a frame as the primary frame for a row.

		Updates:
		1. z and z_err columns with values from selected frame
		2. primary_frame column to track which frame is primary
		"""
		if row < 0 or row >= len(self.estZ):
			return

		z_col = f'z_{frame_name}'
		z_err_col = f'z_err_{frame_name}'

		# Check if columns exist
		if z_col not in self.estZ.columns or z_err_col not in self.estZ.columns:
			print(f"Warning: Frame columns for {frame_name} not found")
			return

		# Update main z and z_err with values from selected frame
		z_value = self.estZ.iloc[row][z_col]
		z_err_value = self.estZ.iloc[row][z_err_col]

		self.estZ.iloc[row, self.estZ.columns.get_loc('z')] = z_value
		self.estZ.iloc[row, self.estZ.columns.get_loc('z_err')] = z_err_value
		self.estZ.iloc[row, self.estZ.columns.get_loc('primary_frame')] = frame_name

		print(f"Row {row}: Set primary frame to {frame_name}, z={z_value}, z_err={z_err_value}")

		# Refresh table display
		self._update_table()


class TableModel(QtCore.QAbstractTableModel):
	def __init__(self, data):
		super().__init__()
		self._data = data

	def data(self, index, role):
		if index.isValid():
			if role == Qt.DisplayRole or role == Qt.EditRole:
				value = self._data.iloc[index.row(), index.column()]

				# Add visual indicator for primary_frame column
				col_name = self._data.columns[index.column()]
				if col_name == 'primary_frame' and role == Qt.DisplayRole:
					# Add star indicator before frame name
					if pd.notna(value) and str(value).strip() and str(value).strip() != 'nan':
						return f"★ {str(value)}"

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
