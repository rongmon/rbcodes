from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, pyqtSignal

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from astropy.convolution import convolve, Box1DKernel
from astropy.modeling import models, fitting
from scipy.interpolate import splrep, splev
import numpy as np


from guess_transition import GuessTransition

matplotlib.use('Qt5Agg')

class MplCanvas(FigureCanvasQTAgg):
	send_message = pyqtSignal(str)
	send_z_est = pyqtSignal(list)
	send_gcenter = pyqtSignal(list)

	def __init__(self, parent=None, width=5, height=3, dpi=100):
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		pad = 0.05
		self.fig.subplots_adjust(left=pad+0.02, bottom=pad*2, right=0.999, top=0.99)
		self.axes = self.fig.add_subplot(111)
		#self.axes.margins(x=0)

		# initialize some plotting parameters
		self.wave, self.flux, self.error = [],[],[]
		self.init_xlims, self.init_ylims = [],[]
		self.gxval, self.gyval = [], []
		self.scale = 1.
		self.lineindex = 0
		self.linelist = [] #pd.DataFrame(columns=['wave', 'name'])
		self.estZ = 0.
		self.estZstd = 0.
		self.guess_ion = 0.
		self.guess_gcenter = []
		self.gauss_num = 1
		self.gauss_profiles = []

		super().__init__(self.fig)

		# connect funcitons to events
		self.fig.canvas.setFocusPolicy(Qt.ClickFocus)
		self.fig.canvas.setFocus()
		self.cid_k = self.fig.canvas.mpl_connect('key_press_event', self.ontype)
		self_cid_m = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
		#self.fig.canvas.mpl_connect('pick_event', self.onpick)

		#self.fig.tight_layout()
		#self.fig.subplots_adjust(0.1, 0.1, 0.9, 0.9)

	def plot_spec(self, wave, flux, error, filename):
		self.axes.cla()
		self.axes.plot(wave, error, color='red')# label='Error')
		self.axes.plot(wave, flux, color='black')#, label='Flux')
		#self.axes.legend(loc='upper right')
		self.axes.set_ylim([-np.min(flux)*0.01, np.median(flux)*3])
		self.axes.set_xlabel('Wavelength')
		self.axes.set_ylabel('Flux')
		self.axes.set_title(filename)
		self.draw()

		self.send_message.emit(f'You are currently working on {filename} spectrum.')

		# update intiial plotting parameters
		self.wave, self.flux, self.error = wave, flux, error
		self.init_xlims = self.axes.get_xlim()
		self.init_ylims = self.axes.get_ylim()

	def ontype(self, event):
		'''Interactivae keyboard events
		Note:
			Always Update to help mannual for new keyboard events
		'''
		#print(event.key)
		if event.key == 'r':
			# reset flux/err and clear all lines
			del self.axes.lines[2:]	# delete everything except flux/err
			self.line_xlim, self.line_ylim = [], []
			self.gxval, self.gyval = [], []
			
			#self.axes.texts = []
			self.axes.set_ylim(self.init_ylims)
			self.axes.set_xlim(self.init_xlims)
			self.draw()
			self.send_message.emit('You reset the canvas!!!')

		elif event.key == 't':
			# set y axis max value
			ylim = self.axes.get_ylim()
			self.axes.set_ylim([ylim[0], event.ydata])
			self.draw()

		elif event.key == 'b':
			# set y axis min value
			ylim = self.axes.get_ylim()
			self.axes.set_ylim([event.ydata, ylim[-1]])
			self.draw()

		elif event.key == 'X':
			# set x axis max value
			xlim = self.axes.get_xlim()
			self.axes.set_xlim([xlim[0], event.xdata])
			self.draw()

		elif event.key == 'x':
			# set x axis min value
			xlim = self.axes.get_xlim()
			self.axes.set_xlim([event.xdata, xlim[-1]])
			self.draw()

		elif event.key == '[':
			# pan window to the left
			xlim = self.axes.get_xlim()
			delx = (xlim[-1] - xlim[0])
			self.axes.set_xlim([xlim[0] - delx, xlim[0]])
			self.draw()

		elif event.key == ']':
			# pan window to the right
			xlim = self.axes.get_xlim()
			delx = (xlim[-1] - xlim[0])
			self.axes.set_xlim([xlim[1], xlim[1] + delx])
			self.draw()

		elif event.key == 'S':
			# smooth ydata
			self.scale += 2
			self.new_spec = convolve(self.flux, Box1DKernel(self.scale))
			self.new_err = convolve(self.error, Box1DKernel(self.scale))
			self.replot(self.new_spec, self.new_err)
			self.draw()

		elif event.key == 'U':
			# unsmooth ydata
			self.scale -= 2
			if self.scale < 0:
				self.scale = 1
			self.new_spec = convolve(self.flux, Box1DKernel(self.scale))
			self.new_err = convolve(self.error, Box1DKernel(self.scale))
			self.replot(self.new_spec, self.new_err)
			self.draw()

		elif event.key == 'Y':
			# set y-axis limits with precise values
			ylimdialog = CustomLimDialog(axis='y')
			if ylimdialog.exec_():
				ylim = ylimdialog._getlim()
				self.axes.set_ylim(ylim)
				self.draw()

		elif event.key == 'W':
			# set y-axis limits with precise values
			xlimdialog = CustomLimDialog(axis='x')
			if xlimdialog.exec_():
				xlim = xlimdialog._getlim()
				self.axes.set_xlim(xlim)
				self.draw()


		elif event.key == 'G':
			# fit a Gaussian profile
			self.gxval = np.append(self.gxval, event.xdata)
			self.gyval = np.append(self.gyval, event.ydata)

			if self.gauss_num > 1:
				dgxval, dgyval = [],[]

			fclick = len(self.gxval)
			self.axes.plot(self.gxval[-1], self.gyval[-1], 'rs', ms=5)
			self.draw()

			if fclick == 1:
				message = 'You need 3 points to model a Gaussian. Please click 2 more points.'
				self.send_message.emit(message)
			elif fclick == 2:
				message = 'Please click 1 more point to model a Gaussian.'
				self.send_message.emit(message)
			elif fclick == 3:
				# sort xdata before fitting
				x_sort = np.argsort(self.gxval)
				gxval = self.gxval[x_sort]
				gyval = self.gyval[x_sort]
				# fit a Gaussian with 3 data points
				
				g_init = models.Gaussian1D(amplitude=gyval[1],
										   mean=gxval[1],
										   stddev=0.5*(gxval[2]-gxval[0]))
				fit_g = fitting.LevMarLSQFitter(calc_uncertainties=True)
				# 1. fit a local continuum
				c_range = np.where((self.wave>=gxval[0]) & (self.wave <= gxval[2]))
				g_wave = self.wave[c_range]
				g_flux = self.flux[c_range]
				g_error = self.error[c_range]
				spline = splrep([gxval[0], gxval[-1]], 
								[gyval[0], gyval[-1]], 
								k=1)
				cont = splev(g_wave, spline)
				errdata = g_error / cont
				# 2. check if it is an absorption or emission line
				if ((gyval[1] < gyval[0]) & (gyval[1] < gyval[2])):
					ydata = 1. - (g_flux / cont)
					g = fit_g(g_init, g_wave, ydata, weights=1.0/errdata)
					g_final = (1. - g(g_wave)) * cont
				else:
					ydata = g_flux / cont
					g = fit_g(g_init, g_wave, ydata, weights=1.0/errdata)
					g_final = g(g_wave) * cont

				model_fit = self.axes.plot(g_wave, g_final, 'r-')
				
				self.draw()
				message = (f'A Gaussian model you fit has the following parameters:\n'
						   f'Amplitude: {g.parameters[0]:.3f}\n'
						   f"Mean: {g.parameters[1]:.3f} with std={g.stds['mean']:.3f}\n"
						   f'Sigma: {g.parameters[2]:.3f}')

				if self.guess_gcenter:
					self.guess_gcenter[0] = g.parameters[1]
					self.guess_gcenter[1] = g.stds['mean']
				else:
					self.guess_gcenter.append(g.parameters[1])
					self.guess_gcenter.append(g.stds['mean'])

				self.gxval, self.gyval = [], []
				self.send_message.emit(message)
				self.send_gcenter.emit(self.guess_gcenter)


				if self.gauss_num == 2:
					self.gauss_profiles = np.append(self.gauss_profiles, g.parameters, axis=0)
					dgxval = np.append(dgxval, gxval, axis=0)
					dgyval = np.append(dgyval, gyval, axis=0)
					x_sort = np.argsort(dgxval)
					dgxval = dgxval[x_sort]
					dgyval = dgyval[x_sort]
					c_range = np.where((self.wave>=dgxval[0]) & (self.wave <= dgxval[-1]))
					dg_wave = self.wave[c_range]
					dg_flux = self.flux[c_range]

					print(self.gauss_profiles)

				else:
					self.gauss_profiles = []

		elif event.key == 'D':
			# delete previous unwanted points for Gaussian profile fitting
			fclick = len(self.gxval)

			if (fclick > 0) & (fclick <3):
				self.gxval = np.delete(self.gxval, -1)
				self.gyval = np.delete(self.gyval, -1)
				del self.axes.lines[-1]
			else:
				del self.axes.lines[-4:]			
			self.draw()

			self.gauss_profiles = []

		#elif event.key == 'J':
		#	self.xdata = event.xdata
		#	self.guess_ion = GuessTransition(self.linelist, event.xdata)
		#	self.guess_ion.show()

	def onclick(self, event):
		'''Mouse click
			Left == 1; Right == 3
		'''
		if event.button == 3:
			#For single Gaussian
			if self.gauss_num < 2:
				self.send_message.emit(f'Currently, we need {self.gauss_num} Gaussian to guess the line position.')
				
				if self.guess_gcenter:
					self.guess_ion = GuessTransition(self.linelist, self.guess_gcenter[0], self.guess_gcenter[-1])
					self.guess_ion.show()
					self.guess_ion.send_z_cal.connect(self._on_estZ_changed)
				else:
					self.send_message.emit(f'Please fit a Gaussian profile FIRST\n'
											f'to locate the line CENTER!!!')
			#print(self.estZ)

			#For double Gaussian
			# placeholder for now...
			else:
				self.send_message.emit(f'Currently, we need {self.gauss_num} Gaussians to guess the line positions.')
				if self.guess_gcenter < 0:
					self.send_message.emit(f'Please fit {self.gauss_num} Gaussian profiles FIRST\n'
											f'to locate the line CENTERS!!!')
				else:
					pass

	def replot(self, new_spec, new_err):
		'''Re-plot smoothed/unsmoothed spectrum
		'''
		axes = self.figure.gca()
		xlim = axes.get_xlim()
		ylim = axes.get_ylim()
		self.axes.lines[0] = axes.plot(self.wave, self.new_err, color='red')# label='Error')
		self.axes.lines[1] = axes.plot(self.wave, self.new_spec, color='black')#, label='Flux')
		#self.axes.legend(loc='upper right')
		del self.axes.lines[-2:]


	def _compute_distance(self, gxval, gyval, event):
		'''Compute the distance between the event xydata and selected point for Gaussian fitting
		'''
		diffxval = gxval - event.xdata
		diffyval = gyval - event.ydata
		distance = np.array([np.sqrt(diffxval[i]**2 + diffyval[i]**2) for i in range(len(gxval))])
		return distance

	def _clear_plotted_lines(self):
		while self.axes.texts:
			self.axes.texts.pop()
		while self.axes.collections:
			self.axes.collections.pop()
		self.draw()

	def _plot_lines(self, lineindex, estZ=0.):
		axes = self.figure.gca()
		xlim, ylim = axes.get_xlim(), axes.get_ylim()
		self._clear_plotted_lines()

		if lineindex < 0:

			tmp_lines = self.linelist['wave'].to_numpy() * (1+estZ)
			self.axes.vlines(x=tmp_lines,
							 ymin=ylim[0], ymax=ylim[-1], color='blue', linestyle='dashed')
			for row in range(self.linelist.shape[0]):
				self.axes.text(x=tmp_lines[row],
							   y=ylim[-1]*0.6,
							   s=self.linelist.at[row, 'name'],
							   color='blue', fontsize=15, rotation='vertical')
			#print(lineindex)
			#print(self.linelist)
		else:
			self.axes.vlines(x=self.linelist.at[lineindex, 'wave'] * (1+estZ),
							 ymin=ylim[0], ymax=ylim[-1], color='blue', linestyle='dashed')
			self.axes.text(x=self.linelist.at[lineindex, 'wave'] * (1+estZ),
						   y=ylim[-1]*0.6,
						   s=self.linelist.at[lineindex, 'name'],
						   color='blue', fontsize=15, rotation='vertical')

		self.axes.set_xlim(xlim)
		self.axes.set_ylim(ylim)
		self.draw()

		#print(self.axes.texts)
		#print('vlines num: ', len(self.axes.collections))
		
	def on_linelist_slot(self, sent_linelist):
		self.linelist = sent_linelist #.append(sent_linelist, ignore_index=True)
		#self._clear_plotted_lines()
		self._plot_lines(-1, 0.)
		#print(self.linelist)

	def on_lineindex_slot(self, sent_lineindex):
		#print(sent_lineindex == 1)
		if sent_lineindex == 0:
			self._clear_plotted_lines()
		elif sent_lineindex == 1:
			self.lineindex = -1
			self._plot_lines(self.lineindex)
		else:
			self.lineindex = sent_lineindex - 2
			self._plot_lines(self.lineindex)



		

	def _on_estZ_changed(self, newz):
		self.estZ = newz[0]
		self.estZstd = newz[1]
		self._plot_lines(self.lineindex, self.estZ)
		self.send_z_est.emit([self.estZ, self.estZstd])

		#print(self.estZ)

	def _on_sent_gauss_num(self, sent_gauss_num):
		self.gauss_num = int(sent_gauss_num)






class CustomLimDialog(QtWidgets.QDialog):
	def __init__(self, axis='y'):
		super().__init__()

		self.setWindowTitle('Set precise limits on flux range')
		QBtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
		self.buttonbox = QtWidgets.QDialogButtonBox(QBtn)

		self.layout = QtWidgets.QGridLayout()
		if axis == 'x':
			lb_min = QtWidgets.QLabel('X-Axis Min')
			lb_max = QtWidgets.QLabel('X-Axis Max')
		if axis == 'y':
			lb_min = QtWidgets.QLabel('Y-Axis Min')
			lb_max = QtWidgets.QLabel('Y-Axis Max')
		else:
			pass

		self.le_min = QtWidgets.QLineEdit()
		self.le_min.setPlaceholderText('Desired minimum value')
		self.le_max = QtWidgets.QLineEdit()
		self.le_max.setPlaceholderText('Desired maximum value')

		self.layout.addWidget(lb_min, 0, 0)
		self.layout.addWidget(lb_max, 1, 0)
		self.layout.addWidget(self.le_min, 0, 1)
		self.layout.addWidget(self.le_max, 1, 1)
		self.layout.addWidget(self.buttonbox, 2, 1)

		self.setLayout(self.layout)

		self.buttonbox.accepted.connect(self.accept)
		self.buttonbox.rejected.connect(self.reject)
		

	def _getlim(self):
		return [float(self.le_min.text()), float(self.le_max.text())]
