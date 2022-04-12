from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, pyqtSignal

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from astropy.convolution import convolve, Box1DKernel
from astropy.modeling import models, fitting
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
import numpy as np
import copy


from guess_transition import GuessTransition
from spec_hist import FluxHistogram, PixelHistogram
from spec_fit_gauss2d import Gaussfit_2d

matplotlib.use('Qt5Agg')

class MplCanvas(FigureCanvasQTAgg):
	send_message = pyqtSignal(str)
	send_z_est = pyqtSignal(list)
	send_z_manual = pyqtSignal(float)
	send_gcenter = pyqtSignal(list)
	send_scale_limits = pyqtSignal(list)
	send_gauss_num = pyqtSignal(int)
	send_extract1d = pyqtSignal(dict)

	def __init__(self, parent=None, width=5, height=3, dpi=100):
		self.figsize = [width, height]
		self.fig = Figure(figsize=(self.figsize[0], self.figsize[-1]), dpi=dpi)
		pad = 0.05
		self.fig.subplots_adjust(left=pad+0.02, bottom=pad*2, right=0.995, top=0.95)

		#self.axes.margins(x=0)

		# initialize some plotting parameters
		self.wave, self.flux, self.error = [],[],[]
		self.init_xlims, self.init_ylims = [],[]
		self.gxval, self.gyval = [], []
		self.scale = 1.
		self.lineindex = -2
		self.linelist = [] #pd.DataFrame(columns=['wave', 'name'])
		self.estZ = 0.
		self.estZstd = 0.
		self.guess_ion = 0.
		self.guess_gcenter = []
		self.gauss_num = 1
		self.gauss_profiles = []
		self.gauss2d = None
		self.linelists2multiG = []
		self.stamp = None

		num_2ndlist = 6
		self.addtional_linelist = {i:[] for i in range(num_2ndlist)}
		self.addtional_linelist_z = {i:0. for i in range(num_2ndlist)}

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
		self.fig.clf()
		self.axes = self.fig.add_subplot(111)
		self.axes.cla()
		self.axes.lines[0] = self.axes.plot(wave, error, color='red')# label='Error')
		self.axes.lines[1] = self.axes.plot(wave, flux, color='black')#, label='Flux')
		#self.axes.legend(loc='upper right')
		self.axes.set_ylim([-np.nanmin(flux)*0.01, np.nanmedian(flux)*3])
		self.axes.set_xlim([np.min(wave), np.max(wave)])
		self.axes.set_xlabel('Wavelength (Angstrom)')
		self.axes.set_ylabel('Flux')
		self.axes.set_title(filename)
		self.draw()

		self.send_message.emit(f'You are currently working on {filename} spectrum.')

		# update intiial plotting parameters
		self.wave, self.flux, self.error = wave, flux, error
		self.axnum = 1
		self.init_xlims = self.axes.get_xlim()
		self.init_ylims = self.axes.get_ylim()

	

	def replot(self, wave, new_spec, new_err):
		'''Re-plot smoothed/unsmoothed spectrum
		'''
		axes = self.figure.gca()
		xlim = axes.get_xlim()
		ylim = axes.get_ylim()
		self.axes.lines[0] = axes.plot(wave, new_err, color='red')# label='Error')
		self.axes.lines[1] = axes.plot(wave, new_spec, color='black')#, label='Flux')
		self.axes.set_ylim([np.nanmin(new_spec), np.nanmax(new_spec)])
		self.axes.set_xlim([np.min(wave), np.max(wave)])

		del self.axes.lines[2:]

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

	def _select_lines_within_xlim(self, linelist, estZ=0., xbound=0.):
		# select lines within axes.xlim
		axes = self.figure.gca()
		xlim = axes.get_xlim()
		tmp_lines = linelist['wave'].to_numpy() * (1+estZ)
		tmp_names = linelist['name'].to_numpy()
		tmp_names = tmp_names[(tmp_lines > xlim[0]*(1+xbound)) & (tmp_lines < xlim[1]*(1-xbound))]
		tmp_lines = tmp_lines[(tmp_lines > xlim[0]*(1+xbound)) & (tmp_lines < xlim[1]*(1-xbound))]
		return tmp_lines, tmp_names

	def _plot_lines(self, lineindex, estZ=0.):
		axes = self.figure.gca()
		xlim, ylim = axes.get_xlim(), axes.get_ylim()
		xbound = 0.00 # leave non-drawing space at boundary  

		if lineindex < 0:
			tmp_lines, tmp_names = self._select_lines_within_xlim(linelist=self.linelist,
																	estZ=estZ,
																	xbound=xbound)		
			if len(tmp_lines) > 0:
				self.axes.vlines(x=tmp_lines,
								 ymin=ylim[0], ymax=ylim[-1], color='blue', linestyle='dashed')
				for row in range(len(tmp_lines)):
					self.axes.text(x=tmp_lines[row],
								   y=ylim[-1]*0.6,
								   s=tmp_names[row],
								   color='blue', fontsize=15, rotation='vertical')
			else:
				message = 'No ions selected. No lines are within current wavelength range.'
				self.send_message.emit(message)
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

	def _plot_additional_lines(self, linelist_dir, z_guess_dir):
		# secondary linelist plotting function
		# linelist_dir values only contain [] or linelist from secondary linelist widget
		axes = self.figure.gca()
		xlim, ylim = axes.get_xlim(), axes.get_ylim()
		xbound = 0.00
		# double check with the colors in linelist_selection.py
		colors = ['#A52A2A', '#FF7F50', '#40E0D0', '#DAA520', '#008000', '#4B0082']

		for key, linelist in linelist_dir.items():
			if len(linelist) > 0:
				# linelist is not empy
				tmp_lines, tmp_names = self._select_lines_within_xlim(linelist=linelist,
																	estZ=z_guess_dir[key],
																	xbound=xbound)
				if len(tmp_lines) > 0:
					self.axes.vlines(x=tmp_lines,
									ymin=ylim[0], ymax=ylim[-1],
									color=colors[key],
									linestyle='dashed')
					for row in range(len(tmp_lines)):
						self.axes.text(x=tmp_lines[row],
										y=ylim[-1]*0.6,
										s=tmp_names[row],
										color=colors[key],
										fontsize=12,
										rotation='vertical')
		#print(self.axes.collections)
		#print(self.axes.texts)

		self.axes.set_xlim(xlim)
		self.axes.set_ylim(ylim)
		self.draw()

	def _lines_in_current_range(self):
		# wrapper to plot primary and secondary linelists
		self._clear_plotted_lines()
		if len(self.linelist) > 0:
			self._plot_lines(self.lineindex, self.estZ)
		self._plot_additional_lines(self.addtional_linelist,
									self.addtional_linelist_z)
		
	def gauss(self, x, amp, mu, sigma):
		return amp * np.exp(-(x-mu)**2/(2. * sigma**2))

#--------------------- Methods for 2D data-------------------------

	def extract_1d(self, flux2d, error2d):
		'''Extract 1d spectrum and error spectrum from input 2d spectrum
		'''
		flux1d=np.nansum(flux2d, axis=0)
		sub_var = error2d**2.
		err1d = np.sqrt(np.nansum(sub_var, axis=0))

        #Create weights
		#wt=np.nansum(flux2d,axis=1)
		#wt+=np.abs(np.min(wt))
		#wt/=np.sum(wt)

		#tt=np.shape(flux2d)
		#print(tt)
		#print(np.shape(flux2d))
		#temp=np.transpose(np.tile(wt,(tt[1],1)))
		#flux2d *=temp
		#norm = (np.max(flux1d) / np.max(flux2d,axis=(0,1)))
		#fl_sum = flux2d * norm
		#fl1d_opt=np.nansum(fl_sum,axis=0)
		#sig=np.sqrt(np.nansum(sub_var * (temp ** 2), axis=0 ))   # * n_spaxels
		#sig = sig * norm
		#print(np.shape(fl_sum))

		#print(np.shape(sig))


		return flux1d, err1d
		#return fl1d_opt, err1d


	def plot_spec2d(self, wave, flux2d, error2d, filename, scale=0, normalization=0, pre_extraction=None):
		'''Display 2D spec in top panel, 1D extraction in bottom panel
		self.flux2d, self.err2d - 2D spec info
		self.flux, self.error - 1D spec extraction
		'''
		self.fig.clf()
		#self.fig.set_size_inches(w=self.figsize[0], h=self.figsize[-1]*2)
		self.ax2d = self.fig.add_subplot(211)
		self.axes = self.fig.add_subplot(212, sharex = self.ax2d)

		# data processing here... this should be generalized
		# REPLACE NAN WITH ZERO
		self.flux2d = flux2d
		self.err2d = error2d
		self.wave = wave

		# copy of flux2d with nan replaced by 0
		img = np.nan_to_num(self.flux2d, nan=0.)
		
		if pre_extraction is None:
			
			# sum in dispersion direction to do initial selection
			tmp = np.sum(img + np.abs(img.min()), axis=1)
			# make sure tmp has no negative values for initial guess
			tmp_cumsum = np.cumsum(tmp) / np.sum(tmp)
			ylist = np.arange(0, len(tmp_cumsum), 1)
			# flux/error are active and extracted from flux2d/error2d
			self.flux, self.error = self.extract_1d(self.flux2d, self.err2d)
			# keep a frozen copy and only used for reset
			self.flux_fix, self.error_fix = self.flux.copy(), self.error.copy()

			#print(tmp_cumsum)
			lower, upper = 0.16, 0.84
			self.extraction_y = [int(np.interp(lower, tmp_cumsum, ylist)),
								int(np.interp(upper, tmp_cumsum, ylist))]			
			#print(self.extraction_y)
		else:
			self.extraction_y = pre_extraction

		self.flux, self.error = self.extract_1d(self.flux2d[self.extraction_y[0]: self.extraction_y[1], :],
												self.err2d[self.extraction_y[0]: self.extraction_y[1], :])
		self.tmp_extraction_y = []
		
		# plot starting...
		# 1d spec plot... (keep same varname as axes in plot_spec)
		self.axes.plot(wave, self.error, color='red')
		self.axes.plot(wave, self.flux, color='black')
		self.axes.set_xlabel('Wavelength (Angstrom)')
		self.axes.set_ylabel('Flux')
		#xlim_spec1d = self.axes.get_xlim()
		#self.axes.set_xlim(xlim_spec1d)
		self.axes.set_xlim([np.min(wave), np.max(wave)])

		# 2 spec plot...
		# scaling first
		if scale == 0:
			# Linear.. nothing happened
			#scaled2d = self.flux2d.copy()
			scaled2d = img.copy()
		elif scale == 1:
			# log transformation.. scaled = log(1+ img)/log(1+img_max)
			#if self.flux2d.min() < 0:
			#	scaled2d = np.log(-self.flux2d.min() + self.flux2d) / np.log(-self.flux2d.min() + self.flux2d.max())
			if img.min() < 0:
				scaled2d = np.log(-img.min() + img) / np.log(-img.min() + img.max())
			else:
				#scaled2d = np.log(1 + self.flux2d) / np.log(1 + self.flux2d.max())
				scaled2d = np.log(1 + img) / np.log(1 + img.max())
		elif scale == 2:
			# square root transformation.. 
			# pixel values > 0 ==> regular sqrt; pixel values <0 ==> 1.absolute value 2.sqrt 3.add minus sign
			#scaled2d = copy.deepcopy(self.flux2d)
			scaled2d = copy.deepcopy(img)
			scaled2d[scaled2d>=0] = np.sqrt(scaled2d[scaled2d>=0])
			scaled2d[scaled2d<0] = -np.sqrt(-scaled2d[scaled2d<0])
		elif scale == 3:
			#scaled2d = self.flux2d**2
			scaled2d = img**2

		# normalization next
		# send scaling limits back to toolbar
		if type(normalization) is int:
			#print(scaled2d)
			if normalization == 0:
				
				self.send_scale_limits.emit([scaled2d.min(), scaled2d.max()])
			elif normalization == 1:
				# minmax 100% range
				scaled2d = (scaled2d - scaled2d.min()) / (scaled2d.max() - scaled2d.min())
			elif normalization < 10: # this magic num from n_combobox in toolbar
				if normalization == 2: # 99.5%
					low, up = np.percentile(scaled2d, [0.25, 99.75])
				elif normalization == 3: # 99%
					low, up = np.percentile(scaled2d, [0.5, 99.5])
				elif normalization == 4: # 98%
					low, up = np.percentile(scaled2d, [1., 99.])
				elif normalization == 5: # 97%
					low, up = np.percentile(scaled2d, [1.5, 98.5])
				elif normalization == 6: # 96%
					low, up = np.percentile(scaled2d, [2., 98.])
				elif normalization == 7: # 95%
					low, up = np.percentile(scaled2d, [2.5, 97.5])
				elif normalization == 8: # 92.5%
					low, up = np.percentile(scaled2d, [3.75, 96.25])
				elif normalization == 9: # 90%
					low, up = np.percentile(scaled2d, [5., 95.])

				scaled2d = (scaled2d - low) / (up - low)
				self.send_scale_limits.emit([scaled2d.min(), scaled2d.max()])

			elif normalization == 10: # Z-score
				scaled2d = (scaled2d - scaled2d.mean()) / scaled2d.std()
				self.send_scale_limits.emit([np.nan, np.nan])	

		elif type(normalization) == list:
			tmp = (scaled2d - scaled2d.min()) / (scaled2d.max() - scaled2d.min())
			scaled2d = tmp*(normalization[1] - normalization[0]) + normalization[0]
			
		if scale == 1:
			pos_ax2d = self.ax2d.imshow(scaled2d, origin='lower', 
									vmin=scaled2d.min(), vmax=scaled2d.max() * 0.01,
									extent=(self.wave[0], self.wave[-1], 0, len(self.flux2d)))
		else:
			pos_ax2d = self.ax2d.imshow(scaled2d, origin='lower', 
									vmin=scaled2d.min(), vmax=scaled2d.max() * 1.0,
									extent=(self.wave[0], self.wave[-1], 0, len(self.flux2d)))
		

		self.ax2d_cb = self.fig.colorbar(pos_ax2d, ax=self.ax2d, location='top')
		ax2d_xlim = self.ax2d.get_xlim()
		self.ax2d.hlines(self.extraction_y[0], ax2d_xlim[0], ax2d_xlim[1], color='red', linestyle='dashed')
		self.ax2d.hlines(self.extraction_y[1], ax2d_xlim[0], ax2d_xlim[1], color='red', linestyle='dashed')
		self.ax2d.tick_params(labelbottom=False)
		#self.ax2d.set_xlim(xlim_spec1d)
		self.ax2d.set_aspect('auto')

		self.draw()
		self.init_xlims = self.axes.get_xlim()
		self.init_ylims = self.axes.get_ylim()

		#del scaled2d

		# a dummy var to tell if ax2d is available later for event.key='C'
		self.axnum = 2

	def replot2d(self, wave, new_spec, new_err, new_extraction):
		'''Re-plot smoothed/unsmoothed spectrum
		'''
		axes = self.figure.gca()
		xlim = axes.get_xlim()
		ylim = axes.get_ylim()
		self.axes.lines[0] = axes.plot(wave, new_err, color='red')# label='Error')
		self.axes.lines[1] = axes.plot(wave, new_spec, color='black')#, label='Flux')
		del self.axes.lines[2:]

		ax2d_xlim = self.ax2d.get_xlim()
		while self.ax2d.collections:
			self.ax2d.collections.pop()
		while self.ax2d.lines:
			self.ax2d.lines.pop()
		self.ax2d.hlines(new_extraction[0], np.min(wave), np.max(wave), color='red', linestyle='dashed')
		self.ax2d.hlines(new_extraction[1], np.min(wave), np.max(wave), color='red', linestyle='dashed')




#------------------- Keyboards/Mouse Events------------------------
	def ontype(self, event):
		'''Interactivae keyboard events
		Note:
			Always Update to help mannual for new keyboard events
		'''
		#print(event.key)
		if event.key == 'r':
			# reset flux/err and clear all lines
			#del self.axes.lines[2:]	# delete everything except flux/err
			self.line_xlim, self.line_ylim = [], []
			self.gxval, self.gyval = [], []

			if self.axnum == 1:
				# for 1D spec display only
				self.replot(self.wave, self.flux, self.error)
			else:
				# for 2D spec
				if event.inaxes == self.axes:
					# cursor in 1d canvas
					self.replot(self.wave, self.flux, self.error)
				else:
					# cursor in 2d canvas
					# reset active flux/error to fixed flux/error
					self.flux, self.error = self.flux_fix, self.error_fix
					self.replot2d(self.wave, self.flux_fix, self.error_fix, self.extraction_y)
			self.axes.set_ylim([np.nanmin(self.flux), np.nanmax(self.flux)])
			self.axes.set_xlim([np.min(self.wave), np.max(self.wave)])
			
			self._lines_in_current_range()

			self.draw()
			self.send_message.emit('You RESET the canvas!!!')

			#print(event.inaxes == self.axes)
		elif event.key == 't':
			# set y axis max value
			if event.inaxes == self.axes:
				ylim = self.axes.get_ylim()
				self.axes.set_ylim([ylim[0], event.ydata])
				self.draw()
				self.send_message.emit('Flux MAX value in 1D plot changed.')

		elif event.key == 'b':
			# set y axis min value
			if event.inaxes == self.axes:
				ylim = self.axes.get_ylim()
				self.axes.set_ylim([event.ydata, ylim[-1]])
				self.draw()
				self.send_message.emit('Flux MIN value in 1D plot changed.')

		elif event.key == 'X':
			# set x axis max value
			if event.inaxes == self.axes:
				xlim = self.axes.get_xlim()
				self.axes.set_xlim([xlim[0], event.xdata])
				self._lines_in_current_range()
				self.draw()
				self.send_message.emit('Wavelength MAX value in 1D plot changed.')

		elif event.key == 'x':
			# set x axis min value
			if event.inaxes == self.axes:
				xlim = self.axes.get_xlim()
				self.axes.set_xlim([event.xdata, xlim[-1]])
				self._lines_in_current_range()
				self.draw()
				self.send_message.emit('Wavelength MIN value in 1D plot changed.')

		elif event.key == '[':
			# pan window to the left
			xlim = self.axes.get_xlim()
			delx = (xlim[-1] - xlim[0])
			self.axes.set_xlim([xlim[0] - delx, xlim[0]])
			self._lines_in_current_range()
			self.draw()
			self.send_message.emit('Display window in 1D plot pans LEFT.')

		elif event.key == ']':
			# pan window to the right
			xlim = self.axes.get_xlim()
			delx = (xlim[-1] - xlim[0])
			self.axes.set_xlim([xlim[1], xlim[1] + delx])
			self._lines_in_current_range()
			self.draw()
			self.send_message.emit('Display window in 1D plot pans RIGHT.')

		elif event.key == 'S':
			# smooth ydata
			if event.inaxes == self.axes:
				self.scale += 2
				self.new_spec = convolve(self.flux, Box1DKernel(self.scale))
				self.new_err = convolve(self.error, Box1DKernel(self.scale))
				self.replot(self.wave, self.new_spec, self.new_err)
				self.draw()

				self.send_message.emit(f'Convolutional kernel size = {int(self.scale//2)}.')

		elif event.key == 'U':
			# unsmooth ydata
			if event.inaxes == self.axes:
				self.scale -= 2
				if self.scale < 0:
					self.scale = 1
				self.new_spec = convolve(self.flux, Box1DKernel(self.scale))
				self.new_err = convolve(self.error, Box1DKernel(self.scale))
				self.replot(self.wave, self.new_spec, self.new_err)
				self.draw()
				self.send_message.emit(f'Convolutional kernel size = {int(self.scale//2)}.')

		elif event.key == 'Y':
			# set y-axis limits with precise values
			if event.inaxes == self.axes:
				ylimdialog = CustomLimDialog(axis='y')
				if ylimdialog.exec_():
					ylim = ylimdialog._getlim()
					self.axes.set_ylim(ylim)

					self.draw()
				self.send_message.emit('Flux LIMITS in 1D spectrum are changed.')
			elif event.inaxes == self.ax2d:
				ylimdialog2d = CustomLimDialog(axis='y')
				if ylimdialog2d.exec_():
					self.tmp_extraction_y = ylimdialog2d._getlim()
					if len(self.tmp_extraction_y) == 2:
						(ext_min_y, ext_max_y) = (int(np.round(min(self.tmp_extraction_y))),
												int(np.round(max(self.tmp_extraction_y))))
						self.new_spec, self.new_err = self.extract_1d(self.flux2d[ext_min_y:ext_max_y, :], 
																	self.err2d[ext_min_y:ext_max_y, :])
						self.replot2d(self.wave, self.new_spec, self.new_err, [ext_min_y,ext_max_y])
						self.send_extract1d.emit({'WAVELENGTH': self.wave,
												'FLUX': self.flux,
												'ERROR': self.error,
												'YMIN': ext_min_y,
												'YMAX': ext_max_y})
						self.tmp_extraction_y = []
					self.draw()
				self.send_message.emit('Y-Axis LIMITS in 2D spectrum are changed.')

		elif event.key == 'W':
			# set y-axis limits with precise values
			if event.inaxes == self.axes:
				xlimdialog = CustomLimDialog(axis='x')
				if xlimdialog.exec_():
					xlim = xlimdialog._getlim()
					self.axes.set_xlim(xlim)
					self._lines_in_current_range()
					self.draw()
				self.send_message.emit('Wavelength LIMITS in 1D spectrum are changed.')
			elif event.inaxes == self.ax2d:
				xlimdialog2d = CustomLimDialog(axis='x')
				if xlimdialog2d.exec_():
					xlim2d = xlimdialog2d._getlim()
					self.ax2d.set_xlim(xlim2d)
					self.draw()
				self.send_message.emit('Wavelength LIMITS in 2D spectrum are changed.')


		elif event.key == 'G':
			# fit a Gaussian profile
			if event.inaxes == self.axes:
				if self.gauss_num == 0:
					self.send_message.emit('You are in manual mode. No need to fit Gaussian.')

				else:
					self.gxval = np.append(self.gxval, event.xdata)
					self.gyval = np.append(self.gyval, event.ydata)

					fclick = len(self.gxval)
					self.axes.plot(self.gxval[-1], self.gyval[-1], 'rs', ms=5)
					self.draw()


					if fclick == 1:
						message = 'You need 2 points to model a Gaussian. Please click 1 more point.'
						self.send_message.emit(message)

					elif fclick == 2:
						# sort xdata before fitting
						x_sort = np.argsort(self.gxval)
						gxval = self.gxval[x_sort]
						gyval = self.gyval[x_sort]

						c_range = np.where((self.wave>gxval[0]) & (self.wave < gxval[-1]))
						g_wave = self.wave[c_range]
						g_flux = self.flux[c_range]
						g_error = self.error[c_range]

						if self.gauss_num == 1:
							# Single Gaussian Fitting
							# fit a Gaussian with 3 data points	
							# 1. fit a local continuum
							spline = splrep([gxval[0],gxval[-1]], 
											[gyval[0], gyval[-1]], 
											k=1)
							cont = splev(g_wave, spline)

							# 2. check if it is an absorption or emission line
							EW = np.sum(cont - g_flux)
							if EW > 0:
								# absorption line
								sign = -1
							else:
								# emission line
								sign = 1

							Aguess = np.max(g_flux - cont)
							Cguess = np.mean(g_wave)
							sguess = 0.1 * np.abs(gxval[0] - gxval[1])

							# prepare ydata for fit
							ydata = sign * (g_flux - cont)
							errdata = sign * (g_error - cont)
							# start fitting
							popt, pcov = curve_fit(self.gauss, g_wave, ydata,
													p0=[Aguess, Cguess, sguess],
													sigma=errdata)
							g_final = sign * (self.gauss(g_wave, *popt)) + cont

							perr = np.sqrt(np.diag(pcov))
							model_fit = self.axes.plot(g_wave, g_final, 'r--')
							
							self.draw()

							message = ("A Gaussian model you fit has the following parameters:\n"
									   f"Amplitude: {popt[0]:.3f}\n"
									   f"Mean: {popt[1]:.3f} with std={perr[1]:.3f}\n"
									   f"Sigma: {popt[2]:.3f}")

							if self.guess_gcenter:
								self.guess_gcenter[0] = popt[1]
								self.guess_gcenter[1] = perr[1]
							else:
								self.guess_gcenter.append(popt[1])
								self.guess_gcenter.append(perr[1])


							
							self.send_message.emit(message)
							self.send_gcenter.emit(self.guess_gcenter)

						else:
							print('Multiple Gaussian fitting starts.')
							# Double Gaussian Fitting
							self.axes.fill_between(g_wave,
													y1=np.max(g_flux)*1.1,
													y2=np.min(g_flux)*0.9,
													alpha=0.5,
													color='pink')
							self.draw()
							# delete the drawn polygon from collection
							self.axes.collections.pop()

							self.gauss2d = Gaussfit_2d(g_wave, g_flux, g_error, 
														gauss_num=self.gauss_num,
														linelists=self.linelists2multiG)



						# clear out selection
						self.gxval, self.gyval = [], []

		elif event.key == 'D':
			# delete previous unwanted points for Gaussian profile fitting
			fclick = len(self.gxval)
			'''
			if (fclick > 0) & (fclick <2):
				self.gxval = np.delete(self.gxval, -1)
				self.gyval = np.delete(self.gyval, -1)
				del self.axes.lines[-1]
			else:
				del self.axes.lines[-4:]
			self.gauss_profiles = []
			'''
			if len(self.axes.lines) > 2:
				self.axes.lines.pop()			
			self.draw()
			self.send_message.emit('A previous added object is deleted.')

			

		elif event.key == 'C':
			# change the size of the extraction box in 2D spec plot
			if event.inaxes != self.axes:
				if self.axnum > 1:
					boxlines = self.ax2d.plot(event.xdata, event.ydata, 'r+')
					self.tmp_extraction_y = np.append(self.tmp_extraction_y, event.ydata)
					if len(self.tmp_extraction_y) == 2:
						(ext_min_y, ext_max_y) = (int(np.round(min(self.tmp_extraction_y))),
												int(np.round(max(self.tmp_extraction_y))))
						self.flux, self.error = self.extract_1d(self.flux2d[ext_min_y:ext_max_y, :], 
																	self.err2d[ext_min_y:ext_max_y, :])
						self.send_extract1d.emit({'WAVELENGTH': self.wave,
												'FLUX': self.flux,
												'ERROR': self.error,
												'YMIN': ext_min_y,
												'YMAX': ext_max_y})
						self.replot2d(self.wave, self.flux, self.error, [ext_min_y,ext_max_y])
						self.tmp_extraction_y = []
						self.send_message.emit('A new extraction box is created.')
					self.draw()

					# reset convolution kernel size for new extraction
					self.scale = 1.
				else:
					message = "You don't have a 2D Spectrum plot available."
					self.send_message.emit(message)
		elif event.key == 'H':
			if event.inaxes != self.axes:
				#bring up the flux histogram dialog
				fhist = FluxHistogram(self.flux2d)
				self.send_message.emit('Flux distribution window is popped up.')
				fhist.exec_()

		elif event.key == 'P':
			if event.inaxes != self.axes:
				#bring up the pixel histogram dialog
				phist = PixelHistogram(self.flux2d)
				self.send_message.emit('Pixel distribution window is popped up.')
				phist.exec_()

		elif event.key != 'shift':
			self.send_message.emit('Undefined keyboard function. Please click Help for all key functions.')
			
	def onclick(self, event):
		'''Mouse click
			Left == 1; Right == 3
		'''
		if event.button == 3:
			#Manual mode
			if self.gauss_num == 0:
				self.send_message.emit('You are in manual mode now.')
				self.guess_ion = GuessTransition(self.linelist, event.xdata, 0.)
				self.guess_ion.show()
				self.guess_ion.send_z_cal.connect(self._on_estZ_changed_manual)


			#For single Gaussian
			elif self.gauss_num == 1:
				self.send_message.emit(f'Currently, we need {self.gauss_num} Gaussian to guess the line position.')
				
				if self.guess_gcenter:
					self.guess_ion = GuessTransition(self.linelist, self.guess_gcenter[0], self.guess_gcenter[-1])
					self.guess_ion.show()
					self.guess_ion.send_z_cal.connect(self._on_estZ_changed)
				else:
					self.send_message.emit(f'Please fit a Gaussian profile FIRST\n'
											f'to locate the line CENTER!!!')
			#print(self.estZ)

			#For multiple Gaussian
			else:
				self.send_message.emit(f'Currently, we need {self.gauss_num} Gaussians to guess the line positions.')
				if self.gauss2d is None:
					self.send_message.emit('Please select 2 points to define the range you want to work with')
				else:
					self.gauss2d.exec_()

#-------------------- Slots for External Signals ------------------
	def on_linelist_slot(self, sent_linelist):
		# if no linelist selected, a str is passed along
		if type(sent_linelist) is str:
			self._clear_plotted_lines()
		else:
			self.linelist = sent_linelist
			self._clear_plotted_lines()
			self._plot_lines(-1)

	def on_lineindex_slot(self, sent_lineindex):
		#print(sent_lineindex == 1)
		if sent_lineindex == 0:
			self._clear_plotted_lines()
		elif sent_lineindex == 1:
			self.lineindex = -1
			self._clear_plotted_lines()
			self._plot_lines(self.lineindex)
		else:
			self.lineindex = sent_lineindex - 2
			self._clear_plotted_lines()
			self._plot_lines(self.lineindex)

	def on_additional_linelist_slot(self, addtional_linelist_dir):
		#self._plot_additional_lines(addtional_linelist)
		idx = list(addtional_linelist_dir.keys())[0]
		if type(addtional_linelist_dir[idx]) != str:
			self.addtional_linelist[idx] = addtional_linelist_dir[idx]
		else:
			self.addtional_linelist[idx] = []
		#print(self.addtional_linelist)
		#print(type(self.addtional_linelist))
		self._lines_in_current_range()

	def on_additional_linelist_slot_z(self, linelist_dir_and_z):
		idx = list(linelist_dir_and_z[0].keys())[0]
		self.addtional_linelist[idx] = linelist_dir_and_z[0][idx]
		self.addtional_linelist_z[idx] = linelist_dir_and_z[1]

		self._lines_in_current_range()

	def _on_estZ_changed(self, newz):
		self.estZ = newz[0]
		self.estZstd = newz[1]
		self._lines_in_current_range()
		self.send_z_est.emit([self.estZ, self.estZstd])

	def _on_estZ_changed_manual(self, newz):
		self.estZ = newz[0]
		self._lines_in_current_range()
		self.send_z_manual.emit(self.estZ)

	def _on_estZ_return_pressed(self, sent_estZ):
		self.estZ = sent_estZ
		self._lines_in_current_range()

	def _on_sent_gauss_num(self, sent_gauss_num):
		self.gauss_num = int(sent_gauss_num)

	def _update_lines_for_newfile(self, sent_filename):
		if len(self.linelist) > 0:
			# default value of self.linindex = -2
			if self.lineindex > -2:
				self._clear_plotted_lines()
				self._plot_lines(self.lineindex)
		else:
			self._clear_plotted_lines()

	def _on_sent_linelists2multiG(self, l):
		self.linelists2multiG = l


#-------------------- Dialog Box for XY Ranges --------------------

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

#---------------------------------------------------------------
