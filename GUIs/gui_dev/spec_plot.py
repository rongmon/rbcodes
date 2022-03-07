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

matplotlib.use('Qt5Agg')

class MplCanvas(FigureCanvasQTAgg):
	send_message = pyqtSignal(str)
	send_z_est = pyqtSignal(list)
	send_gcenter = pyqtSignal(list)
	send_scale_limits = pyqtSignal(list)

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
		self.addtional_linelist_idx = []
		self.addtional_linelist_vlines = []
		self.addtional_linelist_texts = []

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

	def _plot_additional_lines(self, linelist_dir, z_guess=0.):
		axes = self.figure.gca()
		xlim, ylim = axes.get_xlim(), axes.get_ylim()
		# double check with the colors in linelist_selection.py
		colors = ['#A52A2A', '#FF7F50', '#40E0D0', '#DAA520', '#008000', '#4B0082']

		# this is the secondary linelist index
		idx = list(linelist_dir.keys())[0]

		if type(linelist_dir[idx]) is str:
			# None passed here
			if idx in self.addtional_linelist_idx:
				loc = self.addtional_linelist_idx.index(idx)

				# delete previous vlines for selected linelist
				vline_loc = self.axes.collections.index(self.addtional_linelist_vlines[loc])
				del self.axes.collections[vline_loc]
				del self.addtional_linelist_vlines[loc]
				# delete previous texts for selected linelist
				for line_text in self.addtional_linelist_texts[loc]:
					ti = self.axes.texts.index(line_text)
					del self.axes.texts[ti]
				del self.addtional_linelist_texts[loc]
				# delete previous linelist index
				del self.addtional_linelist_idx[loc]
				
		else:
			# selected a linelist
			tmp_lines = linelist_dir[idx]['wave'].to_numpy() * (1+z_guess)
			#print(idx)
			#print(self.addtional_linelist_idx)
			if idx in self.addtional_linelist_idx:
			# check if linelist i plotted or not
				loc = self.addtional_linelist_idx.index(idx)

				# delete previous vlines for selected linelist
				vline_loc = self.axes.collections.index(self.addtional_linelist_vlines[loc])
				del self.axes.collections[vline_loc]
				del self.addtional_linelist_vlines[loc]
				# delete previous texts for selected linelist
				for line_text in self.addtional_linelist_texts[loc]:
					ti = self.axes.texts.index(line_text)
					del self.axes.texts[ti]
				del self.addtional_linelist_texts[loc]
				# delete previous linelist index
				del self.addtional_linelist_idx[loc]
				

			# add new lines..
			# secondary linelist showed up for the first time
			self.addtional_linelist_idx.append(idx)

			# not plotted before, just add new one
			tmp_vlines = self.axes.vlines(x=tmp_lines,
							 ymin=ylim[0], ymax=ylim[-1], 
							 color=colors[idx], 
							 linestyle='dashed')
			self.addtional_linelist_vlines.append(tmp_vlines)

			tmp_texts = []
			for row in range(linelist_dir[idx].shape[0]):
				tmp_txt = self.axes.text(x=tmp_lines[row],
							   y=ylim[-1]*0.6,
							   s=linelist_dir[idx].at[row, 'name'],
							   color=colors[idx], 
							   fontsize=12, 
							   rotation='vertical')
				tmp_texts.append(tmp_txt)
			#print(tmp_texts)
			self.addtional_linelist_texts.append(tmp_texts)

		#print(self.axes.collections)
		#print(self.axes.texts)

		self.axes.set_xlim(xlim)
		self.axes.set_ylim(ylim)
		self.draw()
		
	def gauss(self, x, amp, mu, sigma):
		return amp * np.exp(-(x-mu)**2/(2. * sigma**2))

#--------------------- Methods for 2D data-------------------------

	def extract_1d(self, flux2d, error2d):
		'''Extract 1d spectrum and error spectrum from input 2d spectrum
		'''
		flux1d=np.nansum(flux2d, axis=0)
		sub_var = error2d**2.
		err1d = np.sqrt(np.sum(sub_var, axis=0))
		return flux1d, err1d

	def plot_spec2d(self, wave, flux, error, filename, scale=0, normalization=0):
		'''Display 2D spec in top panel, 1D extraction in bottom panel
		self.flux2d, self.err2d - 2D spec info
		self.flux, self.error - 1D spec extraction
		'''
		self.fig.clf()
		#self.fig.set_size_inches(w=self.figsize[0], h=self.figsize[-1]*2)
		self.ax2d = self.fig.add_subplot(211)
		self.axes = self.fig.add_subplot(212, sharex = self.ax2d)

		# data processing here... this should be generalized
		self.flux2d = flux
		self.err2d = error
		self.wave = wave

		# sum in dispersion direction to do initial selection
		tmp = np.sum(self.flux2d + np.abs(self.flux2d.min()), axis=1)
		# make sure tmp has no negative values for initial guess
		tmp_cumsum = np.cumsum(tmp) / np.sum(tmp)
		ylist = np.arange(0, len(tmp_cumsum), 1)
		self.flux, self.error = self.extract_1d(self.flux2d, self.err2d)

		lower, upper = 0.16, 0.84
		self.extraction_y = [int(np.interp(lower, tmp_cumsum, ylist)),
							int(np.interp(upper, tmp_cumsum, ylist))]
		self.tmp_extraction_y = []
		#print(self.extraction_y)
		self.flux, self.error = self.extract_1d(self.flux2d[self.extraction_y[0]: self.extraction_y[1], :],
												self.err2d[self.extraction_y[0]: self.extraction_y[1], :])

		# plot starting...
		# 1d spec plot... (keep same varname as axes in plot_spec)
		self.axes.plot(wave, self.error, color='red')
		self.axes.plot(wave, self.flux, color='black')
		self.axes.set_xlabel('Wavelength')
		self.axes.set_ylabel('Flux')
		xlim_spec1d = self.axes.get_xlim()
		self.axes.set_xlim(xlim_spec1d)

		# 2 spec plot...
		# scaling first
		if scale == 0:
			# Linear.. nothing happened
			scaled2d = copy.deepcopy(self.flux2d)
		elif scale == 1:
			# log transformation.. scaled = log(1+ img)/log(1+img_max)
			if self.flux2d.min() < 0:
				scaled2d = np.log(-self.flux2d.min() + self.flux2d) / np.log(-self.flux2d.min() + self.flux2d.max())
			else:
				scaled2d = np.log(1 + self.flux2d) / np.log(1 + self.flux2d.max())
		elif scale == 2:
			# square root transformation.. 
			# pixel values > 0 ==> regular sqrt; pixel values <0 ==> 1.absolute value 2.sqrt 3.add minus sign
			scaled2d = copy.deepcopy(self.flux2d)
			scaled2d[scaled2d>=0] = np.sqrt(scaled2d[scaled2d>=0])
			scaled2d[scaled2d<0] = -np.sqrt(-scaled2d[scaled2d<0])
		elif scale == 3:
			scaled2d = self.flux2d**2

		# normalization next
		# send scaling limits back to toolbar
		if type(normalization) == int:
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
			scaled2d = (scaled2d - normalization[0]) / (normalization[1] - normalization[0])
		



		pos_ax2d = self.ax2d.imshow(scaled2d, origin='lower', vmin=scaled2d.min(), vmax=scaled2d.max() * 0.01,
									extent=(wave[0], wave[-1], 0, len(flux)))
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

		# a dummy var to tell if ax2d is available later for event.key='C'
		self.axnum = 2

	def replot2d(self, wave, new_spec):
		'''Re-plot smoothed/unsmoothed spectrum
		'''
		axes = self.figure.gca()
		xlim = axes.get_xlim()
		ylim = axes.get_ylim()
		self.axes.lines[0] = axes.plot(wave, new_spec, color='black')#, label='Flux')
		#self.axes.legend(loc='upper right')
		del self.axes.lines[-1:]

#------------------- Keyboards/Mouse Events------------------------
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

			print(event.inaxes == self.axes)
		elif event.key == 't':
			# set y axis max value
			if event.inaxes == self.axes:
				ylim = self.axes.get_ylim()
				self.axes.set_ylim([ylim[0], event.ydata])
				self.draw()

		elif event.key == 'b':
			# set y axis min value
			if event.inaxes == self.axes:
				ylim = self.axes.get_ylim()
				self.axes.set_ylim([event.ydata, ylim[-1]])
				self.draw()

		elif event.key == 'X':
			# set x axis max value
			if event.inaxes == self.axes:
				xlim = self.axes.get_xlim()
				self.axes.set_xlim([xlim[0], event.xdata])
				self.draw()

		elif event.key == 'x':
			# set x axis min value
			if event.inaxes == self.axes:
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
			if event.inaxes == self.axes:
				self.scale += 2
				self.new_spec = convolve(self.flux, Box1DKernel(self.scale))
				if len(self.error) > 0:
					self.new_err = convolve(self.error, Box1DKernel(self.scale))
					self.replot(self.wave, self.new_spec, self.new_err)
				else:
					self.replot2d(self.wave, self.new_spec)
				self.draw()

		elif event.key == 'U':
			# unsmooth ydata
			if event.inaxes == self.axes:
				self.scale -= 2
				if self.scale < 0:
					self.scale = 1
				self.new_spec = convolve(self.flux, Box1DKernel(self.scale))
				if len(self.error) > 0:
					self.new_err = convolve(self.error, Box1DKernel(self.scale))
					self.replot(self.wave, self.new_spec, self.new_err)
				else:
					self.replot2d(self.wave, self.new_spec)
				self.draw()

		elif event.key == 'Y':
			# set y-axis limits with precise values
			if event.inaxes == self.axes:
				ylimdialog = CustomLimDialog(axis='y')
				if ylimdialog.exec_():
					ylim = ylimdialog._getlim()
					self.axes.set_ylim(ylim)
					self.draw()
			elif event.inaxes == self.ax2d:
				ylimdialog2d = CustomLimDialog(axis='y')
				if ylimdialog2d.exec_():
					ylim2d = ylimdialog2d._getlim()
					self.ax2d.set_ylim(ylim2d)
					self.draw()

		elif event.key == 'W':
			# set y-axis limits with precise values
			if event.inaxes == self.axes:
				xlimdialog = CustomLimDialog(axis='x')
				if xlimdialog.exec_():
					xlim = xlimdialog._getlim()
					self.axes.set_xlim(xlim)
					self.draw()
			elif event.inaxes == self.ax2d:
				xlimdialog2d = CustomLimDialog(axis='x')
				if xlimdialog2d.exec_():
					xlim2d = xlimdialog2d._getlim()
					self.ax2d.set_xlim(xlim2d)
					self.draw()


		elif event.key == 'G':
			# fit a Gaussian profile
			if event.inaxes == self.axes:
				self.gxval = np.append(self.gxval, event.xdata)
				self.gyval = np.append(self.gyval, event.ydata)

				if self.gauss_num > 1:
					dgxval, dgyval = [],[]

				fclick = len(self.gxval)
				self.axes.plot(self.gxval[-1], self.gyval[-1], 'rs', ms=5)
				self.draw()

				if fclick == 1:
					message = 'You need 2 points to model a Gaussian. Please click 2 more points.'
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

		elif event.key == 'C':
			# change the size of the extraction box in 2D spec plot
			if event.inaxes != self.axes:
				if self.axnum > 1:
					boxlines = self.ax2d.plot(event.xdata, event.ydata, 'r+')
					self.tmp_extraction_y = np.append(self.tmp_extraction_y, event.ydata)
					if len(self.tmp_extraction_y) == 2:
						# delete old hlines
						while self.ax2d.collections:
							self.ax2d.collections.pop()

						(ext_min_y, ext_max_y) = (int(np.round(min(self.tmp_extraction_y))),
												int(np.round(max(self.tmp_extraction_y))))
						self.ax2d_xlim = self.ax2d.get_xlim()
						self.ax2d.hlines(ext_min_y, self.ax2d_xlim[0], self.ax2d_xlim[1], color='red', linestyle='dashed')
						self.ax2d.hlines(ext_max_y, self.ax2d_xlim[0], self.ax2d_xlim[1], color='red', linestyle='dashed')
						flux2d = self.flux2d[ext_min_y:ext_max_y, :]
						err2d = self.err2d[ext_min_y:ext_max_y, :]
						self.new_spec = np.nansum(flux2d, axis=0)
						self.new_err = np.sqrt(np.nansum((err2d**2), axis=0))

						self.replot(self.wave, self.new_spec, self.new_err)
						self.tmp_extraction_y = []
						while self.ax2d.lines:
							self.ax2d.lines.pop()
					self.draw()
				else:
					message = "You don't have a 2D Spectrum plot available."
					self.send_message.emit(message)
		elif event.key == 'H':
			if event.inaxes != self.axes:
				#bring up the flux histogram dialog
				fhist = FluxHistogram(self.flux2d)
				fhist.exec_()

		elif event.key == 'P':
			if event.inaxes != self.axes:
				#bring up the pixel histogram dialog
				phist = PixelHistogram(self.flux2d)
				phist.exec_()
			
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

#-------------------- Slots for External Signals ------------------
	def on_linelist_slot(self, sent_linelist):
		# if no linelist selected, a str is passed along
		if type(sent_linelist) is str:
			self._clear_plotted_lines()
		else:
			self.linelist = sent_linelist
			self._plot_lines(-1)

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

	def on_additional_linelist_slot(self, addtional_linelist):
		self._plot_additional_lines(addtional_linelist)
	def on_additional_linelist_slot_z(self, linelist_and_z):
		addtional_linelist = linelist_and_z[0]
		z_guess = linelist_and_z[1]
		self._plot_additional_lines(linelist_dir=addtional_linelist,
									z_guess=z_guess)

	def _on_estZ_changed(self, newz):
		self.estZ = newz[0]
		self.estZstd = newz[1]
		self._plot_lines(self.lineindex, self.estZ)
		self.send_z_est.emit([self.estZ, self.estZstd])

		#print(self.estZ)

	def _on_sent_gauss_num(self, sent_gauss_num):
		self.gauss_num = int(sent_gauss_num)

	def _update_lines_for_newfile(self, sent_filename):
		if len(self.linelist) > 0:
			# default value of self.linindex = -2
			if self.lineindex > -2:
				self._plot_lines(self.lineindex)
		else:
			self._clear_plotted_lines()


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
