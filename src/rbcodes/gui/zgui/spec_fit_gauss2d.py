from PyQt5.QtWidgets import QWidget, QGridLayout, QVBoxLayout, QComboBox, QLineEdit, QLabel, QPushButton, QDialog, QCheckBox, QFileDialog
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QDoubleValidator
from PyQt5 import QtWidgets

import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit

from rbcodes.igm.rb_setline import read_line_list

class Gaussfit_2d(QWidget):
    send_linelist = pyqtSignal(object)
    send_lineindex = pyqtSignal(int)
    send_waves = pyqtSignal(object)
    send_gfinal = pyqtSignal(list)
    send_ransac = pyqtSignal(object)

    def __init__(self,wave, flux1d,error1d, gauss_num=2, linelists=[]):
        super().__init__()
        self.gauss_num = gauss_num
        self.waves = np.array([-1] * gauss_num)
        self.names = np.array([None] * gauss_num)
        self.linelists = linelists
        self.result = None
        self.wave, self.flux1d, self.error1d = wave, flux1d, error1d
        self.kernel_size = 149
        self.cont = None
        self.c = None # placeholder for ransac cont_fitter object
        


        layout = QVBoxLayout()
        # sublayout for Guess z functions
        lines_layout = QGridLayout()
        lines_layout.addWidget(QLabel('Linelist'), 0,0)
        
        self.line_combo = QComboBox()
        self.line_combo.setFixedWidth(120)
        self.line_combo.addItems(self.linelists)
        lines_layout.addWidget(self.line_combo, 1, 0)

        l_zf = QLabel('Estimated z')
        l_zf.setFixedHeight(20)
        lines_layout.addWidget(l_zf, 0, 1)
        self.zf = QLineEdit()
        self.zf.setPlaceholderText('Redshift')
        self.zf.setFixedWidth(100)
        self.zf.setReadOnly(True)
        lines_layout.addWidget(self.zf, 1,1)

        l_zferr = QLabel('Estimated Error')
        l_zferr.setFixedHeight(20)
        lines_layout.addWidget(l_zferr, 0, 2)
        self.zferr = QLineEdit()
        self.zferr.setPlaceholderText('Error')
        self.zferr.setFixedWidth(100)
        self.zferr.setReadOnly(True)
        lines_layout.addWidget(self.zferr, 1,2)

        pb = QPushButton('Fit')
        pb.setFixedWidth(100)
        pb.clicked.connect(self._button_clicked)
        lines_layout.addWidget(pb, 1,3)
        self.fit_result = QLabel('Ready')
        lines_layout.addWidget(self.fit_result, 0, 3)
        adv_pb = QPushButton('Advanced')
        adv_pb.setFixedWidth(100)
        adv_pb.clicked.connect(self._adv_button_clicked)
        lines_layout.addWidget(adv_pb, 2,3)

        apply_pb = QPushButton('Apply')
        apply_pb.setFixedWidth(100)
        apply_pb.clicked.connect(self._apply_button_clicked)
        lines_layout.addWidget(apply_pb, 3,3)

        export_pb = QPushButton('Export')
        export_pb.setFixedWidth(100)
        export_pb.clicked.connect(self._export_button_clicked)
        lines_layout.addWidget(export_pb, 4,3)

        ion1 = self._create_linelist_widget(lines_layout, 0)
        ion2 = self._create_linelist_widget(lines_layout, 1)
        ion_widgets = [ion1, ion2]
        if self.gauss_num > 2:
            ion3 = self._create_linelist_widget(lines_layout, 2)
            ion_widgets.append(ion3)
        self.line_combo.currentTextChanged.connect(lambda s, iw=ion_widgets: self._linelist_changed(s, iw))
        ion1.currentIndexChanged.connect(lambda idx, iw=ion_widgets: self._auto_populate_ions(idx, iw))

        self.line1d = LineCanvas()
        self.line1d._plot_spec(wave,flux1d, error1d)
        mpl_toolbar = NavigationToolbar(self.line1d, self)
        self.send_waves.connect(self.line1d._on_sent_waves)


        # --- possible implementation ---
        # right now it is unstable if linetools is not installed
        c_checkbox = QCheckBox('RANSAC')
        # apply ransac fitting for continuum if checked
        c_checkbox.stateChanged.connect(self._initialize_ransac)
        lines_layout.addWidget(c_checkbox, 0, 4)
        lines_layout.addWidget(QLabel('Kernel Size'), 1, 4)
        self.kernel_ransac = QLineEdit()
        self.kernel_ransac.setPlaceholderText('Kernel Size')
        self.kernel_ransac.setReadOnly(True)
        self.kernel_ransac.setFixedWidth(100)
        self.kernel_ransac.returnPressed.connect(self._fit_ransac_continuum)
        lines_layout.addWidget(self.kernel_ransac, 2, 4)
        cont_pb = QPushButton('Draw')
        cont_pb.setFixedWidth(100)
        cont_pb.clicked.connect(self._draw_button_clicked)
        lines_layout.addWidget(cont_pb, 3,4)

        # main layout
        layout.addWidget(mpl_toolbar)
        layout.addWidget(self.line1d)
        layout.addLayout(lines_layout)
        self.setLayout(layout)
        self.setFixedSize(1200,800)

    def _create_linelist_widget(self, sublayout, col):
        l_ion = QLabel('Ion {}'.format(col+1))
        l_ion.setFixedHeight(20)
        sublayout.addWidget(l_ion, 2, col)

        ion_i = QComboBox()
        ion_i.setFixedWidth(150)
        ion_i.addItem('NONE')
        ion_i.setCurrentIndex(0)
        ion_i.currentIndexChanged.connect(lambda idx, ion_widget_idx=col: self._ion_i_index_changed(idx, ion_widget_idx))
        sublayout.addWidget(ion_i, 3, col)

        return ion_i

    def _get_linelist_df(self, linelist_name):
        llist = pd.DataFrame(columns=['wave', 'name'])
        tmp = read_line_list(linelist_name)

        #need a line to append wrest to name if it doesn't have one
        if any(map(str.isdigit, tmp[1]['ion'])):
            # if name column has wrest
            for li in tmp:
                newrow = {'wave': li['wrest'], 'name': li['ion']}
                llist = llist.append(newrow, ignore_index=True)
        else:
            # if name column doesn't have wrest, need to append
            for li in tmp:
                newrow = {'wave': li['wrest'], 'name': li['ion']+' '+str(round(li['wrest']))}
                llist = llist.append(newrow, ignore_index=True)

        return llist

    def _on_sent_gauss_num(self, sent_gauss_num):
        self.gauss_num = int(sent_gauss_num)
        print(self.gauss_num)
        ion3 = self._create_linelist_widget(lines_layout, 4)
        ion_widgets.append(ion3)

    def _ion_i_index_changed(self, i, ion_widget_idx): # i is an int
        #self.send_lineindex.emit(i)
        if i > 0:
            self.waves[ion_widget_idx] = self.linelist.at[i-1, 'wave']
            self.names[ion_widget_idx] = self.linelist.at[i-1, 'name']

            #print(self.waves)
            #print(self.names)
            if sum(self.waves > 0) == self.gauss_num:        
            # now waves contain all selected ion rest wavelength
                self.send_waves.emit({ni:wi for ni,wi in zip(self.names, self.waves)})


    def _on_sent_linelists2multiG(self, l):
        self.LINELISTS = l


    def _linelist_changed(self, s, ion_widgets):
        for ion_i in ion_widgets:
            if s in 'NONE':
                self.send_linelist.emit(s)
                ion_i.clear()
                ion_i.addItem('NONE')
                ion_i.setCurrentIndex(0)
            else:
                self.linelist = self._get_linelist_df(s)
                ion_i.clear()
                ion_i.addItems(['ALL'] + self.linelist['name'].tolist())
                self.send_linelist.emit(self.linelist)
                ion_i.setCurrentIndex(0)

                #print(self.linelist)

    def _button_clicked(self, check):
        show_sigfig = 5
        print('Begin fitting multiple Gaussians')
        self.result = self.line1d.fit(self.cont)
        if self.result is not None:
            self.zf.setText(str(self.round_to_sigfig(self.result[0], show_sigfig)))
            self.zferr.setText(str(self.round_to_sigfig(self.result[1], show_sigfig)))
            self.fit_result.setText('Success!')
            self.fit_result.setStyleSheet('QLabel {color: #000000}')
        else:
            self.zf.setText('0')
            self.zferr.setText('0')
            self.fit_result.setText('Failure!')
            self.fit_result.setStyleSheet('QLabel {color: #FF0000}')
            self.result = None

    def _apply_button_clicked(self, check):
        if self.result is not None:
            self.send_gfinal.emit(self.result)

    def round_to_sigfig(self, num=0., sigfig=1):
        if np.isinf(float(num)):
            return np.inf
        else:
            return round(num, sigfig - int(np.floor(np.log10(abs(num)))) - 1)

    def _adv_button_clicked(self, check):
        print('Change parameter bounds')
        if (self.line1d.bounds is None) | (self.line1d.init_guess is None):
            print('Please Press Fit button first')
        else:
            constraint = FittingConstraintDialog(init_guess=self.line1d.init_guess,
                                                bounds=self.line1d.bounds)
            if constraint.exec_():
                (new_guess, new_bd_low, new_bd_up) = constraint._getvals()
                self.line1d.bounds = [new_bd_low, new_bd_up]
                self.line1d.init_guess = new_guess

    def _fit_ransac_continuum(self):
        self.kernel_size = int(self.kernel_ransac.text())
        if self.kernel_size % 2 == 0:
            self.kernel_size += 1
        self.c.fit_continuum(window=self.kernel_size)
        self.line1d.axline.lines.pop()
        self.line1d.axline.plot(self.wave, self.c.cont, 'b')
        self.line1d.draw()
        self.cont = self.c.cont

    def _initialize_ransac(self, s):
        # if the checkbox is checked, initialize ransac
        if s == Qt.Checked:
            self.kernel_ransac.setReadOnly(False)
            self.kernel_ransac.setText(str(self.kernel_size))
            from IGM.ransac_contfit import cont_fitter
            self.c = cont_fitter.from_data(self.wave, self.flux1d, error=self.error1d)
            self.c.fit_continuum(window=self.kernel_size)
            self.line1d.axline.plot(self.wave, self.c.cont, 'b')
            self.line1d.draw()
            self.cont = self.c.cont
        else:
            self.kernel_ransac.setReadOnly(True)
            self.kernel_ransac.clear()
            self.cont = None
            del self.line1d.axline.lines[2:]

    def _draw_button_clicked(self, check):
        if self.cont is not None:
            self.send_ransac.emit([self.wave,self.cont])

    def _export_button_clicked(self, check):
        fpath_params, fcheck = QFileDialog.getSaveFileName(None,
            'Save Multi-Gaussian Fitting Parameters',
            '',
            'Text Files (*.txt)')
        print(fpath_params)
        if fcheck:
            if self.line1d.track_fit is not None:
                with open(fpath_params, 'w') as f:
                    for line in self.line1d.track_fit:
                        line_i = line + '\n'
                        f.write(line_i)


    def _auto_populate_ions(self, i, ion_widgets):
        len_items = ion_widgets[0].count()
        if len(ion_widgets) < 3:
            # double Gaussian
            if i+1 < len_items:
                ion_widgets[-1].setCurrentIndex(i+1)
            else:
                ion_widgets[-1].setCurrentIndex(i)
        else:
            # triple Gaussian
            if i+2 > len_items:
                ion_widgets[1].setCurrentIndex(i)
                ion_widgets[2].setCurrentIndex(i)
            if i+2 == len_items:
                ion_widgets[1].setCurrentIndex(i+1)
                ion_widgets[2].setCurrentIndex(i+1)
            else:
                ion_widgets[1].setCurrentIndex(i+1)
                ion_widgets[2].setCurrentIndex(i+2)




class LineCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=3, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
        super().__init__(self.fig)
        self.fig.canvas.setFocusPolicy(Qt.ClickFocus)
        self.fig.canvas.setFocus()
        self.cid_k = self.fig.canvas.mpl_connect('key_press_event', self.ontype)

        # initialization of selected ion wavelength
        # will receive updated wavelength once all ions have been selected
        self.wavelist = []
        self.names = []
        self.init_guess = None
        self.bounds = None
        self.z_guess = 0.
        self.track_fit = None


    def _plot_spec(self, wave,flux1d,error1d):
        self.axline = self.fig.add_subplot(111)
        self.axline.cla()
        self.axline.plot(wave,flux1d,'k', alpha=0.8)
        self.axline.plot(wave,error1d,'r', alpha=0.8)
        self.g_wave, self.g_flux, self.g_error = wave, flux1d, error1d
        self.axline.set_ylim([np.min(flux1d)*0.8, np.max(flux1d)*1.3])
        self.axline.set_xlabel('Wavelength')
        self.axline.set_title('Fit Gaussians')
        
        self.draw()

    def _plot_line(self, wave_obs):
        ax = self.fig.gca()
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        self._clear_plotted_lines()

        self.z_guess = wave_obs/self.wavelist[0] - 1
        self.waves_guess = self.wavelist * (1 + self.z_guess)
        self.axline.vlines(x=self.waves_guess,
                            ymin=ylim[0], ymax=ylim[-1], color='blue', linestyle='dashed')
        for i in range(len(self.waves_guess)):
            self.axline.text(x=self.waves_guess[i],
                            y=ylim[-1]*0.6,
                            s=self.names[i],
                            color='blue', fontsize=12, rotation='vertical')
        self.axline.set_xlim(xlim)
        self.axline.set_ylim(ylim)
        self.draw()

    def _clear_plotted_lines(self):
        while self.axline.texts:
            self.axline.texts.pop()
        while self.axline.collections:
            self.axline.collections.pop()
        del self.axline.lines[2:]
        self.draw()



    def fit(self, ransac_cont=None):
        print('Start fitting multi-Gaussian profile...')
        fit_log = ['Multi-Gaussian Fitting Log', 
                    '----------------------------------']
        # mimic the single Gaussian fitting process
        # 1. fit a local continum across the entire window
        if ransac_cont is None:
            spline = splrep([self.g_wave[0], self.g_wave[-1]],
                        [self.g_flux[0], self.g_flux[-1]],
                        k=1)
            cont = splev(self.g_wave, spline)
        else:
            cont = ransac_cont

        # 2. only fit emission lines or absorption lines at once
        EW = np.sum(cont - self.g_flux)
        if EW > 0:
            # all lines are absorption lines
            sign = -1
        else:
            # all emission lines
            sign = 1

        # initialize guessed values
        if self.init_guess is None:
            sig_guess = [20] * len(self.waves_guess)

            # intialize amp's from sliced spectrum
            amp_guess = []
            for wi in self.waves_guess:
                amp_guess.append(self.g_flux[self.g_wave < wi][-1])
            #amp_guess = [3] * len(self.waves_guess)
            p_guess = [self.z_guess] + sig_guess + amp_guess
            self.init_guess = p_guess.copy()
        else:
            p_guess = self.init_guess

        # prepare ydata for fit
        ydata = sign * (self.g_flux - cont)
        errdata = sign * (self.g_error - cont)
        # start fitting process
        model_guess = MultiGauss(self.wavelist)
        if self.bounds is None:
            SoL = 299792.458 #km/s
            v_uncer = 1000 # km/s
            beta = v_uncer/SoL
            delz = np.sqrt((1.+beta)/(1.-beta)) - 1.
            self.bd_low = [self.z_guess-delz] + [0] * ((len(p_guess)-1)//2) + (np.array(amp_guess)*0.5).tolist()
            self.bd_up = [self.z_guess+delz] + [100] * ((len(p_guess)-1)//2) + (np.array(amp_guess)*1.5).tolist()
            self.bounds = [self.bd_low, self.bd_up]
        else:
            self.bd_low, self.bd_up = self.bounds[0], self.bounds[-1]

        # if we want to use errdata for estimation
        # No errdata
        try:
            #popt, pcov = curve_fit(model_guess.compile_model, self.g_wave, ydata,
            #                    p0=p_guess, bounds=(bd_low, bd_up))
            # With errdata
            popt, pcov = curve_fit(model_guess.compile_model, self.g_wave, ydata,
                                p0=p_guess, bounds=(self.bd_low, self.bd_up), sigma=errdata)

            # interpolate with the model parameters
            gfinal = sign * model_guess.compile_model(self.g_wave, *popt) + cont
            perr = np.sqrt(np.diag(pcov))
            # do not forget to bring back to observed frame


            del self.axline.lines[2:]
            cont_fit = self.axline.plot(self.g_wave, cont, 'b')
            model_fit = self.axline.plot(self.g_wave, gfinal, 'r--')

            self.draw()
            fit_log.append('\nCurrent multi-Gaussian optimal parameters:')
            fit_log.append('-----------------------------------------')
            fit_log.append(f'z = {popt[0]:.10f}, error = {perr[0]:.12f}')
            
            num_g = int(len(popt)-1) // 2
            for i in range(1, num_g+1):
                fit_log.append(f'Sigma {i} = {popt[i]:.4f}, error = {perr[i]:.4f}')
                fit_log.append(f'Amp {i} = {popt[i+num_g]:.4f}, error = {perr[i+num_g]:.4f}')
            fit_log.append('-----------------------------------------')

            for line in fit_log:
                print(line)
            self.track_fit = fit_log
            return [popt[0], perr[0]]
        except (RuntimeError, ValueError):
            print('Fitting failed...')
            self.track_fit = None
            return None

    #------------------- Keyboards/Mouse Events------------------------
    def ontype(self, event):
        '''Interactivae keyboard events
        Note:
            Always Update to help mannual for new keyboard events
        '''
        #print(event.key)
        if event.key == 'C':
            # Pick the Gaussian Center
            self.axline.plot(event.xdata,event.ydata,'r+')
            self.g1x_init = event.xdata
            self.g1y_init = event.ydata

            print('Observed wavelength for Ion 1: {:.2f}'.format(event.xdata))
            print('Fitted line wavelengths are ', event.xdata + self.delw)

            self._plot_line(event.xdata)

            self.draw()

        # enable a few keyevent for navigation
        elif event.key == 'r':
            axline = self.figure.gca()
            xlim, ylim = axline.get_xlim(), axline.get_ylim()
            self.axline.lines[0] = self.axline.plot(wave,error1d,'r', alpha=0.8)
            self.axline.lines[1] = self.axline.plot(wave,flux1d,'k', alpha=0.8)
            self.axline.set_xlim(xlim)
            self.axline.set_ylim(ylim)
        elif event.key == 't':
            ylim = self.axline.get_ylim()
            self.axline.set_ylim([ylim[0], event.ydata])
            self.draw()
        elif event.key == 'b':
            ylim = self.axline.get_ylim()
            self.axline.set_ylim([event.ydata, ylim[-1]])
            self.draw()
        elif event.key == 'X':
            xlim = self.axline.get_xlim()
            self.axline.set_xlim([xlim[0], event.xdata])
            self.draw()
        elif event.key == 'x':
            xlim = self.axline.get_xlim()
            self.axline.set_xlim([event.xdata, xlim[-1]])
            self.draw()
        elif event.key == '[':
            xlim = self.axline.get_xlim()
            delx = (xlim[-1] - xlim[0])
            self.axline.set_xlim([xlim[0] - delx, xlim[0]])
            self.draw()
        elif event.key == ']':
            xlim = self.axline.get_xlim()
            delx = (xlim[-1] - xlim[0])
            self.axline.set_xlim([xlim[1], xlim[1] + delx])
            self.draw()


    #------------------- Signals/Slots --------------------------------
    def _on_sent_waves(self, dict_waves_names):
        self.wavelist = np.array(list(dict_waves_names.values()))
        self.names = list(dict_waves_names.keys())
        self.delw = self.wavelist - self.wavelist[0]
        #print(self.wavelist)
        #print(self.names)

class MultiGauss():
    def __init__(self, linelist):
        self.linelist_rest = linelist
        #self.siglist = [20] * len(wavelist)
        #self.amplist = [5] * len(wavelist)
    def _gauss1d(self, wave_rest, sig, amp, mu_rest):
        return amp * np.exp(-(wave_rest - mu_rest)**2/(2. * sig**2))
    
    def compile_model(self, wave_obs, *params):
        # curve_fit can only take params as 1D array...
        # params includes z, sigs, amps
        # params format: z, sig1, sig2, ..., amp1, amp2, ...
        # siglist have unit Angstrom in rest frame
        # z is always saved at index 0
        z = params[0]
        # number of Gaussians needed
        num_g = int((len(params) - 1) // 2)
        # intializing sigma's and amp's
        self.siglist, self.amplist = [], []
        for i in range(1, num_g+1):
            self.siglist.append(params[i])
            self.amplist.append(params[i + num_g])

        wave_rest = wave_obs/(1. + z)

        ind_model=np.zeros((len(self.linelist_rest), len(wave_rest)))
        for i in range(len(self.linelist_rest)):
            # make sure the model is defined in rest frame
            ind_model[i]=self._gauss1d(wave_rest, self.siglist[i], self.amplist[i], self.linelist_rest[i])
        self.final_model=np.sum(ind_model,axis=0)

        return self.final_model

class FittingConstraintDialog(QDialog):
    def __init__(self, init_guess, bounds):
        super().__init__()
        self.onlyFloat = QDoubleValidator()
        self.setWindowTitle('Modify fitting parameter bounds')
        QBtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        self.buttonbox = QtWidgets.QDialogButtonBox(QBtn)
        self.layout = QGridLayout()
        self.layout.addWidget(QLabel('Names:'), 0,0)
        self.layout.addWidget(QLabel('Guess:'), 0,1)
        self.layout.addWidget(QLabel('Lower Limit:'), 0,2)
        self.layout.addWidget(QLabel('Upper Limit:'), 0,3)

        self.ion_num = (len(init_guess)-1) //2
        bd_low, bd_up = bounds[0], bounds[1]

        self.layout.addWidget(QLabel('z_guess'), 1,0)
        self.le_z = QLineEdit(str(round(init_guess[0],5)))
        self.layout.addWidget(self.le_z, 1,1)
        self.le_z_low = QLineEdit(str(round(bd_low[0],5)))
        self.layout.addWidget(self.le_z_low, 1,2)
        self.le_z_up = QLineEdit(str(round(bd_up[0],5)))
        self.layout.addWidget(self.le_z_up, 1,3)

        self.layout.addWidget(QLabel('Sigma 1 (AA)'), 2,0)
        self.le_s1 = QLineEdit(str(round(init_guess[1],5)))
        self.layout.addWidget(self.le_s1, 2,1)
        self.le_s1_low = QLineEdit(str(round(bd_low[1],5)))
        self.layout.addWidget(self.le_s1_low, 2,2)
        self.le_s1_up = QLineEdit(str(round(bd_up[1],5)))
        self.layout.addWidget(self.le_s1_up, 2,3)

        self.layout.addWidget(QLabel('Sigma 2 (AA)'), 3,0)
        self.le_s2 = QLineEdit(str(round(init_guess[2],5)))
        self.layout.addWidget(self.le_s2, 3,1)
        self.le_s2_low = QLineEdit(str(round(bd_low[2],5)))
        self.layout.addWidget(self.le_s2_low, 3,2)
        self.le_s2_up = QLineEdit(str(round(bd_up[2],5)))
        self.layout.addWidget(self.le_s2_up, 3,3)

        if self.ion_num > 2:
            self.layout.addWidget(QLabel('Sigma 3 (AA)'), 4,0)
            self.le_s3 = QLineEdit(str(round(init_guess[3],5)))
            self.layout.addWidget(self.le_s3, 4,1)
            self.le_s3_low = QLineEdit(str(round(bd_low[3],5)))
            self.layout.addWidget(self.le_s3_low, 4,2)
            self.le_s3_up = QLineEdit(str(round(bd_up[3],5)))
            self.layout.addWidget(self.le_s3_up, 4,3)

            self.layout.addWidget(QLabel('Amp 1'), 5,0)
            self.le_a1 = QLineEdit(str(round(init_guess[4],5)))
            self.layout.addWidget(self.le_a1, 5,1)
            self.le_a1_low = QLineEdit(str(round(bd_low[4],5)))
            self.layout.addWidget(self.le_a1_low, 5,2)
            self.le_a1_up = QLineEdit(str(round(bd_up[4],5)))
            self.layout.addWidget(self.le_a1_up, 5,3)

            self.layout.addWidget(QLabel('Amp 2'), 6,0)
            self.le_a2 = QLineEdit(str(round(init_guess[5],5)))
            self.layout.addWidget(self.le_a2, 6,1)
            self.le_a2_low = QLineEdit(str(round(bd_low[5],5)))
            self.layout.addWidget(self.le_a2_low, 6,2)
            self.le_a2_up = QLineEdit(str(round(bd_up[5],5)))
            self.layout.addWidget(self.le_a2_up, 6,3)

            self.layout.addWidget(QLabel('Amp 3'), 7,0)
            self.le_a3 = QLineEdit(str(round(init_guess[6],5)))
            self.layout.addWidget(self.le_a3, 7,1)
            self.le_a3_low = QLineEdit(str(round(bd_low[6],5)))
            self.layout.addWidget(self.le_a3_low, 7,2)
            self.le_a3_up = QLineEdit(str(round(bd_up[6],5)))
            self.layout.addWidget(self.le_a3_up, 7,3)
            self.layout.addWidget(self.buttonbox, 8, 1)
        else:
            self.layout.addWidget(QLabel('Amp 1'), 4,0)
            self.le_a1 = QLineEdit(str(round(init_guess[3],5)))
            self.layout.addWidget(self.le_a1, 4,1)
            self.le_a1_low = QLineEdit(str(round(bd_low[3],5)))
            self.layout.addWidget(self.le_a1_low, 4,2)
            self.le_a1_up = QLineEdit(str(round(bd_up[3],5)))
            self.layout.addWidget(self.le_a1_up, 4,3)

            self.layout.addWidget(QLabel('Amp 2'), 5,0)
            self.le_a2 = QLineEdit(str(round(init_guess[4],5)))
            self.layout.addWidget(self.le_a2, 5,1)
            self.le_a2_low = QLineEdit(str(round(bd_low[4],5)))
            self.layout.addWidget(self.le_a2_low, 5,2)
            self.le_a2_up = QLineEdit(str(round(bd_up[4],5)))
            self.layout.addWidget(self.le_a2_up, 5,3)
            self.layout.addWidget(self.buttonbox, 6, 1)

        self.setLayout(self.layout)

        self.buttonbox.accepted.connect(self.accept)
        self.buttonbox.rejected.connect(self.reject)

    def _getvals(self):
        if self.ion_num > 2:
            new_guess = [float(self.le_z.text()), 
                        float(self.le_s1.text()), float(self.le_s2.text()), float(self.le_s3.text()),
                        float(self.le_a1.text()), float(self.le_a2.text()), float(self.le_a3.text())]
            new_bd_low = [float(self.le_z_low.text()),
                        float(self.le_s1_low.text()), float(self.le_s2_low.text()), float(self.le_s3_low.text()),
                        float(self.le_a1_low.text()), float(self.le_a2_low.text()), float(self.le_a3_low.text())]
            new_bd_up = [float(self.le_z_up.text()),
                        float(self.le_s1_up.text()), float(self.le_s2_up.text()),float(self.le_s3_up.text()),
                        float(self.le_a1_up.text()), float(self.le_a2_up.text()),float(self.le_a3_up.text())]
        else:
            new_guess = [float(self.le_z.text()), 
                        float(self.le_s1.text()), float(self.le_s2.text()),
                        float(self.le_a1.text()), float(self.le_a2.text())]
            new_bd_low = [float(self.le_z_low.text()),
                        float(self.le_s1_low.text()), float(self.le_s2_low.text()),
                        float(self.le_a1_low.text()), float(self.le_a2_low.text())]
            new_bd_up = [float(self.le_z_up.text()),
                        float(self.le_s1_up.text()), float(self.le_s2_up.text()),
                        float(self.le_a1_up.text()), float(self.le_a2_up.text())]
        return [new_guess,new_bd_low,new_bd_up]