from PyQt5.QtWidgets import QWidget, QGridLayout, QVBoxLayout, QComboBox, QLineEdit, QLabel, QPushButton
from PyQt5.QtCore import Qt, pyqtSignal

import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit

from IGM.rb_setline import read_line_list

class Gaussfit_2d(QWidget):
    send_linelist = pyqtSignal(object)
    send_lineindex = pyqtSignal(int)
    send_waves = pyqtSignal(object)
    send_gfinal = pyqtSignal(list)

    def __init__(self,wave, flux1d,error1d, gauss_num=2, linelists=[]):
        super().__init__()
        self.gauss_num = gauss_num
        self.waves = np.array([-1] * gauss_num)
        self.names = np.array([None] * gauss_num)
        self.linelists = linelists
        


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

        ion1 = self._create_linelist_widget(lines_layout, 0)
        ion2 = self._create_linelist_widget(lines_layout, 1)
        ion_widgets = [ion1, ion2]
        if self.gauss_num > 2:
            ion3 = self._create_linelist_widget(lines_layout, 2)
            ion_widgets.append(ion3)
        self.line_combo.currentTextChanged.connect(lambda s, iw=ion_widgets: self._linelist_changed(s, iw))

        self.line1d = LineCanvas()
        self.line1d._plot_spec(wave,flux1d, error1d)
        mpl_toolbar = NavigationToolbar(self.line1d, self)
        self.send_waves.connect(self.line1d._on_sent_waves)

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
        print('begin fitting multiple gaussians')
        result = self.line1d.fit()
        if result is not None:
            self.zf.setText(str(self.round_to_sigfig(result[0], show_sigfig)))
            self.zferr.setText(str(self.round_to_sigfig(result[1], show_sigfig)))
            self.send_gfinal.emit(result)
            self.fit_result.setText('Success!')
            self.fit_result.setStyleSheet('QLabel {color: #000000}')
        else:
            self.zf.setText('0')
            self.zferr.setText('0')
            self.fit_result.setText('Failure!')
            self.fit_result.setStyleSheet('QLabel {color: #FF0000}')

    def round_to_sigfig(self, num=0., sigfig=1):
        return round(num, sigfig - int(np.floor(np.log10(abs(num)))) - 1)

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
        self.z_guess = 0.


    def _plot_spec(self, wave,flux1d,error1d):
        self.axline = self.fig.add_subplot(111)
        self.axline.cla()
        self.axline.plot(wave,flux1d,'k')
        self.axline.plot(wave,error1d,'r')
        self.g_wave, self.g_flux, self.g_error = wave, flux1d, error1d


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



    def fit(self):
        print('Start fitting multi-Gaussian profile...')
        
        # mimic the single Gaussian fitting process
        # 1. fit a local continum across the entire window
        spline = splrep([self.g_wave[0], self.g_wave[-1]],
                        [self.g_flux[0], self.g_flux[-1]],
                        k=1)
        cont = splev(self.g_wave, spline)

        # 2. only fit emission lines or absorption lines at once
        EW = np.sum(cont - self.g_flux)
        if EW > 0:
            # all lines are absorption lines
            sign = -1
        else:
            # all emission lines
            sign = 1

        # initialize guessed values
        sig_guess = [5] * len(self.waves_guess)
        amp_guess = [3] * len(self.waves_guess)
        p_guess = [self.z_guess] + sig_guess + amp_guess

        # prepare ydata for fit
        ydata = sign * (self.g_flux - cont)
        errdata = sign * (self.g_error - cont)
        # start fitting process
        model_guess = MultiGauss(self.wavelist)
        bd_low = [self.z_guess*0.84] + [0] * (len(p_guess)-1)
        bd_up = [self.z_guess*1.16] + [np.inf] * (len(p_guess)-1)

        # if we want to use errdata for estimation
        # No errdata
        try:
            #popt, pcov = curve_fit(model_guess.compile_model, self.g_wave, ydata,
            #                    p0=p_guess, bounds=(bd_low, bd_up))
            # With errdata
            popt, pcov = curve_fit(model_guess.compile_model, self.g_wave, ydata,
                                p0=p_guess, bounds=(bd_low, bd_up), sigma=errdata)

            gfinal = sign * model_guess.compile_model(self.g_wave, *popt)
            perr = np.sqrt(np.diag(pcov))
            model_fit = self.axline.plot(self.g_wave, gfinal, 'r--')

            self.draw()
            return [popt[0], perr[0]]
        except RuntimeError:
            print('Fitting failed...')
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

    #------------------- Signals/Slots --------------------------------
    def _on_sent_waves(self, dict_waves_names):
        self.wavelist = np.array(list(dict_waves_names.values()))
        self.names = list(dict_waves_names.keys())
        self.delw = self.wavelist - self.wavelist[0]
        print(self.wavelist)
        print(self.names)

class MultiGauss():
    def __init__(self, wavelist):
        self.wave_rest=wavelist
        #self.siglist = [20] * len(wavelist)
        #self.amplist = [5] * len(wavelist)
    def _gauss1d(self,x,sig,amp,mu):
        return amp * np.exp(-(x-mu)**2/(2. * sig**2))
    
    def compile_model(self, wave, *params):
        # curve_fit can only take params as 1D array...
        # params includes z, sigs, amps
        # params format: z, sig1, sig2, ..., amp1, amp2, ...
        z = params[0]
        num_g = int((len(params) - 1) // 2)
        self.siglist, self.amplist = [], []
        self.wave_obs=self.wave_rest*(1. + z)
        for i in range(1, num_g+1):
            self.siglist.append(params[i])
            self.amplist.append(params[i + num_g])

        ind_model=np.zeros((len(self.wave_rest),len(wave)))
        for i in range(len(self.wave_rest)):
            ind_model[i]=self._gauss1d(wave,self.siglist[i],self.amplist[i],self.wave_obs[i])
        self.final_model=np.sum(ind_model,axis=0)

        return self.final_model

#class FittingConstraint():
