from PyQt5.QtWidgets import QWidget, QGridLayout, QVBoxLayout, QDialog, QComboBox, QLineEdit, QLabel, QPushButton
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

class Gaussfit_2d(QDialog):
    send_linelist = pyqtSignal(object)
    send_lineindex = pyqtSignal(int)
    send_waves = pyqtSignal(object)

    def __init__(self,wave, flux1d,error1d, gauss_num=2):
        super().__init__()
        self.gauss_num = gauss_num
        self.waves = np.array([-1] * gauss_num)
        self.names = np.array([None] * gauss_num)
        


        layout = QVBoxLayout()
        # sublayout for Guess z functions
        lines_layout = QGridLayout()
        lines_layout.addWidget(QLabel('Linelist'), 0,0)
        
        self.line_combo = QComboBox()
        self.line_combo.setFixedWidth(120)
        self.line_combo.addItems(['NONE', 'LBG', 'Gal', 'LLS', 'LLS Small', 'DLA', 'atom'])
        lines_layout.addWidget(self.line_combo, 1, 0)

        l_zf = QLabel('Estimated z')
        l_zf.setFixedHeight(20)
        lines_layout.addWidget(l_zf, 0, 1)
        zf = QLineEdit()
        zf.setPlaceholderText('Redshift')
        zf.setFixedWidth(100)
        zf.setReadOnly(True)
        lines_layout.addWidget(zf, 1,1)

        l_zferr = QLabel('Estimated Error')
        l_zferr.setFixedHeight(20)
        lines_layout.addWidget(l_zferr, 0, 2)
        zferr = QLineEdit()
        zferr.setPlaceholderText('Error')
        zferr.setFixedWidth(100)
        zferr.setReadOnly(True)
        lines_layout.addWidget(zferr, 1,2)

        pb = QPushButton('Fit')
        pb.setFixedWidth(100)
        pb.clicked.connect(self._button_clicked)
        lines_layout.addWidget(pb, 1,3)

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
        print('begin fitting multiple gaussians')
        self.line1d.fit()

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
        self.waves = []
        self.names = []


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

        z_guess = wave_obs/self.waves[0] - 1
        self.waves_guess = self.waves * (1 + z_guess)
        self.axline.vlines(x=self.waves_guess,
                            ymin=ylim[0], ymax=ylim[-1], color='blue', linestyle='dashed')
        for i in range(len(self.delw)):
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

    def gauss(self, x, amp, mu, sigma):
        return amp * np.exp(-(x-mu)**2/(2. * sigma**2))

    def multi_gauss(self, x, *params):
        # params should have 6 parameters for 2 gausses
        # and 9 parameters for 3 gausses
        # Parameter order should be [amp1, mu1, sigma1, amp2, mu2, sigma2, ...]

        # double gaussian
        g1 = self.gauss(x, params[0], params[1], params[2])
        g2 = self.gauss(x, params[3], params[4], params[5])

        g_final = g1 + g2
        if len(params)//3 == 3:
            # triple gaussian
            g3 = self.gauss(x, params[6], params[7], params[8])
            g_final += g3

        return g_final



    def fit(self):
        print('fitting goes here...')
        
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

        Aguess = [self.g_flux.max()]*len(self.waves_guess)
        Cguess = self.waves_guess
        sguess = [5]*len(self.waves_guess)

        # prepare ydata for fit
        ydata = sign * (self.g_flux - cont)
        errdata = sign * (self.g_error - cont)
        # start fitting
        # first ion
        popt1, pcov1 = curve_fit(self.gauss, self.g_wave, ydata,
                                p0=[Aguess[0], Cguess[0], sguess[0]],
                                sigma=errdata)
        g1_final = sign * (self.gauss(self.g_wave, *popt1)) + cont
        perr1 = np.sqrt(np.diag(pcov1))
        model1_fit = self.axline.plot(self.g_wave, g1_final, 'r--')

        # second ion
        popt2, pcov2 = curve_fit(self.gauss, self.g_wave, ydata,
                                p0=[Aguess[1], Cguess[1], sguess[1]],
                                sigma=errdata)
        g2_final = sign * (self.gauss(self.g_wave, *popt2)) + cont
        perr2 = np.sqrt(np.diag(pcov2))
        model2_fit = self.axline.plot(self.g_wave, g2_final, 'g--')

        '''
        # fit two gaussians together
        p0_f = np.array([Aguess, Cguess, sguess]).T.flatten()
        popt_f, pcov_f = curve_fit(self.multi_gauss, self.g_wave, ydata,
                                    p0=p0_f,
                                    sigma=errdata)
        g_final = sign * (self.multi_gauss(self.g_wave, *popt_f)) + cont
        perr_f = np.sqrt(np.diag(pcov_f))
        modelf_fit = self.axline.plot(self.g_wave, g_final, 'y--')
        '''
        if len(self.delw) > 2:
            # optionally, third ion
            popt3, pcov3 = curve_fit(self.gauss, self.g_wave, ydata,
                                p0=[Aguess[-1], Cguess[-1], sguess[-1]],
                                sigma=errdata)
            g3_final = sign * (self.gauss(self.g_wave, *popt3)) + cont
            perr3 = np.sqrt(np.diag(pcov3))
            model3_fit = self.axline.plot(self.g_wave, g3_final, 'b--')

        print('Ion 1 center: ', popt1[1])
        print('Ion 1 error: ', perr1[1])
        print('Ion 2 center: ', popt2[1])
        print('Ion 2 error: ', perr2[1])  

        #print('Final: ', popt_f)      

        self.draw()


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
        self.waves = np.array(list(dict_waves_names.values()))
        self.names = list(dict_waves_names.keys())
        self.delw = self.waves - self.waves[0]
        print(self.delw)
        #print(self.names)
