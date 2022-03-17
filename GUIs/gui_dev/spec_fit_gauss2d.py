from PyQt5.QtWidgets import QWidget, QGridLayout, QVBoxLayout, QDialog, QComboBox, QLineEdit, QLabel
from PyQt5.QtCore import Qt, pyqtSignal

import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import pandas as pd

from IGM.rb_setline import read_line_list

class Gaussfit_2d(QDialog):
    send_linelist = pyqtSignal(object)
    send_lineindex = pyqtSignal(int)
    send_waves = pyqtSignal(object)

    def __init__(self,wave, flux1d,error1d, gauss_num=2):
        super().__init__()
        self.gauss_num = gauss_num
        self.waves = np.array([-1] * gauss_num)
        


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

        ion1 = self._create_linelist_widget(lines_layout, 0)
        ion2 = self._create_linelist_widget(lines_layout, 1)
        ion_widgets = [ion1, ion2]
        if self.gauss_num > 2:
            ion3 = self._create_linelist_widget(lines_layout, 2)
            ion_widgets.append(ion3)
        self.line_combo.currentTextChanged.connect(lambda s, iw=ion_widgets: self._linelist_changed(s, iw))

        line1d = LineCanvas()
        line1d._plot_spec(wave,flux1d, error1d)
        mpl_toolbar = NavigationToolbar(line1d, self)
        self.send_waves.connect(line1d._on_sent_waves)

        # main layout
        layout.addWidget(mpl_toolbar)
        layout.addWidget(line1d)
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

    def gauss(self, x, amp, mu, sigma):
        return amp * np.exp(-(x-mu)**2/(2. * sigma**2))

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

    def _ion_i_index_changed(self, i, ion_widget_idx,): # i is an int
        #self.send_lineindex.emit(i)
        if i > 0:
            self.waves[ion_widget_idx] = self.linelist.at[i-1, 'wave']
            print(self.waves)
            if sum(self.waves > 0) == self.gauss_num:        
            # now waves contain all selected ion rest wavelength
                self.send_waves.emit(self.waves)




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


    def _plot_spec(self, wave,flux1d,error1d):
        self.axline = self.fig.add_subplot(111)
        self.axline.cla()
        self.axline.plot(wave,flux1d,'k')
        self.axline.plot(wave,error1d,'r')


        self.axline.set_xlabel('Wavelength')
        self.axline.set_title('Fit Gaussians')
        
        self.draw()

    def _plot_line(self, ionindex, estZ=0.):
        ax = self.fig.gca()
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        if ionindex > 0:
            i_ion = ionindex - 1
            self.axline.vlines(x=self.linelist.at[i_ion, 'wave'] * (1+estZ),
                                ymin=ylim[0], ymax=ylim[-1], color='blue', linestyle='dashed')
            self.axline.text(x=self.linelist.at[i_ion, 'wave'] * (1+estZ),
                            y=ylim[-1]*0.6,
                            s=self.linelist.at[i_ion, 'name'],
                            color='blue', fontsize=12, rotation='vertical')
            self.axline.set_xlim(xlim)
            self.axline.set_ylim(ylim)
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
            self.draw()

            print('Observed wavelength for Ion 1: {:.2f}'.format(event.xdata))

    #------------------- Signals/Slots --------------------------------
    def _on_sent_waves(self, sent_waves):
        self.waves = sent_waves
        self.delw = self.waves - self.waves[0]
        print(self.delw)
