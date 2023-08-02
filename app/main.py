import streamlit as st
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from astropy.io import fits
import numpy as np
import os

# a workaround/hack for st.file-uploader()
spec_tmp_filename = "data/tmp/spec.fits"

if 'spectra' not in st.session_state:
    st.session_state.spectra = {}
if 'A' not in st.session_state.spectra:
    st.session_state.spectra['A'] = {}
if 'B' not in st.session_state.spectra:
    st.session_state.spectra['B'] = {}
if 'C' not in st.session_state.spectra:
    st.session_state.spectra['C'] = {}


def make_local_copy(upfile, tmpfile):
    # st.write("flux file_details: ", upfile.name, upfile.type, upfile.size)
    if os.path.exists(tmpfile):
        os.remove(tmpfile)
        # st.write("old file removed prior to creation")
    spec = open(tmpfile, "w+b")
    spec.write(upfile.getvalue())
    spec.close()
    # st.write(f"tmp file {tmpfile} written")

def load_kcwi_spectra(filename):
    fits.info(filename)
    specfile = fits.open(filename)

    wave = specfile[2].data
    flux = specfile[0].data
    ivar = None
    return wave, flux, ivar

def load_sdss_spectra(statekey, filename):
    hdu = fits.open(filename)
    slices = hdu[1].header['NAXIS2']
    wavemin = hdu[2].data['WAVEMIN'][0]
    wavemax = hdu[2].data['WAVEMAX'][0]
    diff = wavemax - wavemin
    delta = diff / slices
    # shift the spectra a little to demonstrate...
    if statekey == 'A':
        wavemin = wavemin - diff * .40
        wavemax = wavemin + diff
    elif statekey == 'B':
        wavemin = wavemin + diff * .40
        wavemax = wavemin + diff
    # elif statekey == 'C': do nothing :D

    wave = np.arange(start=wavemin, stop=wavemax, step=delta)
    flux = hdu[1].data['flux']
    ivar = hdu[1].data['ivar']    
    return wave, flux, ivar

def load_fits_spectra1d(statekey, filename):
    # load flux and variance fits file    
    # fits.info(filename)
    # specfile = fits.open(filename)
    # print(f"loading {filename}   statekey: {statekey}")
    wave, flux, _ = load_sdss_spectra(statekey, filename)

    st.session_state.spectra[statekey]["name"] = st.session_state[statekey].name
    st.session_state.spectra[statekey]["wave"] = wave
    st.session_state.spectra[statekey]["flux"] = flux
    # hoping we won't need the flux header (eg. for ["UNITS"])
    # specfile.close()

def harmonize_data(statekey):
    if st.session_state.spectra[statekey] is not None:
        wdata = st.session_state.spectra[statekey]["wave"]
        fdata = st.session_state.spectra[statekey]["flux"]

        fdata[fdata < 0] = 0
        wmin = np.min(wdata)
        wmax = np.max(wdata)
        wrange = wmax - wmin
        print(f"wdata min: {wmin}")
        print(f"wdata max: {wmax}")
        print(f"wdata range: {wrange}")
        fmin = np.min(fdata)
        fmax = np.max(fdata)
        print(f"fdata min: {fmin}")
        print(f"fdata max: {fmax}")
        print(f"fdata range: {fmax-fmin}")



def process_upload(*args):
    ''' this is a callback function for st.file_uploader() '''

    if args[0] == "ProcessA":
        statekey = 'A'
    elif args[0] == "ProcessB":
        statekey = 'B'
    elif args[0] == "ProcessC":
        statekey = 'C'
    else:
        return
    
    if st.session_state[statekey] is None:
        print(f"clearing out session_state.spectra[{statekey}]")
        st.session_state.spectra[statekey] = {}

    if statekey in st.session_state and st.session_state[statekey] is not None:
        # print("processing")
        make_local_copy(st.session_state[statekey], spec_tmp_filename)
        load_fits_spectra1d(statekey, spec_tmp_filename)
        # harmonize_data(statekey)


def main():
    # Plot area
    st.title('FITS Data Plots')
    
    # Define sidebar
    st.sidebar.title('Upload FITS Files')

    st.sidebar.file_uploader('Upload 1D Spectra FITS file for plot A', type=['fits'], on_change=process_upload, args=("ProcessA",), key="A") 
    st.sidebar.file_uploader('Upload 1D Spectra FITS file for plot B', type=['fits'], on_change=process_upload, args=("ProcessB",), key="B") 
    st.sidebar.file_uploader('Upload 1D Spectra FITS file for plot C', type=['fits'], on_change=process_upload, args=("ProcessC",), key="C") 
    
    # for i in sorted(st.session_state):
    #     st.write("session_state: ", i,": ", st.session_state[i])

    # Create subplot
    # graphics_mode = "pyplot"
    # graphics_mode = "plotly"
    graphics_mode = st.selectbox('Graphics mode',('plotly', 'pyplot'))
    sharedx = st.checkbox("Shared x-axis", value=False,  label_visibility="visible")
    if graphics_mode == "pyplot":
        drawme = False
        fig, axs = plt.subplots(3, 1, sharex=sharedx)
        # {'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan'} which are the Tableau Colors from the 'tab10' categorical palette (which is the default color cycle);


        if "flux" in st.session_state.spectra['A']:
            axs[0].plot(st.session_state.spectra['A']["wave"], 
                        st.session_state.spectra['A']["flux"], 
                        color='tab:blue', 
                        linewidth=1.0)
            drawme = True
        if "flux" in st.session_state.spectra['B']:
            axs[1].plot(st.session_state.spectra['B']["wave"], 
                        st.session_state.spectra['B']["flux"], 
                        color='tab:cyan', 
                        linewidth=1.0)
            drawme = True            
        if "flux" in st.session_state.spectra['C']:
            axs[2].plot(st.session_state.spectra['C']["wave"], 
                        st.session_state.spectra['C']["flux"], 
                        color='tab:red', 
                        linewidth=1.0)
            drawme = True            

        if drawme: st.pyplot(fig, use_container_width=True)

    if graphics_mode == "plotly":

        fig = make_subplots(rows=3, cols=1, shared_xaxes=sharedx)
        if "flux" in st.session_state.spectra['A']:
            fig.add_trace(go.Scatter(x=st.session_state.spectra['A']["wave"], 
                                     y=st.session_state.spectra['A']["flux"], 
                                     mode='lines', 
                                     name=st.session_state['A'].name), row=1, col=1)

        if "flux" in st.session_state.spectra['B']:
            fig.add_trace(go.Scatter(x=st.session_state.spectra['B']["wave"], 
                                     y=st.session_state.spectra['B']["flux"], 
                                     mode='lines', 
                                     name=st.session_state['B'].name), row=2, col=1)
            
        if "flux" in st.session_state.spectra['C']:
            fig.add_trace(go.Scatter(x=st.session_state.spectra['C']["wave"], 
                                     y=st.session_state.spectra['C']["flux"], 
                                     mode='lines', 
                                     name=st.session_state['C'].name), row=3, col=1)
            

        fig.update_layout(height=600, width=800, title_text="SDSS Spectra")
        st.plotly_chart(fig, use_container_width=True)

if __name__ == "__main__":
    main()
