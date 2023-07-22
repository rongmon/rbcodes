import streamlit as st
from matplotlib import pyplot as plt
# import plotly.figure_factory as ff
# import plotly.express as px
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
    spec = open(tmpfile, "w+b")
    spec.write(upfile.getvalue())
    spec.close()

def load_fits_spectra1d(statekey, filename):
    # load flux and variance fits file    
    fits.info(filename)
    specfile = fits.open(filename)
    # FLUX, ERROR, WAVELENGTH
    st.session_state.spectra[statekey]["flux"] = specfile[0].data
    st.session_state.spectra[statekey]["wave"] = specfile[2].data
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
    ''' we'll use this to  '''

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
        print("processing")
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
    graphics_mode = "plotly"
    if graphics_mode == "pyplot":
        drawme = False
        fig, axs = plt.subplots(3, 1, sharex=True)
        if "flux" in st.session_state.spectra['A']:
            axs[0].plot(st.session_state.spectra['A']["wave"], st.session_state.spectra['A']["flux"], color='black', linewidth=1.0)
            drawme = True
        if "flux" in st.session_state.spectra['B']:
            axs[1].plot(st.session_state.spectra['B']["wave"], st.session_state.spectra['B']["flux"], color='black', linewidth=1.0)
            drawme = True            
        if "flux" in st.session_state.spectra['C']:
            axs[2].plot(st.session_state.spectra['C']["wave"], st.session_state.spectra['C']["flux"], color='black', linewidth=1.0)
            drawme = True            

        if drawme: st.pyplot(fig, use_container_width=True)

    if graphics_mode == "plotly":
        drawme = False
        fig = make_subplots(rows=3, cols=1, shared_xaxes=True)
        if "flux" in st.session_state.spectra['A']:
            fig.add_trace(go.Scatter(x=st.session_state.spectra['A']["wave"], y=st.session_state.spectra['A']["flux"], mode='lines', name='A'), row=1, col=1)
            drawme = True            
        if "flux" in st.session_state.spectra['B']:
            fig.add_trace(go.Scatter(x=st.session_state.spectra['B']["wave"], y=st.session_state.spectra['B']["flux"], mode='lines', name='B'), row=2, col=1)
            drawme = True            
        if "flux" in st.session_state.spectra['C']:
            fig.add_trace(go.Scatter(x=st.session_state.spectra['C']["wave"], y=st.session_state.spectra['C']["flux"], mode='lines', name='C'), row=3, col=1)
            drawme = True            

        fig.update_layout(height=600, width=800, title_text="Spectra")
        # if drawme: 
        st.plotly_chart(fig, use_container_width=True)

if __name__ == "__main__":
    main()
