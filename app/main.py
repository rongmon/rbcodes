import streamlit as st
from matplotlib import pyplot as plt
from data_processing.fits_processing import process_fits_A, process_fits_B, process_fits_C

def main():
    # Define sidebar
    st.sidebar.title('Upload FITS Files')

    file_A = st.sidebar.file_uploader('Upload FITS file for plot A', type=['fits'])
    file_B = st.sidebar.file_uploader('Upload FITS file for plot B', type=['fits'])
    file_C = st.sidebar.file_uploader('Upload FITS file for plot C', type=['fits'])

    # Plot area
    st.title('FITS Data Plots')

    # Create subplot
    fig, axs = plt.subplots(3, 1, sharex=True)

    if file_A is not None:
        fig_A = process_fits_A(file_A)
        axs[0].imshow(fig_A)

    if file_B is not None:
        fig_B = process_fits_B(file_B)
        axs[1].imshow(fig_B)

    if file_C is not None:
        fig_C = process_fits_C(file_C)
        axs[2].imshow(fig_C)

    st.pyplot(fig)

if __name__ == "__main__":
    main()
