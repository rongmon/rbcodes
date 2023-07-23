import streamlit as st

print("main-test.py")
if 'specA' not in st.session_state:
    st.session_state['specA'] = None

for i in st.session_state:
    st.write("sessions_state: ", i,": ", st.session_state[i])

st.sidebar.write("specA: ", st.session_state.specA)
def process_upload(*args):
    ''' we'll use this to  '''
    for i in args:
        st.write("process_upload: ", i)
        st.write("process_upload - A: ", st.session_state.A)

print("process_upload: ")
def main():
    # Plot area
    st.title('FITS Data Plots')
    
    # Define sidebar
    st.sidebar.title('Upload FITS Files')

    st.session_state['specA'] = st.sidebar.file_uploader('Upload file', type=['fits'], on_change=process_upload, args=("ProcessA",), key="A") 

print("main-test.py::start")
if __name__ == "__main__":
    main()
