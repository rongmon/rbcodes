import PySimpleGUI as sg
import pandas as pd
#sg.theme('Topanga')      # Add some color to the window
sg.ChangeLookAndFeel('Dark')      

# Very basic window.  Return values using auto numbered keys
d={'zabs':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],'List':['None','None','None','None','None','None','None','None','None','None'], 'color':['k','k','k','k','k','k','k','k','k','k']} 
df=pd.DataFrame(data=d)    

col1=[ [sg.Text('1. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('2. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('3. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('4. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('5. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('6. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('7. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('8. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('9. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))],
       [sg.Text('10. zabs', size=(5, 1)), sg.In(default_text='0.0',  size=(5, 1))]]

col2= [[sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)],
       [sg.Text('LineList', size=(5, 1)),sg.Drop(values=('LLS', 'LLS Small', 'DLA'), auto_size_text=True)]]



col3=  [[sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))],
        [sg.Text('color', size=(5, 1)), sg.In(default_text='None' ,size=(5, 1))]]
    



layout = [[sg.Column(col1),sg.Column(col2),sg.Column(col3)], [sg.Submit(), sg.Exit()]]


window = sg.Window('Update Selected Absorbers', layout, font=("Helvetica", 12))

while True:
    event, values = window.read()
    #update database
    for i in range(0,10):
        df.at[i, 'zabs'] = values[i]
        df.at[i, 'List'] = values[i+10]
        df.at[i, 'color'] = values[i+20]   
    print(values)    # the input data looks like a simple list when auto numbered
    print(df)
    if event == sg.WIN_CLOSED or event == 'Exit': 
        break 
window.close()
