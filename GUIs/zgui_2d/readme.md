# GUI Development specifically for 2D

Integrate with Plotly Dash or Bokeh in PyQt5.

## Installation
Besides the packages in `requirements_simple.txt`, the following packages are required:
* `plotly`: `conda install -c plotly plotly=5.17`
* `bokeh` : `conda install bokeh=2.4.3`

## `Plotly` demo
For demo purpose, the `plotly_example.py` worked well. It uses `spec2d_coadd_QSO_J0100_sID010242.fits` as an example. To run the demo, simply type
```bash
python <path to plotly_example.py>/plotly_example.py
```

## `rbcodes` installation
* Window users: Click Control Panel > System and Security > System > Advanced System Settings > Environment Variables > System variables > Path > Edit > New > Paste the path of `rbcodes` folder > OK > OK > OK.