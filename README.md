<img src="https://raw.githubusercontent.com/ahotovec/REDPyAlpha/master/img/logo.png" width=800 alt="REDPy Logo" />

## Overview
REDPy (Repeating Earthquake Detector in Python) is a tool for automated detection and analysis of repeating earthquakes in continuous data. It works without any previous assumptions of what repeating seismicity looks like (that is, does not require a template event). Repeating earthquakes are clustered into "families" via the OPTICS ([Ordering Points To Identify the Clustering Structure](https://en.wikipedia.org/wiki/OPTICS_algorithm)) algorithm, with distance defined by cross-correlation. All data, including waveforms, are stored in an HDF5 table using [PyTables](http://www.pytables.org/).

## Installation
REDPy runs on Python 2.7, and currently has the following major package dependencies:  
[numpy](http://www.numpy.org/) | [scipy](http://www.scipy.org/) | [matplotlib](http://www.matplotlib.org/) | [obspy](http://www.obspy.org/) | [pytables](http://www.pytables.org/) | [pandas](http://pandas.pydata.org/) | [bokeh](http://bokeh.pydata.org/) 

Most of these dependencies can be easily installed via [Anaconda](https://www.continuum.io/) on the command line:  
`>> conda install -c obspy obspy`  
`>> conda install pytables`  
`>> conda install pandas`  
`>> conda install bokeh`

However, due to a 'feature' in numpy 1.10.0, processing runs _extremely_ slow (by at least an order of magnitude, if not two or more). Until a fix is released for this issue in numpy, we recommend downgrading numpy, scipy, matplotlib, and pytables to the following previously stable versions:  
`>> conda install numpy=1.9.3 scipy=0.15.0 matplotlib=1.4.2 pytables=3.2.0 pandas=0.16.0 bokeh=0.9.3`

## Usage
Once dependencies are installed and REDPy is downloaded, REDPy can be run out of the box with the following commands to test if the code is working on your computer:  
`>> python initialize.py`  
`>> python catfill.py mshcat.csv`  
`>> python backfill.py -s 2004-09-15 -e 2004-09-24`

Check out the Wiki for more detailed usage!
