<img src="https://raw.githubusercontent.com/ahotovec/REDPy/master/img/logo.png" width=800 alt="REDPy Logo" />

## Overview
REDPy (Repeating Earthquake Detector in Python) is a tool for automated detection and analysis of repeating earthquakes in continuous data. It works without any previous assumptions of what repeating seismicity looks like (that is, does not require a template event). Repeating earthquakes are clustered into "families" via the OPTICS ([Ordering Points To Identify the Clustering Structure](https://en.wikipedia.org/wiki/OPTICS_algorithm)) algorithm, with distance defined by cross-correlation. All data, including waveforms, are stored in an HDF5 table using [PyTables](http://www.pytables.org/).

## Installation
REDPy runs on Python 2.7, and currently has the following major package dependencies:  
[numpy](http://www.numpy.org/) | [scipy](http://www.scipy.org/) | [matplotlib with basemap add-on](http://www.matplotlib.org/) | [obspy](http://www.obspy.org/) | [pytables](http://www.pytables.org/) | [pandas](http://pandas.pydata.org/) | [bokeh](http://bokeh.pydata.org/) 

These dependencies can be easily installed via [Anaconda](https://www.continuum.io/) on the command line:
```
>> conda install -c obspy obspy pytables pandas basemap bokeh=0.9.3 mock nose PIL
```

## Usage
Once dependencies are installed and REDPy is downloaded, REDPy can be run out of the box with the following commands to test if the code is working on your computer:  
```
>> python initialize.py
>> python catfill.py mshcat.csv
>> python backfill.py -s 2004-09-15 -e 2004-09-24
```

Check out the [Wiki](https://github.com/ahotovec/REDPy/wiki) for more detailed usage!
