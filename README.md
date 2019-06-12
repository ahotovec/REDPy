<img src="https://raw.githubusercontent.com/ahotovec/REDPy/master/img/logo.png" width=800 alt="REDPy Logo" />

## Overview
REDPy (Repeating Earthquake Detector in Python) is a tool for automated detection and analysis of repeating earthquakes in continuous data. It works without any previous assumptions of what repeating seismicity looks like (that is, does not require a template event). Repeating earthquakes are clustered into "families" based on cross-correlation across multiple stations. All data, including waveforms, are stored in an HDF5 table using [PyTables](http://www.pytables.org/).

## Installation
REDPy runs on Python 2.7 and Python 3.5, with the following major package dependencies:  
[numpy](http://www.numpy.org/) | [scipy](http://www.scipy.org/) | [matplotlib](http://www.matplotlib.org/) | [obspy](http://www.obspy.org/) | [pytables](http://www.pytables.org/) | [pandas](http://pandas.pydata.org/) | [bokeh](http://bokeh.pydata.org/) | [cartopy](http://scitools.org.uk/cartopy/)

These dependencies can be easily installed via [Anaconda](https://www.continuum.io/) on the command line. I *highly* recommend using a virtual environment so that your REDPy environment does not conflict with any other Python packages you may be using. This can be done with the following commands:
```
>> conda create -n redpy python=3.5
>> conda activate redpy
>> conda install -c obspy obspy
>> conda install -c conda-forge bokeh cartopy
>> conda install pytables pandas
```
You may either use `python=3.5` or `python=2.7`, but other versions are not supported. Whenever you intend to run REDPy, be sure to `conda activate redpy` and then `conda deactivate` when you are done.

There are also environment files that should create a stable environment with either version of Python; to use `python=3.5`:
```
>>> conda env create -f redpy35.yml
```
or `python=2.7`:
```
>>> conda env create -f redpy27.yml
```

## Usage
Once dependencies are installed and REDPy is downloaded, REDPy can be run out of the box with the following commands to test if the code is working on your computer. If it completes without error, it will produce files in a folder named `default` after several minutes.
```
>> python initialize.py
>> python catfill.py -v mshcat.csv
>> python backfill.py -v -s 2004-09-15 -e 2004-09-24
```

Check out the [Wiki](https://github.com/ahotovec/REDPy/wiki) for more detailed usage!

## Reference

If you would like to reference REDPy in your paper, please cite the following abstract until I finish writing the *Electronic Seismologist* paper for it:

> Hotovec-Ellis, A.J., and Jeffries, C., 2016. Near Real-time Detection, Clustering, and Analysis of Repeating Earthquakes: Application to Mount St. Helens and Redoubt Volcanoes – *Invited*, presented at Seismological Society of America Annual Meeting, Reno, Nevada, 20 Apr.
