<img src="https://raw.githubusercontent.com/ahotovec/REDPy/master/img/logo.png" width=800 alt="REDPy Logo" />

## Overview
REDPy (Repeating Earthquake Detector in Python) is a tool for automated detection and analysis of repeating earthquakes in continuous data. It works without any previous assumptions of what repeating seismicity looks like (that is, does not require a template event). Repeating earthquakes are clustered into "families" based on cross-correlation across multiple stations. All data, including waveforms, are stored in an HDF5 table using [PyTables](http://www.pytables.org/).

## Installation
Download the [zip file](https://github.com/ahotovec/REDPy/archive/master.zip) or use `git` to clone the entire repository to a working directory (e.g., mine is `/Users/ahotovecellis/REDPy/`). All scripts will be run from this directory, and all new files will be generated here.

REDPy runs on Python 3.6+, with the following major package dependencies:  
[numpy](http://www.numpy.org/) | [scipy](http://www.scipy.org/) | [matplotlib](http://www.matplotlib.org/) | [obspy](http://www.obspy.org/) | [pytables](http://www.pytables.org/) | [pandas](http://pandas.pydata.org/) | [bokeh](http://bokeh.pydata.org/) | [cartopy](http://scitools.org.uk/cartopy/)

These dependencies can be easily installed via [Anaconda](https://www.continuum.io/) on the command line. I *highly* recommend using a virtual environment so that your REDPy environment does not conflict with any other Python packages you may be using. This can be done with the following commands:
```
>> conda config --add channels conda-forge
>> conda create -n redpy python=3.7 bokeh cartopy obspy pandas pytables
```
You may either use `python=3.6` or `python=3.7`, but other versions are no longer supported.

I have included an environment file that should create a stable environment with Python `python=3.7`:
```
>>> conda env create -f redpy37.yml
```
or at least provide a guide for if the code breaks and the output of your `conda list` does not match.

Whenever you intend to run REDPy, be sure to `conda activate redpy` and then `conda deactivate` when you are done.


## Usage
Once dependencies are installed, REDPy is downloaded, and you are in the `redpy` environment, REDPy can be run out of the box with the following commands to test if the code is working on your computer. If it completes without error, it will produce files in a folder named `default` after several minutes.
```
>> python initialize.py
>> python catfill.py -v mshcat.csv
>> python backfill.py -v -s 2004-09-15 -e 2004-09-24
```

Check out the [Wiki](https://github.com/ahotovec/REDPy/wiki) for more detailed usage!

## Reference

If you would like to reference REDPy in your paper, please cite the following abstract until I finish writing the *Electronic Seismologist* paper for it:

> Hotovec-Ellis, A.J., and Jeffries, C., 2016. Near Real-time Detection, Clustering, and Analysis of Repeating Earthquakes: Application to Mount St. Helens and Redoubt Volcanoes – *Invited*, presented at Seismological Society of America Annual Meeting, Reno, Nevada, 20 Apr.
