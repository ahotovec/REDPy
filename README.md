## Overview
REDPy (Repeating Earthquake Detector in Python) is a tool for automated detection and analysis of repeating earthquakes in continuous data. It works without any previous assumptions of what repeating seismicity looks like (that is, does not require a template event). Repeating earthquakes are clustered into "families" via the OPTICS ([Ordering Points To Identify the Clustering Structure](https://en.wikipedia.org/wiki/OPTICS_algorithm)) algorithm, with distance defined by cross-correlation. All data, including waveforms, are stored in an HDF5 table using [PyTables](http://www.pytables.org/).

## Motivation
Monitoring earthquakes is a challenge because signals from hundreds of seismograms record up to several hundred samples a second. Automatic event detection is required, yet important signals are too buried in noise and variable in character to identify with simple schemes. We aim to systematically mine the numerous events that are near-repetitions.

If two or more earthquakes occur in the same location and the same source, they will have highly similar waveforms. These "repeating" earthquakes are a common occurrence in glacial, tectonic, and volcanic areas worldwide, and are a useful indicator of seismic activity. Repeating earthquakes can be used to track fluid movement (e.g., water, magma), slip on faults, volcanic activity, or subtle changes in structure, and therefore are a topic of active research in the seismological community. In practice, repeating earthquakes are occasionally noticed visually hours after they occur, but usually go undetected by routine seismic monitoring.

We seek to develop an automated system for detecting and cataloging repeating earthquakes at any subset of seismic stations, which could be used in near real-time for monitoring or on archived data for research. 

## Installation
REDPy runs on Python 2.7, and currently has the following major package dependencies:  
[numpy](http://www.numpy.org/) | [scipy](http://www.scipy.org/) | [matplotlib](http://www.matplotlib.org/) | [obspy](http://www.obspy.org/) | [pytables](http://www.pytables.org/) | [pandas](http://pandas.pydata.org/)

All of these dependencies can be easily installed via [Anaconda](https://www.continuum.io/) on the command line:
`>> conda install -c obspy obspy`
`>> conda install pytables`
`>> conda install pandas`

However, due to a 'feature' in numpy 1.10.0, processing runs _extremely_ slow (by at least an order of magnitude, if not two or more). Until a fix is released for this issue in numpy, we recommend downgrading numpy, scipy, matplotlib, and pytables to the following previously stable versions:
`>> conda install numpy=1.9.3 scipy=0.15.0 matplotlib=1.4.2 pytables=3.2.0`

## Usage
Once dependencies are installed and REDPy is downloaded, REDPy can be run out of the box with the following commands to test if the code is working on your computer:  
`>> python initialize.py`  
`>> python backfill.py -s 2004-09-22 -e 2004-09-24`

This will download two days of data during the beginning of the 2004 eruption of Mount St. Helens, where several repeating earthquakes are sure to be found. `initialize.py` sets up the hdf5 pytable where data are stored, and creates all the necessary folders. `backfill.py` fills the table with data between a given start and end date, if provided. If no dates are provided, it will attempt to fill up to the current date and assumes a starting time based on whether or not data exist in the table. If data exist, it will pick up from the time of the last trigger in the table, if not it will download only the last N seconds (defined in configuration file). The table can also be populated from a CSV catalog using `catfill.py`; an example catalog `mshcat.csv` is included that will work with the default configuration.

REDPy is configured using a .cfg file. The default file is `settings.cfg`, but the `-c` flag can be used to specify a different .cfg file.

Images are output to a folder with the same name as the 'group name' in the configuration file.

Typing `-h` after each script will bring up a help usage with options.


## Development
This code is still under major development and not ready for release.

Please contact Alicia Hotovec-Ellis (University of Washington) to be notified on release of this package. Persons interested in collaborating to enhance the functionality of the code are also welcome to contact Alicia.
