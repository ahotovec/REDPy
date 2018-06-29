# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2018  Alicia Hotovec-Ellis (ahotovec@gmail.com)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import argparse
import redpy
import numpy as np
import obspy
from obspy import UTCDateTime
import time
import datetime as dt
import pandas as pd
from matplotlib.dates import num2date, date2num

"""
Run this script to compare table with a specific catalog of events for agreement.
 
usage: compareCatalog.py [-h] [-v] [-c CONFIGFILE] csvfile

positional arguments:
  csvfile               catalog csv file with required column 'Time UTC'

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
"""

parser = argparse.ArgumentParser(description=
    "Compares REDPy catalog with csv catalog")
parser.add_argument("csvfile",
    help="catalog csv file with required column 'Time UTC'")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-c", "--configfile",
    help="use configuration file named CONFIGFILE instead of default settings.cfg")
args = parser.parse_args()

if args.configfile:
    opt = redpy.config.Options(args.configfile)
    if args.verbose: print("Using config file: {0}".format(args.configfile))
else:
    opt = redpy.config.Options("settings.cfg")
    if args.verbose: print("Using config file: settings.cfg")

if args.verbose: print("Opening hdf5 table: {0}".format(opt.filename))
h5file, rtable, otable, ttable, ctable, jtable, dtable, ftable = redpy.table.openTable(opt)

# Read in csv file using pandas
df = pd.read_csv(args.csvfile)

# First off, a huge assumption here is that the event times are simply within some number
# of seconds of the trigger times. Unless otherwise necessary, I'm going to start off with
# the assumption that they're within the length of time used for the correlation window.
terr = opt.winlen/opt.samprate

# I'll append the best candidate family and the number of seconds the event times differ,
# so that the user can identify questionable matches. I may add the ability to assume a
# location and do some simple ray-tracing like in checkComCat(), but at this point I'm
# not convinced it's necessary. Column dt here is the REDPy trigger time - csv catalog
# time; if it is negative, REDPy either triggered early or it may not be a match.
# Amplitude is amplitude on printed station.
df['Cluster'] = ''
df['dt'] = ''
df['Amplitude'] = ''

# Set up times to compare
evtimes = date2num(np.array(pd.to_datetime(df['Time UTC']).tolist()))
rtimes = rtable.cols.startTimeMPL[:]+rtable.cols.windowStart[:]/86400.0/opt.samprate
ramps = rtable.cols.windowAmp[:][:,opt.printsta]
otimes = otable.cols.startTimeMPL[:]+otable.cols.windowStart[:]/86400.0/opt.samprate
oamps = otable.cols.windowAmp[:][:,opt.printsta]
jtimes = date2num(np.array([dt.datetime.strptime(jtable.cols.startTime[i].decode('utf-8'),
    '%Y-%m-%dT%H:%M:%S.%f')+dt.timedelta(
    seconds=jtable.cols.windowStart[i]/opt.samprate) for i in range(len(jtable))]))
ttimes = ttable.cols.startTimeMPL[:]

# Flatten families to list
famlist = np.zeros((len(rtimes),)).astype(int)
for fnum in range(len(ftable)):
    members = np.fromstring(ftable[fnum]['members'], dtype=int, sep=' ')
    famlist[members] = fnum

for i in range(len(df)):
    
    if i%1000 == 0 and i>0:
        print('{:3.2f}% complete'.format(100.0*i/len(df)))
    
    # See if there's junk that matches
    if len(jtimes) > 0:
        dtimesj = jtimes-evtimes[i]
        bestjunk = dtimesj[np.argmin(np.abs(dtimesj))]
        if np.abs(bestjunk) < terr/86400:
            df['Cluster'][i] = 'junk'
            df['dt'][i] = bestjunk*86400
            df['Amplitude'][i] = 'NaN'
            
    # See if there are any expired orphans that match
    if len(ttimes) > 0:
        dtimest = np.array(ttimes)-evtimes[i]
        besttrig = dtimest[np.argmin(np.abs(dtimest))]
        if np.abs(besttrig) < terr/86400:
            df['Cluster'][i] = 'expired'
            df['dt'][i] = besttrig*86400
            df['Amplitude'][i] = 'NaN'
    
    # See if there's an orphan that matches
    if len(otimes) > 0:
        dtimeso = otimes-evtimes[i]
        bestorph = dtimeso[np.argmin(np.abs(dtimeso))]
        if np.abs(bestorph) < terr/86400:
            df['Cluster'][i] = 'orphan'
            df['dt'][i] = bestorph*86400
            df['Amplitude'][i] = oamps[np.argmin(np.abs(dtimeso))]
    
    # See if there's a repeater that matches
    if len(rtimes) > 0:
        dtimesr = rtimes-evtimes[i]
        bestr = dtimesr[np.argmin(np.abs(dtimesr))]
        if np.abs(bestr) < terr/86400:
            df['Cluster'][i] = famlist[np.argmin(np.abs(dtimesr))]
            df['dt'][i] = bestr*86400
            df['Amplitude'][i] = ramps[np.argmin(np.abs(dtimesr))]
        
# Write to matches.csv
if args.verbose: print("Saving to matches_{}.csv".format(opt.groupName))
df.to_csv('matches_{}.csv'.format(opt.groupName), index=False)   


if args.verbose: print("Closing table...")
h5file.close()

if args.verbose: print("Done")