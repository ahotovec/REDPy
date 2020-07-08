# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import redpy.config
import redpy.table
import argparse
import shutil
import os
import numpy as np

"""
Run this script to create space for additional stations while preserving data in an
existing table or to change the directory name for a run. Additional stations should
always be included at the end of the station list; reordering that list is currently not
supported. Running this script will overwrite any existing table with the same name
defined by filename in the new .cfg file. If the table names in both .cfg files are the
same, the original table will be renamed and then deleted. All output files are also
remade to reflect the additional station, unless flagged otherwise.

usage: extendTable.py [-h] [-v] [-n] CONFIGFILE_FROM CONFIGFILE_TO

positional arguments:
  CONFIGFILE_FROM       old .cfg file corresponding to table to be copied from
  CONFIGFILE_TO         new .cfg file corresponding to table to be copied to

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -n, --noplot          do not re-render plots after extending
"""

parser = argparse.ArgumentParser(description=
    "Create space for additional stations based on an existing table")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-n", "--noplot", action="count", default=0,
    help="do not re-render plots after extending")
parser.add_argument('cfgfrom', metavar='CONFIGFILE_FROM', type=str, nargs=1,
    help="old .cfg file corresponding to table to be copied from")
parser.add_argument('cfgto', metavar='CONFIGFILE_TO', type=str, nargs=1,
    help="new .cfg file corresponding to table to be copied to")

args = parser.parse_args()


if args.verbose: print("Using old config file: {0}".format(args.cfgfrom[0]))
optfrom = redpy.config.Options(args.cfgfrom)

if args.verbose: print("Using new config file: {0}".format(args.cfgto[0]))
optto = redpy.config.Options(args.cfgto)

if args.verbose: print("Making working copy of old hdf5 table...")
shutil.copy(optfrom.filename,'{}.old'.format(optfrom.filename))

# Change filename in optfrom to point to the .old version
optfrom.filename = '{}.old'.format(optfrom.filename)
    
if args.verbose: print("Creating empty hdf5 table: {0}".format(optto.filename))
redpy.table.initializeTable(optto)

if args.verbose: print("Opening hdf5 table: {0}".format(optfrom.filename))
h5filefrom, rtablefrom, otablefrom, ttablefrom, ctablefrom, jtablefrom, dtablefrom, ftablefrom = redpy.table.openTable(optfrom)

if args.verbose: print("Opening hdf5 table: {0}".format(optto.filename))
h5fileto, rtableto, otableto, ttableto, ctableto, jtableto, dtableto, ftableto = redpy.table.openTable(optto)

# Define how many stations need to be added
dsta = optto.nsta - optfrom.nsta

# DO ALL THE COPYING
for rfrom in rtablefrom.iterrows():
    rto = rtableto.row
    # These stay the same
    rto['id'] = rfrom['id']
    rto['startTime'] = rfrom['startTime']
    rto['startTimeMPL'] = rfrom['startTimeMPL']
    rto['windowStart'] = rfrom['windowStart']
    # These must be extended
    rto['windowAmp'] = np.append(rfrom['windowAmp'],np.zeros(dsta))
    rto['windowCoeff'] = np.append(rfrom['windowCoeff'],np.zeros(dsta))
    rto['FI'] = np.append(rfrom['FI'],np.empty(dsta)*np.nan)
    rto['waveform'] = np.append(rfrom['waveform'],np.zeros(dsta*optto.wshape))
    rto['windowFFT'] = np.append(rfrom['windowFFT'],np.zeros(dsta*optto.winlen))
    rto.append()
rtableto.attrs.ptime = rtablefrom.attrs.ptime
rtableto.attrs.previd = rtablefrom.attrs.previd
rtableto.flush()

for ofrom in otablefrom.iterrows():    
    oto = otableto.row
    # These stay the same
    oto['id'] = ofrom['id']
    oto['startTime'] = ofrom['startTime']
    oto['startTimeMPL'] = ofrom['startTimeMPL']
    oto['windowStart'] = ofrom['windowStart']
    oto['expires'] = ofrom['expires']
    # These must be extended
    oto['windowAmp'] = np.append(ofrom['windowAmp'],np.zeros(dsta))
    oto['windowCoeff'] = np.append(ofrom['windowCoeff'],np.zeros(dsta))
    oto['FI'] = np.append(ofrom['FI'],np.empty(dsta)*np.nan)
    oto['waveform'] = np.append(ofrom['waveform'],np.zeros(dsta*optto.wshape))
    oto['windowFFT'] = np.append(ofrom['windowFFT'],np.zeros(dsta*optto.winlen))
    oto.append()
otableto.flush()

for tfrom in ttablefrom.iterrows():
    tto = ttableto.row
    # This stays the same
    tto['startTimeMPL'] = tfrom['startTimeMPL']
    tto.append()
ttableto.flush()

for dfrom in dtablefrom.iterrows():
    dto = dtableto.row
    # These stay the same
    dto['id'] = dfrom['id']
    dto['startTime'] = dfrom['startTime']
    dto['startTimeMPL'] = dfrom['startTimeMPL']
    dto['windowStart'] = dfrom['windowStart']
    # These must be extended
    dto['windowAmp'] = np.append(dfrom['windowAmp'],np.zeros(dsta))
    dto['windowCoeff'] = np.append(dfrom['windowCoeff'],np.zeros(dsta))
    dto['FI'] = np.append(dfrom['FI'],np.empty(dsta)*np.nan)
    dto['waveform'] = np.append(dfrom['waveform'],np.zeros(dsta*optto.wshape))
    dto['windowFFT'] = np.append(dfrom['windowFFT'],np.zeros(dsta*optto.winlen))
    dto.append()
dtableto.flush()

for jfrom in jtablefrom.iterrows():
    jto = jtableto.row
    # These stay the same
    jto['startTime'] = jfrom['startTime']
    jto['windowStart'] = jfrom['windowStart']
    jto['isjunk'] = jfrom['isjunk']
    # This must be extended
    jto['waveform'] = np.append(jfrom['waveform'],np.zeros(dsta*optto.wshape))
    jto.append()
jtableto.flush()

for cfrom in ctablefrom.iterrows():
    cto = ctableto.row
    # All stay the same
    cto['id1'] = cfrom['id1']
    cto['id2'] = cfrom['id2']
    cto['ccc'] = cfrom['ccc']
    cto.append()
ctableto.flush()

for ffrom in ftablefrom.iterrows():
    fto = ftableto.row
    # All stay the same, but printme == 1
    fto['members'] = ffrom['members']
    fto['core'] = ffrom['core']
    fto['startTime'] = ffrom['startTime']
    fto['longevity'] = ffrom['longevity']
    fto['lastprint'] = ffrom['lastprint']
    fto['printme'] = 1
    fto.append()
ftableto.attrs.nClust = ftablefrom.attrs.nClust
ftableto.flush()

if args.noplot:
    if args.verbose: print("Skipping plots...")
else:
    if args.verbose: print("Creating plots...")
    redpy.plotting.createPlots(rtableto, ftableto, ttableto, ctableto, otableto, optto)

if args.verbose: print("Closing tables...")
h5filefrom.close()
h5fileto.close()

if args.verbose: print("Deleting working copy of old hdf5 table...")
os.remove(optfrom.filename)

if args.verbose: print("Done")