# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import argparse
import redpy
import numpy as np
import obspy
from obspy import UTCDateTime
import time

"""
Run this script to fill the table with data from the past. If a start time is not
specified, it will check the attributes of the repeater table to pick up where it left
off. Additionally, if this is the first run and a start time is not specified, it will
assume one time chunk prior to the end time. If an end time is not specified, "now" is
assumed. The end time updates at the end of each time chunk processed (default: by hour,
set in configuration). This script can be run as a cron job that will pick up where it
left off if a chunk is missed. Use -n if you are backfilling with a large amount of time;
it will consume less time downloading the data in small chunks if NSEC is an hour or a day
instead of a few minutes, but at the cost of keeping orphans for longer.

usage: backfill.py [-h] [-v] [-t] [-s STARTTIME] [-e ENDTIME] [-c CONFIGFILE] [-n NSEC]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -t, --troubleshoot    run in troubleshoot mode (without try/except)
  -s STARTTIME, --starttime STARTTIME
                        optional start time to begin filling (YYYY-MM-DDTHH:MM:SS)
  -e ENDTIME, --endtime ENDTIME
                        optional end time to end filling (YYYY-MM-DDTHH:MM:SS)
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
  -n NSEC, --nsec NSEC  overwrite opt.nsec from configuration file with NSEC this run only
"""

t = time.time()

parser = argparse.ArgumentParser(description=
    "Backfills table with data from the past")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-t", "--troubleshoot", action="count", default=0,
    help="run in troubleshoot mode (without try/except)")
parser.add_argument("-s", "--starttime",
    help="optional start time to begin filling (YYYY-MM-DDTHH:MM:SS)")
parser.add_argument("-e", "--endtime",
    help="optional end time to end filling (YYYY-MM-DDTHH:MM:SS)")
parser.add_argument("-c", "--configfile",
    help="use configuration file named CONFIGFILE instead of default settings.cfg")
parser.add_argument("-n", "--nsec", type=int,
    help="overwrite opt.nsec from configuration file with NSEC this run only")
args = parser.parse_args()

if args.configfile:
    opt = redpy.config.Options(args.configfile)
    if args.verbose: print("Using config file: {0}".format(args.configfile))
else:
    opt = redpy.config.Options("settings.cfg")
    if args.verbose: print("Using config file: settings.cfg")

if args.nsec:
    opt.nsec = args.nsec

if args.verbose: print("Opening hdf5 table: {0}".format(opt.filename))
h5file, rtable, otable, ttable, ctable, jtable, dtable, ftable = redpy.table.openTable(opt)

# Check for MPL version mismatch
redpy.table.checkMPL(rtable, ftable, ttable, otable, dtable, opt)

if args.endtime:
    tend = UTCDateTime(args.endtime)
else:
    tend = UTCDateTime()

if args.starttime:
    tstart = UTCDateTime(args.starttime)
    if rtable.attrs.ptime:
        rtable.attrs.ptime = UTCDateTime(tstart)
else:
    if rtable.attrs.ptime:
        tstart = UTCDateTime(rtable.attrs.ptime)
    else:
        tstart = tend-opt.nsec

if len(ttable) > 0:
    ttimes = ttable.cols.startTimeMPL[:]
else:
    ttimes = 0

n = 0
rlen = len(rtable)
while tstart+n*opt.nsec < tend:

    ti = time.time()
    print(tstart+n*opt.nsec)

    # Download and trigger
    if args.troubleshoot:
        endtime = tstart+(n+1)*opt.nsec+opt.atrig
        if endtime > tend:
            endtime = tend
        st, stC = redpy.trigger.getData(tstart+n*opt.nsec-opt.atrig, endtime, opt)
        alltrigs = redpy.trigger.trigger(st, stC, rtable, opt)
    else:
        try:
            endtime = tstart+(n+1)*opt.nsec+opt.atrig
            if endtime > tend:
                endtime = tend
            st, stC = redpy.trigger.getData(tstart+n*opt.nsec-opt.atrig, endtime, opt)
            alltrigs = redpy.trigger.trigger(st, stC, rtable, opt)
        except (TypeError, obspy.clients.fdsn.header.FDSNException, Exception):
            print('Could not download or trigger data... moving on')
            alltrigs = []

	# Clean out data spikes etc.
    trigs, junk, junkFI, junkKurt = redpy.trigger.dataClean(alltrigs, opt, flag=1)

	# Save junk triggers in separate table for quality checking purposes
    for i in range(len(junk)):
        redpy.table.populateJunk(jtable, junk[i], 2, opt) # Both types of junk
    for i in range(len(junkKurt)):
        redpy.table.populateJunk(jtable, junkKurt[i], 1, opt) # Just kurtosis junk
    for i in range(len(junkFI)):
        redpy.table.populateJunk(jtable, junkFI[i], 0, opt) # Just 'teleseisms'

    # Append times of triggers to ttable to compare total seismicity later
    redpy.table.populateTriggers(ttable, trigs, ttimes, opt)

    # Check triggers against deleted events
    if len(dtable) > 0:
        trigs = redpy.correlation.compareDeleted(trigs, dtable, opt)

    if len(trigs) > 0:
        id = rtable.attrs.previd
        if len(trigs) == 1:
            ostart = 0
            if len(otable) == 0:
                # First trigger goes to orphans table
                redpy.table.populateOrphan(otable, 0, trigs[0], opt)
                ostart = 1
            else:
                id = id + 1
                redpy.correlation.runCorrelation(rtable, otable, ctable, ftable, ttimes,
                    trigs[0], id, opt)
        else:
            ostart = 0
            if len(otable) == 0:
                # First trigger goes to orphans table
                redpy.table.populateOrphan(otable, 0, trigs[0], opt)
                ostart = 1
            # Loop through remaining triggers
            for i in range(ostart,len(trigs)):
                id = id + 1
                redpy.correlation.runCorrelation(rtable, otable, ctable, ftable, ttimes,
                    trigs[i], id, opt)
        rtable.attrs.previd = id

    redpy.table.clearExpiredOrphans(otable, opt, tstart+(n+1)*opt.nsec)

    # Print some stats
    if args.verbose:
        print("Length of Orphan table: {}".format(len(otable)))
        if len(rtable) > 1:
            print("Number of repeaters: {}".format(len(rtable)))
            print("Number of clusters: {}".format(ftable.attrs.nClust))

    # Update tend if an end date is not specified so this will run until it is fully
    # caught up, instead of running to when the script was originally run.
    if not args.endtime:
        tend = UTCDateTime()

    n = n+1

    if args.verbose: print("Time spent this iteration: {} minutes".format(
        (time.time()-ti)/60))

print("Caught up to: {}".format(endtime-opt.atrig))

if args.verbose: print("Updating plots...")
redpy.plotting.createPlots(rtable, ftable, ttable, ctable, otable, opt)

if args.verbose: print("Closing table...")
h5file.close()

print("Total time spent: {} minutes".format((time.time()-t)/60))
if args.verbose: print("Done")
