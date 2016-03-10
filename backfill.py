import argparse
import redpy
import numpy as np
import obspy
from obspy import UTCDateTime
import time

# Added this to remove the slew of warnings obspy/numpy was throwing at me
import warnings
warnings.filterwarnings("ignore")

"""
Run this script to fill the table with data from the past. If a start time is not
specified, it will check the attributes of the repeater table to pick up where it left
off. Additionally, if this is the first run and a start time is not specified, it will
assume one time chunk prior to the end time. If an end time is not specified, "now" is
assumed. The end time updates at the end of each time chunk processed (default: by hour,
set in configuration). This script can be run as a cron job that will pick up where it
left off if a chunk is missed, but will not run until a full chunk of time has elapsed
since the last trigger. Use -n if you are backfilling with a large amount of time; it will
consume less time downloading the data in small chunks if NSEC is an hour or a day instead
of a few minutes.

WARNING: Does not currently check for duplicates, do not run with any overlap of previous
runs!
 
usage: backfill.py [-h] [-v] [-s STARTTIME] [-e ENDTIME] [-c CONFIGFILE] [-n NSEC]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -s STARTTIME, --starttime STARTTIME
                        optional start time to begin filling (YYYY-MM-DDTHH:MM:SS)
  -e ENDTIME, --endtime ENDTIME
                        optional end time to end filling (YYYY-MM-DDTHH:MM:SS)
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
  -n NSEC, --nsec NSEC  overwrite opt.nsec from configuration file with NSEC this run only
"""

parser = argparse.ArgumentParser(description=
    "Backfills table with data from the past")
parser.add_argument("-s", "--starttime",
    help="optional start time to begin filling (YYYY-MM-DDTHH:MM:SS)")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
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
h5file, rtable, otable, ctable, jtable, dtable, ftable = redpy.table.openTable(opt)
    
if args.endtime:
    tend = UTCDateTime(args.endtime)
else:
    tend = UTCDateTime()
    
if args.starttime:
    tstart = UTCDateTime(args.starttime)
else:
    if rtable.attrs.ptime:
        tstart = UTCDateTime(rtable.attrs.ptime)
    else:
        tstart = tend-opt.nsec

if len(otable) > 0:
    otimes = otable.cols.startTimeMPL[:]
else:
    otimes = 0
if len(rtable) > 0:
    rtimes = rtable.cols.startTimeMPL[:]
else:
    rtimes = 0

t = time.time()
n = 0
rlen = len(rtable)
while tstart+n*opt.nsec <= tend-opt.nsec:
    
    ti = time.time()
    if args.verbose: print(tstart+n*opt.nsec)
    
    # Download and trigger
    try:
        st, stC = redpy.trigger.getData(tstart+n*opt.nsec, opt)
        alltrigs = redpy.trigger.trigger(st, stC, rtable, opt)
    except (TypeError, obspy.fdsn.header.FDSNException, Exception):
        print('Could not download or trigger data... moving on')
        alltrigs = []
    
	# Clean out data spikes etc.
    trigs, junk = redpy.trigger.dataclean(alltrigs, opt, flag=1)
        
	# Save junk triggers in separate table for quality checking purposes
    for i in range(len(junk)):
        redpy.table.populateJunk(jtable,junk[i],0,opt)
            
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
                redpy.correlation.runCorrelation(rtable, otable, ctable, ftable, otimes, rtimes,
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
                redpy.correlation.runCorrelation(rtable, otable, ctable, ftable, otimes, rtimes,
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

print("Caught up to: {}".format(tstart+n*opt.nsec))
print("End time now: {}".format(tend))

if len(rtable) > rlen:
    if args.verbose: print("Creating plots...")
    redpy.plotting.createBokehTimelineFigure(rtable, ctable, ftable, opt)
else:
    if args.verbose: print("No new repeaters to plot.")

if args.verbose: print("Closing table...")
h5file.close()

if args.verbose: print("Total time spent: {} minutes".format((time.time()-t)/60))
if args.verbose: print("Done")