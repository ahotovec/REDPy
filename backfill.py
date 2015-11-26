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
h5file, rtable, otable, ctable, jtable = redpy.table.openTable(opt)
    
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

t = time.time()
n = 0
rlen = len(rtable)
while tstart+n*opt.nsec <= tend-opt.nsec:
    
    if args.verbose: print(tstart+n*opt.nsec)
    
    # Download and trigger
    try:
        st = redpy.trigger.getData(tstart+n*opt.nsec, opt)
        alltrigs = redpy.trigger.trigger(st, rtable, opt)
        
		# Clean out data spikes etc.
        trigs, junk = redpy.trigger.dataclean(alltrigs, opt, flag=1)
        
		# Save junk triggers in separate table for quality checking purposes
        for i in range(len(junk)):
            redpy.table.populateJunk(jtable,junk[i],0,opt)
            
    except (TypeError, obspy.fdsn.header.FDSNException):
	    print('Could not download or trigger data... moving on')
	    trigs = []
    
    # Check number of clusters before adding new triggers
    if len(rtable) > 1:
        maxclust = max(rtable.cols.clusterNumber[:])
    else:
        maxclust = -1
    
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
                redpy.correlation.runCorrelation(rtable, otable, ctable, trigs[0], id, opt)        
        else:
            ostart = 0
            if len(otable) == 0:
                # First trigger goes to orphans table
                redpy.table.populateOrphan(otable, 0, trigs[0], opt)
                ostart = 1        
            # Loop through remaining triggers
            for i in range(ostart,len(trigs)):  
                id = id + 1
                redpy.correlation.runCorrelation(rtable, otable, ctable, trigs[i], id, opt)            
        rtable.attrs.previd = id
    
    redpy.table.clearExpiredOrphans(otable, opt, tstart+(n+1)*opt.nsec)
    
    # Attempt to merge families if a new cluster is born
    if len(rtable) > 1:
        if max(rtable.cols.clusterNumber[:]) > maxclust:
            redpy.cluster.mergeFamilies(rtable, ctable, opt)
            redpy.cluster.runFullOPTICS(rtable, ctable, opt)
    
    # Deal with leftovers (currently thrown away...)
    leftovers = rtable.get_where_list('clusterNumber == -1')
    if leftovers.any():
        leftovers[::-1].sort()
        print("Leftovers in clustering: {0}".format(len(leftovers)))
#         for l in leftovers:
#             rtable.remove_row(l)
    
    # Print some stats
    if args.verbose:
        print("Length of Orphan table: {}".format(len(otable)))
        if len(rtable) > 1:
            print("Number of repeaters: {}".format(len(rtable)))
            print("Number of clusters: {}".format(max(rtable.cols.clusterNumber[:])+1))
    
    # Update tend if an end date is not specified so this will run until it is fully 
    # caught up, instead of running to when the script was originally run.
    if not args.endtime:
        tend = UTCDateTime()
        
    n = n+1

print("Caught up to: {}".format(tstart+n*opt.nsec))
print("End time now: {}".format(tend))

print("Time spent: {} minutes".format((time.time()-t)/60))

if len(rtable) > rlen:
    if args.verbose: print("Creating plots...")
    redpy.cluster.runFullOPTICS(rtable, ctable, opt)
    redpy.plotting.createBokehTimelineFigure(rtable, ctable, opt)
else:
    print("No new repeaters to plot.")

if args.verbose: print("Closing table...")
h5file.close()
if args.verbose: print("Done")