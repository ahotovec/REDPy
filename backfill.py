import argparse
import redpy
from obspy import UTCDateTime

# Added this to remove the slew of warnings obspy/numpy was throwing at me
import warnings
warnings.filterwarnings("ignore")

"""
Run this script to fill the table with data from the past. Requires a start time in the
format "YYYY-MM-DD" or "YYYY-MM-DDTHH:MM:SS". If an end time is not specified, "now" is
assumed, which updates at the end of each time chunk processed (default: by hour, set in
configuration).

WARNING: Does not currently check for duplicates, do not run with any overlap of previous
runs!
 
usage: backfill.py [-h] [-v] [-e ENDTIME] [-c CONFIGFILE] starttime

positional arguments:
  starttime             start time to begin filling (YYYY-MM-DDTHH:MM:SS)

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -e ENDTIME, --endtime ENDTIME
                        optional end time to end filling (YYYY-MM-DDTHH:MM:SS)
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
"""

parser = argparse.ArgumentParser(description=
    "Backfills table with data from the past")
parser.add_argument("starttime", help="start time to begin filling (YYYY-MM-DDTHH:MM:SS)")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-e", "--endtime",
    help="optional end time to end filling (YYYY-MM-DDTHH:MM:SS)")
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
h5file, rtable, otable, ctable, jtable = redpy.table.openTable(opt)

tstart = UTCDateTime(args.starttime)
if args.endtime:
    tend = UTCDateTime(args.endtime)
else:
    tend = UTCDateTime()

n = 0
while tstart+n*opt.nsec <= tend-opt.nsec:

    if args.verbose: print(tstart+n*opt.nsec)
    
    # Download and trigger
    st = redpy.trigger.getData(tstart+n*opt.nsec, opt)
    alltrigs = redpy.trigger.trigger(st, rtable, opt)
    
    # Clean out data spikes etc.
    trigs, junk = redpy.trigger.dataclean(alltrigs, opt, flag=1)
    
    # Save junk triggers in separate table for quality checking purposes
    for i in range(len(junk)):
        redpy.table.populateJunk(jtable,junk[i],0,opt)
    
    if len(trigs) > 0:
        
        ostart = 0
        if len(otable) == 0:
            # First trigger goes to orphans table
            redpy.table.populateOrphan(otable, 0, trigs[0], opt)
            ostart = 1
        
        # Loop through remaining triggers
        for i in range(ostart,len(trigs)):  
            id = rtable.attrs.previd + i
            redpy.correlation.runCorrelation(rtable, otable, ctable, trigs[i], id, opt)
        rtable.attrs.previd = id + 1
    
    redpy.table.clearExpiredOrphans(otable, opt, tstart+(n+1)*opt.nsec)
    if args.verbose: print("Length of Orphan table: {0}".format(len(otable)))
    
    # Update tend if an end date is not specified so this will run until it is fully 
    # caught up, instead of running to when the script was originally run.
    if not args.endtime:
        tend = UTCDateTime()
        
    n = n+1

print("Caught up to: {0}".format(tstart+n*opt.nsec))

if len(rtable) > 0:
    redpy.cluster.runFullOPTICS(rtable, ctable, opt)

if args.verbose: print("Closing table...")
h5file.close()
if args.verbose: print("Done")