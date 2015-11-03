import argparse
import redpy
from obspy import UTCDateTime

# Added this to remove the slew of warnings obspy/numpy was throwing at me
import warnings
warnings.filterwarnings("ignore")

tend = UTCDateTime()

"""
Run this script to fill the table with data from the last N seconds (defined in config as
nsec). Intended to be run every N seconds as a cron job to constantly update the table in
near real-time.
 
usage: single.py [-h] [-v] [-c CONFIGFILE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
"""

parser = argparse.ArgumentParser(description=
    "Fills table with data from the last nsec seconds")
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
h5file, rtable, otable, ctable, jtable = redpy.table.openTable(opt)

tstart = tend-opt.nsec

if args.verbose: print(tstart)

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

redpy.table.clearExpiredOrphans(otable, opt, tend)
if args.verbose: print("Length of Orphan table: {0}".format(len(otable)))

# This probably doesn't need to be run every time, but is included for consistency
#     with backfill.py
if len(rtable) > 0:
    redpy.cluster.runFullOPTICS(rtable, ctable)

if args.verbose: print("Closing table...")
h5file.close()
if args.verbose: print("Done")