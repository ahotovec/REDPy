# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import redpy.config
import redpy.table
import argparse
import os
import numpy as np
import matplotlib.dates

"""
Run this script to manually produce an editable PDF version of the overview page

usage: createPDFOverview.py [-h] [-v] [-c CONFIGFILE] [-s STARTTIME] [-e ENDTIME]
    [-b BINSIZE] [-u] [-m MINMEMBERS] [-o OCCURHEIGHT] [-f FORMAT]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
  -s STARTTIME, --starttime STARTTIME
                        earliest time to plot, defaults to first trigger
  -e ENDTIME, --endtime ENDTIME
                        latest time to plot, defaults to last trigger
  -b BINSIZE, --binsize BINSIZE
                        custom time bin size in days, defaults to overview.html's
                        binsize
  -u, --usehrs          use hours instead of days for definition of binsize
  -m MINMEMBERS, --minmembers MINMEMBERS
                        minimum number of members to plot, defaults to
                        overview.html's minmembers
  -o OCCURHEIGHT, --occurheight OCCURHEIGHT
                        integer multiplier for how much taller the occurrence
                        plot should be compared to other plots, defaults to 3
  -f FORMAT, --format FORMAT
                        comma separated list of plots to be rendered
"""

parser = argparse.ArgumentParser(description=
    "Run this script to manually produce an editable PDF version of the overview page")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-c", "--configfile",
    help="use configuration file named CONFIGFILE instead of default settings.cfg")
parser.add_argument("-s", "--starttime",
    help="earliest time to plot, defaults to first trigger")
parser.add_argument("-e", "--endtime",
    help="latest time to plot, defaults to last trigger")
parser.add_argument("-b", "--binsize",
    help="custom time bin size, defaults to overview.html's binsize")
parser.add_argument("-u", "--usehrs", action="count", default=0,
    help="use hours instead of days for definition of binsize")
parser.add_argument("-m", "--minmembers",
    help="minimum number of members to plot, defaults to overview.html's minmembers")
parser.add_argument("-o", "--occurheight",
    help="integer multiplier for how much taller the occurrence plot should be compared" +
        " to other plots, defaults to 3")
parser.add_argument("-f", "--format",
    help="comma separated list of plots to be rendered")
args = parser.parse_args()

if args.configfile:
    opt = redpy.config.Options(args.configfile)
    if args.verbose: print("Using config file: {0}".format(args.configfile))
else:
    opt = redpy.config.Options("settings.cfg")
    if args.verbose: print("Using config file: settings.cfg")

if args.verbose: print("Opening hdf5 table: {0}".format(opt.filename))
h5file, rtable, otable, ttable, ctable, jtable, dtable, ftable = redpy.table.openTable(opt)

# Check for MPL version mismatch
redpy.table.checkMPL(rtable, ftable, ttable, otable, dtable, opt)

# Process arguments
if args.starttime:
    tmin = matplotlib.dates.date2num(np.datetime64(args.starttime))
else:
    tmin = 0

if args.endtime:
    tmax = matplotlib.dates.date2num(np.datetime64(args.endtime))
else:
    tmax = 0

if args.binsize:
    if args.usehrs:
        binsize = float(args.binsize)/24
    else:
        binsize = float(args.binsize)
else:
    binsize = opt.dybin

if args.minmembers:
    minmembers = int(args.minmembers)
else:
    minmembers = opt.minplot

if args.occurheight:
    occurheight = int(args.occurheight)
else:
    occurheight = 3

if args.format:
    plotformat = args.format
else:
    plotformat = 'eqrate,fi,occurrence,longevity'

if args.verbose: print("Creating overview.pdf in main output directory...")
redpy.plotting.customPDFoverview(rtable, ftable, ttable, tmin, tmax, binsize, minmembers,
    occurheight, plotformat, opt)

if args.verbose: print("Closing table...")
h5file.close()
if args.verbose: print("Done")
