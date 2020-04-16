# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import redpy.config
import redpy.table
import argparse

"""
Run this script to manually produce a more detailed 'report' page for a given family
(or families)

usage: createReport.py [-h] [-v] [-o] [-c CONFIGFILE] N [N ...]

positional arguments:
  N                     family number(s) to be reported on

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -o, --ordered         order plots by OPTICS
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
"""

parser = argparse.ArgumentParser(description=
    "Run this script to manually produce a more detailed 'report' page for a given " +
    "family (or families)")
parser.add_argument('famnum', metavar='N', type=int, nargs='+',
    help="family number(s) to be reported on")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-o", "--ordered", action="count", default=0,
    help="order plots by OPTICS")
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

for fnum in args.famnum:
    if args.verbose: print("Creating report for family {}...".format(fnum))
    redpy.plotting.plotReport(rtable, ftable, ctable, fnum, args.ordered, opt)

if args.verbose: print("Closing table...")
h5file.close()
if args.verbose: print("Done")
