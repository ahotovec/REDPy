# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2018  Alicia Hotovec-Ellis (ahotovec@gmail.com)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import redpy.config
import redpy.table
import argparse
import os

"""
Run this script first to initialize the hdf5 table where everything will be stored.
Warning: Running this script will overwrite an existing table with the same name defined
    by filename in the .cfg file.

usage: initialize.py [-h] [-v] [-c CONFIGFILE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
"""

parser = argparse.ArgumentParser(description=
    "Initialize hdf5 table using configuration, overwrites existing table defined in config")
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

if args.verbose: print("Writing hdf5 table: {0}".format(opt.filename))

redpy.table.initializeTable(opt)

if args.verbose: print("Creating folder to store images named '{}'".format(opt.groupName))
try:
    os.mkdir(opt.groupName)
except OSError:
    print("Folder exists.")
    
if args.verbose: print("Creating folder to store core images named '{}'/clusters".format(
    opt.groupName))
try:
    os.mkdir('{}/clusters'.format(opt.groupName))
except OSError:
    print("Folder exists.")

if args.verbose: print("Done")