# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import argparse
import numpy as np
import os

"""
Run this script to generate 'meta.html' in a specified directory and with a list of runs.
This page gathers the 'meta_recent.html' tabbed overviews within the output directories
into a single page.

usage: makeMeta.py [-h] [-v] [-p PATH] [-r RUNS]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -p PATH, --path PATH  relative path to where meta.html should be created (e.g., './')
  -r RUNS, --runs RUNS  comma-separated list of runs to include, which should match their
                        groupName
"""

parser = argparse.ArgumentParser(description=
    "Run this script to output the contents of the junk table for troubleshooting.")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-p", "--path",
    help="relative path to where meta.html should be created", default='./')
parser.add_argument("-r", "--runs",
    help="comma-separated list of runs to include, which should match their groupName")

args = parser.parse_args()

filename = '{}meta.html'.format(args.path)
if args.verbose: print("Creating {}...".format(filename))

if args.runs:
    if args.verbose: print("Looping over runs: {}".format(args.runs))
    runs = args.runs
else:
    print("No runs supplied, assuming 'default' only")
    runs = 'default'
    
with open(filename, 'w') as f:
    f.write('<html><head><title>REDPy Meta Overview</title></head>')
    f.write('<body style="padding:0;margin:0">')

    for run in runs.split(','):
        f.write("""
        <iframe src="{0}{1}/meta_recent.html" title="{1}"
                style="height:350px;width:1300px;border:none;"></iframe>
                """.format(args.path,run))
    f.write('</body></html>')


if args.verbose: print("Done")