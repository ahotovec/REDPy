# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import argparse
import redpy.config
import glob, os, itertools
import numpy as np

"""
Run this script to print out the median location of each cluster. This currently only
parses the existing .html files rather than querying ComCat directly! Default behavior
only uses locations for local matched earthquakes.

usage: clusterLocs.py [-h] [-v] [-c CONFIGFILE] [-d]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements, including table of matches
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
  -d, --distant         include distant (regional, teleseismic) matches in addition to
                        local seismicity
  -r, --regional        include regional matches in addition to local seismicity
"""

parser = argparse.ArgumentParser(description=
    "Finds families with regional/teleseismic matches by parsing their .html files")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements, including table of matches")
parser.add_argument("-c", "--configfile",
    help="use configuration file named CONFIGFILE instead of default settings.cfg")
parser.add_argument("-d", "--distant", action="count", default=0,
    help="include distant (regional, teleseismic) matches in addition to local seismicity")
parser.add_argument("-r", "--regional", action="count", default=0,
    help="include regional matches in addition to local seismicity")
args = parser.parse_args()

if args.configfile:
    opt = redpy.config.Options(args.configfile)
    if args.verbose: print("Using config file: {0}".format(args.configfile))
else:
    opt = redpy.config.Options("settings.cfg")
    if args.verbose: print("Using config file: settings.cfg")

flist = np.array(list(itertools.chain.from_iterable(glob.iglob(os.path.join(
                root,'*.html')) for root, dirs, files in os.walk(
                '{}{}/clusters/'.format(opt.outputPath,opt.groupName)))))


### Open output file to write to here
with open('{}{}/clusterlocs.txt'.format(opt.outputPath,opt.groupName), 'w') as outfile:

    outfile.write('cnum latitude longitude depth\n')
    
    # Sort by family number (the list is in a strange order)
    fnums = []
    for f in flist:
        fnum = int(f.split("/")[-1][:-5])
        fnums.append(fnum)

    # Parse each file, counting the number of times a word/phrase is matched
    for f in flist[np.argsort(fnums)]:

        file = open(f, "r")
        fnum = f.split("/")[-1][:-5]
        data = file.readlines()
    
        lats = np.array([])
        lons = np.array([])
        deps = np.array([])
    
        lines = data[20].split('>')
        for line in lines:
            if args.distant:
                if line.count("teleseismic") or line.count("regional") or line.count(
                    "Potential local match:"):
                    lats = np.append(lats,float(line.split(' ')[4].strip('(,')))
                    lons = np.append(lons,float(line.split(' ')[5].strip(')')))
                    deps = np.append(deps,float(line.split(' ')[6].strip('km')))
            elif args.regional:
                if line.count("regional") or line.count("Potential local match:"):
                    lats = np.append(lats,float(line.split(' ')[4].strip('(,')))
                    lons = np.append(lons,float(line.split(' ')[5].strip(')')))
                    deps = np.append(deps,float(line.split(' ')[6].strip('km')))
            else: # Default behavior
                if line.count("Potential local match:"):
                    lats = np.append(lats,float(line.split(' ')[4].strip('(,')))
                    lons = np.append(lons,float(line.split(' ')[5].strip(')')))
                    deps = np.append(deps,float(line.split(' ')[6].strip('km')))
    
        if len(lats)>0:
            outfile.write('{} {:6.4f} {:7.4f} {:3.2f}\n'.format(
                fnum,np.median(lats),np.median(lons),np.median(deps)))
        else:
            outfile.write('{}   \n'.format(fnum))  

outfile.close()   
if args.verbose: print('Done writing to {}{}/clusterlocs.txt'.format(
    opt.outputPath,opt.groupName))
