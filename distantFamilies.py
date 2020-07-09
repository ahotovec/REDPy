# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import argparse
import redpy.config
import glob, os, itertools
import numpy as np

"""
Run this script to print out the families with a minimum percentage of regional and/or
teleseismic matches from ComCat that can then be copy/pasted into removeFamily.py. An
optional table is printed that summarizes matches of each type. This only parses the
existing .html files rather than querying ComCat directly!

usage: distantFamilies.py [-h] [-v] [-c CONFIGFILE] [-e ETC] [-p PERCENT]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements, including table of matches
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
  -e ETC, --etc ETC     phrase to explicitly search for, e.g., name of a specific area
  -p PERCENT, --percent PERCENT 
                        minimum percentage of regional/teleseismic/etc matches, default 90
"""

parser = argparse.ArgumentParser(description=
    "Finds families with regional/teleseismic matches by parsing their .html files")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements, including table of matches")
parser.add_argument("-c", "--configfile",
    help="use configuration file named CONFIGFILE instead of default settings.cfg")
parser.add_argument("-e", "--etc",
    help="phrase to explicitly search for, e.g., name of a specific area")
parser.add_argument("-p", "--percent", type=float,
    help="minimum percentage of regional/teleseismic matches, default 90", default=90.0)
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

fnums = []
removeNums = ''
removeNumsReg = ''
removeNumsTele = ''
removeNumsETC = ''

# Sort by family number (the list is in a strange order)
for f in flist:
    fnum = int(f.split("/")[-1][:-5])
    fnums.append(fnum)

# Parse each file, counting the number of times a word/phrase is matched
for f in flist[np.argsort(fnums)]:

    file = open(f, "r")
    fnum = f.split("/")[-1][:-5]
    data = file.read()

    reg = data.count("regional")
    tele = data.count("teleseismic")
    local = data.count("Potential local match:") # excludes the last two lines!
    if args.etc:
        etc = data.count(args.etc)
        local = local-etc
    else:
        etc = 0
    
    if reg+tele+etc > 0:
    
        if args.verbose:
            if args.etc:
                print("Fam {:4} : L {:2} | R {:2} | T {:2} | E {:2} | Distant {:5.1f}% | \
Etc {:5.1f}%".format(
                    fnum, local, reg, tele, etc, 100*(reg+tele+etc)/(reg+tele+local+etc),
                    100*(etc)/(reg+tele+local+etc)))
            else:
                print("Fam {:4} : L {:2} | R {:2} | T {:2} | Distant {:5.1f}%".format(
                    fnum, local, reg, tele, 100*(reg+tele+etc)/(reg+tele+local+etc)))
            
        if 100*(reg+tele)/(reg+tele+local+etc) >= args.percent:
            removeNums+=' {}'.format(fnum)
        
        if 100*(reg)/(reg+tele+local+etc) >= args.percent:
            removeNumsReg+=' {}'.format(fnum)
        
        if 100*(tele)/(reg+tele+local+etc) >= args.percent:
            removeNumsTele+=' {}'.format(fnum)
        
        if args.etc:
            if 100*(etc)/(reg+tele+local+etc) >= args.percent:
                removeNumsETC+=' {}'.format(fnum)
        

print('\n{}%+ Regional+Teleseismic:\n{}\n'.format(args.percent,removeNums))

print('\n{}%+ Regional:\n{}\n'.format(args.percent,removeNumsReg))

print('\n{}%+ Teleseismic:\n{}\n'.format(args.percent,removeNumsTele))

if args.etc:
    print('{}%+ containing {}:\n{}\n'.format(args.percent,args.etc,removeNumsETC))

