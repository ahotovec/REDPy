import argparse
import matplotlib.pyplot as plt
import redpy

"""
Run this script to make plots

usage: makeplots.py [-h] [-v] [-c CONFIGFILE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
"""

parser = argparse.ArgumentParser(description=
    "Make plots using data in hdf5 table noted in configuration")
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

h5file, rtable, otable, ctable, jtable = redpy.table.openTable(opt)

if len(rtable) > 1:
    redpy.cluster.alignAllDeep(rtable, ctable, opt)
    redpy.cluster.runFullOPTICS(rtable, ctable, opt)
    redpy.plotting.createCMatrixFigure(rtable, ctable)
    redpy.plotting.createOrderedWaveformFigure(rtable, opt)
else:
    print("Nothing in the repeater table to plot!")

if args.verbose:
    print("Orphans saved: {0}".format(len(otable)))
    print("Number of junk triggers: {0}".format(len(jtable)))   

plt.show()

h5file.close()