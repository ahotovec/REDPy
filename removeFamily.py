import redpy.config
import redpy.table
import argparse

# Added this to remove the slew of warnings obspy/numpy was throwing at me
import warnings
warnings.filterwarnings("ignore")

"""
Run this script to manually remove families/clusters (e.g., correlated noise that made it
past the 'junk' detector). Reclusters and remakes images when done.

usage: removeFamily.py [-h] [-v] [-c CONFIGFILE] N [N ...]

positional arguments:
  N                     family number(s) to be moved and deleted

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
"""

parser = argparse.ArgumentParser(description=
    "Run this script to manually remove families/clusters")
parser.add_argument('famnum', metavar='N', type=int, nargs='+',
    help="family number(s) to be moved and deleted")
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
h5file, rtable, otable, ctable, jtable, dtable, ftable = redpy.table.openTable(opt)

redpy.table.removeFamilies(rtable, ctable, dtable, ftable, args.famnum, opt)

if args.verbose: print("Creating plots...")
redpy.plotting.createBokehTimelineFigure(rtable, ctable, ftable, opt)

if args.verbose: print("Closing table...")
h5file.close()
if args.verbose: print("Done")