import matplotlib
matplotlib.use('TkAgg')
# ^ Not sure why this needs to be here, but apparently Tk hates me otherwise...

import tkinter as tk
import redpy.config
import redpy.table
import argparse
from PIL import Image
import os
import glob

# Added this to remove the slew of warnings obspy/numpy was throwing at me
import warnings
warnings.filterwarnings("ignore")

"""
Run this script to manually remove families/clusters (e.g., correlated noise that made it
past the 'junk' detector) using a GUI interface. Reclusters and remakes images when done.
Note: using large NCOLS may make the window too wide for your monitor, and the GUI does
not currently support side scrolling...

usage: removeFamilyGUI.py [-h] [-v] [-c CONFIGFILE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
  -n NCOLS, --ncols NCOLS
                        adjust number of columns in layout (default 3)
  -m MINCLUST, --minclust MINCLUST
                        only look at clusters with numbers at or above MINCLUST
"""


# Define some functions specific to this GUI
def remove(*args):
    """
    Run the removal script using checked boxes
    """
    print('\nYou have selected the following families to remove:')
    removethese = []
    for n in range(len(var)):
        if var[n].get() > 0:
            removethese.append(n+m)
    print(removethese)
    root.destroy() # Close the window

    redpy.table.removeFamilies(rtable, ctable, dtable, ftable, removethese, opt)
        
    if len(removethese) > 0:
        print("Creating plots...")
        redpy.plotting.createPlots(rtable, ftable, ttable, opt)

def close(*args):
    """
    Close the window and the table
    """
    root.destroy()

def onFrameConfigure(canvas):
    """
    Reset the scroll region to encompass the inner frame
    """
    canvas.configure(scrollregion=canvas.bbox("all"))

def mouse_wheel(event):
    """
    Mousewheel scrolling is a bit squiffy
    (only scrolls down, but better than nothing?)
    """
    canvas.yview_scroll(-1*(event.delta/120), "units")


parser = argparse.ArgumentParser(description=
    "Run this script to manually remove families/clusters using a GUI")
parser.add_argument("-v", "--verbose", action="count", default=0,
    help="increase written print statements")
parser.add_argument("-c", "--configfile",
    help="use configuration file named CONFIGFILE instead of default settings.cfg")
parser.add_argument("-n", "--ncols", default=3, type=int,
    help="adjust number of columns in layout (default 3)")
parser.add_argument("-m", "--minclust", default=0, type=int,
    help="only look at clusters with numbers at or above MINCLUST")
args = parser.parse_args()

if args.configfile:
    opt = redpy.config.Options(args.configfile)
    if args.verbose: print("Using config file: {0}".format(args.configfile))
else:
    opt = redpy.config.Options("settings.cfg")
    if args.verbose: print("Using config file: settings.cfg")

if args.verbose: print("Opening hdf5 table: {0}".format(opt.filename))
h5file, rtable, otable, ttable, ctable, jtable, dtable, ftable = redpy.table.openTable(opt)

if args.minclust:
    m = args.minclust
else:
    m = 0

# Create GUI window
root = tk.Tk()
root.title("REDPy - Check Families to Permanently Remove")
canvas = tk.Canvas(root, borderwidth=0, width=560*args.ncols, height=1500, background="#ffffff")
frame = tk.Frame(canvas, background="#ffffff")
vsb = tk.Scrollbar(root, orient="vertical", command=canvas.yview)
canvas.configure(yscrollcommand=vsb.set)
vsb.pack(side="right", fill="y")
canvas.pack(side="left", fill="both", expand=True)
canvas.create_window((4,4), window=frame, anchor="nw")
frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

# Build grid of families
fams = range(len(ftable))
r = 1
c = 1
imgobj = []
check = []
var = []
for n in fams:
    if n >= args.minclust:
        im = Image.open('{0}/clusters/{1}.png'.format(opt.groupName,n))
        im.save('{0}/clusters/{1}.gif'.format(opt.groupName,n))
        imgobj.append(tk.PhotoImage(file='{0}/clusters/{1}.gif'.format(opt.groupName,n)))
        var.append(tk.IntVar())
        check.append(tk.Checkbutton(frame, image=imgobj[n-m],
            variable = var[n-m]).grid(column=c, row=r, sticky='N'))
        c = c+1
        if c == args.ncols+1:
            c = 1
            r = r+1
            if r > 255:
                print("Ran out of rows. Use -n or -m flags to view more...")

print("\nIgnore these warning things:")

# Add buttons
tk.Button(frame, text="Remove Checked", background="#ffffff", command=remove).grid(
    column=1, row=r+1, columnspan=args.ncols, sticky='N')
tk.Button(frame, text="Cancel", background="#ffffff", command=close).grid(
    column=1, row=r+2, columnspan=args.ncols, sticky='S')

# Bind MouseWheel, Return, Escape keys to be more useful
root.bind_all("<MouseWheel>", mouse_wheel)
root.bind('<Return>', remove)
root.bind('<Escape>', close)

# Add some padding
for child in frame.winfo_children(): child.grid_configure(padx=15, pady=15)

# Go!
root.mainloop()

# Clean up
print("\nCleaning up .gif files...")
dlist = glob.glob('./{0}/clusters/*.gif'.format(opt.groupName))
for tmp in dlist:
    os.remove(tmp) 
print("Closing table...")
h5file.close()
print("Done")