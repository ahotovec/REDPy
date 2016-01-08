import tkinter as tk
import redpy.config
import redpy.table
import argparse

# Added this to remove the slew of warnings obspy/numpy was throwing at me
import warnings
warnings.filterwarnings("ignore")

"""
Run this script to manually remove families/clusters (e.g., correlated noise that made it
past the 'junk' detector) using a GUI interface. Reclusters and remakes images when done.

usage: removeFamilyGUI.py [-h] [-v] [-c CONFIGFILE]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         increase written print statements
  -c CONFIGFILE, --configfile CONFIGFILE
                        use configuration file named CONFIGFILE instead of
                        default settings.cfg
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
            removethese.append(n)
    print(removethese)
    root.destroy() # Close the window
    
    for f in removethese:
        print("Removing family {}...".format(f))
        redpy.table.removeFamily(rtable, ctable, dtable, f, opt)
        
    if len(removethese) > 0:
        print("Creating plots...")
        redpy.plotting.createBokehTimelineFigure(rtable, ctable, opt)

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
args = parser.parse_args()

if args.configfile:
    opt = redpy.config.Options(args.configfile)
    if args.verbose: print("Using config file: {0}".format(args.configfile))
else:
    opt = redpy.config.Options("settings.cfg")
    if args.verbose: print("Using config file: settings.cfg")

if args.verbose: print("Opening hdf5 table: {0}".format(opt.filename))
h5file, rtable, otable, ctable, jtable, dtable = redpy.table.openTable(opt)


print('\nIgnore these warning things:')

# Create GUI window
root = tk.Tk()
root.title("REDPy - Check Families to Permanently Remove")
canvas = tk.Canvas(root, borderwidth=0, width=1675, height=1000, background="#ffffff")
frame = tk.Frame(canvas, background="#ffffff")
vsb = tk.Scrollbar(root, orient="vertical", command=canvas.yview)
canvas.configure(yscrollcommand=vsb.set)
vsb.pack(side="right", fill="y")
canvas.pack(side="left", fill="both", expand=True)
canvas.create_window((4,4), window=frame, anchor="nw")
frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

# Build grid of families
fams = range(len(rtable.attrs.cores))
r = 1
c = 1
imgobj = []
check = []
var = []
for n in fams:
    imgobj.append(tk.PhotoImage(file='default/clusters/{}.gif'.format(n)))
    var.append(tk.IntVar())
    check.append(tk.Checkbutton(frame, image=imgobj[n], variable = var[n]).grid(
        column=c, row=r, sticky='N'))
    c = c+1
    if c == 4:
        c = 1
        r = r+1

# Add buttons
tk.Button(frame, text="Remove Checked", background="#ffffff", command=remove).grid(
    column=2, row=r+1, sticky='N')
tk.Button(frame, text="Cancel", background="#ffffff", command=close).grid(
    column=2, row=r+2, sticky='S')

# Bind MouseWheel, Return, Escape keys to be more useful
root.bind_all("<MouseWheel>", mouse_wheel)
root.bind('<Return>', remove)
root.bind('<Escape>', close)

# Add some padding
for child in frame.winfo_children(): child.grid_configure(padx=15, pady=15)

# Go!
root.mainloop()
print("\nClosing table...")
h5file.close()
print("Done")