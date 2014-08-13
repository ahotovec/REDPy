import numpy as np
import matplotlib.pyplot as plt
import time
import redpy

from tables import *
from redpy.optics import *
from obspy import UTCDateTime


"""
Control script to test scaling

Currently has a decent amount of printing to check on things as it's running
"""

opt = redpy.config.Options(filename="test.h5", cmin=0.7)

redpy.table.initializeTable(opt) 

# This time range has ~3400 triggers, a small handful of repeating earthquakes, and a 
# TON of noise spikes
t = UTCDateTime("2014-08-01")
st = redpy.trigger.getIRIS(t, opt, nsec=13*86400)

trigs = redpy.trigger.trigger(st, opt)

h5file = open_file(opt.filename, "a")
rtable = eval('h5file.root.'+ opt.groupName + '.repeaters')
otable = eval('h5file.root.'+ opt.groupName + '.orphans')
ctable = eval('h5file.root.'+ opt.groupName + '.correlation')

# First trigger goes to orphans table
redpy.table.populateOrphan(otable, 0, trigs[0], 'Never', opt)

t = time.time()
# Loop through remaining triggers
for i in range(1,len(trigs)):  
    print("{0} of {1}".format(i,len(trigs)-1))  
    redpy.correlation.runCorrelation(rtable, otable, ctable, trigs[i], i, opt)

print("Correlation done in: {:03.2f} seconds".format(time.time()-t))

# Run clustering one more time before plotting
redpy.cluster.runFullOPTICS(rtable, ctable)

# Print some information about the run
print("Repeaters found: {0}".format(len(rtable)))
print("Orphans abandoned: {0}".format(len(otable)))
print("Number of clusters: {0}".format(max(rtable.cols.clusterNumber[:])))
# Leftovers are 'orphans' in the repeater table - OPTICS didn't think they belonged to a family
# May need to alter code to better deal with these guys, currently if there is a leftover,
# it is likely to never be re-associated with a family again...
print("Number of leftovers in clustering: {0}".format(len(rtable.get_where_list('clusterNumber == -1'))))

# Plot ordered waveforms and correlation matrix (unordered and ordered)
redpy.plotting.createOrderedWaveformFigure(rtable)
redpy.plotting.createCMatrixFigure(rtable, ctable)

plt.show()

h5file.close()
