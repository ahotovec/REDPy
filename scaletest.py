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

opt = redpy.config.Options(filename="testsep.h5", cmin=0.7,
    station="SEP", network="CC", samprate=50.0, winlen=256)

redpy.table.initializeTable(opt) 

h5file = open_file(opt.filename, "a")
rtable = eval('h5file.root.'+ opt.groupName + '.repeaters')
otable = eval('h5file.root.'+ opt.groupName + '.orphans')
ctable = eval('h5file.root.'+ opt.groupName + '.correlation')
jtable = eval('h5file.root.'+ opt.groupName + '.junk')

ttimer = time.time()
tstart = UTCDateTime('2014-08-12')
nhour = int(1*24)

previd = 0
for hour in range(nhour):

    t = tstart+hour*3600
    print(t)
    st = redpy.trigger.getIRIS(t, opt, nsec=3600)
    alltrigs, maxratio = redpy.trigger.trigger(st, opt)

    # Clean out data spikes etc.
    trigs, junk = redpy.trigger.dataclean(alltrigs,opt,flag=1)

    #save junk triggers in separate table for quality checking purposes
    for i in range(len(junk)):
        redpy.table.populateJunk(jtable,junk[i],0,opt)

    if len(trigs) > 0:
        
        ostart = 0
        if len(otable) == 0:
            # First trigger goes to orphans table
            redpy.table.populateOrphan(otable, 0, trigs[0], 'Never', opt)
            ostart = 1
        
        # Loop through remaining triggers
        for i in range(ostart,len(trigs)):  
            id = previd + i
            redpy.correlation.runCorrelation(rtable, otable, ctable, trigs[i], id, opt)
        previd = id + 1



print("Correlation done in: {:03.2f} seconds".format(time.time()-ttimer))

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
print("Number of junk triggers: {0}".format(len(jtable)))

# Plot ordered waveforms and correlation matrix (unordered and ordered)
redpy.plotting.createOrderedWaveformFigure(rtable, opt)
redpy.plotting.createCMatrixFigure(rtable, ctable)
# Plot junk events
redpy.plotting.createWigglePlot(jtable, opt)
plt.show()

h5file.close()
