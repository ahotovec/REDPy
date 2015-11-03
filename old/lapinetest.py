import numpy as np
import matplotlib.pyplot as plt
import redpy
import datetime

from tables import *
from redpy.optics import *
from obspy import UTCDateTime

import warnings
warnings.filterwarnings("ignore")

# Changed options to read a configuration file
opt = redpy.config.Options('lapine.cfg')

redpy.table.initializeTable(opt) 

h5file, rtable, otable, ctable, jtable = redpy.table.openTable(opt)

tstart = UTCDateTime('2015-10-10 00:00')
nhour = 24*21

for hour in range(nhour):

    t = tstart+hour*opt.nsec
    print(t)
    st = redpy.trigger.getData(t, opt)
    alltrigs = redpy.trigger.trigger(st, rtable, opt)

    # Clean out data spikes etc.
    trigs, junk = redpy.trigger.dataclean(alltrigs, opt, flag=1)

    #save junk triggers in separate table for quality checking purposes
    for i in range(len(junk)):
        redpy.table.populateJunk(jtable,junk[i],0,opt)

    if len(trigs) > 0:
        
        ostart = 0
        if len(otable) == 0:
            # First trigger goes to orphans table
            redpy.table.populateOrphan(otable, 0, trigs[0], opt)
            ostart = 1
        
        # Loop through remaining triggers
        for i in range(ostart,len(trigs)):  
            id = rtable.attrs.previd + i
            redpy.correlation.runCorrelation(rtable, otable, ctable, trigs[i], id, opt)
        rtable.attrs.previd = id + 1
    
    redpy.table.clearExpiredOrphans(otable, opt, tstart+hour*3600)
    print("Length of Orphan table: {0}".format(len(otable)))
 
if len(rtable) > 0:
    # Run clustering one more time before plotting
    redpy.cluster.runFullOPTICS(rtable, ctable)
        
    # Print some information about the run
    print("Repeaters found: {0}".format(len(rtable)))
    print("Orphans saved: {0}".format(len(otable)))
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
    if len(jtable) > 0:
        redpy.plotting.createWigglePlot(jtable, opt)
    
    plt.show()
else:
    print('Nothing to cluster! Exiting...')

h5file.close()