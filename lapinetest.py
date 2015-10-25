import numpy as np
import matplotlib.pyplot as plt
import time
import redpy
import datetime

from tables import *
from redpy.optics import *
from obspy import UTCDateTime

# Current issues: alignment, holes in cmatrix (from when two clusters combine?)
# Something's also up with needing to run OPTICS before plotting....?

opt = redpy.config.Options(title="La Pine Test", filename="LPtest.h5", groupName="svic",
        groupDesc="LaPine: SVIC", station="SVIC", network="CC", channel="EHZ",
        winlen=1024, ptrig=20.0, atrig=40.0)

redpy.table.initializeTable(opt) 

h5file = open_file(opt.filename, "a")
rtable = eval('h5file.root.'+ opt.groupName + '.repeaters')
otable = eval('h5file.root.'+ opt.groupName + '.orphans')
ctable = eval('h5file.root.'+ opt.groupName + '.correlation')
jtable = eval('h5file.root.'+ opt.groupName + '.junk')

ttimer = time.time()
tstart = UTCDateTime('2015-10-22 00:00')
nhour = 24

previd = 0
ptime = -2000
lasttable = 0
for hour in range(nhour):

    t = tstart+hour*3600
    print(t)
    st = redpy.trigger.getIRIS(t, opt, nsec=3600)
    alltrigs, ptime = redpy.trigger.trigger(st, ptime, opt)

    # Clean out data spikes etc.
    trigs, junk = redpy.trigger.dataclean(alltrigs,opt,flag=1)

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
            id = previd + i
            redpy.correlation.runCorrelation(rtable, otable, ctable, trigs[i], id, opt)
        previd = id + 1
        
        
        if len(rtable) > lasttable:
            redpy.cluster.runFullOPTICS(rtable, ctable)
            redpy.plotting.createCMatrixFigure(rtable, ctable)
            lasttable = len(rtable)

print("Correlation done in: {:03.2f} seconds".format(time.time()-ttimer))

if len(rtable) > 0:
    # Run clustering one more time before plotting
    redpy.cluster.runFullOPTICS(rtable, ctable)
    
    # Clear out expired orphans
    #redpy.table.clearExpiredOrphans(otable,opt,UTCDateTime('2014-08-25')) # For testing orphan removal, type date manually here (make it at least opt.minorph days later than tstart. For real time running, use commented line below
    #redpy.table.clearExpiredOrphans(otable,opt,tstart+nhour*3600)
    
    # Print some information about the run
    print("Repeaters found: {0}".format(len(rtable)))
    print("Orphans saved: {0}".format(len(otable)))
    print("Number of clusters: {0}".format(max(rtable.cols.clusterNumber[:])))
    # Leftovers are 'orphans' in the repeater table - OPTICS didn't think they belonged to a family
    # May need to alter code to better deal with these guys, currently if there is a leftover,
    # it is likely to never be re-associated with a family again...
    print("Number of leftovers in clustering: {0}".format(len(rtable.get_where_list('clusterNumber == -1'))))
    print("Number of junk triggers: {0}".format(len(jtable)))
    
    print(previd)
    print(ptime)
    
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
