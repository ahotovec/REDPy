import numpy as np
from tables import *
from redpy.optics import *
from redpy.correlation import *
from redpy.table import *
import time
import copy

def setClusters(rtable, ftable, order, reach, dist, opt):

    """
    Cuts the clustering order into flat clusters, defines orphans as -1

    rtable: Repeater table, with cluster ordering in columns 6 - 8
    opt: Options object describing station/run parameters

    Sets cluster numbers in column 9 of the repeater table, ordered by first event in
    each cluster
    """
    
    cutoff = opt.cmin
    oreach = reach[order]
    odist = dist[order]
    cluster_id = -1

    oclust = np.zeros((len(oreach),))
    for x in range(len(oreach)):
        if oreach[x] > 1 - cutoff:
            if odist[x] <= 1 - cutoff:
                cluster_id += 1
                oclust[x] = cluster_id
            else:
                oclust[x] = -1 # orphan
        else:
            oclust[x] = cluster_id

    cnum = np.zeros((len(oreach),), dtype=np.int32)
    cnum[order] = oclust
    
    # Figure out earliest member in each family
    dt = rtable.cols.startTimeMPL[:]
    mindt = np.zeros((max(cnum)+1,))
    for clustNum in range(int(max(cnum)+1)):
        mindt[clustNum] = min(dt[cnum==clustNum])
    
    n = 0
    clust = np.zeros((len(oreach),), dtype=np.int32)
    for clustNum in np.argsort(mindt):
        clust[cnum==clustNum] = n
        n = n+1
        
    rtable.cols.plotClust[:] = clust
    rtable.cols.clusterNumber[:] = cnum
    rtable.flush()
    
    # Write to ftable
    
    # First check print options so all members are printed
    np.set_printoptions(threshold=np.nan)
    np.set_printoptions(linewidth=np.nan)
    
    # Check length of ftable compared to max(cnum)
    if len(ftable) <= max(cnum):
        f = ftable.row
        f['pnum'] = -1
        f['members'] = ''
        f.append()
        ftable.flush()
        
    # Populate ftable
    n = 0
    for c in np.argsort(mindt):
        ftable.cols.pnum[c] = n
        ftable.cols.members[c] = np.array2string(np.where(cnum==c)[0])[1:-1]
        n = n+1
    ftable.flush()
    
    return cnum

    
def setCenters(rtable, ftable, order, reach, clust, opt):

    """
    Finds the "center" of each cluster (including orphans, if they exist)
    
    rtable: Repeater table, with clustering order in columns 6 - 8
    opt: Options object describing station/run parameters

    Sets column 10 of the repeater table to 0 if it is not a core, 1 if it is a core,
    and -1 if it is an orphan. Orphans may exist in the repeater table if the cutoff in
    opt is higher than what was used to originally consider them a repeater.
    """
    
    cutoff = opt.cmin
    oo = np.sort(order) # Unsorted row position
    
    cluster_id = np.max(clust).astype(int)
    centers = np.zeros((cluster_id + 1,)).astype(int)
    for clusternum in range(cluster_id + 1):
        clustermembers = oo[clust == clusternum]
        centers[clusternum] = clustermembers[np.argmin(reach[clustermembers])]    
    orphans = order[clust == -1]  

    cores = np.zeros((len(order),)).astype(int) # Reset everything to 0
    cores[centers] = np.ones((len(centers),)).astype(int)
    cores[orphans] = -1*np.ones((len(orphans),)).astype(int)
    
    rtable.cols.isCore[:] = cores
    ftable.attrs.cores = sorted(centers) # Sort is somehow important?
    rtable.flush()
              
   
def alignAll(rtable, ctable, ftable, clusterNumber, opt):
    """
    Aligns events in the table that were misaligned. Uses the column 'alignedTo' to guide
    which events are misaligned and skips events which have already been aligned
    """

    cores = ftable.attrs.cores
    id = rtable.cols.id[:]
    alignedTo = rtable.cols.alignedTo[:]
    
    cid1 = ctable.cols.id1[:]
    cid2 = ctable.cols.id2[:]
    C = np.eye(len(rtable))
    rtable_ids = rtable.cols.id[:]
    r = np.zeros((max(rtable_ids)+1,)).astype('int')
    r[rtable_ids] = range(len(rtable_ids))
    C[r[cid1], r[cid2]] = ctable.cols.ccc[:]
    C[r[cid2], r[cid1]] = ctable.cols.ccc[:]
    
    for core in cores:

        members = np.where(clusterNumber == clusterNumber[core])
        
        # Change alignedTo of the core to itself, as well as any events that share its
        # former alignedTo value
        calignedTo = copy.copy(alignedTo[core])
        alignedTo[np.intersect1d(members, np.where(alignedTo == calignedTo))] = id[core]
            
        # Check for any members with other alignedTo values
        fam = np.intersect1d(members, np.where(alignedTo != id[core]))
        if fam.any():
            # Clean up alignedTo values that reference an id that is not itself
            for u in np.unique(alignedTo[fam]):
                unum = np.where(id == u)[0][0]
                alignedTo[np.intersect1d(members, np.where(
                    alignedTo == u))] = rtable[unum]['alignedTo']
             
            # Loop over remaining unique values
            fftj = rtable[core]['windowFFT']
            coeffj = rtable[core]['windowCoeff']
            
            for u in np.unique(alignedTo[fam]):
                # Align to core, apply same lag to other members
                unum = np.where(id == u)[0][0]
                cor, lag = redpy.correlation.xcorr1x1(fftj, rtable[unum]['windowFFT'],
                    coeffj, rtable[unum]['windowCoeff'])
                        
                # If doesn't correlate well, try a better event
                if cor < opt.cmin + 0.05:
                    tmp = np.intersect1d(members, np.where(alignedTo != u))
                    utmp = np.intersect1d(members, np.where(alignedTo == u))
                    
                    Cslice = C[tmp,:]
                    Cslice = Cslice[:,utmp]
                    [t,f] = np.unravel_index(np.argmax(Cslice), np.shape(Cslice))
                    id1 = utmp[f]
                    id2 = tmp[t]
                    
                    cor1, lag1 = redpy.correlation.xcorr1x1(
                        rtable[id2]['windowFFT'], rtable[id1]['windowFFT'],
                        rtable[id2]['windowCoeff'], rtable[id1]['windowCoeff'])
                    cor2, lag2 = redpy.correlation.xcorr1x1(
                        rtable[core]['windowFFT'], rtable[id2]['windowFFT'],
                        rtable[core]['windowCoeff'], rtable[id2]['windowCoeff'])
                    lag = lag1 + lag2
                
                for f in np.intersect1d(members, np.where(alignedTo == u)):
                    alignedTo[f] = id[core]
                    rtable.cols.windowStart[f] = rtable.cols.windowStart[f] - lag
                    rtable.cols.windowCoeff[f], rtable.cols.windowFFT[f] = redpy.correlation.calcWindow(
                        rtable.cols.waveform[f], rtable.cols.windowStart[f], opt)
                rtable.flush()
    rtable.cols.alignedTo[:] = alignedTo
    rtable.flush()

    
def runFullOPTICS(rtable, ctable, ftable, opt):
    
    """
    Runs a full, brute-force OPTICS clustering using the correlation values in ctable
    
    rtable: Repeater table
    ctable: Correlation matrix table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Sets the order column in rtable
    """
     
    t = time.time()       
    C = np.ones((len(rtable),len(rtable)))
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]
    ccc = 1-ctable.cols.ccc[:]
    
    # Convert id to row
    rtable_ids = rtable.cols.id[:]
    r = np.zeros((max(rtable_ids)+1,)).astype('int')
    r[rtable_ids] = range(len(rtable_ids))
    C[r[id1], r[id2]] = ccc
    C[r[id2], r[id1]] = ccc
    C[range(len(rtable)),range(len(rtable))] = 0
    
    # Cluster with OPTICS
    ttree = setOfObjects(C)
    prep_optics(ttree,1)
    build_optics(ttree,1)
    order = np.array(ttree._ordered_list)

    # Save the ordering to the repeater table
    rtable.cols.order[:] = order
    rtable.flush()
        
    # Update the clusters and cores, too!
    cnum = setClusters(rtable, ftable, order, ttree._reachability, ttree._core_dist, opt)
    setCenters(rtable, ftable, order, ttree._reachability, cnum, opt)
    alignAll(rtable, ctable, ftable, cnum, opt)

    
