import numpy as np
from tables import *
from redpy.optics import *
from redpy.correlation import *
from redpy.table import *
import time
import copy

def setClusters(rtable, ftable, order, reach, dist, opt):

    """
    Cuts the clustering order into flat clusters

    rtable: Repeater table
    ftable: Families table
    order: Ordering from OPTICS
    reach: Reachability from OPTICS
    dist: Correlation distance from OPTICS
    opt: Options object describing station/run parameters

    Sets cluster numbers and members in ftable, row ordered by first event in
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

    # First check print options so all members are printed
    np.set_printoptions(threshold=np.nan)
    np.set_printoptions(linewidth=np.nan)
    
    # Check length of ftable compared to max(cnum)
    if len(ftable) <= max(clust):
        f = ftable.row
        f['members'] = ''
        f['printme'] = 0
        f.append()
        ftable.flush()
        
    # Populate ftable
    n = 0
    for c in np.argsort(mindt):
        if ftable.cols.members[n] != np.array2string(order[oclust==c])[1:-1]:
            ftable.cols.members[n] = np.array2string(order[oclust==c])[1:-1]
            ftable.cols.printme[n] = 1
        n = n+1
    ftable.flush()
    
    return clust

    
def setCenters(rtable, ftable, order, reach, clust, opt):

    """
    Finds the "center" of each cluster (including orphans, if they exist)
    
    rtable: Repeater table
    ftable: Families table
    order: Ordering from OPTICS
    reach: Reachability from OPTICS
    clust: Cluster number from setClusters
    opt: Options object describing station/run parameters

    Sets cluster attributes of ftable to locations of cores
    """
    
    cutoff = opt.cmin
    oo = np.sort(order) # Unsorted row position
    
    cluster_id = np.max(clust).astype('uint32')
    centers = np.zeros((cluster_id + 1,)).astype('uint32')
    for clusternum in range(cluster_id + 1):
        clustermembers = oo[clust == clusternum]
        centers[clusternum] = clustermembers[np.argmin(reach[clustermembers])]    

    ftable.attrs.prevcores = np.array(copy.copy(ftable.attrs.cores)).astype('uint32')
    ftable.attrs.cores = np.array(centers).astype('uint32')
    ftable.attrs.nClust = len(centers)
    ftable.flush()
              
   
def alignAll(rtable, ctable, ftable, clusterNumber, opt):

    """
    Aligns events in the table that were misaligned.
    
    rtable: Repeater table
    ctable: Correlation table
    ftable: Families table
    clusterNumber: Cluster number from OPTICS
    opt: Options object describing station/run parameters
    
    Uses the column 'alignedTo' to guide  which events are misaligned and skips events
    which have already been aligned
    """

    cores = ftable.attrs.cores
    id = rtable.cols.id[:]
    alignedTo = rtable.cols.alignedTo[:]
    alignedToF = copy.copy(alignedTo)
    
    for core in cores:
        members = np.where(clusterNumber == clusterNumber[core])
        alignedToF[members[0]] = id[core]
        if alignedTo[core] != id[core]:
            calignedTo = copy.copy(alignedTo[core])
            alignedTo[np.intersect1d(members, np.where(
                alignedTo == calignedTo))] = id[core]
    
    notAligned = np.where(alignedTo != alignedToF)[0]
    if notAligned.any():
        for n in notAligned:
            alignedTo[n] = alignedTo[np.where(id == alignedTo[n])[0]]
        
        notAligned = np.where(alignedTo != alignedToF)[0]
        if notAligned.any():
            cid1 = ctable.cols.id1[:]
            cid2 = ctable.cols.id2[:]
            ccc = ctable.cols.ccc[:]
            C = np.eye(len(rtable))
            rtable_ids = id
            r = np.zeros((max(rtable_ids)+1,)).astype('uint32')
            r[rtable_ids] = range(len(rtable_ids))
            C[r[cid1], r[cid2]] = ccc
            C[r[cid2], r[cid1]] = ccc
        
            clustNA = clusterNumber[notAligned]
            idNA = id[notAligned]
            alignedNA = alignedTo[notAligned]
        
            for c in np.unique(clustNA):
        
                core = np.intersect1d(np.where(clusterNumber == c), cores)[0]
                fftj = rtable[core]['windowFFT']
                coeffj = rtable[core]['windowCoeff']
            
                for u in np.unique(alignedNA[np.where(clustNA == c)[0]]):
                
                    members = np.where(clusterNumber == c)
                    utmp = np.intersect1d(members, np.where(alignedTo == u))
                    unum = notAligned[np.where(idNA == u)[0]][0]
                    cor, lag = redpy.correlation.xcorr1x1(fftj, rtable[unum]['windowFFT'],
                        coeffj, rtable[unum]['windowCoeff'])
                
                    # If doesn't correlate well, try a better event
                    if cor < opt.cmin + 0.05:
                        
                        tmp = np.intersect1d(members, np.where(alignedTo == id[core]))
                    
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
                
                    for f in utmp:
                        rtable.cols.windowStart[f] = rtable.cols.windowStart[f] - lag
                        rtable.cols.windowCoeff[f], rtable.cols.windowFFT[f] = redpy.correlation.calcWindow(
                            rtable.cols.waveform[f], rtable.cols.windowStart[f], opt)
                    rtable.flush()
    
            rtable.cols.alignedTo[:] = alignedToF
            rtable.flush()


def mergeCores(rtable, ctable, ftable, opt):

    """
    Attempts to merge families by ensuring their cores are compared to each other
    
    rtable: Repeater table
    ctable: Correlation matrix table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Compares current cores together and appends any cores that are good matches to
    ctable prior to reclustering.
    """
    
    # Exclude reprocessing core pairs from the last run
    newcores = np.array(np.intersect1d(ftable.attrs.cores,
        np.setxor1d(ftable.attrs.prevcores, ftable.attrs.cores))).astype('uint32')
    cores = np.array(ftable.attrs.cores).astype('uint32')
    for n in newcores:
        for m in cores[np.where(cores>n)[0]]:
            cor, lag = redpy.correlation.xcorr1x1(rtable[n]['windowFFT'],
                rtable[m]['windowFFT'], rtable[n]['windowCoeff'],
                rtable[m]['windowCoeff'])
            if cor >= opt.cmin:
                redpy.table.appendCorrelation(ctable, rtable[n]['id'],
                    rtable[m]['id'], cor, opt)
    
    
def runFullOPTICS(rtable, ctable, ftable, opt):
    
    """
    Runs a full, brute-force OPTICS clustering using the correlation values in ctable
    
    rtable: Repeater table
    ctable: Correlation matrix table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Also updates clusters and cores, then runs alignAll
    """
    
    mergeCores(rtable, ctable, ftable, opt)
         
    C = np.ones((len(rtable),len(rtable)))
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]
    ccc = 1-ctable.cols.ccc[:]
    
    # Convert id to row
    rtable_ids = rtable.cols.id[:]
    r = np.zeros((max(rtable_ids)+1,)).astype('uint32')
    r[rtable_ids] = range(len(rtable_ids))
    C[r[id1], r[id2]] = ccc
    C[r[id2], r[id1]] = ccc
    C[range(len(rtable)),range(len(rtable))] = 0
    
    # Cluster with OPTICS
    ttree = setOfObjects(C)
    prep_optics(ttree,1)
    build_optics(ttree,1)
    order = np.array(ttree._ordered_list)
        
    # Update the clusters and cores, too!
    cnum = setClusters(rtable, ftable, order, ttree._reachability, ttree._core_dist, opt)
    setCenters(rtable, ftable, order, ttree._reachability, cnum, opt)
    alignAll(rtable, ctable, ftable, cnum, opt)

    
