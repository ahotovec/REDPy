import numpy as np
from tables import *
from redpy.optics import *
from redpy.correlation import *
import time

def setClusters(rtable, cutoff=0.7):

    """
    Cuts the clustering order into flat clusters, defines orphans as -1

    rtable: Repeater table, with cluster ordering in columns 6 - 8
    cutoff: Minimum coefficient to cut the clusters

    Sets cluster numbers in column 9 of the repeater table
    """

    order = rtable.cols.order[:] # Ordering
    oreach = rtable.cols.reachability[:] # Ordered reachability
    oreach = oreach[order]
    odist = rtable.cols.coreDistance[:]
    odist = odist[order] # Ordered core distance
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

    clust = np.zeros((len(oreach),))
    clust[order] = oclust

    rtable.cols.clusterNumber[:] = clust
    rtable.flush()

    
def setCenters(rtable, cutoff=0.7):

    """
    Finds the "center" of each cluster (including orphans, if they exist)
    
    rtable: Repeater table, with clustering order in columns 6 - 8
    cutoff: Minimum coefficient to cut the clusters

    Sets column 10 of the repeater table to 0 if it is not a core, 1 if it is a core,
    and -1 if it is an orphan. Orphans may exist in the repeater table if the cutoff is
    higher than what was used to originally consider them a repeater.
    """

    order = rtable.cols.order[:]
    oo = np.sort(order) # Unsorted row position
    
    reach = rtable.cols.reachability[:]
    clust = rtable.cols.clusterNumber[:]

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
    rtable.flush()
    

def checkCores(rtable, ctable, opt):
    """
    Checks to make sure cores are correlated with each other
    
    rtable: Repeater table
    ctable: Correlation matrix table
    
    Sets appropriate correlation values in ctable. I think this may
    be very computationally expensive to constantly check to make sure
    the values are filled...
    """
    
    cores = rtable.get_where_list('(isCore==1)')
    n = 0
    if cores.any():
        for core1 in cores[0:-2]:
            cid1 = rtable.cols.id[core1]
            n = n+1
            for core2 in cores[n:-1]:
                cid2 = rtable.cols.id[core2]
                
                # This line is super slow...
                clist = ctable.get_where_list('(id1 == {0}) & (id2 == {1})'.format(
                    np.min([cid1,cid2]),np.max([cid1,cid2])))
                    
                if not clist.any():
                    cor, lag = redpy.correlation.xcorr1x1(rtable.cols.windowFFT[core1],
                        rtable.cols.windowFFT[core2], rtable.cols.windowCoeff[core1],
                        rtable.cols.windowCoeff[core2])
                    redpy.table.appendCoreCorrelation(ctable, cid1, cid2, cor, opt)

    
def runFullOPTICS(rtable, ctable, opt):
    
    """
    Runs a full, brute-force OPTICS clustering using the correlation values in ctable
    
    rtable: Repeater table
    ctable: Correlation matrix table
    
    Sets the order, reachability, and coreDistance columns in rtable
    """
    t = time.time()
    
    # May move this over to a 'cleanup' function so it isn't run every run of OPTICS    
    checkCores(rtable, ctable, opt)
    print("Time spent checking cores: {:03.2f} seconds".format(time.time()-t))
    
    C = np.zeros((len(rtable),len(rtable)))
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]
    
    # Convert id to row
    rtable_ids = rtable.cols.id[:]
    r = np.zeros((max(rtable_ids)+1,)).astype('int')
    r[rtable_ids] = range(len(rtable_ids))
    C[r[id1], r[id2]] = ctable.cols.ccc[:]
            
    C = C + C.T + np.eye(len(C))
    
    # Cluster with OPTICS
    ttree = setOfObjects(1-C)
    prep_optics(ttree,1)
    build_optics(ttree,1)
    order = ttree._ordered_list

    # Save the ordering to the repeater table
    rtable.cols.order[:] = order
    rtable.cols.reachability[:] = ttree._reachability
    rtable.cols.coreDistance[:] = ttree._core_dist
    rtable.flush()
        
    # Update the clusters and cores, too!
    setClusters(rtable)
    setCenters(rtable)
    
    print("Total time spent clustering: {:03.2f} seconds".format(time.time()-t))
    
    if len(rtable) > 0:
        print("Repeaters found: {0}".format(len(rtable)))
        print("Number of clusters: {0}".format(max(rtable.cols.clusterNumber[:])+1))
        bigfam = 0
        for n in range(max(rtable.cols.clusterNumber[:])+1):
            tmp = len(rtable.get_where_list('(isCore == 0) & (clusterNumber == {})'.format(n)))
            if tmp > bigfam:
               bigfam = tmp
        print("Members in largest cluster: {0}".format(bigfam+1))
        print("Number of leftovers in clustering: {0}".format(len(rtable.get_where_list('clusterNumber == -1'))))
