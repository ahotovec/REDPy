import numpy as np
from tables import *
from redpy.optics import *

def setClusters(rtable, cutoff=0.7):

    """
    Cuts the clustering order into flat clusters, defines orphans as -1

    rtable: Repeater catalog table, with cluster ordering in columns 6 - 8
    cutoff: Minimum coefficient to cut the clusters

    Sets cluster numbers in column 9 of the repeater table
    """

    order = rtable.cols.order[:] # Ordering
    oreach = rtable.cols.reachability[:] # Ordered reachability
    oreach = oreach[order]
    odist = rtable.cols.coreDistance[:]
    odist = odist[order] # Ordered core distance
    cluster_id = -1
    
    oo = order[order] # This is equivalent to the unsorted rows, needed to convert back
                      # to real position in the saved table ordering

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

    rtable.cols.clusterNumber[:] = oclust[oo[order]]
    rtable.flush()

    
def setCenters(rtable, cutoff=0.7):

    """
    Finds the "center" of each cluster (including orphans, if they exist)
    
    rtable: Repeater catalog table, with clustering order in columns 6 - 8
    cutoff: Minimum coefficient to cut the clusters

    Sets column 10 of the repeater table to 0 if it is not a core, 1 if it is a core,
    and -1 if it is an orphan. Orphans may exist in the repeater table if the cutoff is
    higher than what was used to originally consider them a repeater.
    """

    order = rtable.cols.order[:]
    oreach = rtable.cols.reachability[:]
    oreach = oreach[order]
    oclust = rtable.cols.clusterNumber[:]
    oclust = oclust[order]

    cluster_id = np.max(oclust).astype(int)
    ocenters = np.zeros((cluster_id + 1,)).astype(int)
    for clusternum in range(cluster_id + 1):
        clustermembers = order[oclust == clusternum]
        ocenters[clusternum] = clustermembers[np.argmin(oreach[oclust == clusternum])]
        
    oorphans = order[oclust == -1]
    
    oo = order[order] # This is equivalent to the unsorted rows, needed to convert back
                      # to real position in the saved table ordering

    cores = np.zeros((len(order),)).astype(int) # Reset everything to 0
    cores[oo[ocenters]] = np.ones((len(ocenters),)).astype(int)
    cores[oo[oorphans]] = -1*np.ones((len(oorphans),)).astype(int)
    
    rtable.cols.isCore[:] = cores
    rtable.flush()
    
    
def runFullOPTICS(rtable, ctable):
    
    """
    Runs a full, brute-force OPTICS clustering using the correlation values in ctable
    
    rtable: Repeater catalog table
    ctable: Correlation matrix table
    
    Sets the order, reachability, and coreDistance columns in rtable
    """
    
    C = np.zeros((len(rtable),len(rtable)))
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]

    # Convert id to row
    rtable_ids = rtable.cols.id[:]
    for i in range(len(id1)):
        C[np.where(rtable_ids == id1[i])[0][0],
            np.where(rtable_ids == id2[i])[0][0]] = ctable.cols.ccc[i]
            
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
    
    # Update the clusters and cores, too!
    setClusters(rtable)
    setCenters(rtable)
    rtable.flush()