import numpy as np
from tables import *

def getClusters(rtable, cutoff=0.7):

    """
    Cuts the clustering order into clusters, defines orphans as -1

    rtable: Repeater catalog table, with cluster ordering in columns 6 - 8
    cutoff: Minimum coefficient to cut the clusters

    Returns cluster numbers wrt. ordered list.
    """

    order = rtable.cols.order[:] # Ordering
    oreach = rtable.cols.reachability[order] # Ordered reachability
    odist = rtable.cols.coreDistance[order] # Ordered core distance
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

    return oclust

    
def getCenters(rtable, cutoff=0.7):

    """
    Finds the "center" of each cluster (including orphans)
    
    rtable: Repeater catalog table, with clustering order in columns 6 - 8
    cutoff: Minimum coefficient to cut the clusters

    Returns the id numbers of cluster centers and orphans
    
    Orphans may exist even in the repeater table if the cutoff is higher than what
    was used to originally consider them a repeater
    """

    order = rtable.cols.order[:]
    oreach = rtable.cols.reachability[order] # Ordered reachability
    oclust = getClusters(rtable, cutoff) # Clusters

    cluster_id = np.max(oclust).astype(int)
    o = np.array(order)
    centers = np.zeros((cluster_id + 1,)).astype(int)
    for clusternum in range(cluster_id + 1):
        oo = o[oclust == clusternum]
        centers[clusternum] = oo[np.argmin(oreach[oclust == clusternum])]

    orphans = o[oclust == -1]

    return centers, orphans