import numpy as np
from tables import *

def getClusters(ctable, cutoff=0.7):

    """
    Cuts the clustering order into clusters, defines orphans as -1

    ctable: Correlation table, with clustering order in attributes
    cutoff: Minimum coefficient to cut the clusters

    Returns cluster numbers wrt. ordered list.
    """

    order = np.array(ctable.attrs.order)
    oreach = ctable.attrs.reachability[order]
    odist = ctable.attrs.coredist[order]
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

    
def getCenters(ctable, cutoff=0.7):

    """
    Finds the "center" of each cluster (including orphans)
    
    ctable: Correlation table, with clustering order in attributes
    cutoff: Minimum coefficient to cut the clusters

    Returns the id numbers of cluster centers and orphans
    """

    order = np.array(ctable.attrs.order)
    oreach = ctable.attrs.reachability[order]
    oclust = getClusters(ctable, cutoff)

    cluster_id = np.max(oclust).astype(int)
    o = np.array(order)
    centers = np.zeros((cluster_id + 1,)).astype(int)
    for clusternum in range(cluster_id + 1):
        oo = o[oclust == clusternum]
        centers[clusternum] = oo[np.argmin(oreach[oclust == clusternum])]

    orphans = o[oclust == -1]

    return centers, orphans