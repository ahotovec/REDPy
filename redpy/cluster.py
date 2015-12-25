import numpy as np
from tables import *
from redpy.optics import *
from redpy.correlation import *
from redpy.table import *
import time
import copy

def setClusters(rtable, opt):

    """
    Cuts the clustering order into flat clusters, defines orphans as -1

    rtable: Repeater table, with cluster ordering in columns 6 - 8
    opt: Options object describing station/run parameters

    Sets cluster numbers in column 9 of the repeater table, ordered by first event in
    each cluster
    """
    
    cutoff = opt.cmin

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

    cnum = np.zeros((len(oreach),))
    cnum[order] = oclust
    
    # Figure out earliest member in each family
    dt = rtable.cols.startTimeMPL[:]
    mindt = np.zeros((max(cnum)+1,))
    for clustNum in range(int(max(cnum)+1)):
        mindt[clustNum] = min(dt[cnum==clustNum])
    
    n = 0
    clust = np.zeros((len(oreach),))
    for clustNum in np.argsort(mindt):
        clust[cnum==clustNum] = n
        n = n+1
        
    rtable.cols.plotClust[:] = clust
    rtable.cols.clusterNumber[:] = cnum
    rtable.flush()

    
def setCenters(rtable, opt):

    """
    Finds the "center" of each cluster (including orphans, if they exist)
    
    rtable: Repeater table, with clustering order in columns 6 - 8
    opt: Options object describing station/run parameters

    Sets column 10 of the repeater table to 0 if it is not a core, 1 if it is a core,
    and -1 if it is an orphan. Orphans may exist in the repeater table if the cutoff in
    opt is higher than what was used to originally consider them a repeater.
    """
    
    cutoff = opt.cmin

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


def alignAll(rtable, ctable, opt):
    """
    Aligns events in the table that were misaligned. Uses the column 'alignedTo' to guide
    which events are misaligned and skips events which have already been aligned
    """
    
    cores = rtable.get_where_list('(isCore==1)')
    for core in cores:
        
        # Change alignedTo of the core to itself
        calignedTo = copy.copy(rtable[core]['alignedTo'])
        rtable.cols.alignedTo[core] = rtable[core]['id']
            
        # As well as any events that share its alignedTo value
        for f in rtable.get_where_list(
            '(clusterNumber=={0}) & (isCore==0) & (alignedTo=={1})'.format(
            rtable[core]['clusterNumber'], calignedTo)):
            rtable.cols.alignedTo[f] = rtable[core]['id']
        rtable.flush()
        
        # Check for any members with other alignedTo values
        fam = rtable.get_where_list(
            '(clusterNumber=={0}) & (isCore==0) & (alignedTo!={1})'.format(
            rtable[core]['clusterNumber'], rtable[core]['id']))
        if fam.any():
            # Clean up alignedTo values that reference an id that is not itself 
            for u in np.unique(rtable[fam]['alignedTo']):
                unum = rtable.get_where_list('(id=={})'.format(u))[0]
                for f in rtable.get_where_list(
                    '(clusterNumber=={0}) & (isCore==0) & (alignedTo=={1})'.format(
                    rtable[core]['clusterNumber'], u)):
                    rtable.cols.alignedTo[f] = rtable[unum]['alignedTo']
                rtable.flush()
            
            # Loop over remaining unique values
            fftj = rtable[core]['windowFFT']
            coeffj = rtable[core]['windowCoeff']
            for u in np.unique(rtable[fam]['alignedTo']):

                # Align to core, apply same lag to other members
                unum = rtable.get_where_list('(id=={})'.format(u))[0]
                cor, lag = redpy.correlation.xcorr1x1(fftj,
                    rtable.cols.windowFFT[unum], coeffj,
                    rtable.cols.windowCoeff[unum])
                        
                # If doesn't correlate well, try a better event
                if cor < opt.cmin + 0.05:
                        
                    tmp = rtable.get_where_list(
                        '(clusterNumber=={0}) & (alignedTo!={1})'.format(
                        rtable[core]['clusterNumber'], u))
                    utmp = rtable.get_where_list(
                        '(clusterNumber=={0}) & (alignedTo=={1})'.format(
                        rtable[core]['clusterNumber'], u))
                        
                    cormax = 0
                    id1 = -1
                    for t in tmp:
                        for f in utmp:
                            clist = ctable.get_where_list(
                                '(id1=={0}) & (id2=={1})'.format(
                                min(rtable[f]['id'],rtable[t]['id']),
                                max(rtable[f]['id'],rtable[t]['id'])))
                            if clist.any():
                                if ctable[clist[0]]['ccc'] > cormax:
                                    cormax = ctable[clist[0]]['ccc']
                                    id1 = f
                                    id2 = t
                    if id1 != -1:
                        cor1, lag1 = redpy.correlation.xcorr1x1(
                            rtable.cols.windowFFT[id2], rtable.cols.windowFFT[id1],
                            rtable.cols.windowCoeff[id2],
                            rtable.cols.windowCoeff[id1])
                        cor2, lag2 = redpy.correlation.xcorr1x1(
                            rtable.cols.windowFFT[core], rtable.cols.windowFFT[id2],
                            rtable.cols.windowCoeff[core],
                            rtable.cols.windowCoeff[id2])
                        lag = lag1 + lag2
                    else:
                        print("Couldn't find good event to realign with...")
                    
                for f in rtable.get_where_list(
                    '(clusterNumber=={0}) & (isCore==0) & (alignedTo=={1})'.format(
                    rtable[core]['clusterNumber'], u)):
                    rtable.cols.alignedTo[f] = rtable[core]['id']
                    rtable.cols.windowStart[f] = rtable.cols.windowStart[f] - lag
                    rtable.cols.windowCoeff[f], rtable.cols.windowFFT[f] = redpy.correlation.calcWindow(
                        rtable.cols.waveform[f], rtable.cols.windowStart[f], opt)
                rtable.flush()
    
    
def mergeFamilies(rtable, ctable, opt):
    """
    Cross-correlates the cores of each family and attempts to realign and merge them.
    """
    
    cores = rtable.get_where_list('(isCore==1)')    
    c = 0
    for core1 in cores[0:-1]:
        c = c+1
        if (rtable[core1]['isCore'] == 1):
            coeffj = rtable[core1]['windowCoeff']
            fftj = rtable[core1]['windowFFT']
            for core2 in cores[c:]:
                if (rtable[core2]['isCore'] == 1):
                    # Try comparing cores with no lag
                    cor, lag = redpy.correlation.xcorr1x1(fftj, rtable[core2]['windowFFT'],
                        coeffj, rtable[core2]['windowCoeff'])
                    cormax = cor
                    lagmax = lag
                    if cor <= opt.cmin - 0.05:
                        waveform = rtable[core2]['waveform']
                        windowStart = rtable[core2]['windowStart']
                        # Try comparing cores with window moved forward half window
                        coeffi, ffti = redpy.correlation.calcWindow(waveform,
                            windowStart + opt.winlen/2, opt)
                        cor, lag = redpy.correlation.xcorr1x1(fftj, ffti, coeffj,
                            coeffi)
                        lag = lag + opt.winlen/2
                        if cor > cormax:
                            cormax = cor
                            lagmax = lag
                        if cor <= opt.cmin - 0.05 and windowStart > opt.winlen/2:
                            # Try with window moved back half window
                            coeffi, ffti = redpy.correlation.calcWindow(waveform,
                                windowStart - opt.winlen/2, opt)
                            cor, lag = redpy.correlation.xcorr1x1(fftj, ffti, coeffj,
                                coeffi)
                            lag = lag - opt.winlen/2
                            if cor > cormax:
                                cormax = cor
                                lagmax = lag
                    cor = cormax
                    lag = lagmax
                    if cor > opt.cmin - 0.05:
                        f1 = rtable.get_where_list(
                            '(clusterNumber=={})'.format(rtable[core1]['clusterNumber']))
                        f2 = rtable.get_where_list(
                            '(clusterNumber=={})'.format(rtable[core2]['clusterNumber']))
                        # Determine smaller family to adjust
                        if len(f1) >= len(f2):
                            f = f2
                            ff = rtable[f1]
                            fnum = rtable[core1]['clusterNumber']
                            falign = rtable[core1]['id']
                        else:
                            f = f1
                            ff = rtable[f2]
                            lag = -lag    
                            fnum = rtable[core2]['clusterNumber']   
                            falign = rtable[core2]['id'] 
                        # Shove the smaller family over, remove isCore, set clusterNumber
                        # Problem when it tries to move it over and there isn't enough padding
                        for n in f:
                            coeffn, fftn = redpy.correlation.calcWindow(
                                rtable.cols.waveform[n], rtable.cols.windowStart[n] - lag,
                                opt)
                            rtable.cols.windowCoeff[n] = coeffn
                            rtable.cols.windowFFT[n] = fftn
                            rtable.cols.windowStart[n] = rtable.cols.windowStart[n] - lag
                            rtable.cols.clusterNumber[n] = fnum
                            rtable.cols.isCore[n] = 0
                            rtable.cols.alignedTo[n] = falign
                            rtable.flush()
                            corn, lagn = redpy.correlation.xcorr1xtable(coeffn, fftn, ff,
                                opt)
                            for j in range(len(corn)):
                                redpy.table.appendCorrelation(ctable, rtable[n]['id'],
                                    ff[j]['id'], corn[j], opt)
    

def deepClean(rtable, ctable, opt):

    """
    Completely recalculates the correlation table after realignment.
    
    rtable: Repeater table
    ctable: Correlation matrix table
    opt: Options object describing station/run parameters
    
    May take a VERY long time to complete! Run with extreme caution.
    """
    
    alignAll(rtable, ctable, opt)
    mergeFamilies(rtable, ctable, opt)

    cores = rtable.where('(isCore==1)')
    for core in cores:
        fam = rtable.get_where_list('(clusterNumber=={})'.format(core['clusterNumber']))
        if len(fam) > 2:
            
            # Destroy any lines in ctable relating to the members of this family
            cstring = []
            cstring.append('(id1 == {0}) | (id2 == {0})'.format(rtable[fam[0]]['id']))
            for f in range(1,len(fam)):
                cstring.append(' | (id1 == {0}) | (id2 == {0})'.format(
                    rtable[fam[f]]['id']))
            cindex = ctable.get_where_list(''.join(cstring))
            for n in reversed(cindex):
                ctable.remove_row(n)
            
            c = 0
            # Refill the correlation table
            # There is a LOT of overhead associated with this step!
            print('Cross-correlating {} events...'.format(len(fam)))
            for n in fam[0:-1]:                
                c = c+1
                ffti = rtable[n]['windowFFT']
                coeffi = rtable[n]['windowCoeff']
                id = rtable[n]['id']
                # Correlate
                for m in fam[c:]:
                    tx = time.time()
                    cor, lag = redpy.correlation.xcorr1x1(ffti, rtable[m]['windowFFT'],
                        coeffi, rtable[m]['windowCoeff'])
                    txcorr = txcorr + time.time() - tx
                    ts = time.time()
                    redpy.table.appendCorrelation(ctable, id, rtable[m]['id'], cor, opt)
                    tstore = tstore + time.time() - ts

    
def runFullOPTICS(rtable, ctable, opt):
    
    """
    Runs a full, brute-force OPTICS clustering using the correlation values in ctable
    
    rtable: Repeater table
    ctable: Correlation matrix table
    opt: Options object describing station/run parameters
    
    Sets the order, reachability, and coreDistance columns in rtable
    """
            
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
    setClusters(rtable, opt)
    setCenters(rtable, opt)
    alignAll(rtable, ctable, opt)

    
