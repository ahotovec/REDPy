import numpy as np
from tables import *
from redpy.optics import *
    
def runFamOPTICS(rtable, ctable, ftable, fnum, opt):
    
    """
    Runs OPTICS ordering within a single family
    
    rtable: Repeater table
    ctable: Correlation matrix table
    ftable: Families table
    fnum: Family number to run
    opt: Options object describing station/run parameters
    
    Returns slightly different ordering than full version, but probably better
    """
    
    fam = np.fromstring(ftable[fnum]['members'], dtype=int, sep=' ')
    
    if len(fam) in (3, 4, 5, 6, 10, 15, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000,
        25000, 50000, 100000, 250000, 500000):
    
        # Could be sped up if these three don't have to be called every time
        id1 = ctable.cols.id1[:]
        id2 = ctable.cols.id2[:]
        ccc = 1-ctable.cols.ccc[:]
        
        # Create distance matrix
        ids = rtable[fam]['id']
        ix = np.where(np.in1d(id2,ids))
        r = np.zeros((max(ids)+1,)).astype('int')
        r[ids] = range(len(ids))
        D = np.ones((len(ids),len(ids)))
        D[r[id2[ix]],r[id1[ix]]] = ccc[ix]
        D[r[id1[ix]],r[id2[ix]]] = ccc[ix]
        D[range(len(ids)),range(len(ids))] = 0
    
        # Sort so most connected event is always considered for core
        s = np.argsort(sum(D))[::-1]
        D = D[s,:]
        D = D[:,s]
        fam = fam[s]
    
        # Run OPTICS
        ttree = setOfObjects(D)
        prep_optics(ttree,1)
        build_optics(ttree,1)
        order = np.array(ttree._ordered_list)
        core = fam[np.argmin(ttree._reachability)]
    
        # Write to ftable
        np.set_printoptions(threshold=np.nan)
        np.set_printoptions(linewidth=np.nan)
        ftable.cols.members[fnum] = np.array2string(fam[order])[1:-1]
        ftable.cols.core[fnum] = core
        ftable.cols.startTime[fnum] = np.min(rtable[fam]['startTimeMPL'])
        ftable.cols.printme[fnum] = 1
        ftable.flush()
        