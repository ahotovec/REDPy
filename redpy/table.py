from tables import *
from obspy.core.trace import Trace
from obspy import UTCDateTime
import datetime
import numpy as np
import redpy.correlation
import datetime
import matplotlib

def Repeaters(opt):

    """
    Defines the columns in the 'Repeater Catalog' table based on the Options in opt
    
    id: unique ID number for the event (integer)
    startTime: UTC time of start of the waveform (string)
    startTimeMPL: matplotlib number associated with time (float)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    windowCoeff: amplitude scaling for cross-correlation (float)
    windowFFT: Fourier transform of window (complex ndarray)
    windowAmp: amplitude in first half of window (float)
    
    Returns a dictionary defining the table
    """
    
    dict = {
        "id"            : Int32Col(shape=(), pos=0),
        "startTime"     : StringCol(itemsize=32, pos=1),
        "startTimeMPL"  : Float64Col(shape=(), pos=2),
        "waveform"      : Float64Col(shape=(opt.wshape*opt.nsta,), pos=3),
        "windowStart"   : Int32Col(shape=(), pos=4),
        "windowCoeff"   : Float64Col(shape=(opt.nsta,), pos=5),
        "windowFFT"     : ComplexCol(shape=(opt.winlen*opt.nsta,), itemsize=16, pos=6),
        "windowAmp"     : Float64Col(shape=(opt.nsta,), pos=7),
        }
    
    return dict


def Orphans(opt):

    """
    Defines the columns in the 'Orphans' table based on the Options in opt

    id: unique ID number for the event (integer)
    startTime: UTC time of start of the waveform (string)
    startTimeMPL: matplotlib number associated with time (float)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    windowCoeff: amplitude scaling for cross-correlation (float)
    windowFFT: Fourier transform of window (complex ndarray)
    windowAmp: amplitude in first half of window (float)
    expires: UTC time of when orphan should no longer be considered (string)
    
    Returns a dictionary defining the table
    """
    
    dict = {
        "id"          : Int32Col(shape=(), pos=0),
        "startTime"   : StringCol(itemsize=32, pos=1),
        "startTimeMPL": Float64Col(shape=(), pos=2),
        "waveform"    : Float64Col(shape=(opt.wshape*opt.nsta,), pos=3),
        "windowStart" : Int32Col(shape=(), pos=4),
        "windowCoeff" : Float64Col(shape=(opt.nsta,), pos=5),
        "windowFFT"   : ComplexCol(shape=(opt.winlen*opt.nsta,), itemsize=16, pos=6),
        "windowAmp"   : Float64Col(shape=(opt.nsta,), pos=7),
        "expires"     : StringCol(itemsize=32, pos=8)
        }

    return dict


def Deleted(opt):

    """
    Defines the columns in the 'Deleted' table based on the Options in opt

    id: unique ID number for the event (integer)
    startTime: UTC time of start of the waveform (string)
    startTimeMPL: matplotlib number associated with time (float)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    windowCoeff: amplitude scaling for cross-correlation (float)
    windowFFT: Fourier transform of window (complex ndarray)
    windowAmp: amplitude in first half of window (float)
    
    Returns a dictionary defining the table
    """
    
    dict = {
        "id"          : Int32Col(shape=(), pos=0),
        "startTime"   : StringCol(itemsize=32, pos=1),
        "startTimeMPL": Float64Col(shape=(), pos=2),
        "waveform"    : Float64Col(shape=(opt.wshape*opt.nsta,), pos=3),
        "windowStart" : Int32Col(shape=(), pos=4),
        "windowCoeff" : Float64Col(shape=(opt.nsta,), pos=5),
        "windowFFT"   : ComplexCol(shape=(opt.winlen*opt.nsta,), itemsize=16, pos=6),
        "windowAmp"   : Float64Col(shape=(opt.nsta,), pos=7)
        }

    return dict


def Junk(opt):
    
    """
    Defines the columns in the 'Junk' table, a holding tank for testing suspect events
    
    startTime: UTC time of start of the waveform (string)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    isjunk: Logic holder (integer)
    
    Returns a dictionary defining the table
    """
    
    dict = {
        "startTime"   : StringCol(itemsize=32, pos=1),
        "waveform"    : Float64Col(shape=(opt.wshape*opt.nsta,), pos=2),
        "windowStart" : Int32Col(shape=(), pos=3),
        "isjunk"      : Int32Col(shape=(), pos=0)
        }
        
    return dict


def Correlation(opt):

    """
    Defines the columns in the 'Correlation' table

    id1: unique ID number for the first event (integer)
    id2: unique ID number for the second event (integer)
    ccc: cross-correlation coefficient between those two events (float)
    
    Returns a dictionary defining the table
    """
    
    dict = {
        "id1" : Int32Col(shape=(), pos=0),
        "id2" : Int32Col(shape=(), pos=1),
        "ccc" : Float64Col(shape=(), pos=2)
    }

    return dict


def Families(opt):

    """
    Defines the columns in the 'Families' table

    members: rows in rtable that contain members of the family as ordered string (string)
    core: row in rtable that corresponds to current core event (int)
    startTime: MPL datetime denoting the time of the first event in the family (float)
    printme: describes whether the family has been updated since last printing (int)
    
    The members column is a string so that it can be of completely arbitrary length.
    Well, not completely arbitrary. The itemsize set here (1000000) is a guess at how
    big the string might get for really large families. One hopes to never encounter a
    family this big... (100000+ members?)
        
    Returns a dictionary defining the table
    """
    # COULD ADD MORE METRICS HERE SO IT DOESN'T HAVE TO CALCULATE EACH TIME?
    
    dict = {
        "members"   : StringCol(itemsize=1000000, shape=(), pos=0),
        "core"      : Int32Col(shape=(), pos=1),
        "startTime" : Float64Col(shape=(), pos=2),
        "printme"   : Int32Col(shape=(), pos=3)
    }

    return dict

    
def initializeTable(opt):

    """
    Initializes the hdf5 file with 'Repeater Catalog', 'Orphans', 'Junk', and 'Correlation
    Matrix' tables in a group related to the station where the data come from. This is
    defined via the redpy.config.Options class.
    
    opt: Options object describing the station/run parameters

    Saves table to file and closes it.
    Will likely need extensive editing when more tables get added...
    """

    h5file = open_file(opt.filename, mode="w", title=opt.title)
    group = h5file.create_group("/", opt.groupName, opt.groupDesc)

    rtable = h5file.create_table(group, "repeaters", Repeaters(opt),
        "Repeater Catalog")
    rtable.attrs.scnl = [opt.station, opt.channel, opt.network, opt.location]
    rtable.attrs.samprate = opt.samprate
    rtable.attrs.windowLength = opt.winlen
    rtable.attrs.ptrig = opt.ptrig
    rtable.attrs.atrig = opt.atrig
    rtable.attrs.fmin = opt.fmin
    rtable.attrs.fmax = opt.fmax
    rtable.attrs.previd = 0
    rtable.attrs.ptime = 0
    rtable.flush()
    
    otable = h5file.create_table(group, "orphans", Orphans(opt),
        "Orphan Catalog")
    otable.flush()
    
    jtable = h5file.create_table(group, "junk", Junk(opt), "Junk Catalog")
    jtable.flush()
    
    dtable = h5file.create_table(group, "deleted", Deleted(opt),
        "Manually Deleted Events")
    dtable.flush()

    ctable = h5file.create_table(group, "correlation", Correlation(opt),
        "Correlation Matrix")
    ctable.flush()
    
    ftable = h5file.create_table(group, "families", Families(opt), "Families Table")
    ftable.attrs.nClust = 0
    ftable.flush()

    h5file.close()


def openTable(opt):

    """
    Convenience function to open the catalog and access the tables in it.
    
    opt: Options object describint station/run parameters
    
    Returns handles to h5file, rtable, otable, ctable, and jtable
    """

    h5file = open_file(opt.filename, "a")
    
    rtable = eval('h5file.root.'+ opt.groupName + '.repeaters')
    otable = eval('h5file.root.'+ opt.groupName + '.orphans')
    ctable = eval('h5file.root.'+ opt.groupName + '.correlation')
    jtable = eval('h5file.root.'+ opt.groupName + '.junk')
    dtable = eval('h5file.root.'+ opt.groupName + '.deleted')
    ftable = eval('h5file.root.'+ opt.groupName + '.families')
    
    return h5file, rtable, otable, ctable, jtable, dtable, ftable


def calcAmps(data, windowStart, opt):
    
    amps = []
    for n in range(opt.nsta):
        amps.append(max(abs(data[(n*opt.wshape)+windowStart:int(
            (n*opt.wshape)+windowStart+opt.winlen/2)])))
    
    return amps
    
def populateRepeater(rtable, ftable, id, trig, opt, windowStart=-1):

    """
    Initially populates a new row in the 'Repeater Catalog' table.
    
    rtable: object pointing to the repeater table to populate
        (e.g., h5file.root.groupName.repeaters)
    id: integer id number given to this trigger, should be unique
    trig: ObsPy trace from triggering function
    opt: Options object describing station/run parameters
    windowStart: triggering time (defaults to opt.ptrig seconds)

    Appends this row to Repeaters table, but does not update the clustering parameters
        (sets them to 0)
    """
    
    trigger = rtable.row
    
    if windowStart == -1:
        windowStart = int(opt.ptrig*opt.samprate)
    
    trigger['id'] = id
    trigger['startTime'] = trig.stats.starttime.isoformat()
    try:
        trigger['startTimeMPL'] = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S.%f'))
    except ValueError:
        trigger['startTimeMPL'] = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S'))
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['windowCoeff'], trigger['windowFFT'] = redpy.correlation.calcWindow(
        trig.data, windowStart, opt)
    trigger['windowAmp'] = calcAmps(trig.data, windowStart, opt)
    
    trigger.append()  
    rtable.flush()  

    
def populateOrphan(otable, id, trig, opt):

    """
    Initially populates a new row in the 'Orphans' table.
    
    otable: object pointing to the table to populate
        (e.g., h5file.root.groupName.orphans)
    id: integer id number given to this trigger, should be unique
    trig: ObsPy trace from triggering function
    opt: Options object describing station/run parameters

    Appends this row to Orphans table, adding an expiration date 
    """
    
    trigger = otable.row
    
    windowStart = int(opt.ptrig*opt.samprate)
    
    trigger['id'] = id
    trigger['startTime'] = trig.stats.starttime.isoformat()
    try:
        trigger['startTimeMPL'] = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S.%f'))
    except ValueError:
        trigger['startTimeMPL'] = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S'))
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['windowCoeff'], trigger['windowFFT'] = redpy.correlation.calcWindow(
        trig.data, windowStart, opt)
    trigger['windowAmp'] = calcAmps(trig.data, windowStart, opt)

    adddays = ((opt.maxorph-opt.minorph)/7.)*(trig.stats.maxratio-opt.trigon)+opt.minorph
    trigger['expires'] = (trig.stats.starttime+adddays*86400).isoformat()
    trigger.append()
    otable.flush()


def populateJunk(jtable, trig, isjunk, opt):
    
    """
    Initially populates a new row in the 'Junk' table.
    
    jtable: object pointing to the table to populate
        (e.g., h5file.root.groupName.junk)
    trig: ObsPy trace from triggering function
    isjunk: Integer flag, 0=junk, 1=expired orphan
    opt: Options object describing station/run parameters
    
    """
    
    trigger = jtable.row
    
    windowStart = int(opt.ptrig*opt.samprate)
    
    trigger['startTime'] = trig.stats.starttime.isoformat()
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['isjunk'] = isjunk
    trigger.append()
    jtable.flush()


def moveOrphan(rtable, otable, ftable, oindex, opt):
    
    """
    Moves a row from the 'Orphans' table to the 'Repeater Catalog' table.
    
    rtable: Repeater table
    otable: Orphan table
    ftable: Families table
    oindex: Row in otable to move
    opt: Options object describing station/run parameters
    
    """
    
    trigger = rtable.row
    orow = otable[oindex]
    
    trigger['id'] = orow['id']
    trigger['startTime'] = orow['startTime']
    trigger['startTimeMPL'] = orow['startTimeMPL']
    trigger['waveform'] = orow['waveform']
    trigger['windowStart'] = orow['windowStart']
    if len(otable) > 1:
        trigger['windowCoeff'] = orow['windowCoeff']
        trigger['windowFFT'] = orow['windowFFT']
    else:
        coeff, fft = redpy.correlation.calcWindow(orow['waveform'], orow['windowStart'],
            opt)
        trigger['windowCoeff'] = coeff
        trigger['windowFFT'] = fft
        otable.cols.windowCoeff[oindex] = 0
        otable.cols.expires[oindex] = (UTCDateTime(orow['startTime'])-86400).isoformat()
    trigger['windowAmp'] = orow['windowAmp']
    
    trigger.append()
        
    if len(otable) > 1:
        otable.remove_row(oindex)
    
    rtable.flush()   
    otable.flush()  
    

def removeFamilies(rtable, ctable, dtable, ftable, cnums, opt):

    """
    Moves the core of the families into the dtable, deletes the rest of the members.
    
    rtable: Repeater table
    ctable: Correlation matrix table
    dtable: Deleted table
    ftable: Families table
    cnums: Families to remove
    opt: Options object describing station/run parameters
    """

    cnums = np.sort(cnums)
    old = range(len(rtable))
    transform = np.zeros((len(rtable),)).astype(int)    
    ids = rtable.cols.id[:]
    members = np.array([])
    oldcores = ftable.cols.core[:]
    for cnum in cnums[::-1]:
        members = np.append(members, np.fromstring(ftable[cnum]['members'], dtype=int, sep=' '))
        ftable.remove_row(cnum)
    ftable.attrs.nClust-=len(cnums)
    members = np.sort(members).astype('uint32')
    
    cores = rtable[np.intersect1d(members, oldcores)]
    for core in cores:
        trigger = dtable.row
        trigger['id'] = core['id']
        trigger['startTime'] = core['startTime'][0]
        trigger['startTimeMPL'] = core['startTimeMPL']
        trigger['waveform'] = core['waveform']
        trigger['windowStart'] = core['windowStart']
        trigger['windowCoeff'] = core['windowCoeff']
        trigger['windowFFT'] = core['windowFFT']
        trigger['windowAmp'] = core['windowAmp']
        trigger.append()
        
    ids = ids[members]
    id2 = ctable.cols.id2[:]
    idxc = np.where(np.in1d(id2,ids))[0]
    for c in idxc[::-1]:
        ctable.remove_row(c)

    for m in members[::-1]:
        rtable.remove_row(m)
        old.remove(m)
    
    new = range(len(rtable))
    transform[old] = new
    
    for n in range(len(ftable)):
        members = np.fromstring(ftable[n]['members'], dtype=int, sep=' ')
        core = ftable[n]['core']
        ftable.cols.members[n] = np.array2string(transform[members])[1:-1]
        ftable.cols.core[n] = transform[core]
        if n>=min(cnums):
            ftable.cols.printme[n] = 1
        ftable.flush()
        
    rtable.flush()
    dtable.flush()


def clearExpiredOrphans(otable, opt, tend):

    """
    Deletes orphans that have passed their expiration date
    
    otable: object pointing to the table to populate
        (e.g., h5file.root.groupName.orphans)
    opt: Options object describing station/run parameters
    tend: Time to remove orphans older than, corresponds usually to end of run increment
    
    Removes orphans from table, prints how many were removed. Checks to make sure there
    is always at least one orphan in the table.
    """
    
    index = np.where(otable.cols.expires[:] < tend.isoformat())
    if len(index) != len(otable):
        for n in range(len(index[0])-1,-1,-1):
            otable.remove_row(index[0][n])        
    else:
        print('Warning: All orphans expired...')
        for n in range(len(index[0])-1,0,-1):
            otable.remove_row(index[0][n])
    otable.flush()
    


def appendCorrelation(ctable, id1, id2, ccc, opt):

    """
    Appends a new value to the 'Correlation Matrix' table.

    corr: object pointing to the row in the table to populate
        (e.g., h5file.root.hsr.correlation.row)
    id1: unique id number of first trigger
    id2: unique id number of second trigger
    ccc: cross-correlation between the two triggers in the window
    opt: Options object describing station/run parameters

    Appends this row to the table, and automatically puts the smaller of
    the two id numbers first
    
    Only appends if the value is greater than the minimum defined in opt
    """
    
    if (ccc >= opt.cmin) and (id1!=id2):
        corr = ctable.row
        corr['id1'] = min(id1, id2)
        corr['id2'] = max(id1, id2)
        corr['ccc'] = ccc
        corr.append()        
        ctable.flush()


def createNewFamily(rtable, ftable, members, core, opt):
    
    """
    Creates new family from two or more orphans
    """
    
    f = ftable.row
    f['members'] = np.array2string(members)[1:-1]
    f['core'] = core
    f['startTime'] = np.min(rtable[members]['startTimeMPL'])
    f['printme'] = 1
    f.append()
    ftable.attrs.nClust+=1
    ftable.flush()
    
    if len(ftable)>1:
        reorderFamilies(ftable, opt)
    
    
def reorderFamilies(ftable, opt):

    """
    Cleanup function to ensure families are ordered by start time
    """
    
    startTimes = ftable.cols.startTime[:]
    order = np.argsort(startTimes)
    x = np.arange(len(ftable))
        
    if (x!=order).any():
        members = ftable.cols.members[:]
        cores = ftable.cols.core[:]
        printme = ftable.cols.printme[:]
        for n in np.where(x!=order)[0]:
            ftable.cols.startTime[n] = startTimes[order[n]]
            ftable.cols.members[n] = members[order[n]]
            ftable.cols.core[n] = cores[order[n]]
            ftable.cols.printme[n] = printme[order[n]]
        ftable.flush()
    
    
def mergeFamilies(rtable, ctable, ftable, wfam, wlag, opt):

    """
    Combines families that have been merged by adding a new event
    """
    
    # Determine which family is largest, use that as base family
    nmem = []
    for n in range(len(wfam)):
        nmem.append(len(np.fromstring(ftable[wfam[n]]['members'], dtype=int, sep=' ')))
    wlag = np.array(wlag)
    wlag = wlag - wlag[np.argmax(nmem)]
    
    # Adjust lags
    for n in range(len(wfam)):
        if wlag[n]!=0:
            members = np.fromstring(ftable[wfam[n]]['members'], dtype=int, sep=' ')
            for m in members:
                rtable.cols.windowStart[m] = rtable.cols.windowStart[m] - wlag[n]
                rtable.cols.windowCoeff[m], rtable.cols.windowFFT[m] = redpy.correlation.calcWindow(
                    rtable.cols.waveform[m], rtable.cols.windowStart[m], opt)
            rtable.flush()
    
    # Perform merge, cluster        
    f1 = min(wfam)
    for n in range(len(wfam))[::-1]:
        f2 = np.sort(wfam)[n]
        if f2!=f1:            
            ftable.cols.members[f1]+=' '+ftable[f2]['members']
            ftable.remove_row(f2)
            ftable.attrs.nClust-=1
    ftable.cols.printme[f1:len(ftable)] = np.ones((len(ftable)-f1,))            
    reorderFamilies(ftable, opt)
    redpy.cluster.runFamOPTICS(rtable, ctable, ftable, f1, opt)
    