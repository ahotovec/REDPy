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
    order: Order in the cluster ordering (integer)
    reachability: Reachability in the cluster ordering (float)
    coreDistance: Core distance in the cluster ordering (float)
    clusterNumber: ID of flat cluster (integer)
    isCore: 1 if core, else 0 (integer)
    alignedTo: ID of event this one is aligned to (integer)
    plotClust: ID of flat cluster ordered by time (integer)
    lastClust: ID of flat cluster when last plotted (integer)
    
    Returns a dictionary defining the table
    """
    
    dict = {
        "id"            : Int32Col(shape=(), pos=0),
        "startTime"     : StringCol(itemsize=32, pos=1),
        "startTimeMPL"  : Float64Col(shape=(), pos=2),
        "waveform"      : Float64Col(shape=(opt.wshape,), pos=3),
        "windowStart"   : Int32Col(shape=(), pos=4),
        "windowCoeff"   : Float64Col(shape=(), pos=5),
        "windowFFT"     : ComplexCol(shape=(opt.winlen,), itemsize=16, pos=6),
        "windowAmp"     : Float64Col(shape=(), pos=7),
        "order"         : Int32Col(shape=(), pos=8),
        "reachability"  : Float64Col(shape=(), pos=9),
        "coreDistance"  : Float64Col(shape=(), pos=10),
        "clusterNumber" : Int32Col(shape=(), pos=11),
        "isCore"        : Int32Col(shape=(), pos=12),
        "alignedTo"     : Int32Col(shape=(), pos=13),
        "plotClust"     : Int32Col(shape=(), pos=14),
        "lastClust"     : Int32Col(shape=(), pos=15)
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
        "waveform"    : Float64Col(shape=(opt.wshape,), pos=3),
        "windowStart" : Int32Col(shape=(), pos=4),
        "windowCoeff" : Float64Col(shape=(), pos=5),
        "windowFFT"   : ComplexCol(shape=(opt.winlen,), itemsize=16, pos=6),
        "windowAmp"   : Float64Col(shape=(), pos=7),
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
        "waveform"    : Float64Col(shape=(opt.wshape,), pos=3),
        "windowStart" : Int32Col(shape=(), pos=4),
        "windowCoeff" : Float64Col(shape=(), pos=5),
        "windowFFT"   : ComplexCol(shape=(opt.winlen,), itemsize=16, pos=6),
        "windowAmp"   : Float64Col(shape=(), pos=7)
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
        "waveform"    : Float64Col(shape=(opt.wshape,), pos=2),
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
    rtable.attrs.cores = []
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
    
    return h5file, rtable, otable, ctable, jtable, dtable

    
def populateRepeater(rtable, id, trig, opt, alignedTo, windowStart=-1):

    """
    Initially populates a new row in the 'Repeater Catalog' table.
    
    rtable: object pointing to the repeater table to populate
        (e.g., h5file.root.groupName.repeaters)
    id: integer id number given to this trigger, should be unique
    trig: ObsPy trace from triggering function
    opt: Options object describing station/run parameters
    alignedTo: id number of repeater this one is aligned to (can be itself)
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
    trigger['windowAmp'] = max(abs(trig.data[windowStart:int(windowStart+opt.winlen/2)]))
    trigger['order'] = -1
    trigger['reachability'] = -1.0
    trigger['coreDistance'] = -1.0
    trigger['clusterNumber'] = -1
    trigger['plotClust'] = -1
    trigger['lastClust'] = -1
    trigger['isCore'] = 0 # Set to zero to avoid being counted erroneously as a core
    trigger['alignedTo'] = alignedTo
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
    trigger['windowAmp'] = max(abs(trig.data[windowStart:int(windowStart+opt.winlen/2)]))

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


def moveOrphan(rtable, otable, oindex, alignedTo, opt):
    
    """
    Moves a row from the 'Orphans' table to the 'Repeater Catalog' table.
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
    trigger['order'] = -1
    trigger['reachability'] = -1.0
    trigger['coreDistance'] = -1.0
    trigger['clusterNumber'] = -1
    trigger['plotClust'] = -1
    trigger['lastClust'] = -1
    trigger['isCore'] = 0 # Set to zero to avoid being counted erroneously as a core
    trigger['alignedTo'] = alignedTo
    trigger.append()
        
    if len(otable) > 1:
        otable.remove_row(oindex)
    
    rtable.flush()   
    otable.flush()  
    


def removeFamily(rtable, ctable, dtable, cnum, opt):

    """
    Moves the core of a family into the dtable, deletes the rest of the members.
    """
    
    trigger = dtable.row
    
    members = rtable.get_where_list('(plotClust=={})'.format(cnum))
    idx = rtable.get_where_list('(plotClust=={}) & (isCore==1)'.format(cnum))
    core = rtable[idx]
    
    trigger['id'] = core['id']
    trigger['startTime'] = core['startTime'][0]
    trigger['startTimeMPL'] = core['startTimeMPL']
    trigger['waveform'] = core['waveform']
    trigger['windowStart'] = core['windowStart']
    trigger['windowCoeff'] = core['windowCoeff']
    trigger['windowFFT'] = core['windowFFT']
    trigger['windowAmp'] = core['windowAmp']
    trigger.append()
    
    for m in members[::-1]:
        id = rtable.cols.id[m]
        idxc = ctable.get_where_list('(id1=={0}) | (id2=={0})'.format(id))
        for c in idxc[::-1]:
            ctable.remove_row(c)
        rtable.remove_row(m)
    
    rtable.attrs.cores = rtable.get_where_list('(isCore == 1)')
    
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
        print '%i Orphans aged out of the system' % len(index[0])
    else:
        for n in range(len(index[0])-1,0,-1):
            otable.remove_row(index[0][n])
        print '%i Orphans aged out of the system, 1 left' % len(index[0])-1
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
