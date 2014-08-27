from tables import *
from obspy.core.trace import Trace
from obspy import UTCDateTime
import datetime
import numpy as np
import redpy.correlation

def Repeaters(opt):

    """
    Defines the columns in the 'Repeater Catalog' table based on the Options in opt
    
    id: unique ID number for the event (integer)
    startTime: UTC time of start of the waveform (string)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    windowCoeff: amplitude scaling for cross-correlation (float)
    windowFFT: Fourier transform of window (complex ndarray)
    order: Order in the cluster ordering (integer)
    reachability: Reachability in the cluster ordering (float)
    coreDistance: Core distance in the cluster ordering (float)

    Returns a dictionary defining the table
    """
    
    dict = {
        "id"            : Int32Col(shape=(), pos=0),
        "startTime"     : StringCol(itemsize=32, pos=1),
        "waveform"      : Float64Col(shape=(opt.wshape,), pos=2),
        "windowStart"   : Int32Col(shape=(), pos=3),
        "windowCoeff"   : Float64Col(shape=(), pos=4),
        "windowFFT"     : ComplexCol(shape=(opt.winlen,), itemsize=16, pos=5),
        "order"         : Int32Col(shape=(), pos=6),
        "reachability"  : Float64Col(shape=(), pos=7),
        "coreDistance"  : Float64Col(shape=(), pos=8),
        "clusterNumber" : Int32Col(shape=(), pos=9),
        "isCore"        : Int32Col(shape=(), pos=10),
        }
    
    return dict


def Orphans(opt):

    """
    Defines the columns in the 'Orphans' table based on the Options in opt

    id: unique ID number for the event (integer)
    startTime: UTC time of start of the waveform (string)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    windowCoeff: amplitude scaling for cross-correlation (float)
    windowFFT: Fourier transform of window (complex ndarray)
    expires: UTC time of when orphan should no longer be considered (string)
    
    Returns a dictionary defining the table
    """
    
    dict = {
        "id"          : Int32Col(shape=(), pos=0),
        "startTime"   : StringCol(itemsize=32, pos=1),
        "waveform"    : Float64Col(shape=(opt.wshape,), pos=2),
        "windowStart" : Int32Col(shape=(), pos=3),
        "windowCoeff" : Float64Col(shape=(), pos=4),
        "windowFFT"   : ComplexCol(shape=(opt.winlen,), itemsize=16, pos=5),
        "expires"     : StringCol(itemsize=32, pos=6)
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
    rtable.flush()
    
    otable = h5file.create_table(group, "orphans", Orphans(opt),
        "Orphan Catalog")
    otable.flush()
    
    jtable = h5file.create_table(group, "junk", Junk(opt), "Junk Catalog")
    jtable.flush()

    ctable = h5file.create_table(group, "correlation", Correlation(opt),
        "Correlation Matrix")
    ctable.flush()

    h5file.close()

    
def populateRepeater(rtable, id, trig, opt, windowStart=-1):

    """
    Initially populates a new row in the 'Repeater Catalog' table.
    
    rtable: object pointing to the repeater table to populate
        (e.g., h5file.root.groupName.repeaters)
    id: integer id number given to this trigger, should be unique
    trig: ObsPy trace from triggering function
    opt: Options object describing station/run parameters
    windowStart: triggering time (defaults to opt.ptrig seconds)

    Appends this row to table, but does not update the clustering parameters (sets them
    to 0)
    """
    
    trigger = rtable.row
    
    if windowStart == -1:
        windowStart = int(opt.ptrig*opt.samprate)
    
    trigger['id'] = id
    trigger['startTime'] = trig.stats.starttime.isoformat()
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['windowCoeff'], trigger['windowFFT'] = redpy.correlation.calcWindow(
        trig.data, windowStart, opt)
    trigger['order'] = -1
    trigger['reachability'] = -1.0
    trigger['coreDistance'] = -1.0
    trigger['clusterNumber'] = -1
    trigger['isCore'] = 0 # Set to zero to avoid being counted erroneously as a core
    trigger.append()  
    rtable.flush()  

    
def populateOrphan(otable, id, trig, opt):

    """
    Initially populates a new row in the 'Orphans' table.
    
    otable: object pointing to the table to populate
        (e.g., h5file.root.groupName.orphans)
    id: integer id number given to this trigger, should be unique
    trig: ObsPy trace from triggering function
    expires: Expiration date (as a string) for orphan
    opt: Options object describing station/run parameters

    Appends this row to table, 
    """
    
    trigger = otable.row
    
    windowStart = int(opt.ptrig*opt.samprate)
    
    trigger['id'] = id
    trigger['startTime'] = trig.stats.starttime.isoformat()
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['windowCoeff'], trigger['windowFFT'] = redpy.correlation.calcWindow(
        trig.data, windowStart, opt)
    #equation of a straight line between (opt.trigon,opt.minorph) and (opt.trigon+7,opt.maxorph) - finds y value (in days from trig starttime) that goes with maxratio of current trigger
    adddays = ((opt.maxorph-opt.minorph)/7.)*(trig.stats.maxratio-opt.trigon)+opt.minorph
    trigger['expires'] = (trig.stats.starttime+datetime.timedelta(days=adddays)).isoformat
    #print trigger['expires']
    trigger.append()
    otable.flush()
    #print("Orphans abandoned: {0}".format(len(otable)))


def populateJunk(jtable, trig, isjunk, opt):
    
    """
    Initially populates a new row in the 'Orphans' table.
    
    jtable: object pointing to the table to populate
    (e.g., h5file.root.groupName.junk)
    trig: ObsPy trace from triggering function
    isjunk: Integer flag, 0=junk, 1=expired orphan
    opt: Options object describing station/run parameters
    
    Appends this row to table,
    """
    
    trigger = jtable.row
    
    windowStart = int(opt.ptrig*opt.samprate)
    
    trigger['startTime'] = trig.stats.starttime.isoformat()
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['isjunk'] = isjunk
    trigger.append()
    jtable.flush()

    # print("Junk events: {0}".format(len(jtable)))


def moveOrphan(rtable, otable, oindex, opt):
    
    """
    Moves a row from the 'Orphans' table to the 'Repeater Catalog' table.
    """
    
    trigger = rtable.row
    orow = otable[oindex]
    
    trigger['id'] = orow['id']
    trigger['startTime'] = orow['startTime']
    trigger['waveform'] = orow['waveform']
    trigger['windowStart'] = orow['windowStart']
    trigger['windowCoeff'] = orow['windowCoeff']
    trigger['windowFFT'] = orow['windowFFT']
    trigger['order'] = -1
    trigger['reachability'] = -1.0
    trigger['coreDistance'] = -1.0
    trigger['clusterNumber'] = -1
    trigger['isCore'] = 0 # Set to zero to avoid being counted erroneously as a core
    trigger.append()
    
    otable.remove_row(oindex)
    
    otable.flush()  
    rtable.flush()

def clearExpiredOrphans(otable,opt,tend):
    """
    Deletes orphans that have passed their expiration date (relative to the end time of the current clustering run)
    """

    #index = otable.get_where_list(eval('expires > ' + repr(float(tend))))
    #print len(otable.cols.expires[:])
    print tend
    print otable.cols.expires[0]
    index = np.where(otable.cols.expires[:] > tend.isoformat)
    #print index[0]
    for n in range(len(index[0])):
        otable.remove_row(index[0][n])
    otable.flush()
    print '%i Orphans aged out of the system' % len(index[0])

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
    
    if ccc >= opt.cmin:
        corr = ctable.row
        corr['id1'] = min(id1, id2)
        corr['id2'] = max(id1, id2)
        corr['ccc'] = ccc
        corr.append()        
        ctable.flush()