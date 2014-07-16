from tables import *
from obspy.core.trace import Trace
import redpy.correlation

class Triggers(IsDescription):

    """
    Defines the columns in the "Repeater Catalog" table

    id: unique ID number for the event (integer)
    startTime: UTC time of start of the waveform (string)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    windowCoeff: amplitude scaling for cross-correlation (float)
    windowFFT: Fourier transform of window (complex ndarray)

    Needs work to figure out how to adjust the shape of the waveform and
    windowFFT columns when the window length and padding around the triggers
    are not the same from station to station
    """

    id = Int32Col(shape=(), pos=0)
    startTime = StringCol(itemsize=32, pos=1)
    waveform = Float64Col(shape=(3001,), pos=2)
    windowStart = Int32Col(shape=(), pos=3)
    windowCoeff = Float64Col(shape=(), pos=4)
    windowFFT = ComplexCol(shape=(512,), itemsize=16, pos=5)

    
class Correlation(IsDescription):

    """
    Defines the columns in the "Correlation" table

    id1: unique ID number for the first event (integer)
    id2: unique ID number for the second event (integer)
    ccc: cross-correlation coefficient between those two events (float)    
    """
    
    id1 = Int32Col(shape=(), pos=0)
    id2 = Int32Col(shape=(), pos=1)
    ccc = Float64Col(shape=(), pos=2)

    
def initializeTable(opt):

    """
    Initializes the hdf5 file with a "Repeater Catalog" table and a "Correlation Matrix"
    table in a group related to the station where the data come from. This is defined
    via the redpy.config.Options class.
    
    opt: an Options object describing the table/run

    Saves table to file and closes it.
    Will need extensive editing when more tables get added...
    """

    h5file = open_file(opt.filename, mode="w", title=opt.title)
    group = h5file.create_group("/", opt.groupName, opt.groupDesc)

    rtable = h5file.create_table(group, "repeaters", Triggers,
        "Repeater Catalog")
    rtable.attrs.scnl = [opt.station, opt.channel, opt.network, opt.location]
    rtable.attrs.samprate = opt.samprate
    rtable.attrs.windowLength = opt.winlen
    rtable.attrs.ptrig = opt.ptrig
    rtable.attrs.atrig = opt.atrig
    rtable.attrs.fmin = opt.fmin
    rtable.attrs.fmax = opt.fmax
    rtable.flush()

    ctable = h5file.create_table(group, "correlation", Correlation,
        "Correlation Matrix")
    ctable.attrs.order = 0
    ctable.attrs.reachability = 0
    ctable.attrs.coredist = 0
    ctable.flush()

    h5file.close()

    
def populateTrigger(trigger, id, trig, opt):

    """
    Initially populates the trigger row in the 'Repeater Catalog' table.
    
    trigger: object pointing to the row in the table to populate
        (e.g., h5file.root.hsr.repeaters.row)
    id: integer id number given to this trigger, should be unique
    trig: ObsPy trace from triggering function
    opt: Options object describing station/run parameters

    Appends this row to table
    """
    
    windowStart = int(opt.ptrig*opt.samprate)
    
    trigger['id'] = id
    trigger['startTime'] = trig.stats.starttime.isoformat()
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['windowCoeff'], trigger['windowFFT'] = redpy.correlation.calcWindow(trig.data,
        windowStart)
    trigger.append()    


def appendCorrelation(corr, id1, id2, ccc):

    """
    Appends a new value to the 'Correlation Matrix' table.

    corr: object pointing to the row in the table to populate
        (e.g., h5file.root.hsr.correlation.row)
    id1: unique id number of first trigger
    id2: unique id number of second trigger
    ccc: cross-correlation between the two triggers in the window

    Appends this row to the table, and automatically puts the smaller of
    the two id numbers first
    """

    corr['id1'] = min(id1, id2)
    corr['id2'] = max(id1, id2)
    corr['ccc'] = ccc
    corr.append()


def getCell(table, id, column):
    
    """
    Shorthand way of getting data from the PyTable.

    table: PyTable you're querying
    id: unique id of row you want
    column: column you want (e.g., 'windowFFT' OR position as integer)

    Returns data inside that cell

    While simple, can be time consuming if called a lot!
    """
   
    c = '(id == {})'.format(id)
    t = table.where(c)
    for r in t: data = r[column]

    return data