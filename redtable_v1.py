from tables import *
import numpy as np
import redpy as red
from scipy.fftpack import fft
from obspy import UTCDateTime

# Build the table
class Triggers(IsDescription):
    id = Int32Col(shape=(), pos=0)
    startTime = StringCol(itemsize=32, pos=1)
    waveform = Float64Col(shape=(3001,), pos=2)
    windowStart = Int32Col(shape=(), pos=3)
    windowCoeff = Float64Col(shape=(), pos=4)
    windowFFT = ComplexCol(shape=(512,), itemsize=16, pos=5)

# Initialize
h5file = open_file("redtable_v1.h5", mode="w", title="REDPy Catalog")
group = h5file.create_group("/", 'hsr', 'MSH: HSR-EHZ-UW')
table = h5file.create_table(group, 'repeaters',
    Triggers, "Repeater Catalog")
table.attrs.scnl = ["HSR", "EHZ", "UW", "--"]
table.attrs.samprate = 100
table.attrs.windowLength = 512

# Grab some data to put in it
t = UTCDateTime("2004-09-24")
st = red.getIRIS(t, "HSR", "EHZ", "UW", nsec=3600)
trigs = red.trigger(st)

# Populate table
trigger = table.row
for i in range(len(trigs)):
    trigger['id'] = i
    trigger['startTime'] = trigs[i].stats.starttime.isoformat()
    trigger['waveform'] = trigs[i].data
    trigger['windowStart'] = 500
    trigger['windowCoeff'] = 1/np.sqrt(sum(trigs[i].data[500:1012] * 
        trigs[i].data[500:1012]))
    trigger['windowFFT'] = np.reshape(fft(trigs[i].data[500:1012]), (512,))
    trigger.append()

# Write to disk
table.flush()

# Close file when done
h5file.close()
