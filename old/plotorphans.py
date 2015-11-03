from tables import *
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import datetime

"""
This is a very crude code for browsing orphans - change line 18 to edit which ones are shown
"""

h5file = open_file("testsep.h5", "a")
otable = h5file.root.hsr.orphans
samprate = h5file.root.hsr.repeaters.attrs.samprate

#waveform wiggle plot of input waveforms
fig = plt.figure(figsize=(20, 15))
n=0.
for r in otable.iterrows():
    #if r['expires'] > (datetime.datetime.utcnow()+datetime.timedelta(days=20.)).isoformat():
    if r['expires'] > UTCDateTime('2014-08-24').isoformat():
        dat=r['waveform']
        dat=dat/max(dat)
        tvec = np.arange(0,len(dat)*1/samprate,1/samprate)
        plt.plot(tvec,np.add(dat,n))
        n = n+3
plt.ylabel('index')
plt.xlabel('time(s)')
#plt.title('Junk triggers')
plt.autoscale(tight=True)
plt.show()