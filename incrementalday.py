import numpy as np
import redpy as red
import matplotlib.pyplot as plt
from obspy import UTCDateTime
import time

'''
This script is a testing ground for realigning and clustering based on
cross-correlations using a brute-force approach.

- First it downloads the data
- Triggers on it using STA/LTA
- Then, it goes through each trigger: (THIS IS TIME CONSUMING!)
        Realigns it to the best match
        Finds the optimal cross-correlation coefficient
- Finally, it uses OPTICS to cluster and order the triggers

Plots the correlation and lag matrices of the triggers ordered in time vs.
ordered in clustering order, and also plots the ordered waveforms for visual
inspection.
'''

# Currently does not scale well... so take the first 250

t = UTCDateTime("2004-09-27")
st = red.getIRIS(t, "HSR", "EHZ", "UW", nsec=86400)
trigs = red.trigger(st)
trigs = trigs[0:250]

# Initialize with first two triggers
# Need to define the windows better...
ttrim = trigs.copy()
for n in range(len(ttrim)):
    ttrim[n] = ttrim[n].slice(ttrim[n].stats.starttime + 10, ttrim[n].stats.starttime + 15.11)

print("Correlating!")
C, L, X, wcoeff = red.xcorr_all_fft(ttrim[0:2])
if C[0, 1] > 0.6:
    # Check if first two need to be realigned
    ttrim[1] = trigs[1].slice(trigs[1].stats.starttime + 10 - L[0, 1], trigs[1].stats.starttime + 15.11 - L[0, 1])
    C, L, X, wcoeff = red.xcorr_all_fft(ttrim[0:2])

# Now go through next triggers!
t = time.time()
for n in range(2,len(trigs)):
    Ctmp, Ltmp, Xtmp, wcoefftmp = red.xcorr_addrow_fft(ttrim[0:n], ttrim[n], C, L, X, wcoeff)
    if max(Ctmp[n, 0:n]) > 0.6:
        indnn = np.argmax(Ctmp[n, 0:n])
        ttrim[n] = trigs[n].slice(trigs[n].stats.starttime + 10 - Ltmp[indnn, n], trigs[n].stats.starttime + 15.11 - Ltmp[indnn, n])
        C, L, X, wcoeff = red.xcorr_addrow_fft(ttrim[0:n], ttrim[n], C, L, X, wcoeff)
    else:
        C = Ctmp
        L = Ltmp
        X = Xtmp
        wcoeff = wcoefftmp
    print(100.0*n/len(trigs)) # print % of events done

print("Done in: " + repr(time.time()-t))

maxes = np.max(C - np.eye(len(C)), axis=0)
print("Orphans: " + repr(100.0*len(maxes[maxes<0.7])/len(C)) + "%")

# Cluster with OPTICS
ttree = red.setOfObjects(1-C)
red.prep_optics(ttree,1)
red.build_optics(ttree,1)
order = ttree._ordered_list
Co = C[order, :]
Co = Co[:, order]
Lo = L[order, :]
Lo = Lo[:, order]

# Plot the C matrices (saturated below 0.65)
# Left is ordered by time
# Right is ordered by cluster ordering
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 2, 1)
ax.imshow(C, aspect='auto', vmin=0.65, vmax=1, interpolation='nearest')
ax = fig.add_subplot(1, 2, 2)
ax.imshow(Co, aspect='auto', vmin=0.65, vmax=1, interpolation='nearest')

# Plot the L matrices (same as before, just shows how aligned they are)
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 2, 1)
ax.imshow(L, aspect='auto', vmin=-5, vmax=5, interpolation='nearest', cmap='RdBu')
ax = fig.add_subplot(1, 2, 2)
ax.imshow(Lo, aspect='auto', vmin=-5, vmax=5, interpolation='nearest', cmap='RdBu')

# Plot the ordered waveforms
ttrim = ttrim.normalize()
data = np.zeros((len(ttrim), 1501))
for n in range(len(ttrim)):
    data[n, :] = trigs[n].normalize().slice(ttrim[n].stats.starttime - 5, ttrim[n].stats.starttime + 10).data
datao = data[order, :]
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 1, 1)
ax.imshow(datao, aspect='auto', vmin=-1, vmax=1, interpolation='nearest', cmap='RdBu')

plt.show()

