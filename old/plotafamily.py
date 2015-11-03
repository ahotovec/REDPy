from tables import *
import numpy as np
import matplotlib.pyplot as plt

"""
This is a very crude code for browsing individual families after a run
"""

h5file = open_file("test.h5", "a")
rtable = h5file.root.hsr.repeaters

for fam in range(max(rtable.cols.clusterNumber[:])+1):
    members = rtable[rtable.get_where_list('clusterNumber == {}'.format(fam))]
    
    print('Family {0}'.format(fam))
    
    data = np.zeros((len(members), 2000))
    fig = plt.figure(figsize=(12, 6))
    n=-1
    for r in members:
        n = n+1
        data[n, :] = r['waveform'][r['windowStart']-500:r['windowStart']+1500]
        data[n, :] = data[n, :]/max(data[n, 500:1012])
        
        print(r['startTime'])
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(data, aspect='auto', vmin=-1, vmax=1, interpolation='nearest', cmap='RdBu')
    
    plt.show()

h5file.close()