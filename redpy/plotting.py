from tables import *
import numpy as np
import matplotlib.pyplot as plt
import time
import mpld3

"""
These are very brute force plotting; should be replaced with more sophisticated functions
"""

def createTimelineFigure(rtable, ctable, opt):
    fig = plt.figure(figsize=(15, 10))
    
    dt = rtable.cols.startTimeMPL[:]
    cnum = rtable.cols.clusterNumber[:]
    reach = rtable.cols.reachability[:]
    reach[0] = 1
    reach[reach==1] = 0

    # Figure out earliest member in each family
    mindt = np.zeros((max(cnum)+1,))
    for clustNum in range(max(cnum)+1):
        mindt[clustNum] = min(dt[cnum==clustNum])
    
    n = 0
    for clustNum in np.argsort(mindt):
        if len(dt[cnum==clustNum]) > 2:
            plt.plot_date((min(dt[cnum==clustNum]), max(dt[cnum==clustNum])), (n, n),
                'k:', linewidth=0.5)
            plt.scatter(dt[cnum==clustNum], n*np.ones((len(dt[cnum==clustNum]),)),
                c=(1-reach[cnum==clustNum]), vmin=0.6, vmax=1, cmap='jet', linewidth=0.5)
            n = n+1
    
    plt.margins(0.05)
    
    plt.gcf().autofmt_xdate()
    plt.xlabel('Date')
    plt.ylabel('Cluster by Date (with at least 3 members)')
    
    mpld3.save_html(fig, '{}/timeline.html'.format(opt.groupName))
        
    plt.savefig('{}/timeline.png'.format(opt.groupName))
#     plt.savefig('{0}/timeline_{1}.png'.format(opt.groupName,time.strftime(
#         '%Y%m%dT%H%M%S',time.gmtime())))
        

def createOrderedWaveformFigure(rtable, opt):
    data = np.zeros((len(rtable), int(opt.winlen*2)))
    fig = plt.figure(figsize=(12, 6))
    n=-1
    for r in rtable.iterrows():
        n = n+1
        
        maxwin = max(abs(r['waveform'][r['windowStart']:(r['windowStart']+opt.winlen/2)]))
        
        # Determine padding        
        ppad = int(max(0, opt.ptrig*opt.samprate - r['windowStart']))
        apad = int(max(0, r['windowStart'] - opt.ptrig*opt.samprate - 1))
        
        tmp = r['waveform'][max(0, r['windowStart']-int(
            opt.ptrig*opt.samprate)):min(len(r['waveform']),
            r['windowStart']+int(opt.atrig*opt.samprate))]
            
        tmp = np.hstack((np.zeros(ppad), tmp, np.zeros(apad)))
        data[n, :] = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
            opt.ptrig*opt.samprate + opt.winlen*1.5)]/maxwin
    
    order = rtable.cols.order[:]
    datao = data[order, :]
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(datao, aspect='auto', vmin=-1, vmax=1, interpolation='nearest', cmap='RdBu',
        extent=[-1*opt.winlen*0.5/opt.samprate, opt.winlen*1.5/opt.samprate, n+0.5, -0.5])
    
    clust = rtable.cols.clusterNumber[:]
    clust = clust[order]
    diffclust = np.diff(clust)
    
    for n in range(len(diffclust)):
        if diffclust[n]!=0:
            plt.axhline(y=n+0.5, color='k')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Ordered Event #')
    
    plt.savefig('{}/orderedWaveform.png'.format(opt.groupName))
#     plt.savefig('{0}/orderedWaveform_{1}.png'.format(opt.groupName,time.strftime(
#         '%Y%m%dT%H%M%S',time.gmtime())))


def createCMatrixFigure(rtable, ctable, opt):
    
    order = rtable.cols.order[:]
    
    C = np.zeros((len(rtable),len(rtable)))
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]
    
    # Convert id to row
    rtable_ids = rtable.cols.id[:]
    r = np.zeros((max(rtable_ids)+1,)).astype('int')
    r[rtable_ids] = range(len(rtable_ids))
    C[r[id1], r[id2]] = ctable.cols.ccc[:]
    C = C + C.T + np.eye(len(C))
    Co = C[order, :]
    Co = Co[:, order]

    # Plot the C matrices (saturated below 0.6)
    # Left is ordered by time
    # Right is ordered by cluster ordering
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 2, 1)
    ax.imshow(C, aspect='auto', vmin=0.6, vmax=1, interpolation='nearest')
    ax = fig.add_subplot(1, 2, 2)
    ax.imshow(Co, aspect='auto', vmin=0.6, vmax=1, interpolation='nearest')

    clust = rtable.cols.clusterNumber[:]
    clust = clust[order]
    diffclust = np.diff(clust)
    
    for n in range(len(diffclust)):
        if diffclust[n]!=0:
            plt.axhline(y=n+0.5, color='w')
    
    plt.savefig('{}/cmatrix.png'.format(opt.groupName))
#     plt.savefig('{0}/cmatrix_{1}.png'.format(opt.groupName,time.strftime(
#         '%Y%m%dT%H%M%S',time.gmtime())))


def plotCores(rtable, opt):
    #waveform wiggle plot of input waveforms
    fig = plt.figure(figsize=(15, 10))
    
    cores = rtable.where('isCore==1')
    
    for r in cores:
        dat=r['waveform'][int(
            r['windowStart']-opt.winlen*0.5):int(r['windowStart']+opt.winlen*1.5)]
        maxwin = max(abs(r['waveform'][r['windowStart']:(r['windowStart']+opt.winlen/2)]))
        dat=dat/maxwin
        dat[dat>1] = 1
        dat[dat<-1] = -1
        tvec = np.arange(
            -opt.winlen*0.5/opt.samprate,opt.winlen*1.5/opt.samprate,1/opt.samprate)
        plt.plot(tvec,np.add(dat/2,-1*r['clusterNumber']),'k',linewidth=0.5)
        
    plt.ylabel('Cluster Number')
    plt.xlabel('Time (s)')
    plt.autoscale(tight=True)
    
    mpld3.save_html(fig, '{}/cores.html'.format(opt.groupName))
    
    plt.savefig('{}/cores.png'.format(opt.groupName))
#     plt.savefig('{0}/cores_{1}.png'.format(opt.groupName,time.strftime(
#         '%Y%m%dT%H%M%S',time.gmtime())))
    

def createWigglePlot(jtable, opt):
    #waveform wiggle plot of input waveforms
    fig = plt.figure(figsize=(20, 15))
    n=0.
    for r in jtable.iterrows():
        dat=r['waveform']
        dat=dat/max(dat)
        tvec = np.arange(0,len(dat)*1/opt.samprate,1/opt.samprate)
        plt.plot(tvec,np.add(dat,n))
        n = n+3
    plt.ylabel('index')
    plt.xlabel('time(s)')
    #plt.title('Junk triggers')
    plt.autoscale(tight=True)
    
    mpld3.save_html(fig, '{}/wiggle.html'.format(opt.groupName))
    
    plt.savefig('{}/wiggle.png'.format(opt.groupName))
#     plt.savefig('{0}/wiggle_{1}.png'.format(opt.groupName,time.strftime(
#         '%Y%m%dT%H%M%S',time.gmtime())))
