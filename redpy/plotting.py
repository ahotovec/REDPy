from tables import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import mpld3
from bokeh.plotting import figure, output_file, save, ColumnDataSource
from bokeh.models import HoverTool

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
        if len(dt[cnum==clustNum]) >= opt.minplot:
            plt.plot_date((min(dt[cnum==clustNum]), max(dt[cnum==clustNum])), (n, n),
                'k:', linewidth=0.5)
            plt.scatter(dt[cnum==clustNum], n*np.ones((len(dt[cnum==clustNum]),)),
                c=(1-reach[cnum==clustNum]), vmin=0.6, vmax=1, cmap='jet', linewidth=0.5)
            plt.text(max(dt[cnum==clustNum]), n+0.2, format(len(dt[cnum==clustNum])))
            n = n+1
    
    plt.margins(0.05)
    
    plt.gcf().autofmt_xdate()
    plt.xlabel('Date')
    plt.ylabel('Cluster by Date (with at least {} members)'.format(opt.minplot))
    
    mpld3.save_html(fig, '{}/timeline.html'.format(opt.groupName))
        
    plt.savefig('{}/timeline.png'.format(opt.groupName))
#     plt.savefig('{0}/timeline_{1}.png'.format(opt.groupName,time.strftime(
#         '%Y%m%dT%H%M%S',time.gmtime())))
        
        
def createBokehTimelineFigure(rtable, ctable, opt):

    output_file('{}/timelineBokeh.html'.format(opt.groupName),
        title='{} Timeline'.format(opt.title))

    dt = rtable.cols.startTimeMPL[:]
    cnum = rtable.cols.clusterNumber[:]
#     reach = rtable.cols.reachability[:]
#     
#     # Fix some values of reach
#     reach[0] = 0.3
#     reach[reach==1] = 0
#     reach[reach==-1] = 0.3
#     reach[reach>0.3] = 0.3
    
    amp = np.log10(rtable.cols.windowAmp[:])
    
    # Figure out earliest member in each family
    mindt = np.zeros((max(cnum)+1,))
    for clustNum in range(max(cnum)+1):
        mindt[clustNum] = min(dt[cnum==clustNum])
     
    TOOLS = "pan,box_zoom,reset,resize,save"
    p = figure(tools=TOOLS, plot_width=1250, plot_height=700, x_axis_type='datetime')
    p.title = 'Occurrence Timeline (Color by log10(Amplitude) within Family)'
    p.grid.grid_line_alpha = 0.3
    p.xaxis.axis_label = 'Date'
    p.yaxis.axis_label = 'Cluster by Date - {}+ Members'.format(opt.minplot)
    
    # Steal colormap from matplotlib
    colormap = matplotlib.cm.get_cmap('YlOrRd')
    bokehpalette = [matplotlib.colors.rgb2hex(m) for m in colormap(
        np.arange(colormap.N))]
        
    n = 0
    for clustNum in np.argsort(mindt):
        if len(dt[cnum==clustNum]) >= opt.minplot:
        
            # Date is required as datenum
            p.line((matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                matplotlib.dates.num2date(max(dt[cnum==clustNum]))), (n, n),
                color='black')            
#             ind = [int(255*((1-reach[i])-0.7)/(0.3)) for i in np.where(cnum==clustNum)[0]]

            minamp = min(amp[cnum==clustNum])
            maxamp = max(amp[cnum==clustNum])
            ind = [int(255*((amp[i]-minamp)/(maxamp-minamp))) for i in np.where(
                cnum==clustNum)[0]]
            colors = [bokehpalette[i] for i in ind]
               
            d = matplotlib.dates.num2date(dt[cnum==clustNum])                        
            p.circle(d, n, color=colors, size=8, line_color='black', line_width=0.5,
                fill_alpha=1.0)
            
            # Text doesn't understand datetimes, need to convert to a number and subtract
            # about 8 hours
            p.text(time.mktime(matplotlib.dates.num2date(
                 max(dt[cnum==clustNum])).timetuple())*1000 - 28799000, n,
                 text=['   {}'.format(len(dt[cnum==clustNum]))], text_font_size='9pt',
                 text_baseline='middle')
            
            
            
            n = n+1
        
    save(p)     


def createOrderedWaveformFigure(rtable, opt):
    data = np.zeros((len(rtable), int(opt.winlen*2)))
    fig = plt.figure(figsize=(12, 6))
    n=-1
    for r in rtable.iterrows():
        n = n+1
        
        # Determine padding        
        ppad = int(max(0, opt.ptrig*opt.samprate - r['windowStart']))
        apad = int(max(0, r['windowStart'] - opt.ptrig*opt.samprate - 1))
        
        tmp = r['waveform'][max(0, r['windowStart']-int(
            opt.ptrig*opt.samprate)):min(len(r['waveform']),
            r['windowStart']+int(opt.atrig*opt.samprate))]
            
        tmp = np.hstack((np.zeros(ppad), tmp, np.zeros(apad)))
        data[n, :] = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
            opt.ptrig*opt.samprate + opt.winlen*1.5)]/r['windowAmp']
    
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
        dat=dat/r['windowAmp']
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
