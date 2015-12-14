from tables import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
from obspy import UTCDateTime
from bokeh.plotting import figure, output_file, save
from bokeh.models import HoverTool, ColumnDataSource, OpenURL, TapTool

"""
These are very brute force plotting; should be replaced with more sophisticated functions
"""
        
def createBokehTimelineFigure(rtable, ctable, opt):
    
    # Run plotCores to ensure thumbnails are up to date
    plotCores(rtable, opt)
    plotFamilies(rtable, ctable, opt)

    output_file('{}/timeline.html'.format(opt.groupName),
        title='{} Timeline'.format(opt.title))

    dt = rtable.cols.startTimeMPL[:]
    cnum = rtable.cols.clusterNumber[:]    
    amp = np.log10(rtable.cols.windowAmp[:])
        
    # Build hover to show an image of the cluster core
    hover = HoverTool(
        tooltips="""
        <div>
        <div>
            <img src="./clusters/@famnum.png" style="height: 100px; width: 500px;
                vertical-align: middle;" />
            <span style="font-size: 9px; font-family: Helvetica;">Cluster ID: </span>
            <span style="font-size: 12px; font-family: Helvetica;">@famnum</span>
        </div>
        </div>
        """, names=["patch"])
    
    TOOLS = [hover,'pan,box_zoom,reset,resize,save,tap']
    p = figure(tools=TOOLS, plot_width=1250, plot_height=700, x_axis_type='datetime')
    p.title = 'Occurrence Timeline (Color by log10(Amplitude) within Family)'
    p.grid.grid_line_alpha = 0.3
    p.xaxis.axis_label = 'Date'
    p.yaxis.axis_label = 'Cluster by Date ({}+ Members)'.format(opt.minplot)
    
    # Steal YlOrRd (len=256) colormap from matplotlib
    colormap = matplotlib.cm.get_cmap('YlOrRd')
    bokehpalette = [matplotlib.colors.rgb2hex(m) for m in colormap(
        np.arange(colormap.N))]
    
    # Figure out earliest member in each family
    mindt = np.zeros((max(cnum)+1,))
    for clustNum in range(max(cnum)+1):
        mindt[clustNum] = min(dt[cnum==clustNum])
    
    # Build the lists and dictionaries    
    n = 0
    for clustNum in np.argsort(mindt):
        if len(dt[cnum==clustNum]) >= opt.minplot:
        
            # Date is required as datenum
            p.line((matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                matplotlib.dates.num2date(max(dt[cnum==clustNum]))), (n, n),
                color='black')            
            
            minamp = min(amp[cnum==clustNum])
            maxamp = max(amp[cnum==clustNum])
            
            if minamp==maxamp:
                colors = [bokehpalette[-1] for i in range(len(amp[cnum==clustNum])+1)]
            else:
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
                 
            # Build source for hover patches
            if n == 0:
                xs=[[matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                    matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                    matplotlib.dates.num2date(max(dt[cnum==clustNum])),
                    matplotlib.dates.num2date(max(dt[cnum==clustNum]))]]
                ys=[[n-0.5, n+0.5, n+0.5, n-0.5]]
                famnum=[clustNum]
            else:
                xs.append([matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                           matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                           matplotlib.dates.num2date(max(dt[cnum==clustNum])),
                           matplotlib.dates.num2date(max(dt[cnum==clustNum]))])
                ys.append([n-0.5, n+0.5, n+0.5, n-0.5])
                famnum.append([clustNum])
            
            n = n+1
    
    if n > 0:
        # Patches allow hovering for image of core and cluster number
        source = ColumnDataSource(data=dict(xs=xs, ys=ys, famnum=famnum))
        p.patches(xs=xs, ys=ys, source=source, name="patch", alpha=0)
        
        # Tapping on one of the patches will open a window to a file with more information
        # on the cluster in question. Temporarily, it leads to only the image, but could
        # lead to an HTML file instead. Need to render those files, of course.
        url = "./clusters/@famnum.html"
        renderer = p.select(name="patch")[0]
        renderer.nonselection_glyph=renderer.glyph.clone()
        taptool = p.select(dict(type=TapTool))[0]
        taptool.names.append("patch")
        taptool.callback = OpenURL(url=url)
    
        save(p) 


def plotCores(rtable, opt):
    # Save cores individually in clusters for timeline hover
    cores = rtable.where('isCore==1')
    for r in cores:
        fig = plt.figure(figsize=(5, 1))
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        
        ppad = int(max(0, opt.ptrig*opt.samprate - r['windowStart']))
        apad = int(max(0, r['windowStart'] - opt.ptrig*opt.samprate - 1))
        tmp = r['waveform'][max(0, r['windowStart']-int(
            opt.ptrig*opt.samprate)):min(len(r['waveform']),
            r['windowStart']+int(opt.atrig*opt.samprate))]
        dat = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
            opt.ptrig*opt.samprate + opt.winlen*1.5)]/r['windowAmp']        
        dat[dat>1] = 1
        dat[dat<-1] = -1
        
        ax.plot(dat,'k',linewidth=0.25)
        plt.autoscale(tight=True)
        plt.savefig('{0}/clusters/{1}.png'.format(opt.groupName,r['clusterNumber']))
        plt.close(fig)


def plotFamilies(rtable, ctable, opt):
    # Save ordered waveforms and correlation matrix for each family, as well as a timeline
    
    # Adjust the font face
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['font.size'] = 8.0
    
    # Get full correlation matrix, sort by order
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
    
    # Get waveform data, sort by order
    n=-1
    data = np.zeros((len(rtable), int(opt.winlen*2)))
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

    datao = data[order, :]
    
    q = 0
    cmap = matplotlib.cm.get_cmap('YlOrRd')
    cmap.set_under('k')
    for cnum in range(max(rtable.cols.clusterNumber[:])+1):
        
        fam = rtable.get_where_list('clusterNumber == {}'.format(cnum))
        core = rtable.get_where_list(
            '(clusterNumber == {}) & (isCore == 1)'.format(cnum))[0]
        
        fig = plt.figure(figsize=(10, 11))
        
        # Plot waveforms
        ax1 = fig.add_subplot(3, 3, (1,2))
        if len(fam) > 12:
            ax1.imshow(datao[q:q+len(fam)], aspect='auto', vmin=-1, vmax=1, cmap='RdBu',
                interpolation='nearest', extent=[-1*opt.winlen*0.5/opt.samprate,
                opt.winlen*1.5/opt.samprate, n + 0.5, -0.5])
            ax1.axvline(x=-0.1*opt.winlen/opt.samprate, color='k', ls='dotted')
            ax1.axvline(x=0.9*opt.winlen/opt.samprate, color='k', ls='dotted')
            ax1.get_yaxis().set_visible(False)
        else:
            for o in range(0, len(fam)):
                dat=datao[o+q,:]
                dat[dat>1] = 1
                dat[dat<-1] = -1
                tvec = np.arange(-opt.winlen*0.5/opt.samprate,opt.winlen*1.5/opt.samprate,
                    1/opt.samprate)
                ax1.plot(tvec,dat/2-o,'k',linewidth=0.25)
                ax1.axvline(x=-0.1*opt.winlen/opt.samprate, color='k', ls='dotted')
                ax1.axvline(x=0.9*opt.winlen/opt.samprate, color='k', ls='dotted')
                ax1.get_yaxis().set_visible(False)
                ax1.autoscale(tight=True)
        ax1.set_xlabel('Time Relative to Trigger (sec)')
        
        # Plot correlation
        ax2 = fig.add_subplot(3, 3, 3)
        ax2.imshow(Co[q:q+len(fam), q:q+len(fam)], cmap=cmap, aspect='auto',
            vmin=opt.cmin, vmax=1, interpolation='nearest')
        ax2.get_yaxis().set_visible(False)
        ax2.set_xlabel('Ordered Event')
        if len(fam) < 5:
            ax2.set_xticks(range(0, len(fam)))
        q = q+len(fam)
        
        # Plot amplitude timeline
        ax3 = fig.add_subplot(3, 3, (4,6))
        ax3.plot_date(rtable[fam]['startTimeMPL'], rtable[fam]['windowAmp'],
                'ro', alpha=0.5, markeredgecolor='r', markeredgewidth=0.5)
        myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
        ax3.xaxis.set_major_formatter(myFmt)
        ax3.set_ylim(1, max(rtable.cols.windowAmp[:])+500)
        ax3.margins(0.05)
        ax3.set_ylabel('Amplitude (Counts)')
        ax3.set_yscale('log')
        
        # Prep catalog
        catalog = np.sort(rtable[fam]['startTimeMPL'])
        longevity = catalog[-1] - catalog[0]
        spacing = np.diff(catalog)*24
        utcatalog = [UTCDateTime(rtable[fam]['startTime'][i]) +
            rtable[fam]['windowStart'][i]/opt.samprate for i in
            np.argsort(rtable[fam]['startTimeMPL'])]
        
        # Plot spacing timeline
        ax4 = fig.add_subplot(3, 3, (7,9)) 
        ax4.plot_date(catalog[1:], spacing, 'ro', alpha=0.5, markeredgecolor='r',
            markeredgewidth=0.5)
        myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
        ax4.xaxis.set_major_formatter(myFmt)
        ax4.set_xlim(ax3.get_xlim())
        ax4.set_ylim(1e-3, max(spacing)*2)
        ax4.margins(0.05)
        ax4.set_ylabel('Time since previous event (hours)')
        ax4.set_xlabel('Date')
        ax4.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig('{0}/clusters/fam{1}.png'.format(opt.groupName,cnum))
        plt.close(fig)
        
        
        # Now write a simple HTML file to show image and catalog
        with open('{0}/clusters/{1}.html'.format(opt.groupName, cnum), 'w') as f:
            f.write("""
            <html><head><title>{1} - Cluster {0}</title>
            </head>
            <body><center>
            <span style="font-size: 20px; font-weight: bold; font-family: Helvetica;">
                Cluster {0}</span></br></br>
            <img src="{0}.png"></br></br>
            <span style="font-size: 12px; font-family: Helvetica;">
                Number of events: {2}</br>
                Longevity: {5:.2f} days</br>
                Mean event spacing: {7:.2f} hours</br>
                Median event spacing: {8:.2f} hours</br></br>
                First event: {3}</br>
                Core event: {6}</br>
                Last event: {4}</br>
                </span> 
            <img src="fam{0}.png"></br></br>
            <span style="font-size: 16px; font-weight: bold; font-family: Helvetica;">
            Aligned Catalog Times</br></br></span>
            <span style="font-size: 12px; font-family: Helvetica;">
            """.format(cnum, opt.title, len(fam), utcatalog[0].isoformat(),
                utcatalog[-1].isoformat(), longevity, (UTCDateTime(
                rtable[core]['startTime']) +
                rtable[core]['windowStart']/opt.samprate).isoformat(), np.mean(spacing),
                np.median(spacing)))
            for u in utcatalog:
                f.write("{}</br>".format(u.isoformat()))
            f.write("""
            </br></span></center></body></html>
            """)
        




# These are old...
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
    
    # mpld3.save_html(fig, '{}/wiggle.html'.format(opt.groupName))
    
    plt.savefig('{}/wiggle.png'.format(opt.groupName))


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
    
    # mpld3.save_html(fig, '{}/timeline.html'.format(opt.groupName))
        
    plt.savefig('{}/timeline.png'.format(opt.groupName))
