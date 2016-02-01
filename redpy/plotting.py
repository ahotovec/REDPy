from tables import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import redpy.cluster
from obspy import UTCDateTime
from bokeh.plotting import figure, output_file, save, gridplot
from bokeh.models import HoverTool, ColumnDataSource, OpenURL, TapTool, Range1d

"""
These are very brute force plotting; should be replaced with more sophisticated functions
"""
        
def createBokehTimelineFigure(rtable, ctable, ftable, opt):
    
    # Run OPTICS to ensure everything is up to date
    redpy.cluster.runFullOPTICS(rtable, ctable, ftable, opt)
        
    # Run plotCores to ensure thumbnails are up to date
    plotCores(rtable, ftable, opt)
    plotFamilies(rtable, ctable, ftable, opt)
    printCatalog(rtable, ftable, opt)
    
    # Update lastClust column
    rtable.cols.lastClust[:] = rtable.cols.clusterNumber[:]
    rtable.flush()

    dt = rtable.cols.startTimeMPL[:]
    cnum = rtable.cols.clusterNumber[:]   
    amp = np.log10(rtable.cols.windowAmp[:])
    
    dy = np.arange(np.floor(min(dt)/opt.dybin),np.ceil(max(dt+opt.dybin)/opt.dybin))*opt.dybin
    dyfams = np.zeros((len(dy),))
    dyrept = np.zeros((len(dy),))
    hr = np.arange(np.floor(max(dt-opt.recplot)/(opt.hrbin/24)),np.ceil(
        max(dt+opt.hrbin/24)/(opt.hrbin/24)))*(opt.hrbin/24)
    hrfams = np.zeros((len(hr),))
    hrrept = np.zeros((len(hr),))
        
    # Build hover to show an image of the cluster core
    hover = HoverTool(
        tooltips="""
        <div>
        <div>
            <img src="./clusters/@famnum.gif" style="height: 100px; width: 500px;
                vertical-align: middle;" />
            <span style="font-size: 9px; font-family: Helvetica;">Cluster ID: </span>
            <span style="font-size: 12px; font-family: Helvetica;">@famnum</span>
        </div>
        </div>
        """, names=["patch"])
    
    TOOLS = [hover,'pan,box_zoom,reset,resize,save,tap']
    
    # Build hover to show an image of the cluster core
    hoverr = HoverTool(
        tooltips="""
        <div>
        <div>
            <img src="./clusters/@famnum.gif" style="height: 100px; width: 500px;
                vertical-align: middle;" />
            <span style="font-size: 9px; font-family: Helvetica;">Cluster ID: </span>
            <span style="font-size: 12px; font-family: Helvetica;">@famnum</span>
        </div>
        </div>
        """, names=["patchr"])
    
    TOOLSrec = [hoverr,'pan,box_zoom,reset,resize,save,tap']
    
    p0 = figure(plot_width=1250, plot_height=250, x_axis_type='datetime')
    if opt.dybin>=1:
        p0.title = 'Active Families and Repeaters by {:.1f} Day Bin'.format(opt.dybin)
    else:
        p0.title = 'Active Families and Repeaters by {:.1f} Hour Bin'.format(opt.dybin*24)
    p0.grid.grid_line_alpha = 0.3
    p0.xaxis.axis_label = 'Date'
    p0.yaxis.axis_label = 'Count'
    
    p1 = figure(tools=TOOLS, plot_width=1250, plot_height=500, x_axis_type='datetime',
        x_range=p0.x_range)
    p1.title = 'Occurrence Timeline (Color by log10(Amplitude) within Family)'
    p1.grid.grid_line_alpha = 0.3
    p1.xaxis.axis_label = 'Date'
    p1.yaxis.axis_label = 'Cluster by Date ({}+ Members)'.format(opt.minplot)
    
    r0 = figure(plot_width=1250, plot_height=250, x_axis_type='datetime')
    if opt.hrbin<1:
        r0.title = 'Last {0} Days: Active Families and Repeaters by {1:.1f} Day Bin'.format(
            opt.recplot,opt.hrbin*24)
    else:
        r0.title = 'Last {0} Days: Active Families and Repeaters by {1:.1f} Hour Bin'.format(
            opt.recplot,opt.hrbin)
    r0.grid.grid_line_alpha = 0.3
    r0.xaxis.axis_label = 'Date'
    r0.yaxis.axis_label = 'Count'
    
    r1 = figure(tools=TOOLSrec, plot_width=1250, plot_height=500, x_axis_type='datetime',
        x_range=r0.x_range)
    r1.title = 'Last {} Days: Occurrence Timeline (Color by log10(Amplitude) within Family)'.format(
        opt.recplot)
    r1.grid.grid_line_alpha = 0.3
    r1.xaxis.axis_label = 'Date'
    r1.yaxis.axis_label = 'Cluster by Date ({}+ Members)'.format(opt.minplot)    
    
    # Steal YlOrRd (len=256) colormap from matplotlib
    colormap = matplotlib.cm.get_cmap('YlOrRd')
    bokehpalette = [matplotlib.colors.rgb2hex(m) for m in colormap(
        np.arange(colormap.N))]

    # Build the lists and dictionaries    
    n = 0  
    m = 0  
    for clustNum in np.unique(cnum):
        
        t1 = np.argmin(np.abs(min(dt[cnum==clustNum])-dy))
        t2 = np.argmin(np.abs(max(dt[cnum==clustNum])-dy))+1
        dyfams[t1:t2] = dyfams[t1:t2]+1
        if max(dt[cnum==clustNum])>hr[0]:
            t1 = np.argmin(np.abs(min(dt[cnum==clustNum])-hr))
            t2 = np.argmin(np.abs(max(dt[cnum==clustNum])-hr))+1
            hrfams[t1:t2] = hrfams[t1:t2]+1
        for d in dt[cnum==clustNum]:
            t0 = np.argmin(np.abs(d-dy))
            dyrept[t0] = dyrept[t0]+1
            if d>hr[0]:
                t0 = np.argmin(np.abs(d-hr))
                hrrept[t0] = hrrept[t0]+1
        
        if len(dt[cnum==clustNum]) >= opt.minplot:
        
            # Date is required as datenum
            p1.line((matplotlib.dates.num2date(min(dt[cnum==clustNum])),
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
            p1.circle(d, n, color=colors, size=8, line_color='black', line_width=0.5,
                fill_alpha=1.0)
            
            # Text doesn't understand datetimes, need to convert to a number and subtract
            # about 8 hours
            p1.text(time.mktime(matplotlib.dates.num2date(
                 max(dt[cnum==clustNum])).timetuple())*1000 - 28799000, n,
                 text=['   {}'.format(len(dt[cnum==clustNum]))], text_font_size='9pt',
                 text_baseline='middle')
                 
            # Build source for hover patches
            fnum = clustNum
            if n == 0:
                xs=[[matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                    matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                    matplotlib.dates.num2date(max(dt[cnum==clustNum])),
                    matplotlib.dates.num2date(max(dt[cnum==clustNum]))]]
                ys=[[n-0.5, n+0.5, n+0.5, n-0.5]]
                famnum=[fnum]
            else:
                xs.append([matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                           matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                           matplotlib.dates.num2date(max(dt[cnum==clustNum])),
                           matplotlib.dates.num2date(max(dt[cnum==clustNum]))])
                ys.append([n-0.5, n+0.5, n+0.5, n-0.5])
                famnum.append([fnum])
            
            n = n+1
            
            
            if max(dt[cnum==clustNum])>hr[0]:
            
                if min(dt[cnum==clustNum])<hr[0]:
                    r1.line((matplotlib.dates.num2date(hr[0]-opt.hrbin/6),
                    matplotlib.dates.num2date(max(dt[cnum==clustNum]))), (m, m),
                        color='black')
                    r1.text(time.mktime(matplotlib.dates.num2date(
                        hr[0]-opt.hrbin/6).timetuple())*1000 - 28799000, m, text=['<'], 
                        text_font_size='9pt', text_baseline='middle')
                    dtmp = dt[cnum==clustNum]
                    atmp = amp[cnum==clustNum]
                    amps = atmp[dtmp>hr[0]]
                    d = matplotlib.dates.num2date(dtmp[dtmp>hr[0]])

                    if minamp==maxamp:
                        colors = [bokehpalette[-1] for i in range(len(d)+1)]
                    else:
                        ind = [int(255*((i-minamp)/(maxamp-minamp))) for i in amps]
                        colors = [bokehpalette[i] for i in ind]
                        
                else:
                    r1.line((matplotlib.dates.num2date(min(dt[cnum==clustNum])),
                        matplotlib.dates.num2date(max(dt[cnum==clustNum]))), (m, m),
                        color='black')
                    d = matplotlib.dates.num2date(dt[cnum==clustNum])   
                
                r1.circle(d, m, color=colors, size=8, line_color='black', line_width=0.5,
                    fill_alpha=1.0)                     
                
                # Text doesn't understand datetimes, need to convert to a number and subtract
                # about 8 hours
                r1.text(time.mktime(matplotlib.dates.num2date(
                    max(dt[cnum==clustNum])).timetuple())*1000 - 28799000, m,
                    text=['   {}'.format(len(dt[cnum==clustNum]))], text_font_size='9pt',
                    text_baseline='middle')
                 
                # Build source for hover patches
                fnumr = clustNum
                if m == 0:
                    xsr=[[matplotlib.dates.num2date(max(min(dt[cnum==clustNum]),hr[0])),
                        matplotlib.dates.num2date(max(min(dt[cnum==clustNum]),hr[0])),
                        matplotlib.dates.num2date(max(dt[cnum==clustNum])),
                        matplotlib.dates.num2date(max(dt[cnum==clustNum]))]]
                    ysr=[[m-0.5, m+0.5, m+0.5, m-0.5]]
                    famnumr=[fnumr]
                else:
                    xsr.append([matplotlib.dates.num2date(max(min(dt[cnum==clustNum]),hr[0])),
                               matplotlib.dates.num2date(max(min(dt[cnum==clustNum]),hr[0])),
                               matplotlib.dates.num2date(max(dt[cnum==clustNum])),
                               matplotlib.dates.num2date(max(dt[cnum==clustNum]))])
                    ysr.append([m-0.5, m+0.5, m+0.5, m-0.5])
                    famnumr.append([fnumr])
                m = m+1
            
    if n > 0:
        # Patches allow hovering for image of core and cluster number
        source = ColumnDataSource(data=dict(xs=xs, ys=ys, famnum=famnum))
        p1.patches(xs=xs, ys=ys, source=source, name="patch", alpha=0)
        
        # Tapping on one of the patches will open a window to a file with more information
        # on the cluster in question.
        url = "./clusters/@famnum.html"
        renderer = p1.select(name="patch")[0]
        renderer.nonselection_glyph=renderer.selection_glyph
        taptool = p1.select(dict(type=TapTool))[0]
        taptool.names.append("patch")
        taptool.callback = OpenURL(url=url)
        
        p0.line(matplotlib.dates.num2date(dy), dyfams, color='red', line_width=1.5,
            legend='Families')
        p0.line(matplotlib.dates.num2date(dy), dyrept, color='black', line_width=0.5,
            legend='Repeaters')
        p0.legend.orientation = "top_left"
        
        if n > 30:
            p1.set(plot_height=n*15, y_range=Range1d(-1, n))
        p = gridplot([[p0],[p1]])
        
        output_file('{}/timeline.html'.format(opt.groupName),
            title='{} Timeline'.format(opt.title))
        save(p)
        
        
        sourcer = ColumnDataSource(data=dict(xs=xsr, ys=ysr, famnum=famnumr))
        r1.patches(xs=xsr, ys=ysr, source=sourcer, name="patchr", alpha=0)
        
        renderer = r1.select(name="patchr")[0]
        renderer.nonselection_glyph=renderer.selection_glyph
        taptool = r1.select(dict(type=TapTool))[0]
        taptool.names.append("patchr")
        taptool.callback = OpenURL(url=url)
        
        r0.line(matplotlib.dates.num2date(hr), hrfams, color='red', line_width=1.5,
            legend='Families')
        r0.line(matplotlib.dates.num2date(hr), hrrept, color='black', line_width=0.5,
            legend='Repeaters')
        r0.legend.orientation = "top_left"
        
        if m > 30:
            r1.set(plot_height=m*15, y_range=Range1d(-1, m))
        
        output_file('{}/timeline_recent.html'.format(opt.groupName),
            title='{0} Timeline - Last {1:.1f} Days'.format(opt.title,opt.recplot))
        r = gridplot([[r0],[r1]])
        save(r)
         


def plotCores(rtable, ftable, opt):
    # Save cores individually in clusters for timeline hover
    cores = rtable[ftable.attrs.cores]
    for r in cores:
        if r['lastClust'] != r['clusterNumber']:
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
            plt.savefig('{0}/clusters/{1}.gif'.format(opt.groupName,r['clusterNumber']),
                dpi=100)
            plt.close(fig)


def plotFamilies(rtable, ctable, ftable, opt):
    # Save ordered waveforms and correlation matrix for each family, as well as a timeline
    
    # Adjust the font face
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['font.size'] = 8.0
    
    # Get full correlation matrix, sort by order
    order = rtable.cols.order[:]
    
    C = np.eye(len(rtable))
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]
    
    # Convert id to row
    rtable_ids = rtable.cols.id[:]
    r = np.zeros((max(rtable_ids)+1,)).astype('int')
    r[rtable_ids] = range(len(rtable_ids))
    C[r[id1], r[id2]] = ctable.cols.ccc[:]
    C[r[id2], r[id1]] = ctable.cols.ccc[:]
    
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
    for cnum in np.argsort(ftable.cols.pnum[0:max(rtable.cols.clusterNumber[:])+1]):

        fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')

        if sum(rtable[fam]['clusterNumber'] - rtable[fam]['lastClust']) != 0:
            core = np.intersect1d(fam, ftable.attrs.cores)[0]
        
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
            catalogind = np.argsort(rtable[fam]['startTimeMPL'])
            catalog = rtable[fam]['startTimeMPL'][catalogind]
            longevity = catalog[-1] - catalog[0]
            spacing = np.diff(catalog)*24
            minind = catalogind[0]
            maxind = catalogind[-1]
        
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
            plt.savefig('{0}/clusters/fam{1}.png'.format(opt.groupName,
                rtable.cols.clusterNumber[core]), dpi=100)
            plt.close(fig)
        
        
            # Now write a simple HTML file to show image and catalog
            with open('{0}/clusters/{1}.html'.format(opt.groupName,
                rtable.cols.clusterNumber[core]), 'w') as f:
                f.write("""
                <html><head><title>{1} - Cluster {0}</title>
                </head>
                <body><center>
                <span style="font-size: 20px; font-weight: bold; font-family: Helvetica;">
                    Cluster {0}</span></br></br>
                <img src="{0}.gif"></br></br>
                <span style="font-size: 12px; font-family: Helvetica;">
                    Number of events: {2}</br>
                    Longevity: {5:.2f} days</br>
                    Mean event spacing: {7:.2f} hours</br>
                    Median event spacing: {8:.2f} hours</br></br>
                    First event: {3}</br>
                    Core event: {6}</br>
                    Last event: {4}</br>
                    </span> 
                <img src="fam{0}.png"></br>                
                """.format(rtable.cols.clusterNumber[core], opt.title, len(fam), (UTCDateTime(
                    rtable[minind]['startTime']) +
                    rtable[minind]['windowStart']/opt.samprate).isoformat(), (UTCDateTime(
                    rtable[maxind]['startTime']) +
                    rtable[maxind]['windowStart']/opt.samprate).isoformat(), longevity, (UTCDateTime(
                    rtable[core]['startTime']) +
                    rtable[core]['windowStart']/opt.samprate).isoformat(), np.mean(spacing),
                    np.median(spacing)))
                                
                f.write("""
                </center></body></html>
                """)
        else:
                q = q+len(fam)
        

def printCatalog(rtable, ftable, opt):
    """
    A thoroughly inefficient way of printing out the catalog...
    """

    with open('{}/catalog.txt'.format(opt.groupName), 'w') as f:
        
        startTimes = rtable.cols.startTime[:]
        windowStarts = rtable.cols.windowStart[:]
        
        for cnum in np.unique(rtable.cols.clusterNumber[:]):
            fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
            for i in np.argsort(startTimes[fam]):
                f.write("{0} {1}\n".format(cnum,(UTCDateTime(startTimes[fam][i]) +
                    windowStarts[fam][i]/opt.samprate).isoformat()))
