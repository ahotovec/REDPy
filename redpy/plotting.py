from tables import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import redpy.cluster
import os
import glob
from obspy import UTCDateTime
from bokeh.plotting import figure, output_file, save, gridplot
from bokeh.models import HoverTool, ColumnDataSource, OpenURL, TapTool, Range1d
        
def createBokehTimelineFigure(rtable, ctable, ftable, opt):
    
    """
    Creates all output plots (core images, family plots, and two bokeh .html plots)
    
    rtable: Repeater table
    ctable: Correlation matrix table
    ftable: Families table
    opt: Options object describing station/run parameters
        
    """
    
    dt = rtable.cols.startTimeMPL[:]
    amp = np.log10(rtable.cols.windowAmp[:][:,opt.printsta])
    
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
    for clustNum in range(ftable.attrs.nClust):
        
        members = np.fromstring(ftable[clustNum]['members'], dtype=int, sep=' ')
        t1 = np.argmin(np.abs(min(dt[members])-dy))
        t2 = np.argmin(np.abs(max(dt[members])-dy))+1
        dyfams[t1:t2] = dyfams[t1:t2]+1
        if max(dt[members])>hr[0]:
            t1 = np.argmin(np.abs(min(dt[members])-hr))
            t2 = np.argmin(np.abs(max(dt[members])-hr))+1
            hrfams[t1:t2] = hrfams[t1:t2]+1
        for d in dt[members]:
            t0 = np.argmin(np.abs(d-dy))
            dyrept[t0] = dyrept[t0]+1
            if d>hr[0]:
                t0 = np.argmin(np.abs(d-hr))
                hrrept[t0] = hrrept[t0]+1
        
        if len(dt[members]) >= opt.minplot:
        
            # Date is required as datenum
            p1.line((matplotlib.dates.num2date(min(dt[members])),
                matplotlib.dates.num2date(max(dt[members]))), (n, n),
                color='black')
            
            minamp = min(amp[members])
            maxamp = max(amp[members])
            
            if minamp==maxamp:
                colors = [bokehpalette[-1] for i in range(len(amp[members])+1)]
            else:
                ind = [int(255*((amp[i]-minamp)/(maxamp-minamp))) for i in members]
                colors = [bokehpalette[i] for i in ind]
               
            d = matplotlib.dates.num2date(dt[members])                        
            p1.circle(d, n, color=colors, size=8, line_color='black', line_width=0.5,
                fill_alpha=1.0)
            
            # Text doesn't understand datetimes, need to convert to a number and subtract
            # about 8 hours
            p1.text(time.mktime(matplotlib.dates.num2date(
                 max(dt[members])).timetuple())*1000 - 28799000, n,
                 text=['   {}'.format(len(dt[members]))], text_font_size='9pt',
                 text_baseline='middle')
                 
            # Build source for hover patches
            fnum = clustNum
            if n == 0:
                xs=[[matplotlib.dates.num2date(min(dt[members])),
                    matplotlib.dates.num2date(min(dt[members])),
                    matplotlib.dates.num2date(max(dt[members])),
                    matplotlib.dates.num2date(max(dt[members]))]]
                ys=[[n-0.5, n+0.5, n+0.5, n-0.5]]
                famnum=[fnum]
            else:
                xs.append([matplotlib.dates.num2date(min(dt[members])),
                           matplotlib.dates.num2date(min(dt[members])),
                           matplotlib.dates.num2date(max(dt[members])),
                           matplotlib.dates.num2date(max(dt[members]))])
                ys.append([n-0.5, n+0.5, n+0.5, n-0.5])
                famnum.append([fnum])
            
            n = n+1
            
            
            if max(dt[members])>hr[0]:
            
                if min(dt[members])<hr[0]:
                    r1.line((matplotlib.dates.num2date(hr[0]-opt.hrbin/6),
                    matplotlib.dates.num2date(max(dt[members]))), (m, m),
                        color='black')
                    r1.text(time.mktime(matplotlib.dates.num2date(
                        hr[0]-opt.hrbin/6).timetuple())*1000 - 28799000, m, text=['<'], 
                        text_font_size='9pt', text_baseline='middle')
                    dtmp = dt[members]
                    atmp = amp[members]
                    amps = atmp[dtmp>hr[0]]
                    d = matplotlib.dates.num2date(dtmp[dtmp>hr[0]])

                    if minamp==maxamp:
                        colors = [bokehpalette[-1] for i in range(len(d)+1)]
                    else:
                        ind = [int(255*((i-minamp)/(maxamp-minamp))) for i in amps]
                        colors = [bokehpalette[i] for i in ind]
                        
                else:
                    r1.line((matplotlib.dates.num2date(min(dt[members])),
                        matplotlib.dates.num2date(max(dt[members]))), (m, m),
                        color='black')
                    d = matplotlib.dates.num2date(dt[members])   
                
                r1.circle(d, m, color=colors, size=8, line_color='black', line_width=0.5,
                    fill_alpha=1.0)                     
                
                # Text doesn't understand datetimes, need to convert to a number and subtract
                # about 8 hours
                r1.text(time.mktime(matplotlib.dates.num2date(
                    max(dt[members])).timetuple())*1000 - 28799000, m,
                    text=['   {}'.format(len(dt[members]))], text_font_size='9pt',
                    text_baseline='middle')
                 
                # Build source for hover patches
                fnumr = clustNum
                if m == 0:
                    xsr=[[matplotlib.dates.num2date(max(min(dt[members]),hr[0])),
                        matplotlib.dates.num2date(max(min(dt[members]),hr[0])),
                        matplotlib.dates.num2date(max(dt[members])),
                        matplotlib.dates.num2date(max(dt[members]))]]
                    ysr=[[m-0.5, m+0.5, m+0.5, m-0.5]]
                    famnumr=[fnumr]
                else:
                    xsr.append([matplotlib.dates.num2date(max(min(dt[members]),hr[0])),
                               matplotlib.dates.num2date(max(min(dt[members]),hr[0])),
                               matplotlib.dates.num2date(max(dt[members])),
                               matplotlib.dates.num2date(max(dt[members]))])
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
        renderer.nonselection_glyph=renderer.glyph.clone()
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
        renderer.nonselection_glyph=renderer.glyph.clone()
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

    # Run plotCores to ensure thumbnails are up to date
    printCatalog(rtable, ftable, opt)
    plotCores(rtable, ftable, opt)
    plotFamilies(rtable, ctable, ftable, opt)
    ftable.cols.printme[:] = np.zeros((len(ftable),))
    ftable.cols.lastprint[:] = np.arange(len(ftable))
    
    # Rename .tmp files
    tmplist = glob.glob('./{0}/clusters/*.tmp')
    for tmp in tmplist:
        os.rename(tmp,tmp[0:-4])        


def plotCores(rtable, ftable, opt):

    """
    Plots core waveforms as .gif for hovering in timeline
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    """
    
    for n in range(len(ftable))[::-1]:
        if ftable.cols.lastprint[n] != n and ftable.cols.printme[n] == 0:
            os.rename('{0}/clusters/{1}.gif'.format(opt.groupName,
                ftable.cols.lastprint[n]), '{0}/clusters/{1}.gif.tmp'.format(
                opt.groupName, n))
            os.rename('{0}/clusters/fam{1}.png'.format(opt.groupName,
                ftable.cols.lastprint[n]), '{0}/clusters/fam{1}.png.tmp'.format(
                opt.groupName, n))
    
    cores = rtable[ftable.cols.core[:]]
    n = -1
    for r in cores:
        n = n+1
        if ftable.cols.printme[n] == 1:
            fig = plt.figure(figsize=(5, 1))
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            
            waveform = r['waveform'][opt.printsta*opt.wshape:(opt.printsta+1)*opt.wshape]
            tmp = waveform[max(0, r['windowStart']-int(
                opt.ptrig*opt.samprate)):min(opt.wshape,
                r['windowStart']+int(opt.atrig*opt.samprate))]
            dat = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
                opt.ptrig*opt.samprate + opt.winlen*1.5)]/r['windowAmp'][opt.printsta]        
            dat[dat>1] = 1
            dat[dat<-1] = -1
        
            ax.plot(dat,'k',linewidth=0.25)
            plt.autoscale(tight=True)
            plt.savefig('{0}/clusters/{1}.gif'.format(opt.groupName,n),
                dpi=100)
            plt.close(fig)


def plotFamilies(rtable, ctable, ftable, opt):

    """
    Creates a multi-paneled family plot.
    
    rtable: Repeater table
    ctable: Correlation matrix table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Top row: Ordered waveforms and ordered similarity matrix
    Middle row: Timeline of amplitude
    Bottom row: Timeline of event spacing 
    """
    
    # Adjust the font face
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['font.size'] = 8.0
    
    # Load into memory
    startTimeMPL = rtable.cols.startTimeMPL[:]
    startTime = rtable.cols.startTime[:]
    windowAmp = rtable.cols.windowAmp[:][:,opt.printsta]
    windowStart = rtable.cols.windowStart[:]

    # Get waveform data
    ### THIS CAN PROBABLY BE OPTIMIZED ###
    n=-1
    data = np.zeros((len(rtable), int(opt.winlen*2)))
    for r in rtable.iterrows():
        n = n+1
        
        # Determine padding        
        ppad = int(max(0, opt.ptrig*opt.samprate - windowStart[n]))
        apad = int(max(0, windowStart[n] - opt.ptrig*opt.samprate - 1))
        
        waveform = r['waveform'][opt.printsta*opt.wshape:(opt.printsta+1)*opt.wshape]
        tmp = waveform[max(0, windowStart[n]-int(
            opt.ptrig*opt.samprate)):min(opt.wshape,
            windowStart[n]+int(opt.atrig*opt.samprate))]
            
        tmp = np.hstack((np.zeros(ppad), tmp, np.zeros(apad)))
        data[n, :] = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
            opt.ptrig*opt.samprate + opt.winlen*1.5)]/windowAmp[n]

    cmap = matplotlib.cm.get_cmap('YlOrRd')
    cmap.set_under('k')
    for cnum in range(ftable.attrs.nClust):

        fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')

        if ftable.cols.printme[cnum] != 0:
            core = ftable[cnum]['core']
        
            fig = plt.figure(figsize=(10, 11))
        
            # Plot waveforms
            ax1 = fig.add_subplot(3, 1, 1)
            if len(fam) > 12:
                ax1.imshow(data[fam], aspect='auto', vmin=-1, vmax=1, cmap='RdBu',
                    interpolation='nearest', extent=[-1*opt.winlen*0.5/opt.samprate,
                    opt.winlen*1.5/opt.samprate, n + 0.5, -0.5])
                ax1.axvline(x=-0.1*opt.winlen/opt.samprate, color='k', ls='dotted')
                ax1.axvline(x=0.9*opt.winlen/opt.samprate, color='k', ls='dotted')
                ax1.get_yaxis().set_visible(False)
            else:
                for o in range(len(fam)):
                    dat=data[fam[o],:]
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
        
            # Plot amplitude timeline
            ax3 = fig.add_subplot(3, 1, 2)
            ax3.plot_date(startTimeMPL[fam], windowAmp[fam],
                    'ro', alpha=0.5, markeredgecolor='r', markeredgewidth=0.5)
            myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
            ax3.xaxis.set_major_formatter(myFmt)
            ax3.set_ylim(1, max(rtable.cols.windowAmp[:][:,opt.printsta])+500)
            ax3.margins(0.05)
            ax3.set_ylabel('Amplitude (Counts)')
            ax3.set_yscale('log')
        
            # Prep catalog
            catalogind = np.argsort(startTimeMPL[fam])
            catalog = startTimeMPL[fam][catalogind]
            longevity = catalog[-1] - catalog[0]
            spacing = np.diff(catalog)*24
            minind = fam[catalogind[0]]
            maxind = fam[catalogind[-1]]
        
            # Plot spacing timeline
            ax4 = fig.add_subplot(3, 1, 3) 
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
                cnum), dpi=100)
            plt.close(fig)
        
        # Now write a simple HTML file to show image and catalog
        with open('{0}/clusters/{1}.html'.format(opt.groupName, cnum), 'w') as f:
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
            """.format(cnum, opt.title, len(fam), (UTCDateTime(
                startTime[minind]) + windowStart[minind]/opt.samprate).isoformat(),
                (UTCDateTime(startTime[maxind]) + windowStart[
                maxind]/opt.samprate).isoformat(), longevity, (UTCDateTime(
                startTime[core]) + windowStart[core]/opt.samprate).isoformat(),
                np.mean(spacing), np.median(spacing)))
                            
            f.write("""
            </center></body></html>
            """)
        

def printCatalog(rtable, ftable, opt):
    """
    Prints flat catalog to text file
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Note: Time in text file corresponds to current trigger time by alignment
    """

    with open('{}/catalog.txt'.format(opt.groupName), 'w') as f:
        
        startTimes = rtable.cols.startTime[:]
        windowStarts = rtable.cols.windowStart[:]
        
        for cnum in range(ftable.attrs.nClust):
            fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
            for i in np.argsort(startTimes[fam]):
                f.write("{0} {1}\n".format(cnum,(UTCDateTime(startTimes[fam][i]) +
                    windowStarts[fam][i]/opt.samprate).isoformat()))
