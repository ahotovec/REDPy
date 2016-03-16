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
        
def createPlots(rtable, ftable, ttable, opt):
    
    """
    Creates all output plots (core images, family plots, and two bokeh .html plots)
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
        
    """
    
    plotTimelines(rtable, ftable, ttable, opt)
    printCatalog(rtable, ftable, opt)
    plotCores(rtable, ftable, opt)
    plotFamilies(rtable, ftable, opt)
    ftable.cols.printme[:] = np.zeros((len(ftable),))
    ftable.cols.lastprint[:] = np.arange(len(ftable))
    
    # Rename any .tmp files
    tmplist = glob.glob('./{0}/clusters/*.tmp'.format(opt.groupName))
    for tmp in tmplist:
        os.rename(tmp,tmp[0:-4]) 
    
    
def plotTimelines(rtable, ftable, ttable, opt):
    
    """
    Creates the primary .html Bokeh timelines
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    """
    
    dt = rtable.cols.startTimeMPL[:]
    fi = np.mean(rtable.cols.FI[:], axis=1)
    longevity = ftable.cols.longevity[:]
    famstarts = ftable.cols.startTime[:]
    alltrigs = ttable.cols.startTimeMPL[:]
    
    ### OVERVIEW TIMELINES ###
    
    # Create histogram of events/dybin
    histT, hT = np.histogram(alltrigs, bins=np.arange(min(alltrigs),
        max(alltrigs+opt.dybin), opt.dybin))
    histR, hR = np.histogram(dt, bins=np.arange(min(alltrigs),
        max(alltrigs+opt.dybin), opt.dybin))
        
    # Create histogram of events/hrbin
    histTr, hTr = np.histogram(alltrigs, bins=np.arange(max(alltrigs)-opt.recplot,
        max(alltrigs+opt.hrbin/24), opt.hrbin/24))
    histRr, hRr = np.histogram(dt, bins=np.arange(max(alltrigs)-opt.recplot,
        max(alltrigs+opt.hrbin/24), opt.hrbin/24))
    
    oTOOLS = ['pan,box_zoom,reset,resize,save,tap']
    
    o0 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime')
    if opt.dybin>=1:
        o0.title = 'Triggers vs. Repeaters by {:.1f} Day Bin'.format(opt.dybin)
    else:
        o0.title = 'Triggers vs. Repeaters by {:.1f} Hour Bin'.format(opt.dybin*24)
    o0.grid.grid_line_alpha = 0.3
    o0.xaxis.axis_label = 'Date'
    o0.yaxis.axis_label = 'Events'
    
    o0.line(matplotlib.dates.num2date(hT[0:-1]+opt.dybin/2), histT, color='black',
        legend='All Triggers')
    o0.line(matplotlib.dates.num2date(hR[0:-1]+opt.dybin/2), histR, color='red',
        legend='Repeaters')
    o0.legend.orientation = "top_left"
    
    o0r = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime')
    if opt.hrbin<24:
        o0r.title = 'Triggers vs. Repeaters by {:.1f} Hour Bin'.format(opt.hrbin)
    else:
        o0r.title = 'Triggers vs. Repeaters by {:.1f} Day Bin'.format(opt.hrbin/24)
    o0r.grid.grid_line_alpha = 0.3
    o0r.xaxis.axis_label = 'Date'
    o0r.yaxis.axis_label = 'Events'
    
    o0r.line(matplotlib.dates.num2date(hTr[0:-1]+opt.hrbin/48), histTr, color='black',
        legend='All Triggers')
    o0r.line(matplotlib.dates.num2date(hRr[0:-1]+opt.hrbin/48), histRr, color='red',
        legend='Repeaters')
    o0r.legend.orientation = "top_left"
    
    
    o1 = figure(plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0.x_range)
    o1.title = 'Frequency Index vs. Time'
    o1.grid.grid_line_alpha = 0.3
    o1.xaxis.axis_label = 'Date of Repeater'
    o1.yaxis.axis_label = 'FI'
    o1.circle(matplotlib.dates.num2date(dt), fi, color='red', line_alpha=0,
        size=3, fill_alpha=0.5)
        
    o1r = figure(plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0r.x_range)
    o1r.title = 'Frequency Index vs. Time'
    o1r.grid.grid_line_alpha = 0.3
    o1r.xaxis.axis_label = 'Date of Repeater'
    o1r.yaxis.axis_label = 'FI'
    o1r.circle(matplotlib.dates.num2date(dt[dt>(max(alltrigs)-opt.recplot)]),
        fi[dt>(max(alltrigs)-opt.recplot)], color='red', line_alpha=0,
        size=3, fill_alpha=0.5)
    
    
    o2 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0.x_range)
    o2.title = 'Cluster Longevity'
    o2.grid.grid_line_alpha = 0.3
    o2.xaxis.axis_label = 'Start Date'
    o2.yaxis.axis_label = 'Days'
    o2.circle(matplotlib.dates.num2date(famstarts), longevity, color='red',
        line_alpha=0, size=8, fill_alpha=0.5)
    
    o2r = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0r.x_range)
    o2r.title = 'Cluster Longevity'
    o2r.grid.grid_line_alpha = 0.3
    o2r.xaxis.axis_label = 'Start Date'
    o2r.yaxis.axis_label = 'Days'
    o2r.circle(matplotlib.dates.num2date(famstarts[famstarts>(max(alltrigs)-opt.recplot)]),
        longevity[famstarts>(max(alltrigs)-opt.recplot)], color='red', line_alpha=0,
        size=8, fill_alpha=0.5)
    
        
    o = gridplot([[o0],[o1],[o2]])
    o_recent = gridplot([[o0r],[o1r],[o2r]])
        
    output_file('{}/overview.html'.format(opt.groupName),
        title='{} Overview'.format(opt.title))
    save(o)
    
    output_file('{}/overview_recent.html'.format(opt.groupName),
            title='{0} Overview - Last {1:.1f} Days'.format(opt.title,opt.recplot))
    save(o_recent)
    
    
    ### OCCURRENCE TIMELINES ###
    
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
    p1.title = 'Occurrence Timeline (Color by Events per Hour)'
    p1.grid.grid_line_alpha = 0.3
    p1.xaxis.axis_label = 'Date'
    p1.yaxis.axis_label = 'Cluster by Date ({}+ Members)'.format(opt.minplot)
    
    r0 = figure(plot_width=1250, plot_height=250, x_axis_type='datetime')
    if opt.hrbin>24:
        r0.title = 'Last {0} Days: Active Families and Repeaters by {1:.1f} Day Bin'.format(
            opt.recplot,opt.hrbin/24)
    else:
        r0.title = 'Last {0} Days: Active Families and Repeaters by {1:.1f} Hour Bin'.format(
            opt.recplot,opt.hrbin)
    r0.grid.grid_line_alpha = 0.3
    r0.xaxis.axis_label = 'Date'
    r0.yaxis.axis_label = 'Count'
    
    r1 = figure(tools=TOOLSrec, plot_width=1250, plot_height=500, x_axis_type='datetime',
        x_range=r0.x_range)
    r1.title = 'Last {} Days: Occurrence Timeline (Color by Events per Hour)'.format(
        opt.recplot)
    r1.grid.grid_line_alpha = 0.3
    r1.xaxis.axis_label = 'Date'
    r1.yaxis.axis_label = 'Cluster by Date ({}+ Members)'.format(opt.minplot)    
    
    # Steal YlOrRd (len=256) colormap from matplotlib
    colormap = matplotlib.cm.get_cmap('YlOrRd')
    bokehpalette = [matplotlib.colors.rgb2hex(m) for m in colormap(
        np.arange(colormap.N)[::-1])]

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

            # Create histogram of events/hour
            hist, h = np.histogram(dt[members], bins=np.arange(min(dt[members]),
                max(dt[members]+1.0/24), 1.0/24))
            d1 = matplotlib.dates.num2date(h[hist>0])
            d2 = matplotlib.dates.num2date(h[hist>0]+1.0/24)
            histlog = np.log10(hist[hist>0])
            ind = [int(min(255,255*(i/2))) for i in histlog]
            colors = [bokehpalette[i] for i in ind]
            
            p1.quad(top=n+0.3, bottom=n-0.3, left=d1, right=d2,
                color=colors)
            
            
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

                    idx = np.where(h[hist>0]>hr[0])[0]
                        
                else:
                    r1.line((matplotlib.dates.num2date(min(dt[members])),
                        matplotlib.dates.num2date(max(dt[members]))), (m, m),
                        color='black')
                    idx = np.arange(len(d1))
                
                r1.quad(top=m+0.3, bottom=m-0.3, left=np.array(d1)[idx],
                    right=np.array(d2)[idx], color=np.array(colors)[idx])                   
                
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


def plotFamilies(rtable, ftable, opt):

    """
    Creates a multi-paneled family plot.
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Top row: Ordered waveforms
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
    fi = rtable.cols.FI[:]
    
    for cnum in range(ftable.attrs.nClust):

        fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
        core = ftable[cnum]['core']
        
        # Prep catalog
        catalogind = np.argsort(startTimeMPL[fam])
        catalog = startTimeMPL[fam][catalogind]
        longevity = ftable[cnum]['longevity']
        spacing = np.diff(catalog)*24
        minind = fam[catalogind[0]]
        maxind = fam[catalogind[-1]]

        if ftable.cols.printme[cnum] != 0:
        
            fig = plt.figure(figsize=(10, 11))
        
            # Plot waveforms
            ax1 = fig.add_subplot(3, 3, (1,2))
            
            # If only one station, plot all aligned waveforms
            if opt.nsta==1:
            
                famtable = rtable[fam]
                n=-1
                data = np.zeros((len(fam), int(opt.winlen*2)))
                for r in famtable:
                    n = n+1        
                    waveform = r['waveform'][0:opt.wshape]
                    tmp = waveform[max(0, windowStart[fam[n]]-int(
                        opt.ptrig*opt.samprate)):min(opt.wshape,
                        windowStart[fam[n]]+int(opt.atrig*opt.samprate))]
                    data[n, :] = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
                        opt.ptrig*opt.samprate + opt.winlen*1.5)]/windowAmp[fam[n]]
                if len(fam) > 12:
                    ax1.imshow(data, aspect='auto', vmin=-1, vmax=1, cmap='RdBu',
                        interpolation='nearest', extent=[-1*opt.winlen*0.5/opt.samprate,
                        opt.winlen*1.5/opt.samprate, n + 0.5, -0.5])
                else:
                    tvec = np.arange(
                        -opt.winlen*0.5/opt.samprate,opt.winlen*1.5/opt.samprate,
                        1/opt.samprate)
                    for o in range(len(fam)):
                        dat=data[o,:]
                        dat[dat>1] = 1
                        dat[dat<-1] = -1
                        ax1.plot(tvec,dat/2-o,'k',linewidth=0.25)
            
            # Otherwise, plot cores and stacks from all stations            
            else:
            
                r = rtable[core]
                famtable = rtable[fam]
                tvec = np.arange(-opt.winlen*0.5/opt.samprate,opt.winlen*1.5/opt.samprate,
                    1/opt.samprate)
                for s in range(opt.nsta):
                    
                    dats = np.zeros((int(opt.winlen*2),))
                    waveform = famtable['waveform'][:,s*opt.wshape:(s+1)*opt.wshape]
                    for n in range(len(fam)):
                        tmps = waveform[n, max(0, windowStart[fam[n]]-int(
                            opt.ptrig*opt.samprate)):min(opt.wshape,
                            windowStart[fam[n]]+int(
                            opt.atrig*opt.samprate))]/(famtable['windowAmp'][
                            n,s]+1.0/1000)
                        tmps[tmps>1] = 1
                        tmps[tmps<-1] = -1
                        dats = dats + tmps[int(opt.ptrig*opt.samprate -
                            opt.winlen*0.5):int(opt.ptrig*opt.samprate + opt.winlen*1.5)]
                    dats = dats/(max(dats)+1.0/1000)
                    dats[dats>1] = 1
                    dats[dats<-1] = -1
                    ax1.plot(tvec,dats-1.75*s,'r',linewidth=1)
                    
                    waveformc = r['waveform'][s*opt.wshape:(s+1)*opt.wshape]
                    tmpc = waveformc[max(0, r['windowStart']-int(
                        opt.ptrig*opt.samprate)):min(opt.wshape,
                        r['windowStart']+int(opt.atrig*opt.samprate))]
                    datc = tmpc[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
                        opt.ptrig*opt.samprate + opt.winlen*1.5)]/(
                        r['windowAmp'][s]+1.0/1000)        
                    datc[datc>1] = 1
                    datc[datc<-1] = -1
                    ax1.plot(tvec,datc-1.75*s,'k',linewidth=0.25)
            
            ax1.axvline(x=-0.1*opt.winlen/opt.samprate, color='k', ls='dotted')
            ax1.axvline(x=0.9*opt.winlen/opt.samprate, color='k', ls='dotted')
            ax1.get_yaxis().set_visible(False)
            ax1.autoscale(tight=True)
            ax1.set_xlabel('Time Relative to Trigger (sec)')
            
            # Plot mean FFT
            ax2 = fig.add_subplot(3, 3, 3)
            ax2.set_xlabel('Frequency (Hz)')
            ax2.get_yaxis().set_visible(False)
            r = rtable[core]
            famtable = rtable[fam]
            freq = np.linspace(0,opt.samprate/2,opt.winlen/2)
            fftc = np.zeros((opt.winlen/2,))
            fftm = np.zeros((opt.winlen/2,))
            for s in range(opt.nsta):
                fft = np.abs(np.real(r['windowFFT'][s*opt.winlen:s*opt.winlen+opt.winlen/2]))
                fft = fft/np.amax(fft)
                fftc = fftc+fft
                ffts = np.mean(np.abs(np.real(
                    famtable['windowFFT'][:,s*opt.winlen:s*opt.winlen+opt.winlen/2])),
                    axis=0)
                fftm = fftm + ffts/np.amax(ffts)
            ax2.plot(freq,fftm,'r', linewidth=1)
            ax2.plot(freq,fftc,'k', linewidth=0.25)
            ax2.set_xlim(0,opt.fmax*1.5)
            
            # Plot amplitude timeline
            ax3 = fig.add_subplot(3, 3, (4,6))
            ax3.plot_date(startTimeMPL[fam], windowAmp[fam],
                    'ro', alpha=0.5, markeredgecolor='r', markeredgewidth=0.5)
            myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
            ax3.xaxis.set_major_formatter(myFmt)
            ax3.set_ylim(1, max(rtable.cols.windowAmp[:][:,opt.printsta])+500)
            ax3.margins(0.05)
            ax3.set_ylabel('Amplitude (Counts)')
            ax3.set_yscale('log')
        
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
                cnum), dpi=100)
            plt.close(fig)
        
        if ftable.cols.printme[cnum] != 0 or ftable.cols.lastprint[cnum] != cnum:
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
                    Median event spacing: {8:.2f} hours</br>
                    Mean Frequency Index: {9:.2f}<br></br>
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
                    np.mean(spacing), np.median(spacing), np.mean(np.mean(fi[fam],
                    axis=1))))
                                
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
