# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2018  Alicia Hotovec-Ellis (ahotovec@gmail.com)
# Licensed under GNU GPLv3 (see LICENSE.txt)

from tables import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates
import time
import redpy.cluster
import redpy.correlation
from redpy.optics import *
import os
import shutil
import glob
import urllib
from obspy import UTCDateTime
from obspy.geodetics import locations2degrees
from obspy.taup import TauPyModel
import pandas as pd
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.transforms import offset_copy
from bokeh.plotting import figure, output_file, save, gridplot
from bokeh.models import HoverTool, ColumnDataSource, OpenURL, TapTool, Range1d, Div
from bokeh.models import Arrow, VeeHead, ColorBar, LogColorMapper, LogTicker, LabelSet
from bokeh.layouts import column
try:
    import urllib2
except:
    pass
try:
    import urllib.request
except:
    pass
        
def createPlots(rtable, ftable, ttable, ctable, otable, opt):
    
    """
    Creates all output plots (core images, family plots, and two bokeh .html plots)
    
    rtable: Repeater table
    ftable: Families table
    ttable: Triggers table
    ctable: Correlation table
    otable: Orphan table
    opt: Options object describing station/run parameters
        
    """
    
    printOrphanCatalog(otable, opt)
    if len(rtable)>1:
        plotTimelines(rtable, ftable, ttable, opt)
        if np.sum(ftable.cols.printme[:]):
            if opt.printVerboseCat == True:
                printVerboseCatalog(rtable, ftable, ctable, opt)
            else:
                printCatalog(rtable, ftable, opt)
            printCoresCatalog(rtable, ftable, opt)
            printEventsperDay(rtable, ftable, opt)
            plotCores(rtable, ftable, opt)
            plotFamilies(rtable, ftable, ctable, opt)
            ftable.cols.printme[:] = np.zeros((len(ftable),))
            ftable.cols.lastprint[:] = np.arange(len(ftable))
    else:
        print('Nothing to plot!')
    
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
    fi = np.nanmean(rtable.cols.FI[:], axis=1)
    longevity = ftable.cols.longevity[:]
    famstarts = ftable.cols.startTime[:]
    alltrigs = ttable.cols.startTimeMPL[:]
    
    
    # Create histogram of events/dybin
    histT, hT = np.histogram(alltrigs, bins=np.arange(min(alltrigs),
        max(alltrigs+opt.dybin), opt.dybin))
    histR, hR = np.histogram(dt, bins=np.arange(min(alltrigs),
        max(alltrigs+opt.dybin), opt.dybin))
        
    # Determine padding for hover bars (~1% of window range on each side)
    barpad = (max(alltrigs)-min(alltrigs))*0.01
    barpadr = opt.recplot*0.01
        
    # Create histogram of events/hrbin
    histTr, hTr = np.histogram(alltrigs, bins=np.arange(max(alltrigs)-opt.recplot,
        max(alltrigs+opt.hrbin/24), opt.hrbin/24))
    histRr, hRr = np.histogram(dt, bins=np.arange(max(alltrigs)-opt.recplot,
        max(alltrigs+opt.hrbin/24), opt.hrbin/24))
    
    oTOOLS = ['pan,box_zoom,reset,save,tap']
    
    if opt.dybin>=1:
        o0title = 'Repeaters vs. Orphans by {:.1f} Day Bin'.format(opt.dybin)
    else:
        o0title = 'Repeaters vs. Orphans by {:.1f} Hour Bin'.format(opt.dybin*24)
    
    o0 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title=o0title)
    o0.grid.grid_line_alpha = 0.3
    o0.xaxis.axis_label = 'Date'
    o0.yaxis.axis_label = 'Events'
    
    o0.line(matplotlib.dates.num2date(hT[0:-1]+opt.dybin/2), histT-histR, color='black',
        legend='Orphans')
    o0.line(matplotlib.dates.num2date(hR[0:-1]+opt.dybin/2), histR, color='red',
        legend='Repeaters', line_width=2)
    o0.legend.location = 'top_left'
    
    if opt.hrbin<24:
        o0rtitle = 'Repeaters vs. Orphans by {1:.1f} Hour Bin'.format(
            opt.recplot, opt.hrbin)
    else:
        o0rtitle = 'Repeaters vs. Orphans by {1:.1f} Day Bin'.format(
            opt.recplot, opt.hrbin/24)
    o0r = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title=o0rtitle)
    o0r.grid.grid_line_alpha = 0.3
    o0r.xaxis.axis_label = 'Date'
    o0r.yaxis.axis_label = 'Events'
    
    o0r.line(matplotlib.dates.num2date(hTr[0:-1]+opt.hrbin/48), histTr-histRr,
        color='black', legend='Orphans')
    o0r.line(matplotlib.dates.num2date(hRr[0:-1]+opt.hrbin/48), histRr, color='red',
        legend='Repeaters', line_width=2)
    o0r.legend.location = 'top_left'
        
    o1 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0.x_range, title='Frequency Index')
    o1.grid.grid_line_alpha = 0.3
    o1.xaxis.axis_label = 'Date'
    o1.yaxis.axis_label = 'FI'
    o1.circle(matplotlib.dates.num2date(dt), fi, color='red', line_alpha=0,
        size=3, fill_alpha=0.5)
        
    o1r = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0r.x_range, title='Frequency Index')
    o1r.grid.grid_line_alpha = 0.3
    o1r.xaxis.axis_label = 'Date'
    o1r.yaxis.axis_label = 'FI'
    # Put invisible points in for case that there are no events
    o1r.circle(matplotlib.dates.num2date(hTr[0:2]), [1, 1], line_alpha=0, fill_alpha=0)
    o1r.circle(matplotlib.dates.num2date(dt[dt>(max(alltrigs)-opt.recplot)]),
        fi[dt>(max(alltrigs)-opt.recplot)], color='red', line_alpha=0,
        size=3, fill_alpha=0.5)
        
    o2 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0.x_range, y_axis_type='log', y_range=[0.1,
        np.sort(alltrigs)[-1]-np.sort(alltrigs)[0]], title='Cluster Longevity')
    o2.grid.grid_line_alpha = 0.3
    o2.xaxis.axis_label = 'Date'
    o2.yaxis.axis_label = 'Days'
    for n in range(len(famstarts)):
        o2.line((matplotlib.dates.num2date(famstarts[n]), matplotlib.dates.num2date(
            famstarts[n]+longevity[n])), (longevity[n], longevity[n]), color='red',
            line_alpha=0.5)
        
    o2r = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        x_range=o0r.x_range, y_axis_type='log', y_range=[0.1,
        np.sort(alltrigs)[-1]-np.sort(alltrigs)[0]], title='Cluster Longevity')
    o2r.grid.grid_line_alpha = 0.3
    o2r.xaxis.axis_label = 'Date'
    o2r.yaxis.axis_label = 'Days'
    # Put invisible points in for case that there are no events
    o2r.circle(matplotlib.dates.num2date(hTr[0:2]), [1, 1], line_alpha=0, fill_alpha=0)
    for n in range(len(famstarts)):
        if (max(alltrigs)-opt.recplot)<=famstarts[n]:
            o2r.line((matplotlib.dates.num2date(famstarts[n]), matplotlib.dates.num2date(
                famstarts[n]+longevity[n])), (longevity[n], longevity[n]), color='red',
                line_alpha=0.5)
        elif (max(alltrigs)-opt.recplot)<=famstarts[n]+longevity[n]:
            o2r.add_layout(Arrow(end=VeeHead(size=3, line_color='red', line_alpha=0.5),
                    line_color='red', line_alpha=0.5,
                    x_start=matplotlib.dates.num2date(famstarts[n]+longevity[n]),
                    x_end=matplotlib.dates.num2date(hTr[0]-0.5),
                    y_start=longevity[n], y_end=longevity[n]))
            
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
    
    TOOLS = [hover,'pan,box_zoom,reset,save,tap']
    
    # Build hover to show an image of the cluster core
    hoverr = HoverTool(
        tooltips="""
        <div>
        <div>
            <img src="./clusters/@famnum.png" style="height: 100px; width: 500px;
                vertical-align: middle;" />
            <span style="font-size: 9px; font-family: Helvetica;">Cluster ID: </span>
            <span style="font-size: 12px; font-family: Helvetica;">@famnum</span>
        </div>
        </div>
        """, names=["patchr"])
    
    TOOLSrec = [hoverr,'pan,box_zoom,reset,save,tap']
        
    p1 = figure(tools=TOOLS, plot_width=1250, plot_height=500, x_axis_type='datetime',
        x_range=o0.x_range, title='Occurrence Timeline')
    p1.grid.grid_line_alpha = 0.3
    p1.xaxis.axis_label = 'Date'
    p1.yaxis.axis_label = 'Cluster by Date ({}+ Members)'.format(opt.minplot)
    
    r1 = figure(tools=TOOLSrec, plot_width=1250, plot_height=500, x_axis_type='datetime',
        x_range=o0r.x_range, title = 'Occurrence Timeline')
    r1.grid.grid_line_alpha = 0.3
    r1.xaxis.axis_label = 'Date'
    r1.yaxis.axis_label = 'Cluster by Date'    
    
    # Steal YlOrRd (len=256) colormap from matplotlib
    colormap = matplotlib.cm.get_cmap('YlOrRd')
    bokehpalette = [matplotlib.colors.rgb2hex(m) for m in colormap(
        np.arange(colormap.N)[::-1])]

    # Build the lists and dictionaries    
    n = 0  
    m = 0
    cloc1 = 335
    cloc2 = 335
      
    for clustNum in range(ftable.attrs.nClust):
        
        members = np.fromstring(ftable[clustNum]['members'], dtype=int, sep=' ')
        
        # Create histogram of events/hour
        hist, h = np.histogram(dt[members], bins=np.arange(min(dt[members]),
            max(dt[members]+1.0/24), 1.0/24))
        d1 = matplotlib.dates.num2date(h[np.where(hist>0)])
        d2 = matplotlib.dates.num2date(h[np.where(hist>0)]+1.0/24)
        histlog = np.log10(hist[hist>0])
        ind = [int(min(255,255*(i/2))) for i in histlog]
        colors = [bokehpalette[i] for i in ind]
        
        if len(dt[members]) >= opt.minplot:
        
            # Date is required as datenum
            p1.line((matplotlib.dates.num2date(min(dt[members])),
                matplotlib.dates.num2date(max(dt[members]))), (n, n),
                color='black')
            
            p1.quad(top=n+0.3, bottom=n-0.3, left=d1, right=d2,
                color=colors)            
            
            p1.add_layout(LabelSet(x=matplotlib.dates.num2date(
                max(h[np.where(hist>0)]+1.0/24)),
                y=n, text=['{}'.format(len(dt[members]))], level='glyph',
                x_offset=5, y_offset=0, render_mode='canvas', text_font_size='9pt',
                text_baseline='middle'))
                 
            # Build source for hover patches
            fnum = clustNum
            if n == 0:
                xs=[[matplotlib.dates.num2date(min(dt[members])-barpad),
                    matplotlib.dates.num2date(min(dt[members])-barpad),
                    matplotlib.dates.num2date(max(dt[members])+barpad),
                    matplotlib.dates.num2date(max(dt[members])+barpad)]]
                ys=[[n-0.5, n+0.5, n+0.5, n-0.5]]
                famnum=[[fnum]]
            else:
                xs.append([matplotlib.dates.num2date(min(dt[members])-barpad),
                    matplotlib.dates.num2date(min(dt[members])-barpad),
                    matplotlib.dates.num2date(max(dt[members])+barpad),
                    matplotlib.dates.num2date(max(dt[members])+barpad)])
                ys.append([n-0.5, n+0.5, n+0.5, n-0.5])
                famnum.append([fnum])
            
            n = n+1
            
        if max(dt[members])>hRr[0]:
            
            if min(dt[members])<hRr[0]:
                
                r1.line((matplotlib.dates.num2date(hTr[0]),
                    matplotlib.dates.num2date(max(dt[members]))), (m, m),
                    color='black')
                r1.add_layout(Arrow(end=VeeHead(size=3),
                    x_start=matplotlib.dates.num2date(hTr[0]+0.01),
                    x_end=matplotlib.dates.num2date(hTr[0]-0.5),
                    y_start=m, y_end=m))

                idx = np.where(h[np.where(hist>0)[0]]>hRr[0])[0]
                        
            else:
                r1.line((matplotlib.dates.num2date(min(dt[members])),
                    matplotlib.dates.num2date(max(dt[members]))), (m, m),
                    color='black')
                idx = np.arange(len(d1))
                
            r1.quad(top=m+0.3, bottom=m-0.3, left=np.array(d1)[idx],
                right=np.array(d2)[idx], color=np.array(colors)[idx])                   
                
            r1.add_layout(LabelSet(x=np.array(d2)[-1],
                y=m, text=['{}'.format(len(dt[members]))], level='glyph',
                x_offset=5, y_offset=0, render_mode='canvas', text_font_size='9pt',
                text_baseline='middle'))
            
            # Build source for hover patches
            fnumr = clustNum
            if m == 0:
                xsr=[[matplotlib.dates.num2date(max(min(dt[members]),hRr[0])-barpadr),
                    matplotlib.dates.num2date(max(min(dt[members]),hRr[0])-barpadr),
                    matplotlib.dates.num2date(max(dt[members])+barpadr),
                    matplotlib.dates.num2date(max(dt[members])+barpadr)]]
                ysr=[[m-0.5, m+0.5, m+0.5, m-0.5]]
                famnumr=[[fnumr]]
            else:
                xsr.append([matplotlib.dates.num2date(max(min(dt[members]),hRr[0])-
                    barpadr), matplotlib.dates.num2date(max(min(dt[members]),hRr[0])-
                    barpadr),matplotlib.dates.num2date(max(dt[members])+barpadr),
                    matplotlib.dates.num2date(max(dt[members])+barpadr)])
                ysr.append([m-0.5, m+0.5, m+0.5, m-0.5])
                famnumr.append([fnumr])
            m = m+1
            
    if n > 0:
        # Patches allow hovering for image of core and cluster number
        source = ColumnDataSource(data=dict(xs=xs, ys=ys, famnum=famnum))
        p1.patches('xs', 'ys', source=source, name='patch', alpha=0,
            selection_fill_alpha=0, selection_line_alpha=0, nonselection_fill_alpha=0,
            nonselection_line_alpha=0)
                
        # Tapping on one of the patches will open a window to a file with more information
        # on the cluster in question.        
        url = './clusters/@famnum.html'
        renderer = p1.select(name='patch')
        taptool = p1.select(type=TapTool)[0]
        taptool.names.append('patch')
        taptool.callback = OpenURL(url=url)
        
        if n > 30:
            p1.plot_height = n*15
            p1.y_range = Range1d(-1, n)
            cloc1 = n*15-165
    
    else:
        p1.circle(matplotlib.dates.num2date(hTr[0:2]), [0, 0], line_alpha=0, fill_alpha=0)
        
    if m > 0:
        sourcer = ColumnDataSource(data=dict(xs=xsr, ys=ysr, famnum=famnumr))
        r1.patches('xs', 'ys', source=sourcer, name='patchr', alpha=0,
            selection_fill_alpha=0, selection_line_alpha=0, nonselection_fill_alpha=0,
            nonselection_line_alpha=0)
        
        url = './clusters/@famnum.html'
        renderer = r1.select(name='patchr')
        taptool = r1.select(type=TapTool)[0]
        taptool.names.append('patchr')
        taptool.callback = OpenURL(url=url)
                    
        if m > 30:
            r1.plot_height = m*15
            r1.y_range = Range1d(-1, m)
            cloc2 = m*15-165
        
    else: 
        r1.circle(matplotlib.dates.num2date(hTr[0:2]), [0, 0], line_alpha=0, fill_alpha=0)
    
    color_mapper = LogColorMapper(palette=bokehpalette, low=1, high=100)
    color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
        border_line_color='#eeeeee', location=(7,cloc1), orientation='horizontal',
        width=100, height=15, title='Events per Hour', padding=15,
        major_tick_line_alpha=0)
    color_bar2 = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
        border_line_color='#eeeeee', location=(7,cloc2), orientation='horizontal',
        width=100, height=15, title='Events per Hour', padding=15,
        major_tick_line_alpha=0)
    p1.add_layout(color_bar)
    r1.add_layout(color_bar2)
    
    o = gridplot([[Div(text='<h1>{0}</h1>'.format(
                       opt.title), width=1000)],[o0],[o1],[p1],[o2]])
    o_recent = gridplot([[Div(text='<h1>{0} (Last {1:.1f} Days)</h1>'.format(
                              opt.title,opt.recplot), width=1000)],[o0r],[o1r],[r1],[o2r]])
        
    output_file('{}/overview.html'.format(opt.groupName),
        title='{} Overview'.format(opt.title))
    save(o)
    
    output_file('{}/overview_recent.html'.format(opt.groupName),
            title='{0} Overview - Last {1:.1f} Days'.format(opt.title,opt.recplot))
    save(o_recent)


def plotCores(rtable, ftable, opt):

    """
    Plots core waveforms as .png for hovering in timeline
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    """
    
    for n in range(len(ftable))[::-1]:
        if ftable.cols.lastprint[n] != n and ftable.cols.printme[n] == 0:
            os.rename('{0}/clusters/{1}.png'.format(opt.groupName,
                ftable.cols.lastprint[n]), '{0}/clusters/{1}.png.tmp'.format(
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
            plt.savefig('{0}/clusters/{1}.png'.format(opt.groupName,n),
                dpi=100)
            plt.close(fig)


def plotFamilies(rtable, ftable, ctable, opt):

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
    ids = rtable.cols.id[:]
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]
    ccc = ctable.cols.ccc[:]
    
    # Station names
    stas = opt.station.split(',')
    chas = opt.channel.split(',')
    
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
        
            fig = plt.figure(figsize=(10, 12))
        
            # Plot waveforms
            ax1 = fig.add_subplot(9, 3, (1,8))
            
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
                    tvec = [-1*opt.winlen*0.5/opt.samprate, opt.winlen*1.5/opt.samprate]
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
                        try:
                            dats = dats + tmps[int(opt.ptrig*opt.samprate -
                                opt.winlen*0.5):int(opt.ptrig*opt.samprate +
                                opt.winlen*1.5)]
                        except ValueError:
                           pass
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
                    ax1.text(np.min(tvec)-0.1,-1.75*s,'{0}\n{1}'.format(stas[s],chas[s]),
                        horizontalalignment='right', verticalalignment='center')
            
            ax1.axvline(x=-0.1*opt.winlen/opt.samprate, color='k', ls='dotted')
            ax1.axvline(x=0.9*opt.winlen/opt.samprate, color='k', ls='dotted')
            ax1.get_yaxis().set_visible(False)
            ax1.set_xlim((np.min(tvec),np.max(tvec)))
            if opt.nsta > 1:
                ax1.set_ylim((-1.75*s-1,1))
            ax1.set_xlabel('Time Relative to Trigger (seconds)', style='italic')
            
            # Plot mean FFT
            ax2 = fig.add_subplot(9, 3, (3,9))
            ax2.set_xlabel('Frequency (Hz)', style='italic')
            ax2.get_yaxis().set_visible(False)
            r = rtable[core]
            famtable = rtable[fam]
            freq = np.linspace(0,opt.samprate/2,opt.winlen/2)
            fftc = np.zeros((int(opt.winlen/2),))
            fftm = np.zeros((int(opt.winlen/2),))
            for s in range(opt.nsta):
                fft = np.abs(np.real(r['windowFFT'][int(
                    s*opt.winlen):int(s*opt.winlen+opt.winlen/2)]))
                fft = fft/(np.amax(fft)+1.0/1000)
                fftc = fftc+fft
                ffts = np.mean(np.abs(np.real(famtable['windowFFT'][:,int(
                    s*opt.winlen):int(s*opt.winlen+opt.winlen/2)])),axis=0)
                fftm = fftm + ffts/(np.amax(ffts)+1.0/1000)
            ax2.plot(freq,fftm,'r', linewidth=1)
            ax2.plot(freq,fftc,'k', linewidth=0.25)
            ax2.set_xlim(0,opt.fmax*1.5)
            ax2.legend(['Stack','Core'], loc='upper right', frameon=False)
            
            # Plot amplitude timeline
            ax3 = fig.add_subplot(9, 3, (10,15))
            ax3.plot_date(catalog, windowAmp[fam[catalogind]],
                    'ro', alpha=0.5, markeredgecolor='r', markeredgewidth=0.5,
                    markersize=3)
            myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
            ax3.xaxis.set_major_formatter(myFmt)
            ax3.set_ylim(1, max(rtable.cols.windowAmp[:][:,opt.printsta])+500)
            ax3.margins(0.05)
            ax3.set_ylabel('Amplitude (Counts)', style='italic')
            ax3.set_xlabel('Date', style='italic')
            ax3.set_yscale('log')
        
            # Plot spacing timeline
            ax4 = fig.add_subplot(9, 3, (16,21)) 
            ax4.plot_date(catalog[1:], spacing, 'ro', alpha=0.5, markeredgecolor='r',
                markeredgewidth=0.5, markersize=3)
            myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
            ax4.xaxis.set_major_formatter(myFmt)
            ax4.set_xlim(ax3.get_xlim())
            ax4.set_ylim(1e-3, max(spacing)*2)
            ax4.margins(0.05)
            ax4.set_ylabel('Time since previous event (hours)', style='italic')
            ax4.set_xlabel('Date', style='italic')
            ax4.set_yscale('log')
            
            # Plot correlation timeline
            idf = ids[fam]
            ix = np.where(np.in1d(id2,idf))
            r = np.zeros((max(idf)+1,)).astype('int')
            r[idf] = range(len(idf))
            C = np.zeros((len(idf),len(idf)))
            C[r[id2[ix]],r[id1[ix]]] = ccc[ix]
            C[r[id1[ix]],r[id2[ix]]] = ccc[ix]
            C[range(len(idf)),range(len(idf))] = 1.0
            
            ax5 = fig.add_subplot(9, 3, (22,27))
            ax5.plot_date(catalog, C[np.argmax(np.sum(C,0)),:], 'ro', alpha=0.5,
                markeredgecolor='r', markeredgewidth=0.5, markersize=3)
            ax5.plot_date(catalog, C[np.argmax(np.sum(C,0)),:]+opt.cmin, 'wo', alpha=0.5,
                markeredgecolor='r', markeredgewidth=0.5)
            myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
            ax5.xaxis.set_major_formatter(myFmt)
            ax5.set_xlim(ax3.get_xlim())
            ax5.set_ylim(opt.cmin-0.02, 1.02)
            ax5.margins(0.05)
            ax5.set_ylabel('Cross-correlation coefficient', style='italic')
            ax5.set_xlabel('Date', style='italic')
        
            plt.tight_layout()
            plt.savefig('{0}/clusters/fam{1}.png'.format(opt.groupName,
                cnum), dpi=100)
            plt.close(fig)
        
        if ftable.cols.printme[cnum] != 0 or ftable.cols.lastprint[cnum] != cnum:
            if cnum>0:
                prev = "<a href='{0}.html'>&lt; Cluster {0}</a>".format(cnum-1)
            else:
                prev = " "
            if cnum<len(ftable)-1:
                next = "<a href='{0}.html'>Cluster {0} &gt;</a>".format(cnum+1)
            else:
                next = " "   
            # Now write a simple HTML file to show image and catalog
            with open('{0}/clusters/{1}.html'.format(opt.groupName, cnum), 'w') as f:
                f.write("""
                <html><head><title>{1} - Cluster {0}</title>
                </head><style>
                a {{color:red;}}
                body {{font-family:Helvetica; font-size:12px}}
                h1 {{font-size: 20px;}}
                </style>
                <body><center>
                {10} &nbsp; | &nbsp; {11}</br>             
                <h1>Cluster {0}</h1>                
                <img src="{0}.png" width=500 height=100></br></br>
                    Number of events: {2}</br>
                    Longevity: {5:.2f} days</br>
                    Mean event spacing: {7:.2f} hours</br>
                    Median event spacing: {8:.2f} hours</br>
                    Mean Frequency Index: {9:.2f}<br></br>
                    First event: {3}</br>
                    Core event: {6}</br>
                    Last event: {4}</br>
                <img src="fam{0}.png"></br>                
                """.format(cnum, opt.title, len(fam), (UTCDateTime(
                    startTime[minind]) + windowStart[minind]/opt.samprate).isoformat(),
                    (UTCDateTime(startTime[maxind]) + windowStart[
                    maxind]/opt.samprate).isoformat(), longevity, (UTCDateTime(
                    startTime[core]) + windowStart[core]/opt.samprate).isoformat(),
                    np.mean(spacing), np.median(spacing), np.mean(np.nanmean(fi[fam],
                    axis=1)),prev,next))
                
                if opt.checkComCat:
                    checkComCat(rtable, ftable, cnum, f, startTime, windowStart, opt)
                           
                f.write("""
                </center></body></html>
                """)
        

def checkComCat(rtable, ftable, cnum, f, startTime, windowStart, opt):
    """
    Checks repeater trigger times with projected arrival times from ANSS Comprehensive
    Earthquake Catalog (ComCat) and writes these to HTML and image files. Will also
    check NCEDC catalog if location is near Northern California.
    
    rtable: Repeater table
    ftable: Families table
    cnum: cluster number to check
    f: HTML file to write to
    startTime: startTime column from rtable (convenience)
    windowStart: windowStart column from rtable (convenience)
    opt: Options object describing station/run parameters
    
    Traces through iasp91 global velocity model; checks for local, regional, and
    teleseismic matches for limited set of phase arrivals
    """
    
    pc = ['Potential', 'Conflicting']
    model = TauPyModel(model="iasp91")
    mc = 0
    n = 0
    l = 0
    stalats = np.array(opt.stalats.split(',')).astype(float)
    stalons = np.array(opt.stalons.split(',')).astype(float)
    latc = np.mean(stalats)
    lonc = np.mean(stalons)

    members = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
    order = np.argsort(startTime[members])
    matchstring = ('</br><b>ComCat matches:</b></br>'
        '<div style="overflow-y: auto; height:100px; width:1200px;">')
    
    for m in members[order]:
        t = UTCDateTime(startTime[m])+windowStart[m]/opt.samprate
        cc_url = ('http://earthquake.usgs.gov/fdsnws/event/1/query?'
                  'starttime={}&endtime={}&format=text').format(t-1800,t+30)
        try:
            comcat = pd.read_csv(cc_url,delimiter='|')
            otime = comcat['Time'].tolist()
            lat = comcat['Latitude'].tolist()
            lon = comcat['Longitude'].tolist()
            dep = comcat['Depth/km'].tolist()
            mag = comcat['Magnitude'].tolist()
            place = comcat['EventLocationName'].tolist()
        except (urllib.error.HTTPError, urllib2.HTTPError):
            otime = []
            lat = []
            lon = []
            dep = []
            mag = []
            place = []
        
        # Check if near Northern California, then go to NCEDC for additional events but
        # for shorter time interval
        if latc > 34 and latc < 42 and lonc > -124 and lonc < -116:
            cc_urlnc = ('http://ncedc.org/fdsnws/event/1/query?'
                        'starttime={}&endtime={}&format=text').format((t-60).isoformat(),
                        (t+30).isoformat())
            try:
                ncedc = pd.read_csv(cc_urlnc,delimiter='|')
                otime.extend(ncedc[' Time '].tolist())
                lat.extend(ncedc[' Latitude '].tolist())
                lon.extend(ncedc[' Longitude '].tolist())
                dep.extend(ncedc[' Depth/km '].tolist())
                mag.extend(ncedc[' Magnitude '].tolist())
                place.extend(ncedc[' EventLocationName'].tolist())
            except ValueError:
                pass
        
        n0 = 0
        for c in range(len(otime)):
            deg = locations2degrees(lat[c],lon[c],latc,lonc)
            dt = t-UTCDateTime(otime[c])
        
            if deg <= opt.locdeg:
                mc += 1
                if np.remainder(mc,100) == 0:
                    model = TauPyModel(model="iasp91")
                arrivals = model.get_travel_times(source_depth_in_km=max(0,dep[c]),
                    distance_in_degree=deg, phase_list=['p','s','P','S'])
                if len(arrivals) > 0:
                    pt = np.zeros((len(arrivals),))
                    pname = []
                    for a in range(len(arrivals)):
                        pt[a] = arrivals[a].time - dt
                        pname.append(arrivals[a].name)
                    if np.min(abs(pt)) < opt.serr:
                        amin = np.argmin(abs(pt))
                        matchstring+=('{} local match: {} ({:5.3f}, {:6.3f}) {:3.1f}km '
                            'M{:3.2f} - {} - ({}) {:4.2f} s</br>').format(pc[n0],otime[c],
                            lat[c],lon[c],dep[c],mag[c],place[c],pname[amin],pt[amin])
                        n0 = 1
                        l = l+1
                        if l == 1:
                            llats = np.array(lat[c])
                            llons = np.array(lon[c])
                            ldeps = np.array(dep[c])
                        else:
                            llats = np.append(llats,lat[c])
                            llons = np.append(llons,lon[c])
                            ldeps = np.append(ldeps,dep[c])
            elif deg <= opt.regdeg and mag[c] >= opt.regmag:
                mc += 1
                if np.remainder(mc,100) == 0:
                    model = TauPyModel(model="iasp91")
                arrivals = model.get_travel_times(source_depth_in_km=max(0,dep[c]),
                    distance_in_degree=deg, phase_list=['p','s','P','S','PP','SS'])
                if len(arrivals) > 0:
                    pt = np.zeros((len(arrivals),))
                    pname = []
                    for a in range(len(arrivals)):
                        pt[a] = arrivals[a].time - dt
                        pname.append(arrivals[a].name)
                    if np.min(abs(pt)) < opt.serr:
                        amin = np.argmin(abs(pt))
                        matchstring+=('<div style="color:red">{} regional match: {} '
                            '({:5.3f}, {:6.3f}) {:3.1f}km M{:3.2f} - {} - ({}) {:4.2f} '
                            's</div>').format(pc[n0],otime[c],lat[c],lon[c],dep[c],
                            mag[c],place[c],pname[amin],pt[amin])
                        n0 = 1
            elif deg > opt.regdeg and mag[c] >= opt.telemag:
                mc += 1
                if np.remainder(mc,100) == 0:
                    model = TauPyModel(model="iasp91")
                arrivals = model.get_travel_times(source_depth_in_km=max(0,dep[c]),
                    distance_in_degree=deg, phase_list=['P','S','PP','SS','PcP','ScS',
                        'PKiKP','PKIKP'])
                if len(arrivals) > 0:
                    pt = np.zeros((len(arrivals),))
                    pname = []
                    for a in range(len(arrivals)):
                        pt[a] = arrivals[a].time - dt
                        pname.append(arrivals[a].name)
                    if np.min(abs(pt)) < opt.serr:
                        amin = np.argmin(abs(pt))
                        matchstring+=('<div style="color:red">{} teleseismic match: {} '
                            '({:5.3f}, {:3.1f}) {:4.2f}km M{:3.2f} - {} - ({}) {:4.2f} '
                            's</div>').format(pc[n0],otime[c],lat[c],lon[c],dep[c],
                            mag[c],place[c],pname[amin],pt[amin])
                        n0 = 1
        if n0>1:
            n = n+1
        else:
            n = n+n0
    if n>0:
        matchstring+='</div>'
        matchstring+='Total potential matches: {}</br>'.format(n)
        matchstring+='Potential local matches: {}</br>'.format(l)
        if l>0:
            # Make map centered on seismicity
            stamen_terrain = cimgt.StamenTerrain()
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection=stamen_terrain.crs)
            ax.set_extent([np.median(llons)-opt.locdeg/2,np.median(llons)+opt.locdeg/2,
                np.median(llats)-opt.locdeg/4,np.median(llats)+opt.locdeg/4],
                crs=ccrs.Geodetic())
            # Shaded terrain
            ax.add_image(stamen_terrain, 11, alpha=0.75)
            
            # Set up ticks
            ax.set_xticks(np.arange(np.floor(10*(np.median(llons)-opt.locdeg/2))/10,
                np.ceil(10*(np.median(llons)+opt.locdeg/2))/10,0.1), 
                crs=ccrs.PlateCarree())
            ax.set_yticks(np.arange(np.floor(10*(np.median(llats)-opt.locdeg/4))/10,
                np.ceil(10*(np.median(llats)+opt.locdeg/4))/10, 0.1),
                crs=ccrs.PlateCarree())
            ax.set_extent([np.median(llons)-opt.locdeg/2,np.median(llons)+opt.locdeg/2,
                np.median(llats)-opt.locdeg/4,np.median(llats)+opt.locdeg/4],
                crs=ccrs.Geodetic())
            ax.xaxis.set_major_formatter(LongitudeFormatter())
            ax.yaxis.set_major_formatter(LatitudeFormatter())
            plt.yticks(rotation=90, va='center')
            
            # Seismicity in red (halo of white), stations open black triangles
            ax.scatter(llons, llats, s=20, marker='o', color='white',
                transform=ccrs.Geodetic())
            ax.scatter(llons, llats, s=5, marker='o', color='red',
                transform=ccrs.Geodetic())
            ax.scatter(stalons, stalats, marker='^', color='k', facecolors='None',
                transform=ccrs.Geodetic())
            
            # 10 km scale bar
            sbllon = 0.05*(opt.locdeg)+np.median(llons)-opt.locdeg/2
            sbllat = 0.05*(opt.locdeg/2)+np.median(llats)-opt.locdeg/4
            sbelon = sbllon + np.arctan2(np.sin(np.pi/2)*np.sin(
                10./6378.)*np.cos(sbllat*np.pi/180.), np.cos(10./6378.)-np.sin(
                sbllat*np.pi/180.)*np.sin(sbllat*np.pi/180.))*180./np.pi
            ax.plot((sbllon, sbelon), (sbllat,sbllat), 'k-', transform=ccrs.Geodetic(),
                lw=2)
            geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
            text_transform = offset_copy(geodetic_transform, units='dots', y=5)
            ax.text((sbllon+sbelon)/2., sbllat, '10 km', ha='center',
                transform=text_transform)
                       
            plt.title('{} potential local matches (~{:3.1f} km depth)'.format(l,
                np.mean(ldeps)))
            plt.tight_layout()
            plt.savefig('./{}/clusters/map{}.png'.format(opt.groupName,cnum),
                dpi=100)
            plt.close()
            f.write('<img src="map{}.png"></br>'.format(cnum))            
    else:
        matchstring+='No matches found</br></div>'  
    f.write(matchstring)


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

                    
def printOrphanCatalog(otable, opt):
    """
    Prints flat catalog of current orphans to text file
    
    otable: Orphans table
    opt: Options object describing station/run parameters
    
    Note: Time in text file corresponds to original STA/LTA trigger time
    """

    with open('{}/orphancatalog.txt'.format(opt.groupName), 'w') as f:
        
        startTimes = otable.cols.startTime[:]
        
        for i in np.argsort(startTimes):
            f.write("{0}\n".format((UTCDateTime(startTimes[i])+opt.ptrig).isoformat()))


def createJunkPlots(jtable, opt):
    """
    Description goes here!
    
    """
    
    # Write out times of junk triggers
    printJunk(jtable, opt)
    
    for r in jtable:
        fig = plt.figure(figsize=(15, 0.5))
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        
        
        for s in range(opt.nsta):
            waveformc = r['waveform'][s*opt.wshape:(s+1)*opt.wshape]
            tmpc = waveformc[max(0, r['windowStart']-int(
                opt.ptrig*opt.samprate)):min(opt.wshape,
                r['windowStart']+int(opt.atrig*opt.samprate))]
            datc = tmpc[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
                opt.ptrig*opt.samprate + opt.winlen*1.5)]
            datc = datc/np.max(np.abs(datc)+1.0/1000)
            datc[datc>1] = 1
            datc[datc<-1] = -1
            if s == 0:
                dat = datc
            else:
                dat = np.append(dat,datc)
        
        ax.plot(dat,'k',linewidth=0.25)
        plt.autoscale(tight=True)
        plt.savefig('{0}/junk/{1}-{2}.png'.format(opt.groupName,
            UTCDateTime(r['startTime'])+opt.ptrig,r['isjunk']), dpi=100)
        plt.close(fig)

            
def printJunk(jtable, opt):
    """
    Prints flat catalog of contents of junk table to text file for debugging
    
    jtable: Junk table
    opt: Options object describing station/run parameters
    
    Note: Time in text file corresponds to original STA/LTA trigger time
    """

    with open('{}/junk.txt'.format(opt.groupName), 'w') as f:
        
        startTimes = jtable.cols.startTime[:]
        jtype = jtable.cols.isjunk[:]
        
        for i in np.argsort(startTimes):
            f.write("{0} - {1}\n".format((
                UTCDateTime(startTimes[i])+opt.ptrig).isoformat(),jtype[i]))


def printCoresCatalog(rtable, ftable, opt):
    """
    Prints flat catalog of only core events to text file
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Note: Time in text file corresponds to current trigger time by alignment
    """

    with open('{}/cores.txt'.format(opt.groupName), 'w') as f:
        
        startTimes = rtable.cols.startTime[:]
        windowStarts = rtable.cols.windowStart[:]
        
        for cnum in range(ftable.attrs.nClust):
            fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
            core = ftable[cnum]['core']
            f.write("{0} {1}\n".format(cnum,(UTCDateTime(startTimes[core]) +
                windowStarts[core]/opt.samprate).isoformat()))

                
def printEventsperDay(rtable, ftable, opt):
    """
    Prints daily counts of each family in a tablulated text file
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    Each column (with the exception of first and last) correspond to individual families;
    first column is date and last column is total across all families.
    """
    
    with open('{}/dailycounts.txt'.format(opt.groupName), 'w') as f:
        
        startTimes = rtable.cols.startTimeMPL[:]
        firstDay = np.floor(np.min(startTimes)).astype(int)
        lastDay = np.ceil(np.max(startTimes)).astype(int)
        hists = np.zeros((ftable.attrs.nClust,lastDay-firstDay))
        
        # Calculate histograms
        for cnum in range(ftable.attrs.nClust):
            fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
            hists[cnum,:], edges = np.histogram(startTimes[fam], bins=np.arange(
                firstDay,lastDay+1,1))
        
        # Header
        f.write("      Date\t")
        for cnum in range(ftable.attrs.nClust):
            f.write("{}\t".format(cnum))
        f.write("Total\n")
        
        # Write daily counts
        for day in range(firstDay,lastDay):
            f.write("{}\t".format(matplotlib.dates.num2date(day).strftime('%Y/%m/%d')))
            for cnum in range(ftable.attrs.nClust):
                f.write("{}\t".format(hists[cnum,day-firstDay].astype(int)))
            f.write("{}\n".format(np.sum(hists[:,day-firstDay].astype(int))))

            
def printVerboseCatalog(rtable, ftable, ctable, opt):
    """
    Prints flat catalog to text file with additional columns
    
    rtable: Repeater table
    ftable: Families table
    ctable: Correlation table
    opt: Options object describing station/run parameters
    
    Columns correspond to cluster number, event time, frequency index, amplitude, time
    since last event in hours, and correlation coefficient with respect to the best
    correlated event.
    """
    
    with open('{}/catalog.txt'.format(opt.groupName), 'w') as f:
                
        startTimes = rtable.cols.startTime[:]
        startTimeMPL = rtable.cols.startTimeMPL[:]
        windowStarts = rtable.cols.windowStart[:]
        windowAmps = rtable.cols.windowAmp[:][:,opt.printsta]
        ids = rtable.cols.id[:]
        id1 = ctable.cols.id1[:]
        id2 = ctable.cols.id2[:]
        ccc = ctable.cols.ccc[:]
        fi = np.nanmean(rtable.cols.FI[:], axis=1)
        
        f.write("cnum\tevTime                        \tfi\tamp\tdt\t\txcorr\n")
        for cnum in range(ftable.attrs.nClust):
            fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
            
            catalogind = np.argsort(startTimeMPL[fam])
            catalog = startTimeMPL[fam][catalogind]
            spacing = np.diff(catalog)*24
            
            idf = ids[fam]
            ix = np.where(np.in1d(id2,idf))
            r = np.zeros((max(idf)+1,)).astype('int')
            r[idf] = range(len(idf))
            C = np.zeros((len(idf),len(idf)))
            C[r[id2[ix]],r[id1[ix]]] = ccc[ix]
            C[r[id1[ix]],r[id2[ix]]] = ccc[ix]
            C[range(len(idf)),range(len(idf))] = 1.0
            xcorr = C[np.argmax(np.sum(C,0)),:]
            
            j = -1
            for i in catalogind:
                evTime = (UTCDateTime(startTimes[fam][i]) +
                    windowStarts[fam][i]/opt.samprate)
                amp = windowAmps[fam][i]
                if j == -1:
                    dt = 'NaN         '
                else:
                    dt = spacing[j]
                j += 1
                
                f.write("{0}\t{1}\t{2:4.3f}\t{3:5.2f}\t{4}\t\t{5:2.1f}\n".format(cnum,evTime.isoformat(),
                    fi[fam][i],amp,dt,xcorr[i]))


def plotReport(rtable, ftable, ctable, opt, fnum, ordered):
    
    """
    Creates more detailed output plots for a single family
    
    rtable: Repeater table
    ftable: Families table
    ctable: Correlation table
    opt: Options object describing station/run parameters
    fnum: Family to be inspected
        
    """
    
    # Adjust the font face
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['font.size'] = 8.0    
    
    # Set up variables
    fam = np.fromstring(ftable[fnum]['members'], dtype=int, sep=' ')
    startTimeMPL = rtable.cols.startTimeMPL[:]
    startTime = rtable.cols.startTime[:]
    windowStart = rtable.cols.windowStart[:]
    windowAmp = rtable.cols.windowAmp[:][:,opt.printsta]
    windowAmps = rtable.cols.windowAmp[:]
    fi = rtable.cols.FI[:]
    ids = rtable.cols.id[:]
    id1 = ctable.cols.id1[:]
    id2 = ctable.cols.id2[:]
    ccc = ctable.cols.ccc[:]
    core = ftable[fnum]['core']
    catalogind = np.argsort(startTimeMPL[fam])
    catalog = startTimeMPL[fam][catalogind]
    famcat = fam[catalogind]
    longevity = ftable[fnum]['longevity']
    spacing = np.diff(catalog)*24
    minind = fam[catalogind[0]]
    maxind = fam[catalogind[-1]]

    idf = ids[fam]
    ix = np.where(np.in1d(id2,idf))
    r = np.zeros((max(idf)+1,)).astype('int')
    r[idf] = range(len(idf))
    C = np.zeros((len(idf),len(idf)))
    C[r[id2[ix]],r[id1[ix]]] = ccc[ix]
    C[r[id1[ix]],r[id2[ix]]] = ccc[ix]
    C[range(len(idf)),range(len(idf))] = 1.0
    
    # Copy static preview image in case cluster changes
    shutil.copy('{0}/clusters/{1}.png'.format(opt.groupName, fnum),
                '{0}/clusters/{1}-report.png'.format(opt.groupName, fnum))
    
    # Fill in full correlation matrix
    print('Computing full correlation matrix; this will take time if the family is large!')
    famtable = rtable[famcat]
    Cind = C[catalogind,:]
    Cind = Cind[:,catalogind]
    Cfull = Cind.copy()    
    for i in range(len(famcat)-1):
        for j in range(i+1,len(famcat)):
            if Cfull[i,j]==0:
                # Compute correlation
                cor, lag, nthcor = redpy.correlation.xcorr1x1(famtable['windowFFT'][i],
                    famtable['windowFFT'][j], famtable['windowCoeff'][i],
                    famtable['windowCoeff'][j], opt)                
                Cfull[i,j] = cor
                Cfull[j,i] = cor 
    
    ### BOKEH PLOTS    
    oTOOLS = ['pan,box_zoom,reset,save,tap']    
    
    # Amplitude vs. time
    o0 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title='Amplitude on {} with Time'.format(opt.station.split(',')[opt.printsta]),
        y_axis_type='log', y_range=[1, 2*np.max(windowAmp)])
    o0.grid.grid_line_alpha = 0.3
    o0.xaxis.axis_label = 'Date'
    o0.yaxis.axis_label = 'Counts'
    o0.circle(matplotlib.dates.num2date(startTimeMPL[fam]), windowAmp[fam], color='red',
        line_alpha=0, size=4, fill_alpha=0.5)
    
    # Time since last event
    o1 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title='Time since Previous Event', x_range=o0.x_range, y_axis_type='log',
        y_range=[1e-3, 2*np.max(spacing)])
    o1.grid.grid_line_alpha = 0.3
    o1.xaxis.axis_label = 'Date'
    o1.yaxis.axis_label = 'Interval (hr)'
    o1.circle(matplotlib.dates.num2date(catalog[1:]), spacing, color='red',
        line_alpha=0, size=4, fill_alpha=0.5)
    
    # Cross-correlation wrt. core
    o2 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title='Cross-correlation Coefficient with Core Event', x_range=o0.x_range,
        y_range=[0, 1.02])
    o2.grid.grid_line_alpha = 0.3
    o2.xaxis.axis_label = 'Date'
    o2.yaxis.axis_label = 'CCC'   
    o2.circle(matplotlib.dates.num2date(catalog), Cfull[np.where(famcat==core)[0],:][0],
        color='red', line_alpha=0, size=4, fill_alpha=0.5)
    
    # Combine and save
    o = gridplot([[o0],[o1],[o2]])
    output_file('{0}/clusters/{1}-report-bokeh.html'.format(opt.groupName,fnum),
        title='{0} - Cluster {1} Detailed Report'.format(opt.title,fnum))
    save(o)
    
    ### OPTICS ORDERING (OPTIONAL)
    if ordered:
        # Order by OPTICS rather than by time
        D = 1-Cfull
        s = np.argsort(sum(D))[::-1]
        D = D[s,:]
        D = D[:,s]
        famcat = famcat[s]
        Cind = Cind[s,:]
        Cind = Cind[:,s]
        Cfull = Cfull[s,:]
        Cfull = Cfull[:,s]
        ttree = setOfObjects(D)
        prep_optics(ttree,1)
        build_optics(ttree,1)
        order = np.array(ttree._ordered_list)
        famcat = famcat[order]
        Cind = Cind[order,:]
        Cind = Cind[:,order]
        Cfull = Cfull[order,:]
        Cfull = Cfull[:,order]
    
    ### CORRELATION MATRIX
    fig = plt.figure(figsize=(14,5.4))
    ax1 = fig.add_subplot(1,2,1)
    cax = ax1.imshow(Cind, vmin=opt.cmin-0.05, cmap='Spectral_r')
    cbar = plt.colorbar(cax, ticks=np.arange(opt.cmin-0.05,1.05,0.05))
    tix = cbar.ax.get_yticklabels()
    tix[0] = 'Undefined'
    cbar.ax.set_yticklabels(tix)
    if ordered:
        plt.title('Stored Correlation Matrix (Ordered)', fontweight='bold')
    else:
        plt.title('Stored Correlation Matrix', fontweight='bold')
    ax2 = fig.add_subplot(1,2,2)
    cax2 = ax2.imshow(Cfull, vmin=opt.cmin-0.05, cmap='Spectral_r')
    cbar2 = plt.colorbar(cax2, ticks=np.arange(opt.cmin-0.05,1.05,0.05))
    tix = cbar2.ax.get_yticklabels()
    tix[0] = '< {:1.2f}'.format(opt.cmin-0.05)
    cbar2.ax.set_yticklabels(tix)
    if ordered:
        plt.title('Full Correlation Matrix (Ordered)', fontweight='bold')
    else:
        plt.title('Full Correlation Matrix', fontweight='bold')
    plt.tight_layout()
    plt.savefig('{0}/clusters/{1}-reportcmat.png'.format(opt.groupName,fnum), dpi=100)
    plt.close(fig)
    
    ### WAVEFORM IMAGES
    famtable = rtable[famcat]
    fig2 = plt.figure(figsize=(10, 12))
    
    for sta in range(opt.nsta):
        n = -1
        data = np.zeros((len(fam), int(opt.winlen*2)))
        ax = fig2.add_subplot(np.ceil((opt.nsta)/2), 2, sta+1)
        for r in famtable:
            if ordered:
                plt.title('{0}.{1} (Ordered)'.format(opt.station.split(',')[sta],
                          opt.channel.split(',')[sta]), fontweight='bold')
            else:
                plt.title('{0}.{1}'.format(opt.station.split(',')[sta],
                          opt.channel.split(',')[sta]), fontweight='bold')
            n = n+1        
            waveform = r['waveform'][sta*opt.wshape:(sta+1)*opt.wshape]
            tmp = waveform[max(0, windowStart[famcat[n]]-int(
                opt.ptrig*opt.samprate)):min(opt.wshape,
                windowStart[famcat[n]]+int(opt.atrig*opt.samprate))]
            data[n, :] = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
                opt.ptrig*opt.samprate + opt.winlen*1.5)]/windowAmps[famcat[n]][sta]            
        if len(fam) > 12:
            ax.imshow(data, aspect='auto', vmin=-1, vmax=1, cmap='RdBu',
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
                ax.plot(tvec,dat/2-o*0.75,'k',linewidth=0.5)
            plt.xlim([np.min(tvec),np.max(tvec)])
            plt.ylim([-o*0.75-0.5,0.5])
        ax.yaxis.set_visible(False)
        plt.xlabel('Time Relative to Trigger (seconds)', style='italic')
    plt.tight_layout()
    plt.savefig('{0}/clusters/{1}-reportwaves.png'.format(opt.groupName,
                fnum), dpi=100)
    plt.close(fig2)
    
    ### HTML OUTPUT PAGE
    tstamp = UTCDateTime.now()
    with open('{0}/clusters/{1}-report.html'.format(opt.groupName, fnum), 'w') as f:
        f.write("""
        <html><head><title>{1} - Cluster {0} Detailed Report</title>
        </head><style>
        a {{color:red;}}
        body {{font-family:Helvetica; font-size:12px}}
        h1 {{font-size: 20px;}}
        </style>
        <body><center>
        <em>Last updated: {10}</em></br>
        <h1>Cluster {0} - Detailed Report</h1>                
        <img src="{0}-report.png" width=500 height=100></br></br>
            Number of events: {2}</br>
            Longevity: {5:.2f} days</br>
            Mean event spacing: {7:.2f} hours</br>
            Median event spacing: {8:.2f} hours</br>
            Mean Frequency Index: {9:.2f}<br></br>
            First event: {3}</br>
            Core event: {6}</br>
            Last event: {4}</br>
            
            <img src='{11}-reportwaves.png'></br></br>
            
            <iframe src="{11}-report-bokeh.html" width=1350 height=800
            style="border:none"></iframe>
            
            </br>
            <img src='{11}-reportcmat.png'></br></br></br>
            
        """.format(fnum, opt.title, len(fam), (UTCDateTime(
            startTime[minind]) + windowStart[minind]/opt.samprate).isoformat(),
            (UTCDateTime(startTime[maxind]) + windowStart[
            maxind]/opt.samprate).isoformat(), longevity, (UTCDateTime(
            startTime[core]) + windowStart[core]/opt.samprate).isoformat(),
            np.mean(spacing), np.median(spacing), np.mean(np.nanmean(fi[fam],
            axis=1)),tstamp,fnum))
            
        f.write("""
        </center></body></html>
        """)
