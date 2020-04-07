# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2018  Alicia Hotovec-Ellis (ahotovec@gmail.com)
# Licensed under GNU GPLv3 (see LICENSE.txt)

from tables import *
import numpy as np
import matplotlib
import datetime
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates
import time
import pandas as pd
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
from bokeh.models import HoverTool, ColumnDataSource, OpenURL, TapTool, Range1d, Div, Span
from bokeh.models import Arrow, VeeHead, ColorBar, LogColorMapper, LogTicker, LabelSet
from bokeh.models.glyphs import Line, Quad
from bokeh.layouts import column
from bokeh.palettes import inferno, all_palettes
import urllib.request

        
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
    
    printTriggerCatalog(ttable, opt)
    printOrphanCatalog(otable, opt)
    if len(rtable)>1:
        plotTimelines(rtable, ftable, ttable, opt)
        if np.sum(ftable.cols.printme[:]):
            if opt.printVerboseCat == True:
                printVerboseCatalog(rtable, ftable, ctable, opt)
            else:
                printCatalog(rtable, ftable, opt)
            printSwarmCatalog(rtable, ftable, ttable, opt)
            printCoresCatalog(rtable, ftable, opt)
            printEventsperDay(rtable, ftable, opt)
            plotCores(rtable, ftable, opt)
            plotFamilies(rtable, ftable, ctable, opt)
            ftable.cols.printme[:] = np.zeros((len(ftable),))
            ftable.cols.lastprint[:] = np.arange(len(ftable))
    else:
        print('Nothing to plot!')
    
    # Rename any .tmp files
    tmplist = glob.glob('{}{}/clusters/*.tmp'.format(opt.outputPath, opt.groupName))
    for tmp in tmplist:
        os.rename(tmp,tmp[0:-4]) 
    
  
def plotTimelines(rtable, ftable, ttable, opt):
    
    """
    Creates the primary .html Bokeh timelines
    
    rtable: Repeater table
    ftable: Families table
    ttable: Triggers table
    opt: Options object describing station/run parameters
    
    """
    
    dt = rtable.cols.startTimeMPL[:]
    fi = np.nanmean(rtable.cols.FI[:], axis=1)
    longevity = ftable.cols.longevity[:]
    famstarts = ftable.cols.startTime[:]
    alltrigs = ttable.cols.startTimeMPL[:]
    
    # Read in annotation file (if it exists)
    if opt.anotfile != '':
        df = pd.read_csv(opt.anotfile)
    
          
    # Determine padding for hover bars (~1% of window range on each side)
    barpad = (max(alltrigs)-min(alltrigs))*0.01
    barpadr = opt.recplot*0.01
        
    
    # Initialize list of produced plots (plot_types is a stub variable until the config file offers this option)    
    plot_types = 'eqrate,fi,occurrence,longevity'.split(',')
    overview_plots = []
    recent_plots = []
    
    
    # Create each of the subplots specified in the configuration file
    for p in plot_types:
    
        if p == 'eqrate':
            # Plot EQ Rates (Repeaters and Orphans)
            overview_plots.append(plot_eqrate(alltrigs, dt, opt.dybin))
            recent_plots.append(  plot_eqrate(alltrigs[dt>(max(alltrigs)-opt.recplot)], dt[dt>(max(alltrigs)-opt.recplot)], opt.hrbin/24))
            
        elif p == 'fi':
            # Plot Frequency Index
            overview_plots.append(plot_fi(dt, fi)) # Overview
            recent_plots.append(plot_fi(dt[dt>(max(alltrigs)-opt.recplot)], fi[dt>(max(alltrigs)-opt.recplot)])) # Recent
    
        elif p == 'longevity':
	        # Plot Cluster Longevity — This needs to be further functionalized!
            overview_plots.append(plot_longevity(alltrigs, famstarts, longevity, min(alltrigs),             barpad, opt))
            recent_plots.append(  plot_longevity(alltrigs, famstarts, longevity, max(alltrigs)-opt.recplot, barpadr, opt))
	
        elif p == 'occurrence':
            # Plot family occurrence
            overview_plots.append(plot_family_occurrence(dt, ftable, min(alltrigs),             opt.minplot, opt.occurbin, barpad))
            recent_plots.append(  plot_family_occurrence(dt, ftable, max(alltrigs)-opt.recplot, 0,           opt.recbin,   barpadr))
    
        else:
        	print('{} is not a valid plot type. Moving on.'.format(p))
    
    
    # Set ranges
    for i in overview_plots:
        i.x_range = overview_plots[0].x_range
    for i in recent_plots:
        i.x_range = recent_plots[0].x_range
    
    
    # Add annotations
    if opt.anotfile != '':
        for i in [overview_plots, recent_plots]:
            for j in i:
                for row in df.itertuples():
                    spantime = (datetime.datetime.strptime(row[1]
                        ,'%Y-%m-%dT%H:%M:%S')-datetime.datetime(1970, 1, 1)).total_seconds()
                    j.add_layout(Span(location=spantime*1000, dimension='height',
                        line_color=row[2], line_width=row[3], line_dash=row[4],
                        line_alpha=row[5], level='underlay'))
 

    # Create output and save
    # gridplot_items should look like this: [[Div(text=...)],[panel1],[panel2],...]
    gridplot_items = [[Div(text='<h1>{0}</h1>'.format(opt.title), width=1000)]] + [[el] for el in overview_plots]
    o = gridplot(gridplot_items)
    output_file('{}{}/overview.html'.format(opt.outputPath, opt.groupName), title='{} Overview'.format(opt.title))
    save(o)
    
    gridplot_items = [[Div(text='<h1>{0}</h1>'.format(opt.title), width=1000)]] + [[el] for el in recent_plots]
    r = gridplot(gridplot_items)
    output_file('{}{}/overview_recent.html'.format(opt.outputPath, opt.groupName), title='{0} Overview - Last {1:.1f} Days'.format(opt.title, opt.recplot))
    save(r)

         
def bokeh_figure(**kwargs):
    
    # default values for bokeh_figures
    if 'tools' not in kwargs:
        kwargs['tools'] = ['pan,box_zoom,reset,save,tap']
    if 'plot_width' not in kwargs:
        kwargs['plot_width'] = 1250
    if 'plot_height' not in kwargs:
        kwargs['plot_height'] = 250
    if 'x_axis_type' not in kwargs:
        kwargs['x_axis_type'] = 'datetime'
        
    # Create figure
    fig = figure(**kwargs)
        
    fig.grid.grid_line_alpha = 0.3
    fig.xaxis.axis_label = 'Date'
    fig.yaxis.axis_label = ''
    
    return fig


def plot_eqrate(alltrigs, dt, binsize):

    hr_days = 'Day Bin' if binsize>=1 else 'Hour Bin'
    dt_offset = binsize/2 # used to create the lines
    title = 'Repeaters vs. Orphans by {:.1f} {}'.format(binsize,hr_days)
     
    # Create histogram of events/dybin
    histT, hT = np.histogram(alltrigs, bins=np.arange(min(alltrigs),
        max(alltrigs+binsize), binsize))
    histR, hR = np.histogram(dt, bins=np.arange(min(alltrigs),
        max(alltrigs+binsize), binsize))  
            
    # Plot data
    fig = bokeh_figure(title=title)
    fig.yaxis.axis_label = 'Events'
    fig.line(matplotlib.dates.num2date(hT[0:-1]+dt_offset), histT-histR, color='black',
        legend_label='Orphans')
    fig.line(matplotlib.dates.num2date(hR[0:-1]+dt_offset), histR, color='red',
        legend_label='Repeaters', line_width=2)
    fig.legend.location = 'top_left'
    
    return fig


def plot_fi(dt, fi):
	
    fig = bokeh_figure(title='Frequency Index')
    fig.yaxis.axis_label = 'FI'
    fig.circle(matplotlib.dates.num2date(dt), fi, color='red', line_alpha=0,
        size=3, fill_alpha=0.5)
    
    return fig


def plot_longevity(alltrigs, famstarts, longevity, t1, barpad, opt):

    fig = bokeh_figure(y_axis_type='log',
        y_range=[0.1, np.sort(alltrigs)[-1]-np.sort(alltrigs)[0]],
        title='Cluster Longevity')
    fig.yaxis.axis_label = 'Days'

    add_line = True    # Draw a line for the longevity data (turns off if data don't fall within time window)
    add_arrow = True   # Draw an arrow if longevity line extends beyond the data window

    # Plot Data            
    for n in range(len(famstarts)):
        if (max(alltrigs)-opt.recplot)<=famstarts[n]: # family start is after recplot
            
            x1 = famstarts[n]+longevity[n]
            add_line = True
            add_arrow = False

        elif (max(alltrigs)-opt.recplot)<=famstarts[n]+longevity[n]: # family end is after recplot

            x1 = t1-barpad
            add_line = True
            add_arrow = True

        else: # family start and end are before recplot (ie, entire family is before plot window)
            print('I did get here.')
            add_line = False
            add_arrow = False
            #x1=0# this shouldn't (couldn't?) happen
    
        if add_line:
            source = ColumnDataSource(dict(
                x=np.array(
                (matplotlib.dates.num2date(x1),
                matplotlib.dates.num2date(famstarts[n]+longevity[n]))),
                y=np.array((longevity[n],longevity[n]))))
            fig.add_glyph(source, Line(x="x", y="y", line_color='red',
                line_alpha=0.5))
            if add_arrow:
                fig.add_layout(Arrow(end=VeeHead(size=5, fill_color='red', line_color='red',
                    line_alpha=0.5), line_alpha=0,
                    x_start=matplotlib.dates.num2date(famstarts[n]+longevity[n]),
                    x_end=matplotlib.dates.num2date(t1-barpad),
                    y_start=longevity[n], y_end=longevity[n]))

    return fig


def plot_family_occurrence(dt, ftable, t1, minplot, binsize, barpad):

    t0 = t1
    #print(t1)
    #print(t0)
   
    p1 = figure(tools=[create_hover_tool(),'pan,box_zoom,reset,save,tap'],
        title='Occurrence Timeline', plot_height=500, plot_width=1250)
    p1.yaxis.axis_label = 'Cluster by Date' + (' ({}+ Members)'.format(minplot) if minplot>0 else '')
    
    # Steal YlOrRd (len=256) colormap from matplotlib
    colormap = matplotlib.cm.get_cmap('YlOrRd')
    bokehpalette = [matplotlib.colors.rgb2hex(m) for m in colormap(
        np.arange(colormap.N)[::-1])]

    # Build the lists and dictionaries    
    n = 0  
    cloc1 = 335
    
    legtext = determine_legend_text(binsize)
      
    for clustNum in range(ftable.attrs.nClust):
        
        members = np.fromstring(ftable[clustNum]['members'], dtype=int, sep=' ')
        
        # Create histogram of events/hour
        hist, h = np.histogram(dt[members], bins=np.arange(min(dt[members]),
            max(dt[members]+binsize), binsize))
        d1 = matplotlib.dates.num2date(h[np.where(hist>0)])
        d2 = matplotlib.dates.num2date(h[np.where(hist>0)]+1.0/24)
        histlog = np.log10(hist[hist>0])
        ind = [int(min(255,255*(i/2))) for i in histlog]
        colors = [bokehpalette[i] for i in ind]
                    

        if len(dt[members]) >= minplot:

            if max(dt[members])>t1:
            
                if min(dt[members])<t1:
                    
                    # add line w arrow
                    p1.line((matplotlib.dates.num2date(t0),
                        matplotlib.dates.num2date(max(dt[members]))), (n, n),
                        color='black')
                    p1.add_layout(Arrow(end=VeeHead(size=3),
                        x_start=matplotlib.dates.num2date(t0+0.01),
                        x_end=matplotlib.dates.num2date(t0-0.5),
                        y_start=n, y_end=n))

                    idx = np.where(h[np.where(hist>0)[0]]>t1)[0]
                        
                else:
                
                    # Just add line
                    # Date is required as datenum
                    p1.line((matplotlib.dates.num2date(min(dt[members])),
                        matplotlib.dates.num2date(max(dt[members]))), (n, n),
                        color='black')
                    idx = np.arange(len(d1))
                
                # always add a box
                p1.quad(top=n+0.3, bottom=n-0.3,
                    left=np.array(d1)[idx],
                    right=np.array(d2)[idx],
                    color=np.array(colors)[idx])                   
                
                # Text doesn't understand datetimes, need to convert to a number and subtract
                # about 8 hours
                p1.text(time.mktime(matplotlib.dates.num2date(
                    max(dt[members])).timetuple())*1000 - 28799000, n,
                    text=['   {}'.format(len(dt[members]))], text_font_size='9pt',
                    text_baseline='middle')
                 
                # Build source for hover patches
                fnum = clustNum
                if n == 0:
                    xs=[[matplotlib.dates.num2date(max(min(dt[members]),t1)-barpad),
                        matplotlib.dates.num2date(max(min(dt[members]),t1)-barpad),
                        matplotlib.dates.num2date(max(dt[members])+barpad),
                        matplotlib.dates.num2date(max(dt[members])+barpad)]]
                    ys=[[n-0.5, n+0.5, n+0.5, n-0.5]]
                    famnum=[[fnum]]
                else:
                    xs.append([matplotlib.dates.num2date(max(min(dt[members]),t1)-barpad),
                        matplotlib.dates.num2date(max(min(dt[members]),t1)-barpad),
                        matplotlib.dates.num2date(max(dt[members])+barpad),
                        matplotlib.dates.num2date(max(dt[members])+barpad)])
                    ys.append([n-0.5, n+0.5, n+0.5, n-0.5])
                    famnum.append([fnum])

                n = n+1
         
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
        p1.circle(matplotlib.dates.num2date(t0), 0, line_alpha=0, fill_alpha=0)


    color_bar = ColorBar(color_mapper=determine_color_mapper(binsize), ticker=LogTicker(),
        border_line_color='#eeeeee', location=(7,cloc1), orientation='horizontal',
        width=150, height=15, title='Events per {}'.format(determine_legend_text(binsize)), padding=15,
        major_tick_line_alpha=0)    
#    color_mapper = LogColorMapper(palette=bokehpalette, low=1, high=100)
#    color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
#        border_line_color='#eeeeee', location=(7,cloc1), orientation='horizontal',
#        width=100, height=15, title='Events per Hour', padding=15,
#        major_tick_line_alpha=0)
    p1.add_layout(color_bar)

    return p1


def determine_legend_text(binsize):
    # Determine legend wording
    if binsize == 1/24:
        legtext = 'Hour'
    elif binsize == 1:
        legtext = 'Day'
    elif binsize == 7:
        legtext = 'Week'
    elif binsize < 2:
        legtext = '{} Hours'.format(binsize*24)
    else:
        legtext = '{} Days'.format(binsize)

    return legtext


def determine_color_mapper(binsize):
    # Steal YlOrRd (len=256) colormap from matplotlib
    colormap = matplotlib.cm.get_cmap('YlOrRd')
    bokehpalette = [matplotlib.colors.rgb2hex(m) for m in colormap(
        np.arange(colormap.N)[::-1])]
    if binsize >= 1:
        color_mapper = LogColorMapper(palette=bokehpalette, low=1, high=1000)
    else:
        color_mapper = LogColorMapper(palette=bokehpalette, low=1, high=100)


def create_hover_tool():
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

    return hover


def plotCores(rtable, ftable, opt):

    """
    Plots core waveforms as .png for hovering in timeline
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    """
    
    for n in range(len(ftable))[::-1]:
        if ftable.cols.lastprint[n] != n and ftable.cols.printme[n] == 0:
            os.rename('{}{}/clusters/{}.png'.format(opt.outputPath, opt.groupName,
                ftable.cols.lastprint[n]), '{}{}/clusters/{}.png.tmp'.format(
                opt.outputPath, opt.groupName, n))
            os.rename('{}{}/clusters/fam{}.png'.format(opt.outputPath, opt.groupName,
                ftable.cols.lastprint[n]), '{}{}/clusters/fam{}.png.tmp'.format(
                opt.outputPath, opt.groupName, n))
    
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
            plt.savefig('{}{}/clusters/{}.png'.format(opt.outputPath, opt.groupName, n),
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
                    try:
                        data[n, :] = tmp[int(opt.ptrig*opt.samprate - opt.winlen*0.5):int(
                            opt.ptrig*opt.samprate + opt.winlen*1.5)]/windowAmp[fam[n]]
                    except (ValueError, Exception):
                        print('Error in printing family {}, moving on...'.format(cnum))
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
            
            # Set min/max for plotting
            if opt.amplims == 'family':
                windowAmpFam = windowAmp[fam[catalogind]]
                try:
                    ymin = 0.5*np.min(windowAmpFam[np.nonzero(windowAmpFam)])
                    ymax = 2*np.max(windowAmpFam)
                except ValueError:
                    # Use global if all zeros
                    ymin = 0.5*np.min(windowAmp[np.nonzero(windowAmp)])
                    ymax = 2*np.max(windowAmp)
            else:
                # Use global maximum/minimum
                ymin = 0.5*np.min(windowAmp[np.nonzero(windowAmp)])
                ymax = 2*np.max(windowAmp)
            
            # Plot amplitude timeline
            ax3 = fig.add_subplot(9, 3, (10,15))
            ax3.plot_date(catalog, windowAmp[fam[catalogind]],
                    'ro', alpha=0.5, markeredgecolor='r', markeredgewidth=0.5,
                    markersize=3)
            myFmt = matplotlib.dates.DateFormatter('%Y-%m-%d\n%H:%M')
            ax3.xaxis.set_major_formatter(myFmt)
            ax3.set_ylim(ymin, ymax)
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
            plt.savefig('{}{}/clusters/fam{}.png'.format(opt.outputPath, opt.groupName,
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
            with open('{}{}/clusters/{}.html'.format(opt.outputPath, opt.groupName,
                     cnum), 'w') as f:
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

    if opt.matchMax > 0:
        windowAmp = rtable.cols.windowAmp[:][:,opt.printsta]

    members = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
    if opt.matchMax == 0 or opt.matchMax > len(members):
        order = np.argsort(startTime[members])
        matchstring = ('</br><b>ComCat matches (all events):</b></br>'
            '<div style="overflow-y: auto; height:100px; width:1200px;">')
    else:
        nlargest = np.argsort(windowAmp[members])[::-1][:opt.matchMax]
        members = members[nlargest]
        order = np.argsort(startTime[members])
        matchstring = ('</br><b>ComCat matches ({} largest events):</b></br>'
            '<div style="overflow-y: auto; height:100px; width:1200px;">').format(
            opt.matchMax)
    
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
        except (urllib.error.HTTPError, urllib.error.URLError):
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
            except (ValueError, urllib.error.HTTPError, urllib.error.URLError):
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
            plt.savefig('{}{}/clusters/map{}.png'.format(opt.outputPath, opt.groupName,
                cnum), dpi=100)
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

    with open('{}{}/catalog.txt'.format(opt.outputPath, opt.groupName), 'w') as f:
        
        startTimes = rtable.cols.startTime[:]
        windowStarts = rtable.cols.windowStart[:]
        
        for cnum in range(ftable.attrs.nClust):
            fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
            for i in np.argsort(startTimes[fam]):
                f.write("{0} {1}\n".format(cnum,(UTCDateTime(startTimes[fam][i]) +
                    windowStarts[fam][i]/opt.samprate).isoformat()))
    
    """
    startTimes = rtable.cols.startTime[:]
    windowStarts = rtable.cols.windowStart[:]
    famNums = np.ones(len(windowStarts), dtype=int)
    for cnum in range(ftable.attrs.nClust):
        fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
        famNums[fam] = cnum
    print(famNums)   
    """ 


def printTriggerCatalog(ttable, opt):
    """
    Prints flat catalog of all triggers to text file
    
    ttable: Triggers table
    opt: Options object describing station/run parameters
    
    Note: Time in text file corresponds to original STA/LTA trigger time
    """

    with open('{}{}/triggers.txt'.format(opt.outputPath, opt.groupName), 'w') as f:
        
        startTimes = ttable.cols.startTimeMPL[:]
        
        for i in np.argsort(startTimes):
            f.write("{0}\n".format((UTCDateTime(matplotlib.dates.num2date(
                startTimes[i]))+opt.ptrig).isoformat()))

                    
def printOrphanCatalog(otable, opt):
    """
    Prints flat catalog of current orphans to text file
    
    otable: Orphans table
    opt: Options object describing station/run parameters
    
    Note: Time in text file corresponds to original STA/LTA trigger time
    """

    with open('{}{}/orphancatalog.txt'.format(opt.outputPath, opt.groupName), 'w') as f:
        
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
        plt.savefig('{}{}/junk/{}-{}.png'.format(opt.outputPath, opt.groupName,
            (UTCDateTime(r['startTime'])+opt.ptrig).strftime('%Y%m%d%H%M%S'),
            r['isjunk']), dpi=100)
        plt.close(fig)

            
def printJunk(jtable, opt):
    """
    Prints flat catalog of contents of junk table to text file for debugging
    
    jtable: Junk table
    opt: Options object describing station/run parameters
    
    Note: Time in text file corresponds to original STA/LTA trigger time
    """

    with open('{}{}/junk.txt'.format(opt.outputPath, opt.groupName), 'w') as f:
        
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

    with open('{}{}/cores.txt'.format(opt.outputPath, opt.groupName), 'w') as f:
        
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
    
    with open('{}{}/dailycounts.txt'.format(opt.outputPath, opt.groupName), 'w') as f:
        
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
    
    with open('{}{}/catalog.txt'.format(opt.outputPath, opt.groupName), 'w') as f:
                
        startTimes = rtable.cols.startTime[:]
        startTimeMPL = rtable.cols.startTimeMPL[:]
        windowStarts = rtable.cols.windowStart[:]
        windowAmps = rtable.cols.windowAmp[:]
        ids = rtable.cols.id[:]
        id1 = ctable.cols.id1[:]
        id2 = ctable.cols.id2[:]
        ccc = ctable.cols.ccc[:]
        fi = np.nanmean(rtable.cols.FI[:], axis=1)
        
        f.write("cnum\tevTime                        \tfi\tamps\tdt\t\txcorr\n")
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
                amp = windowAmps[fam[i],:]
                if j == -1:
                    dt = 'NaN         '
                else:
                    dt = spacing[j]
                j += 1
                
                f.write("{0}\t{1}\t{2:4.3f}\t[".format(cnum,evTime.isoformat(),fi[fam][i]))
                for a in amp:
                    f.write(" {} ".format(a))
                f.write("]\t{0}\t\t{1:3.2f}\n".format(dt,xcorr[i]))


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
    
    # Read in annotation file (if it exists)
    if opt.anotfile != '':
        df = pd.read_csv(opt.anotfile)
    
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
    shutil.copy('{}{}/clusters/{}.png'.format(opt.outputPath, opt.groupName, fnum),
                '{}{}/clusters/{}-report.png'.format(opt.outputPath, opt.groupName, fnum))
    
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
    
    # Amplitude vs. time on all stations with interactive show/hide
    # Set min/max for plotting
    if opt.amplims == 'family':
        windowAmpFam = windowAmps[fam[catalogind]][:]
        ymin = 0.25*np.amin(windowAmpFam[np.nonzero(windowAmpFam)])
        ymax = 4*np.amax(windowAmpFam)
    else:
        # Use global maximum
        ymin = 0.25*np.amin(windowAmps[np.nonzero(windowAmps)])
        ymax = 4*np.amax(windowAmps)            
    
    o0 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title='Amplitude with Time (Click name to hide)', y_axis_type='log',
        y_range=[ymin,ymax])
    o0.grid.grid_line_alpha = 0.3
    o0.xaxis.axis_label = 'Date'
    o0.yaxis.axis_label = 'Counts'
    if opt.anotfile != '':
        for row in df.itertuples():
            spantime = (datetime.datetime.strptime(row[1]
                ,'%Y-%m-%dT%H:%M:%S')-datetime.datetime(1970, 1, 1)).total_seconds()
            o0.add_layout(Span(location=spantime*1000, dimension='height',
                line_color=row[2], line_width=row[3], line_dash=row[4],
                line_alpha=row[5]))
    if opt.nsta <= 8:
        palette = all_palettes['YlOrRd'][9]
    else:
        palette = inferno(opt.nsta+1)
    for sta, staname in enumerate(opt.station.split(',')):
        o0.circle(matplotlib.dates.num2date(startTimeMPL[fam]), windowAmps[fam][:,sta],
            color=palette[sta], line_alpha=0, size=4, fill_alpha=0.5,
            legend_label='{}.{}'.format(staname,opt.channel.split(',')[sta]))    
    o0.legend.location='bottom_left'
    o0.legend.orientation='horizontal'
    o0.legend.click_policy='hide'
    
    
    # Time since last event
    o1 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title='Time since Previous Event', x_range=o0.x_range, y_axis_type='log',
        y_range=[1e-3, 2*np.max(spacing)])
    o1.grid.grid_line_alpha = 0.3
    o1.xaxis.axis_label = 'Date'
    o1.yaxis.axis_label = 'Interval (hr)'
    if opt.anotfile != '':
        for row in df.itertuples():
            spantime = (datetime.datetime.strptime(row[1]
                ,'%Y-%m-%dT%H:%M:%S')-datetime.datetime(1970, 1, 1)).total_seconds()
            o1.add_layout(Span(location=spantime*1000, dimension='height',
                line_color=row[2], line_width=row[3], line_dash=row[4],
                line_alpha=row[5]))
    o1.circle(matplotlib.dates.num2date(catalog[1:]), spacing, color='red',
        line_alpha=0, size=4, fill_alpha=0.5)
    
    # Cross-correlation wrt. core
    o2 = figure(tools=oTOOLS, plot_width=1250, plot_height=250, x_axis_type='datetime',
        title='Cross-correlation Coefficient with Core Event', x_range=o0.x_range,
        y_range=[0, 1.02])
    o2.grid.grid_line_alpha = 0.3
    o2.xaxis.axis_label = 'Date'
    o2.yaxis.axis_label = 'CCC'  
    if opt.anotfile != '':
        for row in df.itertuples():
            spantime = (datetime.datetime.strptime(row[1]
                ,'%Y-%m-%dT%H:%M:%S')-datetime.datetime(1970, 1, 1)).total_seconds()
            o2.add_layout(Span(location=spantime*1000, dimension='height',
                line_color=row[2], line_width=row[3], line_dash=row[4],
                line_alpha=row[5])) 
    o2.circle(matplotlib.dates.num2date(catalog), Cfull[np.where(famcat==core)[0],:][0],
        color='red', line_alpha=0, size=4, fill_alpha=0.5)
    
    # Combine and save
    o = gridplot([[o0],[o1],[o2]])
    output_file('{}{}/clusters/{}-report-bokeh.html'.format(opt.outputPath, opt.groupName,
        fnum), title='{} - Cluster {} Detailed Report'.format(opt.title, fnum))
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
        if opt.anotfile!='':
            for anot in range(len(df)):
                hloc = np.interp(matplotlib.dates.date2num(
                    pd.to_datetime(df['Time'][anot])),startTimeMPL[fam][catalogind],
                    np.array(range(len(fam))))
                if hloc!=0:
                    ax1.axhline(np.floor(hloc)+0.5,color='k',
                        linewidth=df['Weight'][anot]/2.,linestyle=df['Line Type'][anot])
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
        if opt.anotfile!='':
            for anot in range(len(df)):
                hloc = np.interp(matplotlib.dates.date2num(
                    pd.to_datetime(df['Time'][anot])),startTimeMPL[fam][catalogind],
                    np.array(range(len(fam))))
                if hloc!=0:
                    ax2.axhline(np.floor(hloc)+0.5,color='k',
                        linewidth=df['Weight'][anot]/2.,linestyle=df['Line Type'][anot])
    plt.tight_layout()
    plt.savefig('{}{}/clusters/{}-reportcmat.png'.format(opt.outputPath, opt.groupName,
                                                         fnum), dpi=100)
    plt.close(fig)
    
    ### WAVEFORM IMAGES
    famtable = rtable[famcat]
    fig2 = plt.figure(figsize=(10, 12))
    
    for sta in range(opt.nsta):
        n = -1
        data = np.zeros((len(fam), int(opt.winlen*2)))
        ax = fig2.add_subplot(np.ceil((opt.nsta)/2.), 2, sta+1)
        for r in famtable:
            if ordered:
                plt.title('{0}.{1} (Ordered)'.format(opt.station.split(',')[sta],
                          opt.channel.split(',')[sta]), fontweight='bold')
            else:
                plt.title('{0}.{1}'.format(opt.station.split(',')[sta],
                          opt.channel.split(',')[sta]), fontweight='bold')
                if opt.anotfile!='':
                    for anot in range(len(df)):
                        hloc = np.interp(matplotlib.dates.date2num(
                            pd.to_datetime(df['Time'][anot])),
                            startTimeMPL[fam][catalogind],np.array(range(len(fam))))
                        if hloc!=0:
                            ax.axhline(np.floor(hloc)+0.5,color='k',
                                linewidth=df['Weight'][anot]/2.,
                                linestyle=df['Line Type'][anot])
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
    plt.savefig('{}{}/clusters/{}-reportwaves.png'.format(opt.outputPath, opt.groupName,
                fnum), dpi=100)
    plt.close(fig2)
    
    ### HTML OUTPUT PAGE
    tstamp = UTCDateTime.now()
    with open('{}{}/clusters/{}-report.html'.format(opt.outputPath, opt.groupName, fnum),
              'w') as f:
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


def printSwarmCatalog(rtable, ftable, ttable, opt):
        
    """
    Writes a .csv file for use in annotating repeating events in Swarm v2.8.5+
    
    rtable: Repeater table
    ftable: Families table
    opt: Options object describing station/run parameters
    
    """
    
    nets = opt.network.split(',')
    stas = opt.station.split(',')
    locs = opt.location.split(',')
    chas = opt.channel.split(',')
    
    with open('{}{}/swarm.csv'.format(opt.outputPath, opt.groupName), 'w') as f:
        
        startTimes = rtable.cols.startTime[:]
        windowStarts = rtable.cols.windowStart[:]
        
        for cnum in range(ftable.attrs.nClust):
            fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
            for i in np.argsort(startTimes[fam]):
                # Format for Swarm is 'Date Time, STA CHA NET LOC, label'
                # The SCNL defaults to whichever station was chosen for the preview,
                # which can be changed by a global search/replace in a text editor.
                # The label name is the same as the folder name (groupName) followed by
                # the family number. Highlighting families of interest in a different
                # color can be done by editing the EventClassifications.config file in
                # the Swarm folder, and adding a line for each cluster of interest
                # followed by a hex code for color, such as:
                # default1, #ffff00
                # to highlight family 1 from the default run in yellow compared to other
                # repeaters in red/orange.
                f.write("{}, {} {} {} {}, {}{}\n".format((UTCDateTime(startTimes[fam][i])+
                    windowStarts[fam][i]/opt.samprate).isoformat(sep=' '),
                    stas[opt.printsta],chas[opt.printsta],nets[opt.printsta],
                    locs[opt.printsta],opt.groupName,cnum))
                    
    with open('{}{}/triggerswarm.csv'.format(opt.outputPath, opt.groupName), 'w') as f:
        
        startTimes = ttable.cols.startTimeMPL[:]
        
        for i in np.argsort(startTimes):
            f.write("{}, {} {} {} {}, trigger\n".format((UTCDateTime(
                matplotlib.dates.num2date(startTimes[i]))+opt.ptrig).isoformat(sep=' '),
                stas[opt.printsta],chas[opt.printsta],nets[opt.printsta],
                    locs[opt.printsta]))
