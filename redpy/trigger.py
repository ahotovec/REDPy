# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2018  Alicia Hotovec-Ellis (ahotovec@gmail.com)
# Licensed under GNU GPLv3 (see LICENSE.txt)

from obspy import UTCDateTime
import obspy
from obspy.clients.fdsn import Client
from obspy.clients.earthworm import Client as EWClient
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.signal.trigger import coincidence_trigger
import numpy as np
from scipy import stats
from scipy.fftpack import fft
import glob, os, itertools

import warnings
warnings.filterwarnings("ignore")

def getData(tstart, tend, opt):

    """
    Download data from files in a folder, from IRIS, or a Earthworm waveserver
    
    A note on SAC/miniSEED files: as this makes no assumptions about the naming scheme of
    your data files, please ensure that your headers contain the correct SCNL information!

    tstart: UTCDateTime of beginning of period of interest
    tend: UTCDateTime of end of period of interest
    opt: Options object describing station/run parameters
    
    Returns ObsPy stream objects, one for cutting and the other for triggering
    """    
    
    nets = opt.network.split(',')
    stas = opt.station.split(',')
    locs = opt.location.split(',')
    chas = opt.channel.split(',')
    
    st = Stream()
    
    if opt.server == 'file':
    
        # Generate list of files
        if opt.server == 'file':
            flist = list(itertools.chain.from_iterable(glob.iglob(os.path.join(
                root,opt.filepattern)) for root, dirs, files in os.walk(opt.searchdir)))
                
        # Determine which subset of files to load based on start and end times and
        # station name; we'll fully deal with stations below
        flist_sub = []
        for f in flist:
            # Load header only
            stmp = obspy.read(f, headonly=True)
            # Check if station is contained in the stas list
            if stmp[0].stats.station in stas:
                # Check if contains either start or end time
                ststart = stmp[0].stats.starttime
                stend = stmp[-1].stats.endtime
                if (ststart<=tstart and tstart<=stend) or (ststart<=tend and
                    tend<=stend) or (tstart<=stend and ststart<=tend):
                    flist_sub.append(f)
        
        # Fully load data from file
        stmp = Stream()
        for f in flist_sub:
            tmp = obspy.read(f, starttime=tstart, endtime=tend+opt.maxdt)
            if len(tmp) > 0:
                stmp = stmp.extend(tmp)
    
        # Filter and merge
        stmp = stmp.filter('bandpass', freqmin=opt.fmin, freqmax=opt.fmax, corners=2,
            zerophase=True)
        stmp = stmp.taper(0.05,type='hann',max_length=opt.mintrig)
        for m in range(len(stmp)):
            if stmp[m].stats.sampling_rate != opt.samprate:
                stmp[m] = stmp[m].resample(opt.samprate)
        stmp = stmp.merge(method=1, fill_value=0)
        
        # Only grab stations/channels that we want and in order
        netlist = []
        stalist = []
        chalist = []
        loclist = []
        for s in stmp:
            stalist.append(s.stats.station)
            chalist.append(s.stats.channel)
            netlist.append(s.stats.network)
            loclist.append(s.stats.location)
            
        # Find match of SCNL in header or fill empty
        for n in range(len(stas)):
            for m in range(len(stalist)):
                if (stas[n] in stalist[m] and chas[n] in chalist[m] and nets[n] in
                    netlist[m] and locs[n] in loclist[m]):
                    st = st.append(stmp[m])
            if len(st) == n:
                print("Couldn't find "+stas[n]+'.'+chas[n]+'.'+nets[n]+'.'+locs[n])
                trtmp = Trace()
                trtmp.stats.sampling_rate = opt.samprate
                trtmp.stats.station = stas[n]
                st = st.append(trtmp.copy())
    
    else:   
     
        if '.' not in opt.server:
            client = Client(opt.server)
        else:
            client = EWClient(opt.server, opt.port)
        
        for n in range(len(stas)):
            try:
                stmp = client.get_waveforms(nets[n], stas[n], locs[n], chas[n],
                        tstart, tend+opt.maxdt)
                for m in range(len(stmp)):
                    stmp[m].data = np.where(stmp[m].data == -2**31, 0, stmp[m].data) # replace -2**31 (Winston NaN token) w 0
                stmp = stmp.filter('bandpass', freqmin=opt.fmin, freqmax=opt.fmax,
                    corners=2, zerophase=True)
                stmp = stmp.taper(0.05,type='hann',max_length=opt.mintrig)
                for m in range(len(stmp)):
                    if stmp[m].stats.sampling_rate != opt.samprate:
                        stmp[m] = stmp[m].resample(opt.samprate)
                stmp = stmp.merge(method=1, fill_value=0)
            except (obspy.clients.fdsn.header.FDSNException):
                try: # try again
                    stmp = client.get_waveforms(nets[n], stas[n], locs[n], chas[n],
                            tstart, tend+opt.maxdt)
                    for m in range(len(stmp)):
                        stmp[m].data = np.where(stmp[m].data == -2**31, 0, stmp[m].data) # replace -2**31 (Winston NaN token) w 0
                    stmp = stmp.filter('bandpass', freqmin=opt.fmin, freqmax=opt.fmax,
                        corners=2, zerophase=True)
                    stmp = stmp.taper(0.05,type='hann',max_length=opt.mintrig)
                    for m in range(len(stmp)):
                        if stmp[m].stats.sampling_rate != opt.samprate:
                            stmp[m] = stmp[m].resample(opt.samprate)
                    stmp = stmp.merge(method=1, fill_value=0)
                except (obspy.clients.fdsn.header.FDSNException):
                    print('No data found for {0}.{1}'.format(stas[n],nets[n]))
                    trtmp = Trace()
                    trtmp.stats.sampling_rate = opt.samprate
                    trtmp.stats.station = stas[n]
                    stmp = Stream().extend([trtmp.copy()])
                                            
            # Last check for length; catches problem with empty waveserver
            if len(stmp) != 1:
                print('No data found for {0}.{1}'.format(stas[n],nets[n]))
                trtmp = Trace()
                trtmp.stats.sampling_rate = opt.samprate
                trtmp.stats.station = stas[n]
                stmp = Stream().extend([trtmp.copy()])
                
            st.extend(stmp.copy()) 
    
    # Edit 'start' time if using offset option
    if opt.maxdt:
        dts = np.fromstring(opt.offset, sep=',')
        for n, tr in enumerate(st):
            tr.stats.starttime = tr.stats.starttime-dts[n]
    
    st = st.trim(starttime=tstart, endtime=tend, pad=True, fill_value=0)
    stC = st.copy()
    
    return st, stC


def trigger(st, stC, rtable, opt):

    """
    Run triggering algorithm on a stream of data.

    st: OBSPy stream of data
    rtable: Repeater table contains reference time of previous trigger in samples
    opt: Options object describing station/run parameters

    Returns triggered traces as OBSPy trace object updates ptime for next run 
    """
    
    tr = st[0]
    t = tr.stats.starttime

    cft = coincidence_trigger(opt.trigalg, opt.trigon, opt.trigoff, stC, opt.nstaC,
        sta=opt.swin, lta=opt.lwin, details=True)
            
    if len(cft) > 0:
        
        ind = 0
        
        # Slice out the data from st and save the maximum STA/LTA ratio value for
        # use in orphan expiration
        
        # Convert ptime from time of last trigger to seconds before start time
        if rtable.attrs.ptime:
            ptime = (UTCDateTime(rtable.attrs.ptime) - t)
        else:
            ptime = -opt.mintrig
                
        for n in range(len(cft)):
                    
            ttime = cft[n]['time'] # This is a UTCDateTime, not samples
            
            if (ttime >= t + opt.atrig) and (ttime >= t + ptime +
                opt.mintrig) and (ttime < t + len(tr.data)/opt.samprate -
                2*opt.atrig):
                
                ptime = ttime - t
                
                # Cut out and append all data to first trace              
                tmp = st.slice(ttime - opt.ptrig, ttime + opt.atrig)
                ttmp = tmp.copy()
                ttmp = ttmp.trim(ttime - opt.ptrig, ttime + opt.atrig + 0.05, pad=True,
                    fill_value=0)
                ttmp[0].data = ttmp[0].data[0:opt.wshape] - np.mean(
                    ttmp[0].data[0:opt.wshape])
                for s in range(1,len(ttmp)):
                    ttmp[0].data = np.append(ttmp[0].data, ttmp[s].data[
                        0:opt.wshape] - np.mean(ttmp[s].data[0:opt.wshape]))
                ttmp[0].stats.maxratio = np.max(cft[n]['cft_peaks'])
                if ind is 0:
                    trigs = Stream(ttmp[0])
                    ind = ind+1
                else:
                    trigs = trigs.append(ttmp[0])
                                                         
        if ind is 0:
            return []
        else:
            rtable.attrs.ptime = (t + ptime).isoformat()
            return trigs
    else:
        return []


def dataClean(alltrigs, opt, flag=1):

    """
    Examine triggers and weed out spikes and calibration pulses using kurtosis and
    outlier ratios
    
    alltrigs: triggers output from triggering
    opt: opt from config
    flag: 1 if defining window to check, 0 if want to check whole waveform for spikes
        (note that different threshold values should be used for different window lengths)
    
    Returns good trigs (trigs) and several junk types (junk, junkFI, junkKurt)
    """
    
    trigs=Stream()
    junkFI=Stream()
    junkKurt=Stream()
    junk=Stream()
    for i in range(len(alltrigs)):
            
        njunk = 0
        ntele = 0
        
        for n in range(opt.nsta):
            
            dat = alltrigs[i].data[n*opt.wshape:(n+1)*opt.wshape]
            if flag == 1:
                datcut=dat[range(int((opt.ptrig-opt.kurtwin/2)*opt.samprate),
                    int((opt.ptrig+opt.kurtwin/2)*opt.samprate))]
            else:
                datcut=dat
            
            if np.sum(np.abs(dat))!=0.0:
                # Calculate kurtosis in window
                k = stats.kurtosis(datcut)
                # Compute kurtosis of frequency amplitude spectrum next
                datf = np.absolute(fft(dat))
                kf = stats.kurtosis(datf)
                # Calculate outlier ratio using z ((data-median)/mad)
                mad = np.nanmedian(np.absolute(dat - np.nanmedian(dat)))
                z = (dat-np.median(dat))/mad
                # Outliers have z > 4.45
                orm = len(z[z>4.45])/np.array(len(z)).astype(float)
            
                if k >= opt.kurtmax or orm >= opt.oratiomax or kf >= opt.kurtfmax:
                    njunk+=1
                
                winstart = int(opt.ptrig*opt.samprate - opt.winlen/10)
                winend = int(opt.ptrig*opt.samprate - opt.winlen/10 + opt.winlen)
                fftwin = np.reshape(fft(dat[winstart:winend]),(opt.winlen,))
                if np.median(np.abs(dat[winstart:winend]))!=0:
                    fi = np.log10(np.mean(np.abs(np.real(
                        fftwin[int(opt.fiupmin*opt.winlen/opt.samprate):int(
                        opt.fiupmax*opt.winlen/opt.samprate)])))/np.mean(np.abs(np.real(
                        fftwin[int(opt.filomin*opt.winlen/opt.samprate):int(
                        opt.filomax*opt.winlen/opt.samprate)]))))
                    if fi<opt.telefi:
                        ntele+=1
        
        # Allow if there are enough good stations to correlate
        if njunk <= (opt.nsta-opt.ncor) and ntele <= opt.teleok:
            trigs.append(alltrigs[i])
        else:
            if njunk > 0:
                if ntele > 0:
                    junk.append(alltrigs[i])
                else:
                    junkKurt.append(alltrigs[i])
            else:
                junkFI.append(alltrigs[i])
                
    return trigs, junk, junkFI, junkKurt


def aicpick(st, initialTrigger, opt):
    
    """
    An autopicker to (hopefully) improve consistency in triggering
    
    st: OBSPy stream of data containing trigger
    initialTrigger: initial guess at trigger time (in number of samples into stream)
    opt: Options object describing station/run parameters
    
    Returns updated trigger time
    
    AIC stands for Akaike Information Criterion. This code is based on the formula in
    Zhang, Thurber, and Rowe [2003] (originally from Maeda [1985]) to calculate AIC
    directly from the waveform. It is a purely statistical picker, and the minimum of
    the AIC corresponds to where one can divide the signal into two different parts (in
    this case, noise followed by signal).
    """
    
    t = st[0].stats.starttime   
    x0 = st.slice(t - opt.ptrig/2 + initialTrigger/opt.samprate,
                  t + opt.ptrig/2 + initialTrigger/opt.samprate)
    x = x0[0].data
    nsamp = int(opt.ptrig*opt.samprate)
    
    AIC = np.zeros([nsamp,1])
        
    for k in range(nsamp):
            
        # Calculate the Akaike Information Criteria
        var1 = np.var(x[0:k+1])
        var2 = np.var(x[k:nsamp+1])
        
        if var1 == 0 or var2 == 0:
            AIC[k] = np.NaN
        else:
            AIC[k] = (k+1)*np.log10(var1) + (nsamp-k)*np.log10(var2)
            
    # Pad 10 samples on either end to prevent encountering edge effects
    picksamp = np.argmin(AIC[10:nsamp-10]) + initialTrigger - nsamp/2
                
    return picksamp
