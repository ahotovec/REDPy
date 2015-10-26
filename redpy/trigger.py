from obspy import UTCDateTime
from obspy.fdsn import Client
from obspy.earthworm import Client as EWClient
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.signal.trigger import classicSTALTA, triggerOnset
import numpy as np
from scipy import stats
from scipy.fftpack import fft

def getIRIS(date, opt, nsec=86400):

    """
    Download data from IRIS with padding and filter it.

    date: UTCDateTime of beginning of period of interest
    opt: Options object describing station/run parameters
    nsec: Number of seconds to download without padding
        (default 86400 s, or 1 day)
    
    Returns ObsPy stream object
    """    

    client = Client("IRIS")
    
    # Download data with padding to account for triggering algorithm
    # Make overlap symmetric
    st = client.get_waveforms(opt.network, opt.station, opt.location, opt.channel,
        date - opt.atrig, date + nsec + opt.atrig)

    st = st.detrend() # can create noise artifacts??
    st = st.merge(method=1, fill_value='interpolate')
    st = st.filter("highpass", freq=opt.fhigh, corners=2,
            zerophase=True)

    return st


def getEarthworm(date, server, port, opt, nsec=86400):
    
    """
    Download data from a (local) Earthworm waveserver.
    
    date: UTCDateTime of beginning of period of interest
    server: String name of server
    port: Port number
    opt: Options object describing station/run parameters
    nsec: Number of seconds to download without padding
        (default 86400 s, or 1 day)
    
    Returns ObsPy stream object
    """
    
    client = EWClient(server, port)
    
    # Download data with padding to account for triggering algorithm
    # Make overlap symmetric
    st = client.getWaveform(opt.network, opt.station, opt.location, opt.channel,
        date - opt.atrig, date + nsec + opt.atrig)

    st = st.detrend() # can create noise artifacts??
    st = st.merge(method=1, fill_value='interpolate')
    st = st.filter("highpass", freq=opt.fhigh, corners=2,
            zerophase=True)

    return st


def trigger(st, ptime, opt):

    """
    Run triggering algorithm on a stream of data.

    st: OBSPy stream of data
    ptime: time or previous trigger in samples
    opt: Options object describing station/run parameters

    Returns triggered traces as OBSPy trace object and ptime for next run 
    """

    # Filter the data for triggering

    st_f = st.copy()
    st_f = st_f.filter("bandpass", freqmin=opt.fmin, freqmax=opt.fmax, corners=2,
               zerophase=True)
    tr = st[0]
    tr_f = st_f[0]
    t = tr.stats.starttime

    cft = classicSTALTA(tr_f.data, opt.swin*opt.samprate, opt.lwin*opt.samprate)
    on_off = triggerOnset(cft, opt.trigon, opt.trigoff)
    
    if len(on_off) > 0:
    
        pick = on_off[:,0]    
        ind = 0
        
        # Slice out the raw data (not filtered except for highpass to reduce long period
        # drift) and save the maximum STA/LTA ratio value for
        # use in orphan expiration

        for n in range(len(pick)):
            
            ttime = pick[n]
            
            if (ttime >= opt.atrig*opt.samprate) and (ttime >= ptime + opt.mintrig*opt.samprate) and (ttime < len(tr.data) - 2*opt.atrig*opt.samprate):
                
                ptime = ttime
                if ind is 0:
                    # Slice and save as first trace              
                    trigs = st.slice(t - opt.ptrig + ttime/opt.samprate,
                                     t + opt.atrig + ttime/opt.samprate)
                    trigs[ind].stats.maxratio = np.amax(cft[on_off[n,0]:on_off[n,1]])
                    ind = ind+1
                else:
                    # Slice and append to previous traces
                    trigs = trigs.append(tr.slice(
                            t - opt.ptrig + ttime/opt.samprate,
                            t + opt.atrig + ttime/opt.samprate))
                    trigs[ind].stats.maxratio = np.amax(cft[on_off[n,0]:on_off[n,1]])
                    ind = ind+1
        
        ptime = ptime - len(tr.data) + opt.ptrig*opt.samprate                                      
        if ind is 0:
            return [], -len(tr.data)
        else:
            return trigs, ptime
    else:
        return [], -len(tr.data)


def dataclean(alltrigs, opt, flag=1):
    """
    Examine triggers and weed out spikes and calibration pulses using kurtosis and outlier ratios
    
    alltrigs: triggers output from triggering
    opt: opt from config
    flag: 1 if defining window to check, 0 if want to check whole waveform for spikes (note that different threshold values should be used for different window lengths)
    
    Returns good trigs (trigs) and junk (junk)
    
    """
    trigs=Stream()
    junk=Stream()
    for i in range(len(alltrigs)):
        #define data
        dat=alltrigs[i].data
        if flag==0:
            datcut=dat
        else:
            datcut=alltrigs[i].data[range(int((opt.ptrig-opt.kurtwin/2)*opt.samprate),int((opt.ptrig+opt.kurtwin/2)*opt.samprate))]
        
        #calculate kurtosis in window
        k = stats.kurtosis(datcut)
        #compute kurtosis of frequency amplitude spectrum next
        
        datf = np.absolute(fft(dat))
        kf = stats.kurtosis(datf)
        #calculate outlier ratio using z ((data-median)/mad), outliers have z>4.45
        mad = np.median(np.absolute(dat - np.median(dat)))
        z=(dat-np.median(dat))/mad
        orm = len(z[z>4.45])/len(z)
        if k<opt.kurtmax and orm<opt.oratiomax and kf<opt.kurtfmax:
            trigs.append(alltrigs[i])
        else:
            junk.append(alltrigs[i])
                
    return trigs, junk


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
