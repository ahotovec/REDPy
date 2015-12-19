from obspy import UTCDateTime
from obspy.fdsn import Client
from obspy.earthworm import Client as EWClient
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.signal.trigger import classicSTALTA, triggerOnset
import numpy as np
from scipy import stats
from scipy.fftpack import fft

def getData(date, opt):

    """
    Download data from IRIS or Earthworm waveserver with padding and filter it.

    date: UTCDateTime of beginning of period of interest
    opt: Options object describing station/run parameters
    
    Returns ObsPy stream object
    """    
    
    # Choose where data are downloaded automatically via options
    # Download data with padding to account for triggering algorithm
    # Make overlap symmetric
    
    if opt.server == "IRIS":
        client = Client("IRIS")
        st = client.get_waveforms(opt.network, opt.station, opt.location, opt.channel,
            date - opt.atrig, date + opt.nsec + opt.atrig)
    else:
        client = EWClient(opt.server, opt.port)
        st = client.getWaveform(opt.network, opt.station, opt.location, opt.channel,
            date - opt.atrig, date + opt.nsec + opt.atrig)

    st = st.detrend() # can create noise artifacts??
    st = st.merge(method=1, fill_value='interpolate')

    return st


def getCatData(date, opt):

    """
    Download data from IRIS or Earthworm waveserver with padding and filter it. This is
    a specialized version getData() for catalog events, pulling a smaller amount of time
    around a known event.

    date: UTCDateTime of known catalog event
    opt: Options object describing station/run parameters
    
    Returns ObsPy stream object
    """    
    
    # Choose where data are downloaded automatically via options
    # Download data with padding to account for triggering algorithm
    # Make overlap symmetric
    
    if opt.server == "IRIS":
        client = Client("IRIS")
        st = client.get_waveforms(opt.network, opt.station, opt.location, opt.channel,
            date - opt.atrig, date + 3*opt.atrig)
    else:
        client = EWClient(opt.server, opt.port)
        st = client.getWaveform(opt.network, opt.station, opt.location, opt.channel,
            date - opt.atrig, date + 3*opt.atrig)

    st = st.detrend() # can create noise artifacts??
    st = st.merge(method=1, fill_value='interpolate')

    return st


def trigger(st, rtable, opt):

    """
    Run triggering algorithm on a stream of data.

    st: OBSPy stream of data
    rtable: Repeater table contains reference time of previous trigger in samples
    opt: Options object describing station/run parameters

    Returns triggered traces as OBSPy trace object updates ptime for next run 
    """
    
    if st[0].stats.sampling_rate != opt.samprate:
        print('Warning: Resampling to match configuration...')
        st = st.resample(opt.samprate)

    # Filter the data for triggering
    st = st.filter("bandpass", freqmin=opt.fmin, freqmax=opt.fmax, corners=2,
               zerophase=True)
    tr = st[0]
    t = tr.stats.starttime

    cft = classicSTALTA(tr.data, opt.swin*opt.samprate, opt.lwin*opt.samprate)
    on_off = triggerOnset(cft, opt.trigon, opt.trigoff)
    
    if len(on_off) > 0:
    
        pick = on_off[:,0]    
        ind = 0
        
        # Slice out the data and save the maximum STA/LTA ratio value for
        # use in orphan expiration
        
        # Convert ptime from time of last trigger to samples before start time 
        if rtable.attrs.ptime:
            ptime = (UTCDateTime(rtable.attrs.ptime) - t)*opt.samprate
        else:
            ptime = -opt.mintrig*opt.samprate
        
        for n in range(len(pick)):
            
            ttime = pick[n]
            
            if (ttime >= opt.atrig*opt.samprate) and (ttime >= ptime +
                opt.mintrig*opt.samprate) and (ttime < len(tr.data) -
                2*opt.atrig*opt.samprate):
                
                ptime = ttime
                if ind is 0:
                    # Slice and save as first trace              
                    trigs = st.slice(t - opt.ptrig + ttime/opt.samprate,
                                     t + opt.atrig + ttime/opt.samprate)
                    trigs[ind].stats.maxratio = np.amax(cft[on_off[n,0]:(on_off[n,1]+1)])
                    ind = ind+1
                else:
                    # Slice and append to previous traces
                    trigs = trigs.append(tr.slice(
                            t - opt.ptrig + ttime/opt.samprate,
                            t + opt.atrig + ttime/opt.samprate))
                    trigs[ind].stats.maxratio = np.amax(cft[on_off[n,0]:(on_off[n,1]+1)])
                    ind = ind+1
                                                         
        if ind is 0:
            rtable.attrs.ptime = (t + len(tr.data)/opt.samprate -
                opt.mintrig*opt.samprate).isoformat()
            return []
        else:
            rtable.attrs.ptime = (t + ptime/opt.samprate).isoformat()
            return trigs
    else:
        rtable.attrs.ptime = (t + len(tr.data)/opt.samprate -
            opt.mintrig*opt.samprate).isoformat()
        return []


def dataclean(alltrigs, opt, flag=1):

    """
    Examine triggers and weed out spikes and calibration pulses using kurtosis and
    outlier ratios
    
    alltrigs: triggers output from triggering
    opt: opt from config
    flag: 1 if defining window to check, 0 if want to check whole waveform for spikes
        (note that different threshold values should be used for different window lengths)
    
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
            datcut=alltrigs[i].data[range(int((opt.ptrig-opt.kurtwin/2)*opt.samprate),
                int((opt.ptrig+opt.kurtwin/2)*opt.samprate))]
        
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
