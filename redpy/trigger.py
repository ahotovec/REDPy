from obspy import UTCDateTime
from obspy.fdsn import Client
from obspy.core.trace import Trace
from obspy.signal.trigger import classicSTALTA, triggerOnset
import numpy as np

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
    st = client.get_waveforms(opt.network, opt.station, opt.location, opt.channel,
        date - opt.ptrig, date + nsec + opt.atrig)

    st = st.detrend() # can create noise artifacts??
    st = st.merge(method=1, fill_value='interpolate')
    st = st.filter("highpass", freq=opt.fhigh, corners=2,
            zerophase=True)
    #st = st.filter("bandpass", freqmin=opt.fmin, freqmax=opt.fmax, corners=2,
    #                 zerophase=True)

    return st


def trigger(st, opt):

    """
    Run triggering algorithm on a stream of data.

    st: OBSPy stream of data
    opt: Options object describing station/run parameters

    Returns triggered traces as OBSPy trace object
    """

    #filter the data for triggering
    st_f = st.filter("bandpass", freqmin=opt.fmin, freqmax=opt.fmax, corners=2,
               zerophase=True)
    tr = st_f[0]
    t = tr.stats.starttime

    cft = classicSTALTA(tr.data, opt.swin*opt.samprate, opt.lwin*opt.samprate)
    on_off = triggerOnset(cft, opt.trigon, opt.trigoff)
    
    pick = np.zeros([len(on_off),1])
    for n in range(len(pick)):
        pick[n] = aicpick(st_f, on_off[n, 0], opt)

    ttime = 0
    
    #slice out the raw data, not filtered except for a lowpass to reduce long period drift
    for n in range(len(on_off)):
        if on_off[n, 0] > ttime + opt.mintrig*opt.samprate:
            if ttime is 0 and pick[n] > ttime + opt.ptrig*opt.samprate:
                ttime = pick[n]
                trigs = st.slice(t - opt.ptrig + ttime/opt.samprate,
                                 t + opt.atrig + ttime/opt.samprate)
            else:
                ttime = pick[n]
                if ttime < len(tr.data) - (opt.atrig + opt.ptrig)*opt.samprate:
                    trigs = trigs.append(tr.slice(
                        t - opt.ptrig + ttime/opt.samprate,
                        t + opt.atrig + ttime/opt.samprate))

    return trigs
    

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
