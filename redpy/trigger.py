from obspy import UTCDateTime
from obspy.fdsn import Client
from obspy.core.trace import Trace
from obspy.signal.trigger import classicSTALTA, triggerOnset

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
    st = st.filter("bandpass", freqmin=opt.fmin, freqmax=opt.fmax, corners=2,
        zerophase=True)

    return st


def trigger(st, opt):

    """
    Run triggering algorithm on a stream of data.

    st: OBSPy stream of data
    opt: Options object describing station/run parameters

    Returns triggered traces as OBSPy trace object
    """

    tr = st[0]
    t = tr.stats.starttime

    cft = classicSTALTA(tr.data, opt.swin*opt.samprate, opt.lwin*opt.samprate)
    on_off = triggerOnset(cft, opt.trigon, opt.trigoff)

    ttime = 0
    for n in range(len(on_off)):
        if on_off[n, 0] > ttime + opt.mintrig*opt.samprate:
            if ttime is 0 and on_off[n, 0] > ttime + opt.ptrig*opt.samprate:
                ttime = on_off[n, 0]
                trigs = st.slice(t - opt.ptrig + ttime/opt.samprate,
                                 t + opt.atrig + ttime/opt.samprate)
            else:
                ttime = on_off[n,0]
                if ttime < len(tr.data) - (opt.atrig + opt.ptrig)*opt.samprate:
                    trigs = trigs.append(tr.slice(
                        t - opt.ptrig + ttime/opt.samprate,
                        t + opt.atrig + ttime/opt.samprate))

    return trigs