from obspy import UTCDateTime
from obspy.fdsn import Client
from obspy.core.trace import Trace
from obspy.signal.trigger import classicSTALTA, triggerOnset

def getIRIS(
    date, sta, chan, net, loc="--", nsec=86400, ptrig=10.0, atrig=20.0,
    fmin=1.0, fmax=10.0):

    """
    Download data from IRIS with padding and filter it.

    date: UTCDateTime of beginning of period of interest
    sta: String of station
    chan: String of channel
    net: String of network
    loc: String of location (default "--")
    nsec: Number of seconds to download without padding
        (default 86400 s, or 1 day)
    ptrig: Length of window to keep prior to trigger (default 10.0 s)
    atrig: Length of window to keep after trigger (default 20.0 s)
    fmin: Lower bound of bandpass filter (default 1.0 Hz)
    fmax: Upper bound of bandpass filter (default 10.0 Hz)

    Returns ObsPy stream object
    """    

    client = Client("IRIS")

    # Download data with padding to account for triggering algorithm
    st = client.get_waveforms(
        net, sta, loc, chan, date - ptrig, date + nsec + atrig)

    st = st.detrend() # can create noise artifacts??
    st = st.merge(method=1, fill_value=0)
    st = st.filter("bandpass", freqmin=fmin, freqmax=fmax,
                   corners=2, zerophase=True)

    return st

def trigger(
    st, lwin=7.0, swin=0.8, trigon=3.0, trigoff=2.0, mintrig=10.0,
    ptrig=10.0, atrig=20.0):

    """
    Run triggering algorithm on a stream of data.

    st: OBSPy stream of data
    lwin: Length of long window for STALTA (default 7.0 s)
    swin: Length of short window for STALTA (default 0.8 s)
    trigon: Cutoff ratio for triggering STALTA (default 3.0)
    trigoff: Cutoff ratio for ending STALTA trigger (default 2.0)
    mintrig: Minimum spacing between triggers (default 10.0 s)
    ptrig: Amount to cut prior to trigger (default 10.0 s)
    atrig: Amount to cut after trigger (default 20.0 s)

    Returns triggered traces as OBSPy trace object
    """

    tr = st[0]
    srate = tr.stats.sampling_rate
    t = tr.stats.starttime

    cft = classicSTALTA(tr.data, swin*srate, lwin*srate)
    on_off = triggerOnset(cft, trigon, trigoff)

    ttime = 0
    for n in range(len(on_off)):
        if on_off[n, 0] > ttime + mintrig*srate:
            if ttime is 0 and on_off[n, 0] > ttime + ptrig*srate:
                ttime = on_off[n, 0]
                trigs = st.slice(t - ptrig + ttime/srate,
                                 t + atrig + ttime/srate)
            else:
                ttime = on_off[n,0]
                if ttime < len(tr.data) - (atrig + ptrig)*srate:
                    trigs = trigs.append(tr.slice(
                        t - ptrig + ttime/srate,
                        t + atrig + ttime/srate))

    return trigs