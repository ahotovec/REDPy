class Options(object):
    
    def __init__(self, title="REDPy Catalog", filename="redtable.h5", groupName="hsr",
        groupDesc="MSH: HSR-EHZ-UW Default", station="HSR", channel="EHZ", network="UW",
        location="--", samprate=100.0, lwin=7.0, swin=0.8, trigon=3.0, trigoff=2.0,
        mintrig=10.0, kurtmax=80., kurtfmax=150., oratiomax=0.06, kurtwin=5., winlen=512,
        ptrig=10.0, atrig=20.0, fhigh=0.25, fmin=1., fmax=10.0, cmin=0.7):
        
        """
        Defines the settings that are often passed to routines and that define the table.
        These are also written to the attributes of the table for posterity.
    
        TABLE DEFINITIONS:
        title: Name of the table (default "REDPy Catalog")
        filename: Filename for the table (default "redtable.h5")
        groupName: Short string describing the name of the station (default "hsr")
        groupDesc: Longer string describing the run (default "MSH: HSR-EHZ-UW Default")
    
        STATION PARAMETERS:
        station: String of station name (default "HSR")
        channel: String of channel of interest, no wildcards supported yet (default "EHZ")
        network: String of network code (default "UW")
        location: String of location code (default "--")
        samprate: Sampling rate of that station (default 100.0 Hz)
        
        TRIGGERING PARAMETERS:
        lwin: Length of long window for STALTA (default 7.0 s)
        swin: Length of short window for STALTA (default 0.8 s)
        trigon: Cutoff ratio for triggering STALTA (default 3.0)
        trigoff: Cutoff ratio for ending STALTA trigger (default 2.0)
        mintrig: Minimum spacing between triggers (default 10.0 s)
        kurtmax: Maximum kurtosis allowed for event window, to eliminate spikes (default 80) ~80-100 is appropriate for 5s window, ~130 for 15s, ~200 for 25s
        kurtfmax: Maximum kurtosis of frequency amplitude spectrum (default 150) to eliminate calibration pulses with unnaturally harmonic signals, be careful not to set too low or you could eliminate real harmonic events
        kurtwin: Length of window to use for kurtosis, in seconds, around the trigger time, will be centered on the trigger time (default 5s)
        oratiomax: Maximum ratio of outliers to total number of datapoints in trace (default 0.06, 6%)

    
        WINDOWING PARAMETERS:
        winlen: Length of window for cross-correlation (default 512 samples, 2^n is best)
        ptrig: Length of time cut prior to trigger (default 10.0 s)
        atrig: Length of time cut after trigger (default 20.0 s)
    
        FILTERING PARAMETERS:
        fhigh: Highpass filter to reduce microseism in "raw" waveforms (default 0.25 Hz)
        fmin: Lower band of bandpass filter for triggering and xcorr (default 1.0 Hz)
        fmax: Upper band of bandpass filter for triggering and xcorr (default 10.0 Hz)
        
        CLUSTERING PARAMETERS:
        cmin: Minimum correlation to be considered a repeater (default 0.7)
    
        I envision that these could eventually be set interactively or by control file.
        This list will likely expand.       
        """
        
        self.title = title
        self.filename = filename
        self.groupName = groupName
        self.groupDesc = groupDesc
        self.station = station
        self.channel = channel
        self.network = network
        self.location = location
        self.samprate = samprate
        self.lwin = lwin
        self.swin = swin
        self.trigon = trigon
        self.trigoff = trigoff
        self.mintrig = mintrig
        self.kurtmax = kurtmax
        self.kurtfmax = kurtfmax
        self.oratiomax = oratiomax
        self.kurtwin = kurtwin
        self.winlen = 512         # NOTE: These are all currently hardcoded until we
        self.ptrig = 10.0         # figure out how to not have redpy.table.Trigger and
        self.atrig = 20.0         # redpy.table.Correlation not hardcoded
        self.fhigh = fhigh
        self.fmin = fmin
        self.fmax = fmax
        self.cmin = cmin