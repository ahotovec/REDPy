import ConfigParser

class Options(object):
    
    def __init__(self, configfile='settings.cfg'):
        
        """
        Defines the settings that are often passed to routines and that define the table.
        These are also written to the attributes of the table for posterity.
        
        Requires a configuration file with section header [Settings] on the first line and
        any of the following configurations below it. Passing a configuration file with
        only the header defaults to a test run laid out with the settings below. Format of
        the file below the header is simply:
        
        name=value
                
        where name is the name of the parameter, and value is either a string (no quotes)
        or number. Comments are allowed on separate lines beginning with a #, and the
        parameters may be in any order desired. An example configuration file called
        'settings.cfg' is included in the distribution that contains all of the default
        settings and may be edited. The name of the configfile used is also stored.       
    
        TABLE DEFINITIONS:
        title: Name of the table (default 'REDPy Catalog')
        filename: Filename for the table, should end in .h5 (default 'redpytable.h5')
        groupName: Short string describing the name of the station, may not contain spaces
            (default 'default')
        groupDesc: Longer string describing the run (default 'Default Test Run')
    
        STATION PARAMETERS:
        station: String of station name (default 'HSR')
        channel: String of channel of interest, no wildcards supported yet (default 'EHZ')
        network: String of network code (default 'UW')
        location: String of location code (default '--')
        samprate: Sampling rate of that station (default 100.0 Hz)
        server: Source of data (default "IRIS", otherwise name of waveserver)
        port: Port number for server (default 16017, not used if using IRIS)
        nsec: Number of seconds to download from server at a time (default 3600 s) 
                
        TRIGGERING PARAMETERS:
        lwin: Length of long window for STALTA (default 7.0 s)
        swin: Length of short window for STALTA (default 0.8 s)
        trigon: Cutoff ratio for triggering STALTA (default 3.0)
        trigoff: Cutoff ratio for ending STALTA trigger (default 2.0)
        mintrig: Minimum spacing between triggers (default 10.0 s)
        kurtmax: Maximum kurtosis allowed for event window, to eliminate spikes, ~80-100
            is appropriate for 5 s window, ~130 for 15 s, ~200 for 25 s (default 80.0)
        kurtfmax: Maximum kurtosis of frequency amplitude spectrum to eliminate
            calibration pulses with unnaturally harmonic signals; be careful not
            to set too low or you could eliminate real harmonic events (default 150.0)
        kurtwin: Length of window to use for kurtosis, in seconds, around the trigger
            time, will be centered on the trigger time (default 5 s)
        oratiomax: Maximum ratio of outliers to total number of datapoints in trace
            (default 0.06 (6%))
    
        WINDOWING PARAMETERS:
        winlen: Length of window for cross-correlation (default 512 samples, 2^n is best)
        ptrig: Length of time cut prior to trigger (default 10.0 s)
        atrig: Length of time cut after trigger (default 20.0 s)
        wshape: A derived value (cannot be explicitly defined) corresponding to the number
            of samples that will be cut based on ptrig and atrig
    
        FILTERING PARAMETERS:
        fhigh: Highpass filter to reduce microseism in "raw" waveforms (default 0.25 Hz)
        fmin: Lower band of bandpass filter for triggering and xcorr (default 1.0 Hz)
        fmax: Upper band of bandpass filter for triggering and xcorr (default 10.0 Hz)
        
        CLUSTERING PARAMETERS:
        cmin: Minimum correlation to be considered a repeater (default 0.7)
        
        ORPHAN EXPIRATION PARAMETERS
        minorph: minimum amount of time (days) to keep the smaller orphans alive
            (corresponds to trigon) (default 7 days)
        maxorph: maximum amount of time (days) to keep the largest orphans alive
            (corresponds to trigon+7) (default 30 days)
    
        This list will likely expand.       
        """
        
        self.configfile = configfile
                
        # Load parameters from config file
        config = ConfigParser.RawConfigParser()
        config.read(self.configfile)
        
        # Set parameters to default if not in config file       
        self.title=config.get('Settings','title') if config.has_option(
            'Settings','title') else 'REDPy Catalog'
        self.filename=config.get('Settings','filename') if config.has_option(
            'Settings','filename') else 'redpytable.h5'
        self.groupName=config.get('Settings','groupName') if config.has_option(
            'Settings','groupName') else 'default'
        self.groupDesc=config.get('Settings','groupDesc') if config.has_option(
            'Settings','groupDesc') else 'Default Test Run'
        self.station=config.get('Settings','station') if config.has_option(
            'Settings','station') else 'HSR'
        self.channel=config.get('Settings','channel') if config.has_option(
            'Settings','channel') else 'EHZ'
        self.network=config.get('Settings','network') if config.has_option(
            'Settings','network') else 'UW'
        self.location=config.get('Settings','location') if config.has_option(
            'Settings','location') else '--'
        self.samprate=config.getfloat('Settings','samprate') if config.has_option(
            'Settings','samprate') else 100.
        self.server=config.get('Settings','server') if config.has_option(
            'Settings','server') else 'IRIS'
        self.port=config.getint('Settings','port') if config.has_option(
            'Settings','port') else 16017
        self.nsec=config.getint('Settings','nsec') if config.has_option(
            'Settings','nsec') else 3600
        self.lwin=config.getfloat('Settings','lwin') if config.has_option(
            'Settings','lwin') else 7.
        self.swin=config.getfloat('Settings','swin') if config.has_option(
            'Settings','swin') else 0.8
        self.trigon=config.getfloat('Settings','trigon') if config.has_option(
            'Settings','trigon') else 3.
        self.trigoff=config.getfloat('Settings','trigoff') if config.has_option(
            'Settings','trigoff') else 2.
        self.mintrig=config.getfloat('Settings','mintrig') if config.has_option(
            'Settings','mintrig') else 10.
        self.kurtmax=config.getfloat('Settings','kurtmax') if config.has_option(
            'Settings','kurtmax') else 80.
        self.kurtfmax=config.getfloat('Settings','kurtfmax') if config.has_option(
            'Settings','kurtfmax') else 150.
        self.oratiomax=config.getfloat('Settings','oratiomax') if config.has_option(
            'Settings','oratiomax') else 0.06
        self.kurtwin=config.getfloat('Settings','kurtwin') if config.has_option(
            'Settings','kurtwin') else 5.
        self.winlen=config.getint('Settings','winlen') if config.has_option(
            'Settings','winlen') else 512
        self.ptrig=config.getfloat('Settings','ptrig') if config.has_option(
            'Settings','ptrig') else 10.
        self.atrig=config.getfloat('Settings','atrig') if config.has_option(
            'Settings','atrig') else 20.
        self.wshape = int((self.ptrig + self.atrig)*self.samprate) + 1
        self.fhigh=config.getfloat('Settings','fhigh') if config.has_option(
            'Settings','fhigh') else 0.25
        self.fmin=config.getfloat('Settings','fmin') if config.has_option(
            'Settings','fmin') else 1.
        self.fmax=config.getfloat('Settings','fmax') if config.has_option(
            'Settings','fmax') else 10.
        self.cmin=config.getfloat('Settings','cmin') if config.has_option(
            'Settings','cmin') else 0.7
        self.minorph=config.getfloat('Settings','minorph') if config.has_option(
            'Settings','minorph') else 7.
        self.maxorph=config.getfloat('Settings','maxorph') if config.has_option(
            'Settings','maxorph') else 30.