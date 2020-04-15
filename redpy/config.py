# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import numpy as np
try:
    # Python 3
    import configparser
except ImportError:
    # Python 2.7
    import ConfigParser as configparser

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
        title: Name of the table, used also in plotting titles (default 'REDPy Catalog')
        filename: Filename/path for the table, should end in .h5 (default 'redpytable.h5')
        outputPath: Absolute or relative path to outputs (defaults to current directory)
        groupName: Short string describing the name of the station, may not contain spaces
            (default 'default')
        
        STATION PARAMETERS:
        nsta: Number of stations (default 8)
        station: String of ordered station names
            (default 'SEP,YEL,HSR,SHW,EDM,STD,JUN,SOS')
        channel: String of channels of interest, no wildcards
            (default 'EHZ,EHZ,EHZ,EHZ,EHZ,EHZ,EHZ,EHZ')
        network: String of network code (default 'UW,UW,UW,UW,UW,UW,UW,UW')
        location: String of location code (default '--,--,--,--,--,--,--,--')
        samprate: Sampling rate of that station (default 100.0 Hz)
        server: Source of data (default "IRIS", otherwise "file" or name of waveserver)
        port: Port number for server (default 16017, not used if using IRIS)
        searchdir: Path to directory with local files ending in / (default './', not used
            if using IRIS or waveserver)
        filepattern: Wildcard for selecting subset of files based on their name
            (default "*")
        nsec: Number of seconds to download from server at a time (default 3600 s) 
        
        WINDOWING PARAMETERS:
        winlen: Length of window for cross-correlation (default 1024 samples, 2^n is best)
        ptrig: Length of time cut prior to trigger (default 10.0 s)
        atrig: Length of time cut after trigger (default 20.0 s)
        wshape: A derived value (cannot be explicitly defined) corresponding to the number
            of samples that will be cut based on ptrig and atrig
                
        TRIGGERING PARAMETERS:
        trigalg: Trigger algorithm to be used for STALTA (default 'classicstalta')
        lwin: Length of long window for STALTA (default 7.0 s)
        swin: Length of short window for STALTA (default 0.8 s)
        trigon: Cutoff ratio for triggering STALTA (default 3.0)
        trigoff: Cutoff ratio for ending STALTA trigger (default 2.0)
        mintrig: A derived value (set to 75% of winlen) for the minimum spacing between
            subsequent triggers
        nstaC: Minimum number of stations a trigger must show up on (default 4)
        offset: Optional time offset to advance waveforms as a list of positive floats
            (default 0.0)
        kurtmax: Maximum kurtosis allowed for event window, to eliminate spikes; ~80-100
            is appropriate for 5 s window, ~130 for 15 s, ~200 for 25 s (default 80.0)
        kurtfmax: Maximum kurtosis of frequency amplitude spectrum to eliminate
            calibration pulses with unnaturally harmonic signals; be careful not
            to set too low or you could eliminate real harmonic events (default 150.0)
        kurtwin: Length of window to use for kurtosis, in seconds, around the trigger
            time, will be centered on the trigger time (default 5 s)
        oratiomax: Maximum ratio of outliers to total number of datapoints in trace
            (default 0.15 (15%))
    
        FILTERING PARAMETERS:
        fmin: Lower band of bandpass filter (default 1.0 Hz)
        fmax: Upper band of bandpass filter (default 10.0 Hz)
        
        FREQUENCY INDEX WINDOWS:
        filomin: Lower bound on low window (default 2.0 Hz)
        filomax: Upper bound on low window (default 4.0 Hz)
        fiupmin: Lower bound on upper window (default 7.0 Hz)
        fiupmax: Upper bound on upper window (default 9.0 Hz)
        
        CLUSTERING PARAMETERS:
        cmin: Minimum correlation to be considered a repeater (default 0.7)
        ncor: Number of stations correlation must be exceeded on (default 4)
        
        ORPHAN EXPIRATION PARAMETERS
        minorph: Minimum amount of time (days) to keep the smaller orphans alive
            (corresponds to trigon) (default 7 days)
        maxorph: Maximum amount of time (days) to keep the largest orphans alive
            (corresponds to trigon+7) (default 30 days)
            
        PLOTTING PARAMETERS
        minplot: Minimum number of members required in order to be plotted to timeline
        dybin: Width of bin in days for full histogram (default 1 day)
        hrbin: Width of bin in hours for recent histogram (default 1 hour)
        occurbin: Width of bin for occurrence plot; specified in .cfg as hours, converted to days in redpy/config (default 1 hr -> 1/24 day)
        recbin: Width of bin for recent occurrence; specified in .cfg as hours, converted to days in redpy/config (default 1 hr -> 1/24 day) 
        recplot: Number of days for 'recent' plot (default 14 days)
        plotsta: Station index in station list to be plotted (default 2)
        verbosecatalog: Add additional columns to the catalog file (default False)
        amplims: Use 'global' or 'family' to define amplitude plot limits (default global)
        
        COMCAT PARAMETERS
        checkComCat: Use ComCat to find located seismicity that might match repeaters
            (default False)
        stalats: List of station latitudes (defaults to MSH network:
            '46.200210,46.209550,46.174280,46.193470,46.197170,46.237610,46.147060,
            46.243860')
        stalons: List of station longitudes (defaults to MSH network:
            '-122.190600,-122.188990,-122.180650,-122.236350,-122.151210,-122.223960,
            -122.152430,-122.137870') 
        serr: Seconds of allowable difference in trigger and projected arrival time
            (default 5.0 s)
        locdeg: Degrees of distance to be considered a local event (default 0.5 degrees)
        regdeg: Degrees of distance to be considered a regional event (default 2.0 degrees)
        regmag: Minimum magnitude for regional events (default M2.5)
        telemag: Minimum magnitude for teleseismic events (default M4.5)
        matchMax: Number of largest events to match (default 0 (all))
          
        This list will likely expand.       
        """
        
        self.configfile = configfile
                
        # Load parameters from config file
        config = configparser.ConfigParser()
        config.read(self.configfile)
        
        # Set parameters to default if not in config file       
        self.title=config.get('Settings','title') if config.has_option(
            'Settings','title') else 'REDPy Catalog'
        self.filename=config.get('Settings','filename') if config.has_option(
            'Settings','filename') else 'redpytable.h5'
        self.outputPath=config.get('Settings','outputPath') if config.has_option(
            'Settings','outputPath') else ''
        self.groupName=config.get('Settings','groupName') if config.has_option(
            'Settings','groupName') else 'default'
        self.groupDesc=config.get('Settings','groupDesc') if config.has_option(
            'Settings','groupDesc') else 'Default Test Run'
        self.nsta=config.getint('Settings','nsta') if config.has_option(
            'Settings','nsta') else 8 
        self.station=config.get('Settings','station') if config.has_option(
            'Settings','station') else 'SEP,YEL,HSR,SHW,EDM,STD,JUN,SOS'
        self.channel=config.get('Settings','channel') if config.has_option(
            'Settings','channel') else 'EHZ,EHZ,EHZ,EHZ,EHZ,EHZ,EHZ,EHZ'
        self.network=config.get('Settings','network') if config.has_option(
            'Settings','network') else 'UW,UW,UW,UW,UW,UW,UW,UW'
        self.location=config.get('Settings','location') if config.has_option(
            'Settings','location') else '--,--,--,--,--,--,--,--'
        self.samprate=config.getfloat('Settings','samprate') if config.has_option(
            'Settings','samprate') else 100.
        self.nstaC=config.getint('Settings','nstaC') if config.has_option(
            'Settings','nstaC') else 5
        self.printsta=config.getint('Settings','printsta') if config.has_option(
            'Settings','printsta') else 2
        self.server=config.get('Settings','server') if config.has_option(
            'Settings','server') else 'IRIS'
        self.port=config.getint('Settings','port') if config.has_option(
            'Settings','port') else 16017
        self.searchdir=config.get('Settings','searchdir') if config.has_option(
            'Settings','searchdir') else './'
        self.filepattern=config.get('Settings','filepattern') if config.has_option(
            'Settings','filepattern') else '*'
        self.nsec=config.getint('Settings','nsec') if config.has_option(
            'Settings','nsec') else 3600
        self.trigalg=config.get('Settings','trigalg') if config.has_option(
            'Settings','trigalg') else 'classicstalta'
        self.lwin=config.getfloat('Settings','lwin') if config.has_option(
            'Settings','lwin') else 7.
        self.swin=config.getfloat('Settings','swin') if config.has_option(
            'Settings','swin') else 0.8
        self.trigon=config.getfloat('Settings','trigon') if config.has_option(
            'Settings','trigon') else 3.
        self.trigoff=config.getfloat('Settings','trigoff') if config.has_option(
            'Settings','trigoff') else 2.
        self.offset=config.get('Settings','offset') if config.has_option(
            'Settings','offset') else '0'
        self.kurtmax=config.getfloat('Settings','kurtmax') if config.has_option(
            'Settings','kurtmax') else 80.
        self.kurtfmax=config.getfloat('Settings','kurtfmax') if config.has_option(
            'Settings','kurtfmax') else 150.
        self.oratiomax=config.getfloat('Settings','oratiomax') if config.has_option(
            'Settings','oratiomax') else 0.15
        self.kurtwin=config.getfloat('Settings','kurtwin') if config.has_option(
            'Settings','kurtwin') else 5.
        self.winlen=config.getint('Settings','winlen') if config.has_option(
            'Settings','winlen') else 1024
        self.fmin=config.getfloat('Settings','fmin') if config.has_option(
            'Settings','fmin') else 1.
        self.fmax=config.getfloat('Settings','fmax') if config.has_option(
            'Settings','fmax') else 10.
        self.filomin=config.getfloat('Settings','filomin') if config.has_option(
            'Settings','filomin') else 1.
        self.filomax=config.getfloat('Settings','filomax') if config.has_option(
            'Settings','filomax') else 2.5
        self.fiupmin=config.getfloat('Settings','fiupmin') if config.has_option(
            'Settings','fiupmin') else 5.
        self.fiupmax=config.getfloat('Settings','fiupmax') if config.has_option(
            'Settings','fiupmax') else 10.
        self.telefi=config.getfloat('Settings','telefi') if config.has_option(
            'Settings','telefi') else -1.
        self.teleok=config.getint('Settings','teleok') if config.has_option(
            'Settings','teleok') else 1    
        self.cmin=config.getfloat('Settings','cmin') if config.has_option(
            'Settings','cmin') else 0.7
        self.ncor=config.getint('Settings','ncor') if config.has_option(
            'Settings','ncor') else 4
        self.minorph=config.getfloat('Settings','minorph') if config.has_option(
            'Settings','minorph') else 0.05
        self.maxorph=config.getfloat('Settings','maxorph') if config.has_option(
            'Settings','maxorph') else 7.
        self.minplot=config.getint('Settings','minplot') if config.has_option(
            'Settings','minplot') else 3
        self.dybin=config.getfloat('Settings','dybin') if config.has_option(
            'Settings','dybin') else 1.
        self.hrbin=config.getfloat('Settings','hrbin') if config.has_option(
            'Settings','hrbin') else 1.
        self.occurbin=config.getfloat('Settings','occurbin')/24 if config.has_option(   # settings.cfg (hours) immediately converted to days
            'Settings','occurbin') else 1/24
        self.recbin=config.getfloat('Settings','recbin')/24 if config.has_option( # settings.cfg (hours) immediately converted to days
            'Settings','recbin') else 1/24
        self.recplot=config.getfloat('Settings','recplot') if config.has_option(
            'Settings','recplot') else 14.
        self.printVerboseCat=config.getboolean('Settings','verbosecatalog') if config.has_option(
            'Settings','verbosecatalog') else False
        self.amplims=config.get('Settings','amplims') if config.has_option(
            'Settings','amplims') else 'global'
        self.anotfile=config.get('Settings','anotfile') if config.has_option(
            'Settings','anotfile') else ''
        self.checkComCat=config.getboolean('Settings','checkComCat') if config.has_option(
            'Settings','checkComCat') else False
        self.matchMax=config.getint('Settings','matchMax') if config.has_option(
            'Settings','matchMax') else 0
        self.stalats=config.get('Settings','stalats') if config.has_option(
            'Settings','stalats') else ('46.200210,46.209550,46.174280,46.193470,'
                '46.197170,46.237610,46.147060,46.243860')
        self.stalons=config.get('Settings','stalons') if config.has_option(
            'Settings','stalons') else ('-122.190600,-122.188990,-122.180650,-122.236350,'
                '-122.151210,-122.223960,-122.152430,-122.137870')
        self.serr=config.getfloat('Settings','serr') if config.has_option(
            'Settings','serr') else 5.
        self.locdeg=config.getfloat('Settings','locdeg') if config.has_option(
            'Settings','locdeg') else 0.5
        self.regdeg=config.getfloat('Settings','regdeg') if config.has_option(
            'Settings','regdeg') else 2.
        self.regmag=config.getfloat('Settings','regmag') if config.has_option(
            'Settings','regmag') else 2.5
        self.telemag=config.getfloat('Settings','telemag') if config.has_option(
            'Settings','telemag') else 4.5
        
        # Derived Settings
        self.ptrig=1.5*self.winlen/self.samprate
        self.atrig=3*self.winlen/self.samprate
        self.mintrig=0.75*self.winlen/self.samprate
        self.wshape = int((self.ptrig + self.atrig)*self.samprate) + 1
        self.maxdt = np.max(np.fromstring(self.offset, sep=','))
        
