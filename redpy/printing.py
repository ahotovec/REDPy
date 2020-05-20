# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import numpy as np
import matplotlib.dates
from obspy import UTCDateTime

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
	since last event in hours, correlation coefficient with respect to the best
	correlated event, and correlation coefficient with respect to the core event.
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
	
		f.write("cnum\tevTime                    \tfi\txcormax\txcorcore\tdt(hr)\tamps\n")
		for cnum in range(ftable.attrs.nClust):
			fam = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
		
			catalogind = np.argsort(startTimeMPL[fam])
			catalog = startTimeMPL[fam][catalogind]
			spacing = np.diff(catalog)*24
			idf = ids[fam]
			ix = np.where(np.in1d(id2,idf))
			C = np.eye(len(idf))
			r1 = [np.where(idf==xx)[0][0] for xx in id1[ix]]
			r2 = [np.where(idf==xx)[0][0] for xx in id2[ix]]
			C[r1,r2] = ccc[ix]
			C[r2,r1] = ccc[ix]
			xcorrmax = C[np.argmax(np.sum(C,0)),:]
			core = ftable[cnum]['core']
			xcorrcore = C[np.where(fam==core)[0][0],:]
		
			j = -1
			for i in catalogind:
				evTime = (UTCDateTime(startTimes[fam][i]) +
					windowStarts[fam][i]/opt.samprate)
				amp = windowAmps[fam[i],:]
				if j == -1:
					dt = np.nan
				else:
					dt = spacing[j]
				j += 1
			
				f.write("{0}\t{1}\t{2: 4.3f}\t{4:3.2f}\t{5:3.2f}\t{3:12.6f}\t[".format(
				    cnum,evTime.isoformat(),fi[fam][i],dt,xcorrmax[i],xcorrcore[i]))
				for a in amp:
					f.write(" {:10.2f} ".format(a))
				f.write("]\n")


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