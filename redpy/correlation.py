import numpy as np
import obspy.core.trace as trace
import redpy.table
import redpy.cluster
import datetime
import matplotlib
from scipy.fftpack import fft, ifft

def calcWindow(waveform, windowStart, opt, winlen=1):

    """
    Calculates the amplitude coefficient and FFT for a window of data.

    waveform: numpy array of waveform data
    windowStart: starting sample of window
    opt: Options object describing station/run parameters
    winlen: Fraction of window to use (optional)

    Returns windowCoeff and windowFFT
    """
    
    # Testing: shift window left by 10% of winlen
    windowStart = windowStart - opt.winlen/10
    
    windowCoeff = 1/np.sqrt(sum(waveform[windowStart:(windowStart + opt.winlen*winlen)] *
        waveform[windowStart:(windowStart + opt.winlen*winlen)]))
    windowFFT = np.reshape(fft(waveform[windowStart:(windowStart + opt.winlen*winlen)]),
        (opt.winlen*winlen,))

    return windowCoeff, windowFFT


def xcorr1x1(windowFFT1, windowFFT2, windowCoeff1, windowCoeff2):

    """
    Calculates the cross-correlation coefficient and lag for two windows.

    windowFFT1: FFT of first window
    windowFFT2: FFT of second window
    windowCoeff1: amplitude coefficient of first window
    windowCoeff2: amplitude coefficient of second window

    Order matters for sign of lag, but not CCC.

    Returns maximum cross-correlation and optimal lag (in samples)
    """

    M = len(windowFFT1)
    coeff = windowCoeff1 * windowCoeff2

    lags = np.roll(np.linspace(-M/2 + 1, M/2, M, endpoint=True),
        M/2 + 1).astype(int)
    cors = np.real(ifft(windowFFT1 * np.conj(windowFFT2))) * coeff

    indx = np.argmax(cors)
    maxcor = cors[indx]
    maxlag = lags[indx]

    return maxcor, maxlag
    

def xcorr1xtable(coeffi, ffti, subtable, opt):

    """
    Correlates a new event with all events in a subtable.
    
    coeffi: amplitude coefficient of new event
    ffti: FFT of new event
    subtable: a table of either repeaters or orphans to compare to the new event
    opt: Options object describing station/run parameters
    
    Returns correlation and lag arrays
    
    The 'subtable' can be a full table (the full orphan table) or a selection of
    rows (cluster centers from repeaters, or a full family)
    
    Contemplating figuring how how to run this in parallel...
    """
    
    cor = np.zeros((len(subtable),))
    lag = np.zeros((len(subtable),))
    
    j = -1
    for rj in subtable:
        fftj = rj['windowFFT']
        coeffj = rj['windowCoeff']
        j = j+1
        cor[j], lag[j] = xcorr1x1(ffti, fftj, coeffi, coeffj)
        
    return cor, lag


def compare2Family(rtable, ctable, rnumber, cnum, opt):

    """
    Correlates a known repeater with all events in a family except the core.
    
    rtable: Repeater table
    ctable: Correlation matrix table
    rnumber: Row of repeater in rtable
    cnum: Cluster or family number
    opt: Options object describing station/run parameters
    
    Writes correlations to ctable
    """

    famtable = rtable[rtable.get_where_list(
        '(isCore == 0) & (clusterNumber == {})'.format(cnum))]
    cor, lag = xcorr1xtable(rtable[rnumber]['windowCoeff'], rtable[rnumber]['windowFFT'],
        famtable, opt)
    for j in range(len(cor)):
        redpy.table.appendCorrelation(ctable, rtable[rnumber]['id'],
            famtable[j]['id'], cor[j], opt)


def compareDeleted(trigs, dtable, opt):

    """
    Compares trigger against deleted events
    
    trigs: Triggers to be checked
    dtable: Deleted table (manually removed from rtable)
    opt: Options object describing station/run parameters
    
    Returns trigs that do not match deleted events
    """
    
    for t in trigs:
    
        coeffi, ffti = calcWindow(t.data, int(opt.ptrig*opt.samprate), opt)  
        cor, lag = xcorr1xtable(coeffi, ffti, dtable, opt)
        
        if np.where(cor >= opt.cmin - 0.05)[0].any():
            trigs.remove(t)
    
    return trigs


def compareGoodOrphans(rtable, otable, ctable, trig, id, coeffi, ffti, cor, lag, opt):

    """
    Goes and finds the matches of the new event in the orphan table, appends them to
    the repeater table, and then compares to cores

    rtable: Repeater table
    otable: Orphan table
    ctable: Correlation matrix table
    trig: New trigger to compare
    id: Unique ID of new trigger
    coeffi: Scaling coefficient for trigger
    ffti: FFT of trigger
    cor: Correlation of trigger to orphans
    lag: Lag ot trigger to orphans
    opt: Options object describing station/run parameters
    """
    
    # Loop through potential matches
    written = 0
    while len(cor[cor >= opt.cmin - 0.05]) > 0:
        
        # If not written to rtable yet, realign new event
        if written == 0:
            lagmax = lag[np.argmax(cor)]
            coeffi2, ffti2 = calcWindow(trig.data, int(opt.ptrig*opt.samprate + lagmax),
                opt)
            coeffj2 = otable[np.argmax(cor)]['windowCoeff']
            fftj2 = otable[np.argmax(cor)]['windowFFT']
        # If written already, realign older orphan to new event
        else:
            coeffj2, fftj2 = calcWindow(otable[np.argmax(cor)]['waveform'],
                int(opt.ptrig*opt.samprate + lagmax - lag[np.argmax(cor)]), opt)
            
        cor2, lag2 = xcorr1x1(ffti2, fftj2, coeffi2, coeffj2)
        
        # If actually matches...
        if cor2 >= opt.cmin:
            # Move both the orphans to the repeater table
            if written == 0:
                redpy.table.populateRepeater(rtable, id, trig, opt,
                    id, int(opt.ptrig*opt.samprate + lagmax))
                redpy.table.moveOrphan(rtable, otable, np.argmax(cor), id, opt)
                written = 2
            # Update the table to reflect the new window, then move it
            else:
                otable.cols.windowFFT[np.argmax(cor)] = fftj2
                otable.cols.windowCoeff[np.argmax(cor)] = coeffj2
                otable.cols.windowStart[np.argmax(cor)] = int(opt.ptrig*opt.samprate +
                    lagmax - lag[np.argmax(cor)])
                redpy.table.moveOrphan(rtable, otable, np.argmax(cor), id, opt)
                written = written+1
                
        lag = np.delete(lag, np.argmax(cor))
        cor = np.delete(cor, np.argmax(cor))
    
    # If there are no actual matches in the orphans, check new event with cores
    if written == 0:
        if len(rtable) > 0:
            compareSingleOrphan2Cores(rtable, otable, ctable, trig, id, coeffi, ffti, opt)
        else:
            redpy.table.populateOrphan(otable, id, trig, opt)
    # If there is a match, check new event and its matches with cores
    else:
        compareMultipleOrphans2Cores(rtable, ctable, written, opt)


def compareMultipleOrphans2Cores(rtable, ctable, written, opt):

    """
    Compares multiple orphans that have already been written to the end of the repeater
    table to the other repeaters
    
    rtable: Repeater table
    ctable: Correlation matrix table
    written: Number of new repeaters written to rtable 
    opt: Options object describing station/run parameters
    
    Note: Currently only runs clustering if there are no matches to cores, and this
    is the ONLY case where full clustering is run
    """
    
    # Compare 'key' orphan to cores
    centers = rtable.get_where_list('isCore != 0')
    cores = rtable[centers]
    coeffi = rtable.cols.windowCoeff[-written]
    ffti = rtable.cols.windowFFT[-written]
    cor, lag = xcorr1xtable(coeffi, ffti, cores, opt)
    
    found = 0
    # Loop through families that match
    while len(cor[cor >= opt.cmin - 0.05]) > 0:    
    
        if found == 0:
            lagmax2 = lag[np.argmax(cor)]
            coeffi2, ffti2 = calcWindow(rtable[-written]['waveform'],
                int(rtable[-written]['windowStart'] + lagmax2), opt)
                        
        cor2, lag2 = xcorr1x1(ffti2, cores[np.argmax(cor)]['windowFFT'], coeffi2,
            cores[np.argmax(cor)]['windowCoeff'])   
            
        if cor2 >= opt.cmin:
            if found == 0:
                found = 1
                # Realign all new events in the repeater catalog to the matched family
                for i in range(-written,0):
                    rtable.cols.windowCoeff[i], rtable.cols.windowFFT[i] = calcWindow(
                        rtable.cols.waveform[i], int(rtable.cols.windowStart[i] +
                        lagmax2), opt)
                    rtable.cols.windowStart[i] = int(rtable.cols.windowStart[i] + lagmax2)
                    rtable.cols.clusterNumber[i] = cores[np.argmax(cor)]['clusterNumber']
                    rtable.cols.alignedTo[i] = cores[np.argmax(cor)]['id']
                    rtable.flush()
            
            # Compare to full family, write to correlation table
            for i in range(-written,0):
                cor3, lag3 = xcorr1x1(rtable[i]['windowFFT'],
                    cores[np.argmax(cor)]['windowFFT'], rtable[i]['windowCoeff'],
                    cores[np.argmax(cor)]['windowCoeff'])
                redpy.table.appendCorrelation(ctable, rtable[i]['id'],
                    cores[np.argmax(cor)]['id'], cor3, opt)
                compare2Family(rtable, ctable, i,
                    cores[np.argmax(cor)]['clusterNumber'], opt)
            
        lag = np.delete(lag, np.argmax(cor))
        cor = np.delete(cor, np.argmax(cor))
    
    # Make sure to save correlation of new events with each other
    for i in range(-written+1,0):
        cor4, lag4 = xcorr1x1(rtable[i]['windowFFT'], rtable[-written]['windowFFT'],
            rtable[i]['windowCoeff'], rtable[-written]['windowCoeff'])
        redpy.table.appendCorrelation(ctable, rtable[-written]['id'], rtable[i]['id'],
            cor4, opt)
        
    # Run clustering if events create a new family
    if found == 0:
        redpy.cluster.runFullOPTICS(rtable, ctable, opt)
        

def compareSingleOrphan2Cores(rtable, otable, ctable, trig, id, coeffi, ffti, opt):

    """
    Compares a single orphan to the cluster cores, adds the orphan to the best cluster
    if it matches, else appends to the orphan table
    
    rtable: Repeater table
    otable: Orphan table
    ctable: Correlation matrix table
    trig: New trigger to compare
    id: Unique ID of new trigger
    coeffi: Scaling coefficient for trigger
    ffti: FFT of trigger
    opt: Options object describing station/run parameters
    """
    
    centers = rtable.get_where_list('isCore != 0')
    cores = rtable[centers]
    cor, lag = xcorr1xtable(coeffi, ffti, cores, opt)
    
    written = 0
    # Loop through potential matching families
    while len(cor[cor >= opt.cmin - 0.05]) > 0:    

        if written == 0:
            lagmax = lag[np.argmax(cor)]
            coeffi2, ffti2 = calcWindow(trig.data, int(opt.ptrig*opt.samprate + lagmax),
                opt)

        cor2, lag2 = xcorr1x1(ffti2, cores[np.argmax(cor)]['windowFFT'], coeffi2,
            cores[np.argmax(cor)]['windowCoeff'])
        
        # If it matches a family...
        if cor2 >= opt.cmin:
            if written == 0:
                # Move the orphan to the repeater table
                redpy.table.populateRepeater(rtable, id, trig, opt,
                    cores[np.argmax(cor)]['id'], int(opt.ptrig*opt.samprate + lagmax))
                rtable.cols.clusterNumber[-1] = cores[np.argmax(cor)]['clusterNumber']
            
            # Correlate with other members of the family
            redpy.table.appendCorrelation(ctable, id,
                cores[np.argmax(cor)]['id'], cor2, opt)
            compare2Family(rtable, ctable, -1, cores[np.argmax(cor)]['clusterNumber'],
                opt)
            
            written = 1
            
        lag = np.delete(lag, np.argmax(cor))
        cor = np.delete(cor, np.argmax(cor))
        
    # If doesn't match anything, append as orphan   
    if written == 0:
        redpy.table.populateOrphan(otable, id, trig, opt)


def runCorrelation(rtable, otable, ctable, trig, id, opt):

    """
    Adds a new trigger to the correct table, runs the correlations and clustering
    
    rtable: Repeater table
    otable: Orphan table
    ctable: Correlation matrix table
    trig: New trigger to compare
    id: Unique ID of new trigger
    opt: Options object describing station/run parameters
    
    This is the top-level logic for processing; detailed logic is within the two compare
    functions.
    """
    
    # Check to ensure this isn't a duplicate in either rtable or otable
    try:
        stime = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S.%f'))
    except ValueError:
        stime = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S'))
    
    if not (len(otable.get_where_list('(startTimeMPL > {0}) & (startTimeMPL < {1})'.format(
        stime - opt.mintrig/86400, stime + opt.mintrig/86400))) or
        len(rtable.get_where_list('(startTimeMPL > {0}) & (startTimeMPL < {1})'.format(
        stime - opt.mintrig/86400, stime + opt.mintrig/86400)))):

        coeffi, ffti = calcWindow(trig.data, int(opt.ptrig*opt.samprate), opt)
        
        # Correlate with the new event with all the orphans
        cor, lag = xcorr1xtable(coeffi, ffti, otable, opt)
        
        try:
            # If there's a match, run the most complex function
            if max(cor) >= opt.cmin - 0.05:
                compareGoodOrphans(rtable, otable, ctable, trig, id, coeffi, ffti, cor,
                    lag, opt)
            else:
                # Compare that orphan to the cores in the repeater table
                if len(rtable) > 0:
                    compareSingleOrphan2Cores(rtable, otable, ctable, trig, id, coeffi,
                        ffti, opt)
                # Populate as an orphan if there are no repeaters yet
                else:
                    redpy.table.populateOrphan(otable, id, trig, opt)
        except ValueError:
            print('Could not properly correlate, moving on...')
            redpy.table.populateOrphan(otable, id, trig, opt)
