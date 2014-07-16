# -*- coding: utf-8 -*-

# REDPy: Repeating Earthquake Detector (Python)
# Codes for online detection and clustering of repeating earthquakes
#
# Written by Alicia Hotovec-Ellis (ahotovec@uw.edu)
# Based on codes by Kate Allstadt, Josh Carmichael, and others
#
# Requires ObsPy for waveform handling (https://github.com/obspy/obspy)
# Requires PyTables for database handling (http://www.pytables.org)

import sys
import scipy
import time
import numpy as np
import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft
from obspy.fdsn import Client
from obspy import UTCDateTime
from obspy.signal.trigger import classicSTALTA, triggerOnset
from obspy.core.trace import Trace
from tables import *

class Triggers(IsDescription):

    """
    Defines the columns in the "Repeater Catalog" table

    id: unique ID number for the event (integer)
    startTime: UTC time of start of the waveform (string)
    waveform: Waveform data (ndarray)
    windowStart: "trigger" time, in samples from start (integer)
    windowCoeff: amplitude scaling for cross-correlation (float)
    windowFFT: Fourier transform of window (complex ndarray)

    Needs work to figure out how to adjust the shape of the waveform and
    windowFFT columns when the window length and padding around the triggers
    are not the same from station to station
    """

    id = Int32Col(shape=(), pos=0)
    startTime = StringCol(itemsize=32, pos=1)
    waveform = Float64Col(shape=(3001,), pos=2)
    windowStart = Int32Col(shape=(), pos=3)
    windowCoeff = Float64Col(shape=(), pos=4)
    windowFFT = ComplexCol(shape=(512,), itemsize=16, pos=5)


class Correlation(IsDescription):

    """
    Defines the columns in the "Correlation" table

    id1: unique ID number for the first event (integer)
    id2: unique ID number for the second event (integer)
    ccc: cross-correlation coefficient between those two events (float)    
    """
    
    id1 = Int32Col(shape=(), pos=0)
    id2 = Int32Col(shape=(), pos=1)
    ccc = Float64Col(shape=(), pos=2)


def initializeTable(groupName, groupDesc, scnl, samprate=100.0, winlen=512,
    ptrig=10.0, atrig=20.0, fmin=1.0, fmax=10.0, title="REDPy Catalog",
    filename="redtable.h5"):

    """
    Initializes the hdf5 file with a "Repeater Catalog" table in a group
    related to the station where the data come from.

    groupName: Short string describing the name of the station, e.g., "hsr"
    groupDesc: Longer string describing the station, "MSH: HSR-EHZ-UW"
    scnl: List fully describing station, ["HSR", "EHZ", "UW", "--"]
    samprate: Sampling rate of that station (default 100.0 Hz)
    winlen: Length of window for cross-correlation (default 512 samples)
    ptrig: Length of time cut prior to trigger (default 10.0 s)
    atrig: Length of time cut after trigger (default 20.0 s)
    fmin: Lower band for bandpass filter (default 1.0 Hz)
    fmax: Upper band for bandpass filter (default 10.0 Hz)
    title: Name of the table (default "REDPy Catalog")
    filename: Filename for the table (default "redtable.h5")

    Saves table to file and closes it.
    Will need extensive editing when more tables get added...
    """

    h5file = open_file(filename, mode="w", title=title)
    group = h5file.create_group("/", groupName, groupDesc)

    rtable = h5file.create_table(group, "repeaters", Triggers,
        "Repeater Catalog")
    rtable.attrs.scnl = scnl
    rtable.attrs.samprate = samprate
    rtable.attrs.windowLength = winlen
    rtable.attrs.fmin = fmin
    rtable.attrs.fmax = fmax
    rtable.flush()

    ctable = h5file.create_table(group, "correlation", Correlation,
        "Correlation Matrix")
    ctable.attrs.order = 0
    ctable.attrs.reachability = 0
    ctable.attrs.coredist = 0
    ctable.flush()

    h5file.close()


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

    st = st.detrend()
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


def calcWindow(waveform, windowStart, winlen=512):

    """
    Calculates the amplitude coefficient and FFT for a window of data.

    waveform: numpy array of waveform data
    windowStart: starting sample of window
    winlen: window length (default: 512 samples)

    Returns windowCoeff and windowFFT
    """

    windowCoeff = 1/np.sqrt(sum(waveform[windowStart:(windowStart + winlen)] *
        waveform[windowStart:(windowStart + winlen)]))
    windowFFT = np.reshape(fft(waveform[windowStart:(windowStart + winlen)]),
        (winlen,))

    return windowCoeff, windowFFT


def populateTrigger(trigger, id, trig, windowStart):

    """
    Initially populates the trigger row in the 'Repeater Catalog' table.
    
    trigger: object pointing to the row in the table to populate
        (e.g., h5file.root.hsr.repeaters.row)
    id: integer id number given to this trigger, should be unique
    trig: ObsPy trace from triggering function
    wstart: starting sample of window

    Appends this row to table
    """
    
    trigger['id'] = id
    trigger['startTime'] = trig.stats.starttime.isoformat()
    trigger['waveform'] = trig.data
    trigger['windowStart'] = windowStart
    trigger['windowCoeff'], trigger['windowFFT'] = calcWindow(trig.data,
        windowStart)
    trigger.append()


def appendCorrelation(corr, id1, id2, ccc):

    """
    Appends a new value to the 'Correlation Matrix' table.

    corr: object pointing to the row in the table to populate
        (e.g., h5file.root.hsr.correlation.row)
    id1: unique id number of first trigger
    id2: unique id number of second trigger
    ccc: cross-correlation between the two triggers in the window

    Appends this row to the table, and automatically puts the smaller of
    the two id numbers first
    """

    corr['id1'] = min(id1, id2)
    corr['id2'] = max(id1, id2)
    corr['ccc'] = ccc
    corr.append()


def getCell(table, id, column):
    
    """
    Shorthand way of getting data from the PyTable.

    table: PyTable you're querying
    id: unique id of row you want
    column: column you want (e.g., 'windowFFT' OR position as integer)

    Returns data inside that cell

    While simple, can be time consuming if called a lot!
    """
   
    c = '(id == {})'.format(id)
    t = table.where(c)
    for r in t: data = r[column]

    return data


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


def getClusters(ctable, cutoff=0.7):

    """
    Cuts the clustering order into clusters, defines orphans as -1

    ctable: Correlation table, with clustering order in attributes
    cutoff: Minimum coefficient to cut the clusters

    Returns cluster numbers wrt. ordered list.
    """

    order = np.array(ctable.attrs.order)
    oreach = ctable.attrs.reachability[order]
    odist = ctable.attrs.coredist[order]
    cluster_id = -1

    oclust = np.zeros((len(oreach),))
    for x in range(len(oreach)):
        if oreach[x] > 1 - cutoff:
            if odist[x] <= 1 - cutoff:
                cluster_id += 1
                oclust[x] = cluster_id
            else:
                oclust[x] = -1 # orphan
        else:
            oclust[x] = cluster_id

    return oclust    


def getCenters(ctable, cutoff=0.7):

    """
    Finds the "center" of each cluster (including orphans)
    
    ctable: Correlation table, with clustering order in attributes
    cutoff: Minimum coefficient to cut the clusters

    Returns the id numbers of cluster centers and orphans
    """

    order = np.array(ctable.attrs.order)
    oreach = ctable.attrs.reachability[order]
    oclust = getClusters(ctable, cutoff)

    cluster_id = np.max(oclust).astype(int)
    o = np.array(order)
    centers = np.zeros((cluster_id + 1,)).astype(int)
    for clusternum in range(cluster_id + 1):
        oo = o[oclust == clusternum]
        centers[clusternum] = oo[np.argmin(oreach[oclust == clusternum])]

    orphans = o[oclust == -1]

    return centers, orphans



# These two need to be more modular and speedy! Currently VERY SLOW
# Contemplating making a class of it...
def xcorr_all_fft(trigs):
    """
    Cross-correlate all triggers against each other.
    Runs fastest when number of samples in triggers is a power of 2.

    trigs: OBSPy traces of triggers to cross-correlate

    Returns correlation and lag matrices and fft of traces 
    """

    samprate = trigs[0].stats.sampling_rate
    N, M = np.shape(trigs)
    data = np.zeros((N, M))
    for n in range(N):
        data[n, :] = trigs[n].data

    # Rotate so events are in columns
    data = data.T

    # Prep some other necessary terms
    lags = np.roll(np.linspace(
        (-M/2+1), (M/2), M, endpoint=True), M/2+1)/samprate
    wcoeff = 1/np.sqrt(sum(data*data));
    C = np.eye(N)
    L = np.zeros((N, N))

    # Get FFT of traces
    X = fft(data, axis=0)
    Xc = np.conjugate(X)

    # Loop through rows of similarity matrix
    for n in range(N):
        cols = range(n, N)
        CC = np.multiply(np.tile(X[:, n],
                         (len(cols),1)).T, Xc[:, cols])
        corr = np.real(ifft(CC, axis=0))
        maxval = np.amax(corr, axis=0)
        indx = np.argmax(corr, axis=0)
        C[n, cols] = np.multiply(maxval, wcoeff[cols])*wcoeff[n]
        L[n, cols] = lags[indx]

    # Fill in the lower triangular part of matrix
    C = C + C.T - np.eye(N)
    L = L - L.T

    return C, L, X, wcoeff


def xcorr_addrow_fft(trigs, newtrig, C, L, X, wcoeff):
    """
    Appends rows to the end of an existing corrleation/lag matrix
    CURRENTLY ONLY WORKS FOR SINGLE NEW TRIGGER AT A TIME

    trigs: OBSPy triggers that have already been correlated
    newtrig: OBSPy trigger that needs to be correlated
    C: Filled correlation matrix corresponding to 'trigs'
    L: Filled lag matrix corresponding to 'trigs'
    X: FFT of trigs

    Returns: new C, new L, new fft
    """

    samprate = trigs[0].stats.sampling_rate
    N, M = np.shape(trigs)
    dnew = newtrig.data
    wcoeffnew = 1/np.sqrt(sum(dnew*dnew))
    Xnew = np.reshape(fft(dnew), (M, 1))
    trigs += newtrig
    Nnew = np.shape(trigs)[0]
    data = np.zeros((Nnew, M))
    for n in range(Nnew):
        data[n, :] = trigs[n].data

    # Rotate so events are in columns
    data = data.T

    # Prep some other necessary terms
    lags = np.roll(np.linspace(
        (-M/2+1), (M/2), M, endpoint=True), M/2+1)/samprate
    wcoeff = np.append(wcoeff, wcoeffnew)
    Cnew = np.eye(Nnew)
    Lnew = np.zeros((Nnew, Nnew))

    # Get FFT of traces
    X = np.hstack((X, Xnew))
    Xc = np.conjugate(X)

    # Loop through rows of similarity matrix
    for n in range(N, Nnew):
        cols = range(0, n)
        CC = np.multiply(np.tile(X[:, n],
                         (len(cols),1)).T, Xc[:, cols])
        corr = np.real(ifft(CC, axis=0))
        maxval = np.amax(corr, axis=0)
        indx = np.argmax(corr, axis=0)
        Cnew[n, cols] = np.multiply(maxval, wcoeff[cols])*wcoeff[n]
        Lnew[n, cols] = lags[indx]

    # Fill in the lower triangular part of matrix
    Cnew = Cnew + Cnew.T - np.eye(Nnew)
    Lnew = Lnew - Lnew.T

    Cnew[0:N, 0:N] = C
    Lnew[0:N, 0:N] = L

    return Cnew, Lnew, X, wcoeff


# Old version, keeping temporarily for some older codes
# Will be removed soon(ish)
def xcorr_all(trigs):
    """
    Cross-correlate all triggers against each other.
    Runs fastest when number of samples in triggers is a power of 2.

    trigs: OBSPy traces of triggers to cross-correlate

    Returns correlation and lag matrices
    """

    samprate = trigs[0].stats.sampling_rate
    N, M = np.shape(trigs)
    data = np.zeros((N, M))
    for n in range(N):
        data[n, :] = trigs[n].data

    # Rotate so events are in columns
    data = data.T

    # Prep some other necessary terms
    lags = np.roll(np.linspace(
        (-M/2+1), (M/2), M, endpoint=True), M/2+1)/samprate
    wcoeff = 1/np.sqrt(sum(data*data));
    C = np.eye(N)
    L = np.zeros((N, N))

    # Get FFT of traces
    X = fft(data, axis=0)
    Xc = np.conjugate(X)

    # Loop through rows of similarity matrix
    for n in range(N):
        cols = range(n, N)
        CC = np.multiply(np.tile(X[:, n],
                         (len(cols),1)).T, Xc[:, cols])
        corr = np.real(ifft(CC, axis=0))
        maxval = np.amax(corr, axis=0)
        indx = np.argmax(corr, axis=0)
        C[n, cols] = np.multiply(maxval, wcoeff[cols])*wcoeff[n]
        L[n, cols] = lags[indx]

    # Fill in the lower triangular part of matrix
    C = C + C.T - np.eye(N)
    L = L - L.T

    return C, L


def xcorr_addrows(trigs, newtrigs, C, L):
    """
    Appends rows to the end of an existing corrleation/lag matrix

    trigs: OBSPy triggers that have already been correlated
    newtrigs: OBSPy triggers that need to be correlated
    C: Filled correlation matrix corresponding to 'trigs'
    L: Filled lag matrix corresponding to 'trigs'

    Returns: new C, and new L
    """

    samprate = trigs[0].stats.sampling_rate
    N, M = np.shape(trigs)
    trigs += newtrigs
    Nnew = np.shape(trigs)[0]
    data = np.zeros((Nnew, M))
    for n in range(Nnew):
        data[n, :] = trigs[n].data

    # Rotate so events are in columns
    data = data.T

    # Prep some other necessary terms
    lags = np.roll(np.linspace(
        (-M/2+1), (M/2), M, endpoint=True), M/2+1)/samprate
    wcoeff = 1/np.sqrt(sum(data*data));
    Cnew = np.eye(Nnew)
    Lnew = np.zeros((Nnew, Nnew))

    # Get FFT of traces
    X = fft(data, axis=0)
    Xc = np.conjugate(X)

    # Loop through rows of similarity matrix
    for n in range(N, Nnew):
        cols = range(0, n)
        CC = np.multiply(np.tile(X[:, n],
                         (len(cols),1)).T, Xc[:, cols])
        corr = np.real(ifft(CC, axis=0))
        maxval = np.amax(corr, axis=0)
        indx = np.argmax(corr, axis=0)
        Cnew[n, cols] = np.multiply(maxval, wcoeff[cols])*wcoeff[n]
        Lnew[n, cols] = lags[indx]

    # Fill in the lower triangular part of matrix
    Cnew = Cnew + Cnew.T - np.eye(Nnew)
    Lnew = Lnew - Lnew.T

    Cnew[0:N, 0:N] = C
    Lnew[0:N, 0:N] = L

    return Cnew, Lnew





# CODES FOR OPTICS CLUSTERING

# These all need to be cleaned up something fierce
# I have at least removed the inheritance of the BallTree structure
class setOfObjects(object):

    """Build data structure with processing index from given data
in preparation for OPTICS Algorithm

Parameters
----------
distance_pairs: array [n_samples, n_samples]"""

    def __init__(self, distance_pairs):

        self.data = distance_pairs
        self._n = len(self.data)
        # Start all points as 'unprocessed' ##
        self._processed = scipy.zeros((self._n, 1), dtype=bool)
        self._reachability = scipy.ones(self._n) * scipy.inf
        self._core_dist = scipy.ones(self._n) * scipy.nan
        # Might be faster to use a list below? ##
        self._index = scipy.array(range(self._n))
        self._nneighbors = scipy.ones(self._n, dtype=int)*self._n
        # Start all points as noise ##
        self._cluster_id = -scipy.ones(self._n, dtype=int)
        self._is_core = scipy.ones(self._n, dtype=bool)
        # Ordering is important below... ###
        self._ordered_list = []

# Paralizeable! #


def prep_optics(SetofObjects, epsilon):
    """Prep data set for main OPTICS loop

Parameters
----------
SetofObjects: Instantiated instance of 'setOfObjects' class
epsilon: float or int
Determines maximum object size that can be extracted.
Smaller epsilons reduce run time

Returns
-------
Modified setOfObjects tree structure"""

    for j in SetofObjects._index:
        # Find smallest nonzero distance
        SetofObjects._core_dist[j] = np.sort(SetofObjects.data[j,:])[1]
    print(
        'Core distances and neighborhoods prepped for ' + str(
        SetofObjects._n) + ' points.')

# Main OPTICS loop #


def build_optics(SetOfObjects, epsilon):
    """Builds OPTICS ordered list of clustering structure

Parameters
---------
SetofObjects: Instantiated and prepped instance of 'setOfObjects' class
epsilon: float or int
Determines maximum object size that can be extracted. Smaller
epsilons reduce run time. This should be equal to epsilon
in 'prep_optics'
Output_file_name: string
Valid path where write access is available.
Stores cluster structure"""

    for point in SetOfObjects._index:
        if not SetOfObjects._processed[point]:
            expandClusterOrder(SetOfObjects, point, epsilon)

# OPTICS helper functions; these should not be public #

# NOT Paralizeable! The order that entries are written to
# the '_ordered_list' is important!

def expandClusterOrder(SetOfObjects, point, epsilon):
    if SetOfObjects._core_dist[point] <= epsilon:
        while not SetOfObjects._processed[point]:
            SetOfObjects._processed[point] = True
            SetOfObjects._ordered_list.append(point)
#            # Comment following two lines to not write to a text file ##
#            with open(Output_file_name, 'a') as file:
#                file.write((str(point) + ', ' + str(
#                    SetOfObjects._reachability[point]) + '\n'))
                # Keep following line! ##
            point = set_reach_dist(SetOfObjects, point, epsilon)
        print('Object Found!')
    else:
        SetOfObjects._processed[point] = True # Probably not needed... #


# As above, NOT paralizable! Paralizing would allow items in
# 'unprocessed' list to switch to 'processed' ###
def set_reach_dist(SetOfObjects, point_index, epsilon):

    # Assumes that the query returns ordered (smallest distance first)
    # entries. This is the case for the balltree query...

#    distances, indices = SetOfObjects.query(SetOfObjects.data[point_index],
#                                            SetOfObjects._nneighbors[point_index])

    row = [SetOfObjects.data[point_index,:]]
    indices = np.argsort(row)
    distances = np.sort(row)

    # Checks to see if there more than one member in the neighborhood ##
    if scipy.iterable(distances):

        # Masking processed values ##
        unprocessed = indices[(SetOfObjects._processed[indices] < 1)[0].T]
        rdistances = scipy.maximum(
            distances[(SetOfObjects._processed[indices] < 1)[0].T],
            SetOfObjects._core_dist[point_index])
        SetOfObjects._reachability[
            unprocessed] = scipy.minimum(
                SetOfObjects._reachability[
                    unprocessed],
                rdistances)

        # Checks to see if everything is already processed;
        # if so, return control to main loop ##
        if unprocessed.size > 0:
            # Define return order based on reachability distance ###
            return sorted(zip(SetOfObjects._reachability[unprocessed], unprocessed), key=lambda reachability: reachability[0])[0][1]
        else:
            return point_index
    else: # Not sure if this else statement is actaully needed... ##
        return point_index

# Extract DBSCAN Equivalent cluster structure ##

# Important: Epsilon prime should be less than epsilon used in OPTICS #


def ExtractDBSCAN(SetOfObjects, epsilon_prime):
    """Performs DBSCAN equivalent extraction for arbitrary epsilon.
Can be run multiple times.

Parameters
----------
SetOfObjects: Prepped and build instance of setOfObjects
epsilon_prime: float or int
Must be less than or equal to what was used for prep and build steps

eturns
-------
Modified setOfObjects with cluster_id and is_core attributes."""

    # Start Cluster_id at zero, incremented to '1' for first cluster
    cluster_id = 0
    for entry in SetOfObjects._ordered_list:
        if SetOfObjects._reachability[entry] > epsilon_prime:
            if SetOfObjects._core_dist[entry] <= epsilon_prime:
                cluster_id += 1
                SetOfObjects._cluster_id[entry] = cluster_id
            else:
                # This is only needed for compatibility for repeated scans.
                # -1 is Noise points
                SetOfObjects._cluster_id[entry] = -1
        else:
            SetOfObjects._cluster_id[entry] = cluster_id
            if SetOfObjects._core_dist[entry] <= epsilon_prime:
                # One (i.e., 'True') for core points #
                SetOfObjects._is_core[entry] = 1
            else:
                # Zero (i.e., 'False') for non-core, non-noise points #
                SetOfObjects._is_core[entry] = 0

# End Algorithm #
