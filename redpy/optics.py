# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016  Alicia Hotovec-Ellis (ahotovec@gmail.com)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import scipy
import numpy as np

# Based on https://github.com/espg/OPTICS

class setOfObjects(object):

    """
    Build data structure with processing index from given data
    in preparation for OPTICS Algorithm
    
    distance_pairs: Distance matrix (array [n_samples, n_samples])
    
    """

    def __init__(self, distance_pairs):

        """
        NOTE HERE: The way this is initialized is throwing the following warning:
        ModuleDeprecationWarning: The oldnumeric module will be dropped in Numpy 1.9
        warnings.warn(_msg, ModuleDeprecationWarning)
        """

        self.data = distance_pairs
        self._n = len(self.data)
        self._processed = scipy.zeros((self._n, 1), dtype=bool)
        self._reachability = scipy.ones(self._n) * scipy.inf
        self._core_dist = scipy.ones(self._n) * scipy.nan
        self._index = scipy.array(range(self._n))
        self._nneighbors = scipy.ones(self._n, dtype=int)*self._n
        self._cluster_id = -scipy.ones(self._n, dtype=int)
        self._is_core = scipy.ones(self._n, dtype=bool)
        self._ordered_list = []
        

def prep_optics(SetofObjects, epsilon):

    """
    Prep data set for main OPTICS loop
    
    SetofObjects: Instantiated instance of 'setOfObjects' class
    epsilon: Determines maximum object size that can be extracted. Smaller epsilons
        reduce run time.
    
    Returns modified setOfObjects tree structure
    
    """

    for j in SetofObjects._index:
        # Find smallest nonzero distance
        SetofObjects._core_dist[j] = np.sort(SetofObjects.data[j,:])[1]


def build_optics(SetOfObjects, epsilon):

    """
    Builds OPTICS ordered list of clustering structure
    
    SetofObjects: Instantiated and prepped instance of 'setOfObjects' class
    epsilon: Determines maximum object size that can be extracted. Smaller epsilons
        reduce run time.

    """

    for point in SetOfObjects._index:
        if not SetOfObjects._processed[point]:
            expandClusterOrder(SetOfObjects, point, epsilon)


def expandClusterOrder(SetOfObjects, point, epsilon):

    """
    Expands OPTICS ordered list of clustering structure
    
    SetofObjects: Instantiated and prepped instance of 'setOfObjects' class
    epsilon: Determines maximum object size that can be extracted. Smaller epsilons
        reduce run time.

    """
    
    if SetOfObjects._core_dist[point] <= epsilon:
        while not SetOfObjects._processed[point]:
            SetOfObjects._processed[point] = True
            SetOfObjects._ordered_list.append(point)
            point = set_reach_dist(SetOfObjects, point, epsilon)
    else:
        SetOfObjects._processed[point] = True


def set_reach_dist(SetOfObjects, point_index, epsilon):

    """
    Sets reachability distance and ordering. This function is the primary workhorse of
    the OPTICS algorithm.
    
    SetofObjects: Instantiated and prepped instance of 'setOfObjects' class
    epsilon: Determines maximum object size that can be extracted. Smaller epsilons
        reduce run time. (float)

    """
    
    row = [SetOfObjects.data[point_index,:]]
    indices = np.argsort(row)
    distances = np.sort(row)

    if scipy.iterable(distances):

        unprocessed = indices[(SetOfObjects._processed[indices] < 1)[0].T]
        rdistances = scipy.maximum(distances[(SetOfObjects._processed[indices] < 1)[0].T],
            SetOfObjects._core_dist[point_index])
        SetOfObjects._reachability[unprocessed] = scipy.minimum(
            SetOfObjects._reachability[unprocessed], rdistances)

        if unprocessed.size > 0:
            return unprocessed[np.argsort(np.array(SetOfObjects._reachability[
                unprocessed]))[0]]
        else:
            return point_index
    else:
        return point_index
