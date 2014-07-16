import scipy
import numpy as np

# This needs a lot of cleaning!
# Based on https://github.com/espg/OPTICS

class setOfObjects(object):

    """Build data structure with processing index from given data
in preparation for OPTICS Algorithm

Parameters
----------
distance_pairs: array [n_samples, n_samples]"""

    def __init__(self, distance_pairs):


# NOTE HERE: The way this is initialized is throwing the following deprecation warning:
# /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/
#   numpy/oldnumeric/__init__.py:11: ModuleDeprecationWarning: The oldnumeric module will
#   be dropped in Numpy 1.9
#   warnings.warn(_msg, ModuleDeprecationWarning)


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


# This will probably go unused...
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