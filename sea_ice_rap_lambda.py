# (C) 2025 Francesco Paparella & Faiq Raees
# This file is part of Sea_Ice_RAP.
#
# Sea_Ice_RAP is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Sea_Ice_RAP is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Sea_Ice_RAP. If not, see <https://www.gnu.org/licenses/>.

from pylab import *
from datetime import datetime, timedelta
from lowess import lowess
from scipy.spatial import cKDTree
from numpy.random import random, seed


# ------------------------------------------------------------
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res

# ------------------------------------------------------------
def read_data_file(file_name):
    """Gets as input the NSIDC's .cvs daily extent file, as downloaded
from
    ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/
or from
    https://masie_web.apps.nsidc.org/pub/DATASETS/NOAA/G02135/north/daily/data/
    
Returns the array of days (starting from zero), extent (in millions
of km^2.) and the date of the first datapoint (as datetime.datetime object).

Note that most of the dataset has daily data (so "days" will just
contain consecutive integers) but initially there is one datapoint
every two days, and there is a hole between December 1987 and January
1988. So the data, overall, is not evenly sampled.

    """
    fdata = open(file_name)
    year = []
    month = []
    day = []
    extent = []
    for line in fdata:
        ll = line.split(',')
        if converts_to_integer(ll[0]):
            year.append(int(ll[0]))
            month.append(int(ll[1]))
            day.append(int(ll[2]))
            extent.append(float(ll[3]))

    fdata.close()
    extent = array(extent)
    days = zeros(len(year), dtype='float')
    FirstDay = datetime(year[0], month[0], day[0])
    for i in range(len(year)):
        curdate = datetime(year[i], month[i], day[i])
        days[i] = (curdate - FirstDay).days
    return days, extent, FirstDay

# ------------------------------------------------------------


def converts_to_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

# ------------------------------------------------------------
# Start of class sea_ice
# ------------------------------------------------------------


DefaultParameters = {
    # parameters of the lowess smoother
    "resampling_interval":  5,  # days
    "window_length": 15,  # days
    "polydegree":  2,  # Keep this default! Use 1 if window_length <= 8
    # parameters of the embedding
    "lag": 36,  # resampling intervals
    "embedding_dimension":  5,
    # parameters for the selection rule
    "lambda":  1.,
}


class sea_ice(object):
    def __init__(self, file_name, parameters=DefaultParameters):
        self._parameters = parameters
        self.rawdays, self.rawextent, self.firstday = read_data_file(file_name)
        self.embedding = None
        self.points = None
        self.idx_pts_with_successor = None
        self.forecast_starts = None
        self.forecast_length = None
        self.start_idx = None
        self.end_idx = None
        self.forecast_t = None
        self.start_raw = None
        self.end_raw = None
        self.raw_t = None
        self.smoothdays = None
        self.smoothextent = None

    # ------------------------------------------------------------
    def smooth(self):
        """Evenly resamples the raw data by applying a low-pass lowess
filter.

Computes:
  self.smoothdays   : Time of the resampled datapoints, in days. Starts at
                      zero, progressses with resampling_interval steps.
  self.smoothextent : Resampled sea ice extent in millions of km^2.
                      Missing data in the 1987/1988 hole are marked by nan.
        """
        days = arange(0., self.rawdays[-1],
                      self._parameters["resampling_interval"])
        gap_start = 1714  # These are true as long as NSIDC desn't
        gap_end = 1715  # change the old data
        idx_early = where(days <= self.rawdays[gap_start])[0]
        idx_late = where(days >= self.rawdays[gap_end])[0]

        extent = days * nan

        wl = self._parameters["window_length"]
        pd = self._parameters["polydegree"]
        extent[:len(idx_early)] = lowess(self.rawdays, self.rawextent,
                                         days[idx_early],
                                         deltax=wl, degree=pd)[:, 0]
        extent[-len(idx_late):] = lowess(self.rawdays, self.rawextent,
                                         days[idx_late],
                                         deltax=wl, degree=pd)[:, 0]
        self.smoothdays = days
        self.smoothextent = extent

    # ------------------------------------------------------------
    def embed(self, holed_extent_data=None):
        """Creates the embedding data structure (a KDTree instance). Note that
the actual embedding dimension is self.parameters["embedding_dimension"]+1,
because the phase trick that takes into account the seasonality requires
one extra dimension.

Adds to the class instance:
  self.points : Array of embedded points. Because of the hole(s), some will
                have one or more components set to nan.
  self.idx_pts_with_successor : self.points[self.idx_pts_with_successor] yields
                                all the points with non-nan components which
                                have a successor with non-nan components.
  self.embedding = embedding : a cKDTree instance, useful for, e.g., k-nearest
                               neighbour searches.

        """
        if holed_extent_data is None:
            if self.smoothextent is None:
                self.smooth()
            holed_extent_data = self.smoothextent
        edim = self._parameters["embedding_dimension"]
        yearlength = 365.25
        points = []
        idx = (edim - 1) * self._parameters["lag"]
        while idx < len(holed_extent_data):
            p = []
            for i in range(edim):
                p.append(holed_extent_data[idx - i * self._parameters["lag"]])
            lastp = p[-1]
            phase = self._parameters["resampling_interval"] * \
                idx * 2 * pi / yearlength
            p[-1] = lastp * sin(phase)
            p.append(lastp * cos(phase))
            points.append(p)
            idx += 1
        points = array(points)
        # Exclude all the points which are invalid or do not have a valid
        # successor in time.
        goodpoints = []
        idxmap = []
        for i in range(len(points[:-1])):
            if not isnan(norm(points[i]) * norm(points[i + 1])):
                goodpoints.append(points[i])
                idxmap.append(i)
        embedding = cKDTree(array(goodpoints))
        self.points = points
        self.idx_pts_with_successor = idxmap
        self.embedding = embedding

    # ------------------------------------------------------------
    def forecast(self, startdate=None, length=None, perturbation=None):
        """Performs a hindcast beginning at a scecified data within the date range covered by the embedded dataset.
        'startdate'    is a datetime.datetime object;
        'length'       is the length of the hindcast in days (an integer);
        'perturbation' is a vector having the same length as the vectors in self.points.
        
The resulting hindcast is sampled in steps having length equal to
'resampling_interval' (as specified by the parameter set).

The initial state of the hindcast is taken as the closest antecedent
to the given start date among all the points in the embedding. The
embedding is performed by discarding all datapoints after the date of the
initial state, up to a number of days equal to 'length'.

If 'startdate' or 'length' are omitted, the last used 'startdate' and
'length' are used for the forecast. The embedding is not
re-computed. An error is raised if the last 'startdate' or 'length'
are not available.
       
If 'perturbation' is specified, it has to be a vector with
'embedding_dimension'+1 components. Then the initial state is computed
by adding this vector to the closest antecedent point to startdate (see above).

Computes:
   self.forecast_starts : the starting YYYY/MM/DD of the last forecast performed
                          (a datetime.datetime instance).
                          NOTE: this is not necessarily equal to 'startdate'.
                          It is the date of the closest antecedent resampled
                          datapoint to 'stardate'

   self.forecast_length : the length in days of the last forecast performed

   self.start_idx       : the index of the resampled extent corresponding to
                          self.forecast_starts

   self.end_idx         : the index of the resampled extent corresponding to
                          the end of the forecast

   self.forecast_t      : time in days of the forecast; zero corresponds to
                          the date self.forecast_starts

   self.start_raw       : same as .start_idx but for the unsmooted NSIDC data

   self.end_raw         : same as .end_idx but for the unsmooted NSIDC data

   self.raw_t           : same as .forecast_t but for the unsmooted NSIDC data

        """
        # This is a wrapper that does some sanity checks. The actual forecast
        # is performed by __forecast()

        if startdate is None and type(self.forecast_starts) != datetime:
            raise TypeError("I don't have a forecast starting date memorized:" +
                            " you should specify one!")
        if length is None and self.forecast_length is None:
            raise TypeError("I don't have a forecast length memorized:" +
                            " you should specify one!")
        if self.smoothextent is None:
            self.smooth()
        hole_data = False
        res = self._parameters["resampling_interval"]
        NumEmbedPoints = (len(self.smoothextent) -
                          (self._parameters["embedding_dimension"] - 1) *
                          self._parameters["lag"])

        if startdate is not None:
            if type(startdate) != datetime:
                raise TypeError(
                    "'startdate' must be of type datetime.datetime")
            lastdataday = self.firstday + timedelta(self.rawdays[-1])
            # moves the starting date to the closest antecedent data point
            self.forecast_starts = self.firstday + timedelta(
                res * floor((startdate - self.firstday).total_seconds() /
                            (res * 24. * 3600.)))
            if self.forecast_starts > lastdataday:
                raise ValueError("     startdate = " +
                                 str(startdate) + "\n"
                                 "implies\n" +
                                 "       forecast start date = " +
                                 str(self.forecast_starts) + "\n" +
                                 "which is later than\n"
                                 "   date of last data point = " +
                                 str(lastdataday))
            self.start_idx = int(floor(
                ((self.forecast_starts - self.firstday).total_seconds() -
                 self.rawdays[-1] * 24. * 3600.) / (res * 24. * 3600.)))
            if NumEmbedPoints + self.start_idx < 0:
                raise ValueError("'startdate' is earlier than the first point" +
                                 " in the embedding\n" +
                                 "            Start date = " +
                                 str(startdate) + "\n" +
                                 "   Date of first point = " +
                                 str(lastdataday -
                                     timedelta(NumEmbedPoints * res)))
            hole_data = True

        if length is not None:
            self.forecast_length = int(float(length) / res)
            hole_data = True

        if hole_data:
            self.holed_extent_data = self.smoothextent.copy()
            ileft = self.start_idx + 1 if self.start_idx < -1 else \
                len(self.smoothextent)
            iright = self.start_idx + self.forecast_length + 1 if \
                self.start_idx + self.forecast_length + 1 < 0 else None
            self.holed_extent_data[ileft:iright] *= nan
            self.embed(self.holed_extent_data)
            self.forecast_t = arange(0, self.forecast_length + 1) * res
            self.end_idx = self.start_idx + self.forecast_length + 1
            self.start_raw = find(self.rawdays <= (
                self.forecast_starts - self.firstday).days)[-1]
            aux = find(self.rawdays >= self.rawdays[self.start_raw] +
                                       self.forecast_length * res)
            self.end_raw = None if len(aux)==0 else aux[0]
            self.raw_t = (self.rawdays[self.start_raw:self.end_raw] -
                          self.rawdays[self.start_raw])

        start_state = self.points[self.start_idx]
        if isnan(norm(start_state)):
            raise ValueError("The starting date " + str(self.forecast_starts) +
                             " corresponds to a hole in the dataset. " +
                             "(Is it Dec 1987 - Jan 1988 ?)")

        if perturbation is not None:
            return self.__forecast(start_state + perturbation)
        else:
            return self.__forecast(start_state)

    # ------------------------------------------------------------
    def __select_one_neighbour(self, point):
        """Right now the probability for a neighbour of being selected is
    inversely proportional to its distance from the given point.

        """
        # nn is the n. of embedding points within a sphere of radius lambda from point.
        # If point is an embedding point it will be returned by query with distance 0
        # thus idx with distance zero must be discarded. If there's less than two points
        # in the sphere, query for 2 embedding points, one of these will be the
        # closest embedding point to point. Keep that even if it's outside of
        # the sphere.
        nn = len(self.embedding.query_ball_point(point,
                                                 self._parameters["lambda"]))
        if nn > 1:
            d, idx = self.embedding.query(point, k=nn)
        else:
            d, idx = self.embedding.query(point, k=2)
            if d[0] > 0:
                d = d[:-1]
                idx = idx[:-1]
        d = array(d) ; idx = array(idx)
        if d[0]==0:
            d = d[1:]
            idx = idx[1:]
        invd = 1./d
        cdf = cumsum(invd/sum(invd))
        #---debug stuff---
        #self.d = d
        #self.cdf = cdf
        #---/debug---
        return idx[where(cdf > random())[0][0]]

    # ------------------------------------------------------------
    def __forecast(self, start_state):
        """Don't use this stand-alone. Run it through .forecast() or
        .predict() which perform prerequisite computations."""
        forecast = [start_state]
        for loop in range(self.forecast_length):
            cp = forecast[-1]
            # Remember: the i-th point in the embedding is the
            # self.idx_pts_with_successor[i]-th point in self.points
            idx = self.idx_pts_with_successor[self.__select_one_neighbour(cp)]
            displacement = (self.points[idx + 1] - self.points[idx])
            forecast.append(cp + displacement)
        return array(forecast)

    # ------------------------------------------------------------

    def predict(self, start_state, length):
        """Produces a forecast beginning from 'start_state', with a duration
of 'length' days. The resulting forecast is sampled in steps having
length equal to 'resampling_interval' days (as specified by the
parameter set).

The vector 'start_state' must have 'embedding_dimension'+1 components.

After a call:

    self.forecast_starts : is a copy of the vector start_state

        """
        self.embed()
        res = self._parameters["resampling_interval"]
        self.forecast_length = int(length / res)
        self.forecast_starts = start_state.copy()
        self.forecast_t = arange(0, self.forecast_length + 1) * res
        self.start_idx = None
        self.end_idx = None
        self.start_raw = None
        self.end_raw = None
        self.raw_t = None
        return self.__forecast(start_state)

    # ------------------------------------------------------------
    def get_parameters(self):
        return self._parameters.copy()

    # ------------------------------------------------------------
    def distance_to_next(self):
        """Returns the average, median, 10th, 90th percenntile of the distance of each vector in the embedding to its successor (next-in-time vector).
If an embedding is available, this uses the current one. Keep in mind that if 'self.forecast' was called before, then the current embedding has a hole starting at the forecast date and as long as the forecast length. If this is undesired, call 'self.embed()' before calling this method.
        """
        if self.embedding is None:
            self.embed()
        vectors = self.embedding.data
        distance_to_next = [norm(vectors[i] - vectors[i+1])
                            for i in range(len(vectors)-1)]
        return (mean(distance_to_next), median(distance_to_next),
                percentile(distance_to_next, 10),
                percentile(distance_to_next, 90))
    # ------------------------------------------------------------
    def distance_to_closest(self):
        """Returns the average, median, 10th, 90th percenntile of the distance of each vector in the embedding to its successor (next-in-time vector).
If an embedding is available, this uses the current one. Keep in mind that if 'self.forecast' was called before, then the current embedding has a hole starting at the forecast date and as long as the forecast length. If this is undesired, call 'self.embed()' before calling this method.
        """
        if self.embedding is None:
            self.embed()
        vectors = self.embedding.data
        distance_to_closest = [self.embedding.query(v, k=2)[0][1] for v in vectors]
        return (mean(distance_to_closest), median(distance_to_closest),
                percentile(distance_to_closest, 10),
                percentile(distance_to_closest, 90))
    # ------------------------------------------------------------
    def __getitem__(self, key):
        if key not in self._parameters.keys():
            raise KeyError("Unknown parameter " + str(key))
        return self._parameters[key]

    # ------------------------------------------------------------
    def __setitem__(self, key, value):
        if key not in self._parameters.keys():
            raise KeyError("Unknown parameter " + str(key))
        self._parameters[key] = value
        if key != "lambda":
            self.embedding = None
            self.forecast_starts = None
            self.forecast_length = None
            self.start_idx = None
            self.end_idx = None
            self.forecast_t = None
            self.start_raw = None
            self.end_raw = None
            self.raw_t = None
        if key in ["polydegree", "resampling_interval", "window_length"]:
            self.smoothextent = None
            self.smoothdays = None


# ------------------------------------------------------------
# End of class sea_ice
# ------------------------------------------------------------
