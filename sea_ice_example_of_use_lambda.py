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

#-----------------------------------------------------------------------

#This also imports matplotlib's pylab, as well as datetime(), timedelta()
#from datetime and random(), seed() from numpy.random
from sea_ice_rap_lambda import *
ion() # this is useful at the end, once in the IPython shell


#---------------------------------------------------------------------
# sea_ice_rap provides a class called "sea_ice" whose instances can
# be used to perform Random Analogue Predictions (rap) of sea ice
# extent.
# You initialize a class instance as follows

q = sea_ice('N_seaice_extent_daily_v4.0.csv')

# the argument 'N_seaice_extent_daily_v2.1.csv' is the csv file containing
# the north emisphere sea ice extent data, as provided by the National Snow
# and Ice Data Center (NSIDC). It must be found in the same folder where the
# program runs.

# After initialization the instance contains the raw NSIDC data and
# the corresponding days. Day zero corresponds to the date of the
# first valid datapoint, which is q.firstday.
# You can plot the data in this way:

plot(q.rawdays, q.rawextent, '.')
xlabel("days", fontsize=14)
ylabel("Millions of square km", fontsize=14)
title("NSIDC data. Starting date: "+str(q.firstday), fontsize=16)


#---------------------------------------------------------------------
# When initialized in this way, the a sea_ice instance comes with a
# default set of parameters. You can get a copy of them as follows:

print()
print(q.get_parameters())
print()

# The default parameters and their meaning are listed below.

DefaultParameters = {

    # parameters of the lowess smoother
    "resampling_interval":  5,  # days
    "window_length"      : 15,  # days
    "polydegree"         :  2,  # Keep this default! Use 1 if window_length <= 8

    # parameters of the embedding
    "lag"                : 36,  # measured in resampling intervals
    "embedding_dimension":  5,

    # parameters for the selection rule
    "lambda"             :  1., #take states within this distance from
                                #the current one. See the 'methods'
                                #section of the paper for a better way
                                #to pick "lambda".

}

# The first three are parameters of the low-pass filter that smooths the
# NSIDC data. I would not change "polydegree".

# Meaningful values for "window_length" range between 10 days and 90 days:
# we want to filter out weather systems (that last for one-two weeks) but
# we should not filter out the seasonal variations.

# Meaningful values for "resampling_interval" should range between 2 and
# 10 days. At the low end we'll end up having that the closest neighbour to
# a given state in phase state is often also the closest in time, but this
# defies the very concept of "analog prediction". At the high end of the range
# we are essentially throwing away data.

# The next two parameters detemine the geometry of the embedding. They
# must be integer.

# IMPORTANT: the "lag" parameter is measured in units of the
# "resampling_interval". Meaningful values of the product
# "lag"*"resampling_interval" should range between 30 and 300 days.

# The "embedding_dimension" should range between 1 and 5. Note that the actual
# dimension of the phase space is "embedding_dimension"+1 because of the
# 'compactification of time' trick that I use to impose to the data a notion
# of seasonality.

# Finally, "lambda" is the max distance such that makes a state a
# neighbour ('analog state') of the current state. If no neighbors are
# found within that radius, then the nearest state is
# deterministically picked. Small values of 'lambda' give little
# randomness, high values of lambda give a lot more variability, but
# the far away neighbours may not be good analogs. A reasonable
# initial choice for lambda is a value double than the average
# first-neighbor distance in the embedded dataset.

#---------------------------------------------------------------------------
# You can also read a single parameter directly from the sea_ice instance,
# as if it were a python dictionary:
print("The default value of 'lambda' is :", q['lambda'])

# You can also change a parameter by assigning to the corresponding key:

q['lambda'] = 0.75
print("The updated 'lambda' is :", q['lambda'])
print()
q['resampling_interval'] = 4.25
print("The updated 'resampling_interval' is :", q['resampling_interval'])
print()

# You can override completely the default parameters at initialization:

MyOwnParameters = {
    "resampling_interval" :  3, # days
    "window_length"       : 30, # days
    "polydegree"          :  2, # Keep this default! Use 1 if window_length <= 8
    "lag"                 : 10, # resampling intervals
    "embedding_dimension" :  3,
    "lambda"              :0.5,
}
q1 = sea_ice('N_seaice_extent_daily_v4.0.csv', MyOwnParameters)
print("Parameters of the instance q1:")
print(q1.get_parameters())
print()

# More on the effect of changing parameters later.

#---------------------------------------------------------------------------
# Here is how to  perform a forecast of a given number of days,
# starting at a given date.

F = q.forecast(datetime(2010, 6, 1), 150)

# The first argument, "datetime(2010, 6, 1)" determines the starting date
# and the second '150' is for how many days you should go on.
# Because the smoothed data are resampled only every "resampling_interval"
# days, the date that you choose may not correspond to a point in the phase
# space. The program selects the closest antecedent. The actual starting date
# of the forecast is stored in "q.forecast_starts".
# IMPORTANT NOTE: when you use .forecast() the datapoints corresponding
# to the specified period (the 150 days following June 1st 2010, in this
# example) are not used in the embedding!
# The result F is an array of points with "embedding_dimension"+1 components.
# The actual forecast is the first component of each of these points:

f = F[:,0]

# Now we can visualize the forecast (black line) vs reality (red line is the
# smoothed dataset, blue line is the actual NSIDC data). Time is in days, and
# zero corresponds to the beginning of the forecats (April 27th 2013).

figure()
plot(q.forecast_t, q.smoothextent[q.start_idx:q.end_idx], 'r-',
     q.raw_t, q.rawextent[q.start_raw:q.end_raw], 'b-',
     q.forecast_t, f, 'k-',
     linewidth=3)
legend(["Smoothed data", "Raw Data", "Forecast"])
xlabel("days", fontsize=18)
ylabel("Millions of square km", fontsize=18)

# q.forecast_t is an array of days spanning the length of the forecast
# [q.start_idx:q.end_idx] is the slice of the smoothed data corresponding
# to the forecast (and this is the part that has been discarded from
# the embedding).
# The corresponding quantities for the NSIDC data are q.raw_t and
# [q.start_raw:q.end_raw].

# If you have run the code, you have surely noticed that
# 'q.forecast(datetime(2013, 6, 1), 150)' is SLOW.
# This is because the first time .forecast() is called the NSIDC data
# is smoothed with the lowess filter (this is very slow), then it is
# embedded (this is almost fast), then finally the forecast is performed.

# If you want to perform another (or many other) forecasts with the
# same parameter set, then just call .forecast without any arguments.
# For example:

figure()
NForecasts = 100
for i in range(NForecasts):
    f = q.forecast()[:,0]
    plot(q.forecast_t, f, 'k-')
plot(q.forecast_t, q.smoothextent[q.start_idx:q.end_idx], 'r-',
     q.raw_t, q.rawextent[q.start_raw:q.end_raw], 'b-',
     linewidth=3)
legend(["Smoothed data", "Raw Data"])
xlabel("days", fontsize=18)
ylabel("Millions of square km", fontsize=18)
title(str(NForecasts) + " forecasts", fontsize=20)


#--------------------------------------------------------------------------
# When you run .forecast by specifying a date or a length, the previous
# embedding is discarded, and a new one is performed. However, the low-pass
# smooting is not performed again, so it's a relatively fast process.

# The same happens when you change 'embedding_dimension' or 'lag': the
# next forecast will compute a new embedding. In fact, you also have to
# specify a new date and length: the change of parameter deletes the
# old ones, because they are now meaningless.

# If you change 'polydegree', 'resampling_interval', or 'window_length'
# then, in addition to discarding the old embedding, also the smoothed data
# are discarded, and the next call to .forecast will trigger a new
# SLOW low-pass filtering.

# Finally, if you change 'number_of_neighbours', nothing needs to be
# re-computed, and thus subsequent forecasts are as fast as they can be.

#--------------------------------------------------------------------------
# The initial state selected by specifying a starting date can be
# perturbed by passing an optional perturbation vector (having
# "embedding_dimension"+1 components). For example, the following generates
# 100 forecast starting from the same randomly perturbe starting state.

figure()
NForecasts = 100
for i in range(NForecasts):
    pert = random(q['embedding_dimension']+1) - 0.5 #each component
                                                    #uniformly
                                                    #distributed in
                                                    #[-0.5, 0.5)
    f = q.forecast(perturbation=pert)[:,0]
    plot(q.forecast_t, f, 'k-')
plot(q.forecast_t, q.smoothextent[q.start_idx:q.end_idx], 'r-',
     q.raw_t, q.rawextent[q.start_raw:q.end_raw], 'b-',
     linewidth=3)
legend(["Smoothed data", "Raw Data"])
xlabel("days", fontsize=18)
ylabel("Millions of square km", fontsize=18)
title(str(NForecasts) + " initially perturbed forecasts", fontsize=20)

#--------------------------------------------------------------------------
# The random number generator can be reset to a repeatable
# initial state by using seed().
# For example, the following generates two identical forecasts:

seed(12345)
f1 = q.forecast()[:,0]
seed(12345)
f2 = q.forecast()[:,0]
print("Difference between the two forecatst:")
print(f1-f2)
print("\n\n")

#--------------------------------------------------------------------------
# Finally, there is a second interface for producing data that mimic
# the (smoothed) real NSIDC data. This is the method .predict(), which
# takes as input a vector with "embedding_dimension"+1 components and
# a length in days. It returns a forecast in steps of
# "resampling_interval". The first component of the forecast vectors
# can be compared with NSIDC data.
# For example:

start = array([10.,  13.,  5.,  14., -5, -3])
if len(start)!= q['embedding_dimension']+1:
    raise ValueError(f"This example assumes that the embedding dimension is {len(start)-1}")
f = q.predict(start, 150)[:,0]

figure()
plot(q.forecast_t, f, 'k-', linewidth=3)
xlabel("days", fontsize=18)
ylabel("Millions of square km", fontsize=18)
title("A prediction starting from a arbitrary state", fontsize=20)
# However, because the starting state is not associated with a date,
# you don't get q.start_idx, q.end_idx, q.start_raw, q.end_raw,
# q.raw_t. In addition, there is no check watsoever that the starting
# state is meaningful. If you are careless, and the starting state is
# far from the region sampled by the points of the embedding, results
# will also be meaningless.

#--------------------------------------------------------------------------
# Useful diagnostics: distances of each embedded vector from its next
# (in time) and from its closest

q.embed() #this insures that the embedding doesn't have any spurious
          #holes due to calls to 'forecast' which drops the data
          #corresponding to the period to be forecasted.
da, dm, d10, d90 = q.distance_to_next()
print(f"Average, median, 10th, 90th percentile of\ndistance to next-in-time embedding vector: {da:.3f}, {dm:.3f}, {d10:.3f}, {d90:.3f} ")
da, dm, d10, d90 = q.distance_to_closest()
print(f"Average, median, 10th, 90th percentile of\ndistance to closest embedding vector: {da:.3f}, {dm:.3f}, {d10:.3f}, {d90:.3f} ")


#--------------------------------------------------------------------------
