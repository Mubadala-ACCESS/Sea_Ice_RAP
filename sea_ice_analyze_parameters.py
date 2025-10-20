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
import pickle

fname = "explore_parameters_Nforecasts_1000_seed_7654321_Apert_0.20_method_MBD.pkl"
with open(fname, "rb") as f:
    results = pickle.load(f)

# The dictionary 'results' should contain keys consisting of tuples
# with 3 numbers (the parameter values) associated to a list of four terms:
# 1. list of next-in-time vector distance statistics (mean, med., 10th, 90th);
# 2. list of closest vector distance statistics (mean, med., 10th, 90th);
# 3. list of years;
# 4. list of band depth scores.
# Additional 'seed', 'method', 'amplitude_of_perturbations',
# 'Nforecasts' keys specify the rng seed used for the simulation, the
# method used for the band depth calculation, and the maximum distance
# of the perturbed initial condition from the observed initial state,
# and the number of hindcasts used for each score entry in the results
# database.

# Extract the parameters as arrays
windows = set([])
embedding_dimensions = set([])
lags = set([])
for k in results.keys():
    if (type(k) is tuple) and (len(k)==3):
        w, e, l = k
        windows.add(w)
        embedding_dimensions.add(e)
        lags.add(l)
windows = sort(array(list(windows)))
embedding_dimensions = sort(array(list(embedding_dimensions)))
lags = sort(array(list(lags)))

# Extract the years-score data
minyear     = zeros((len(windows), len(embedding_dimensions), len(lags)))
minscore    = zeros_like(minyear)
meanscore   = zeros_like(minyear)
medianscore = zeros_like(minyear)

for iw, w in enumerate(windows):
    for ie, e in enumerate(embedding_dimensions):
        for il, l in enumerate(lags):
            #years, scores = results[(w, e, l)]
            _, _, years, scores = results[(w, e, l)]
            years = array(years)
            scores = array(scores)
            mi = amin(scores)
            me = mean(scores)
            md = median(scores)
            minscore[iw, ie, il] = mi
            meanscore[iw, ie, il] = me
            medianscore[iw, ie, il] = md
            minyear[iw, ie, il] = where(scores==mi)[0][0]
            
#------------------------------------------------------

iw, ie, il = where(medianscore == amax(medianscore))
print(f"Highest median score: {amax(medianscore)}")
print(f"attained for W={windows[iw[0]]}; M={embedding_dimensions[ie[0]]}; tau={lags[il[0]]}")
