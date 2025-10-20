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

from sea_ice_rap_lambda import *
from  statsmodels.graphics.functional import banddepth
import pickle

Nforecasts = 1000 #number of hindcasts of each year and parameter set.
method = 'MBD' #alternatives are 'MBD' or 'BD2'
rngseed = 7654321
Apert = 0.2 #amplitude of the perturbation on the initial condition.
seed(rngseed) #seeding the rng makes the simulation repeatable.
q = sea_ice('N_seaice_extent_daily_v4.0.csv')

# Sets the fixed parameters
q['resampling_interval'] = 5
q['polydegree'] = 2

window_length = [15, 30, 45, 60, 75, 90]
embedding_dimension = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
lag = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60]

#---------------------------------------------------------------------------

def compute_scores(q, years=arange(1979, 2025), Apert=0.2, Nforecasts=1000,
                   method='MBD'):
    """Computes band-depth scores of the forecasts obtained using q.
    q is a seaice instance; years is a list of years (integers); the band depth is computed using Nforecasts. For 'method' see the documentation of statsmodels.graphics.functional.banddepth; alternatives are 'MBD' (newer, less sensitive to ties) or 'BD2' (the original López-Pintado & Romo definition)
    Returns a list of years, and a list of matchning scores. The returned years may be a subset of the input list of years, because a data hole may make the forecast impossible for some years.
    Defaults to making 150 days long forecasts starting on June 1st.
    Defaults on using initial perturbations that are 1/2 of the median nearest neighbor distance."""
#    Apert = q.distance_to_closest()[1]/2
    scores = []
    actual_years = []
    for y in years:
        pert = 2*Apert*(random(q['embedding_dimension']+1)-0.5)
        try:
            data = q.forecast(datetime(y, 6, 1), 150, perturbation=pert)[:,0]
        except ValueError:
            continue
        for i in range(Nforecasts-1):
            pert = 2*Apert*(random(q['embedding_dimension']+1)-0.5)
            data = vstack((data, q.forecast(perturbation=pert)[:,0] ))
        data_plus = vstack((data, q.smoothextent[q.start_idx:q.end_idx]))
        bd = banddepth(data_plus, method=method)
        # the '+1' accounts for the fact that python indexes start at
        # zero.  In the score ranking, a ranking '1' is the lowest band
        # depth, and len(bd) is the highest band depth.
        score = (where(sort(bd)==bd[-1])[0][0]+1) / len(bd)
        scores.append(score)
        actual_years.append(y)
    return actual_years, scores

#---------------------------------------------------------------------------
results = {'seed': rngseed,
           'method': method,
           'amplitude_of_perturbation':Apert,
           'Nforecasts': Nforecasts
           }
for wl in window_length:
    q['window_length'] = wl
    for ed in embedding_dimension:
        q['embedding_dimension'] = ed
        for lg in lag:
            q['lag'] = lg
            q.embed()
            q['lambda'] = q.distance_to_closest()[1]*2 #search radius:
                                                       #twice the median
                                                       #n. neig. distance.
            print(f"window_length: {wl}, embedding_dimension: {ed}, lag: {lg}")
            dan, dmn, d10n, d90n = q.distance_to_next()
            print(f"   Average, median, 10th, 90th percentile of\n   distance to next-in-time embedding vector: {dan:.3f}, {dmn:.3f}, {d10n:.3f}, {d90n:.3f} ")
            dac, dmc, d10c, d90c = q.distance_to_closest()
            print(f"   Average, median, 10th, 90th percentile of\n   distance to closest embedding vector: {dac:.3f}, {dmc:.3f}, {d10c:.3f}, {d90c:.3f} ")
            #
            y, s = compute_scores(q, Apert=Apert, method=method)
            results[(wl, ed, lg)] = [(dan, dmn, d10n, d90n), (dac, dmc, d10c, d90c), y.copy(), s.copy()]
            print(f"   Mean score: {mean(s)}")
            print("---------------------------")

            
with open(f"explore_parameters_Nforecasts_{Nforecasts}_seed_{rngseed}_Apert_{Apert:.2f}_method_{method}.pkl", "wb") as f:
    pickle.dump(results, f)
