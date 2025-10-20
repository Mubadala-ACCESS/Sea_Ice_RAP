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

rngseed = 7654321
Nforecasts = 1000
Apert = 0.2 #amplitude of the perturbation on the initial condition.
seed(rngseed) #seeding the rng makes the simulation repeatable.
q = sea_ice('N_seaice_extent_daily_v4.0.csv')
q['resampling_interval'] = 5
q['window_length'] = 15
q['polydegree'] = 2
q['lag'] = 36
q['embedding_dimension'] = 5 
q['lambda'] = q.distance_to_closest()[1]*2
method = 'MBD' #band-depth computation method
years=arange(1979, 2026) #years for which the computation should be
                         #attempted. Note that for some years it can
                         #fail if the start data of June 1st falls
                         #into a hole of the embedding. Those years
                         #will be missing in the final results
iseptember = slice(19,25) #indexes of the forecast days in September,
                          #if starting on June 1st, and with
                          #resampling_interval = 5

actual_years = []
forecasted_september_averages = []
most_central_september_average = []
real_september_average = []
for y in years:
    pert = 2*Apert*(random(q['embedding_dimension']+1)-0.5)
    try:
        data = q.forecast(datetime(y, 6, 1), 150, perturbation=pert)[:,0]
    except ValueError:
        continue
    for i in range(Nforecasts-1):
        pert = 2*Apert*(random(q['embedding_dimension']+1)-0.5)
        data = vstack((data, q.forecast(perturbation=pert)[:,0] ))
    bd = banddepth(data, method=method)
    ibd = where(bd == amax(bd))[0][0]
    most_central_september_average.append(mean(data[ibd,iseptember]))
    actual_years.append(y)
    forecasted_september_averages.append(mean(data[:,iseptember], axis=1))
    obsslice = slice(q.start_idx, q.end_idx if q.end_idx<0 else None)
    real_september_average.append(mean(q.smoothextent[obsslice][iseptember]))

#-----------------------------------------------------------
fmin = [amin(x) for x in forecasted_september_averages]
f10p = [percentile(x, 10) for x in forecasted_september_averages]
fmed = [median(x) for x in forecasted_september_averages]
f90p = [percentile(x, 90) for x in forecasted_september_averages]
fmax = [amax(x) for x in forecasted_september_averages]
    
figure(figsize=(9,6))
fill_between(actual_years, f10p, f90p, color='gray', alpha=0.3,
             label=r"10$^\mathrm{th}$ - 90$^\mathrm{th}$ percentile")
plot(actual_years, real_september_average, '*-', c='tab:orange',
     label="NSIDC SIE Data")
plot(actual_years, most_central_september_average, 'o-', c='blue',
     label="Deepest RAP Hindcast")
plot(actual_years, fmed, '-', c='purple',
     label="Median of $10^3$ RAP hindcasts")
legend(loc="upper right")
ylabel("SIE September Average ($10^6$ km)", fontsize=16)
xlabel("Year", fontsize=16)
xticks(fontsize=14)
yticks(fontsize=14)
tight_layout()
#-----------------------------------------------------------------
print(f"RMSE data-rap: {sqrt(mean((array(real_september_average[:-1])-array(most_central_september_average[:-1]))**2))}")
print(f"RMSE data-persistence: {sqrt(mean((array(real_september_average[1:-1])-array(real_september_average[:-2]))**2))}")
grid(True, linewidth=0.5)
tight_layout()

#-----------------------------------------------------------
#savefig("SIE_September_average.png", dpi=300)
