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
from scipy import stats
import matplotlib.dates as mdates

#-- Load the data saved by sea_ice_calc_statistics.py --
fname = 'Statistics_2025-10-02-23-38.pkl'
with open(fname, 'rb') as f:
    statsdata = pickle.load(f)
for key in statsdata.keys():
    exec(f"{key} = statsdata['{key}']")


#--------------------------------------------------------------
total_bias_array_mean_over_years = mean(bias_array, axis = -1)
total_rmse_array_mean_over_years = mean(rmse_array, axis = -1)
total_acc_array_mean_over_years  = mean(acc_array, axis = -1)

#--------------------------------------------------------------
# See: https://en.wikipedia.org/wiki/Sign_test
def one_sample_sign_test(data, hypothesized_median, alternative='two-sided'):
    differences = data - hypothesized_median
    
    # Exclude zero differences
    non_zero_differences = differences[differences != 0]
    
    positive_signs = sum(non_zero_differences > 0)
    negative_signs = sum(non_zero_differences < 0)
    n = len(non_zero_differences)
    # Determine the test statistic (k) based on the alternative hypothesis
    if alternative == 'two-sided':
        k = np.min(positive_signs, negative_signs)
        p_value = stats.binomtest(k, n, p=0.5, alternative='two-sided').pvalue
    # Testing if median > hypothesized_median (more positive differences)
    elif alternative == 'greater': 
        k = negative_signs
        p_value = stats.binomtest(k, n, p=0.5, alternative='less').pvalue
    # Testing if median < hypothesized_median (more negative differences)
    elif alternative == 'less': 
        k = positive_signs
        p_value = stats.binomtest(k, n, p=0.5, alternative='less').pvalue 
    else:
        raise ValueError("Alternative must be 'two-sided', 'greater', or 'less'.")  
    return p_value
#--------------------------------------------------------------


#-- Plot the Bias --
fig = figure(figsize=(5, 5))
L,M = meshgrid(months, leads)
limit = max([abs(np.min(total_bias_array_mean_over_years)),
             np.max(total_bias_array_mean_over_years)])
map = pcolormesh(months, leads, total_bias_array_mean_over_years.transpose(),
                 shading = 'auto', cmap = cm.RdBu_r)
map.set_clim(-1*limit, limit)
biasbar = colorbar(map,orientation='horizontal')
biasbar.set_label("Bias ($10^6$ km$^2$)", fontsize=12)
ylabel("Lead time in months", fontsize=16)
xlabel("Target Month", fontsize=16)
xticks(months, ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
yticks(leads)
for lead in leads:
    for m in months:
        if stats.ttest_1samp(bias_array[m-1,lead-1, :], 0)[1] < 0.05:
            plot(L[lead-1,m-1], M[lead-1,m-1], 'o', color = "black")
text(-0.85, 9.5, 'A', fontsize=16)
tight_layout()

#-- Plot the RMSE --
fig = figure(figsize=(5, 5))
limit = np.max(total_rmse_array_mean_over_years)
map = pcolormesh(months, leads, total_rmse_array_mean_over_years.transpose(),
                 shading = 'auto', cmap = cm.viridis)
map.set_clim(0, limit)
RMSEbar = colorbar(map,orientation='horizontal')
RMSEbar.set_label("RMSE ($10^6$ km$^2$)", fontsize=12)
ylabel("Lead time in months", fontsize=16)
xlabel("Target Month", fontsize=16)
xticks(months, ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
yticks(leads)
for lead in leads:
    for m in months:
        if one_sample_sign_test(rmse_array[m-1,lead-1, :], Apert,
                                alternative="greater") < 0.05:
            plot(L[lead-1,m-1], M[lead-1,m-1], 'o', color = "black")
text(-0.85, 9.5, 'B', fontsize=16)
tight_layout()

#-- Plot ACC --
fig = figure(figsize=(5, 5))
limit = max([abs(np.min(total_acc_array_mean_over_years)),
             np.max(total_acc_array_mean_over_years)])
map = pcolormesh(months, leads, total_acc_array_mean_over_years.transpose(),
                 shading = 'auto', cmap = cm.RdBu_r)
map.set_clim(-1*limit, limit)
ACCbar = colorbar(map,orientation='horizontal')
ACCbar.set_label("ACC", fontsize=12)
ylabel("Lead time in months", fontsize=16)
xlabel("Target Month", fontsize=16)
xticks(months, ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
yticks(leads)
for lead in leads:
    for m in months:
        if stats.ttest_1samp(acc_array[m-1,lead-1, :], 0)[1] < 0.05:
            plot(L[lead-1,m-1], M[lead-1,m-1], 'o', color = "black")
text(-0.85, 9.5, 'C', fontsize=16)
tight_layout()
