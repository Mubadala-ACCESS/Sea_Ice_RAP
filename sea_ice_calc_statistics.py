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
from dateutil.relativedelta import relativedelta
from statsmodels.graphics.functional import banddepth
from scipy.stats import pearsonr
import matplotlib.dates as mdates
import pickle
seed(7654321)  #seeding the rng makes the simulation repeatable.

#--------------------------------------------------------------
def rmse(x): 
    return np.sqrt(np.mean((x) ** 2)) 

#--------------------------------------------------------------
def bias(x):
    return np.mean(x) 

#--------------------------------------------------------------
def hindcast(q, startdate, length, Nforecasts, Apert, month, method='MBD'):
    pert  = 2*Apert*(random(q['embedding_dimension']+1)-0.5)
    data  = q.forecast(startdate, length, perturbation=pert)[:,0]
    dates = [q.forecast_starts + timedelta(float(d)) for d in q.forecast_t]
    month_idx =  [int(dates.index(d)) for d in dates if (d.month == month) ]
    for i in range(Nforecasts-1):
        pert = 2*Apert*(random(q['embedding_dimension']+1)-0.5)
        data = vstack((data, q.forecast(perturbation=pert)[:,0] ))
    bd = banddepth(data, method=method)
    idxrap = where(bd==amax(bd))[0][0]
    mean_data = mean(data, axis = 0)
    return idxrap, data, mean_data, dates, month_idx
    
#--------------------------------------------------------------
def compute_climate_mean(q, months): 
    smooth = q.smoothextent
    smooth_dates = [q.firstday + timedelta(q.smoothdays[i])
                    for i in range(len(q.smoothdays))]
    bins = []
    day = 1 
    while day <= 28:
        end = min([day + q['resampling_interval']  - 1, 31])
        if end < 28: 
            bins.append((day, int(end)))
        else :
            bins.append((day, 31))
        day += q['resampling_interval'] 
    climate_dates = []
    climate_mean  = []
    for m in months: 
        for i in range(len(bins)):
            binidx  = [smooth_dates.index(d)
                       for d in smooth_dates
                       if (d.month == m and (d.day >= bins[i][0]
                                             and d.day <= bins[i][1])
                           )
                       ]
            bindata = smooth[binidx]
            climate_mean.append(np.nanmean(bindata))
            climate_dates.append(datetime(2000,m,bins[i][0]))
    return climate_dates, array(climate_mean)
    
#--------------------------------------------------------------
#Parameters
NForecasts = 1000
years  = arange(1979,2025)
leads  = arange(1,10)
months = arange(1,13)
Apert  = 0.2

#--------------------------------------------------------------
# q = sea_ice('N_seaice_extent_daily_v3.0.csv')
q = sea_ice('N_seaice_extent_daily_v4.0.csv')
q['resampling_interval'] = 5
q['window_length'] = 15
q['polydegree'] = 2
q['lag'] = 36
q['embedding_dimension'] = 5
q.embed()
q['lambda'] = q.distance_to_closest()[1]*2
q.smooth()
smooth_dates = [q.firstday + timedelta(q.smoothdays[i])
                for i in range(len(q.smoothdays))]
smooth = q.smoothextent

climate_dates, climate_mean = compute_climate_mean(q, months)

#--------------------------------------------------------------
# Plot the climatological average
figure() 
plot(climate_dates, climate_mean, color = "blue", linewidth =5)
gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
gca().xaxis.set_major_locator(mdates.MonthLocator())
plot( [d.replace(year = climate_dates[0].year) for d in smooth_dates if d.year == 2000], smooth[[smooth_dates.index(d) for d in smooth_dates if d.year ==2000]], color = "grey", alpha = 0.3 )
legend(["Climate Mean", "NSIDC {}-{}".format(years[0], years[-1])])
for y in years: 
    plot( [d.replace(year = climate_dates[0].year) for d in smooth_dates if d.year == y] , smooth[[smooth_dates.index(d) for d in smooth_dates if d.year ==y]], color = "grey", alpha = 0.3 )
show()

#--------------------------------------------------------------
# Main loop
#--------------------------------------------------------------
bias_array = zeros((len(months), len(leads), len(years)))
rmse_array = zeros((len(months), len(leads), len(years)))
acc_array  = zeros((len(months), len(leads), len(years)))

for year in years:
    print(year)
    for lead in leads:
        for m in months:
            targetmonth = datetime(year,m,1)
            startdate = targetmonth - relativedelta(months=lead)
            length = (lead+1)*31
            # try/except to handle holes in the NSIDC dataset 
            try:
                (idxrap, hindcasts, meanhindcast,
                 forecast_dates, month_idx) = hindcast(q, startdate, length,
                                                       NForecasts, Apert, m)
            except ValueError:
                continue
            # if the forecast covers the hole in the NSIDC data skip the
            # calculation of statistics
            if isnan(sum(smooth[array(q.start_idx) + array(month_idx)])):
                continue
            climate_dates_month = [d for i,d in enumerate(climate_dates)
                                   if d.month == array(forecast_dates
                                                       )[month_idx][0].month
                                   ] 
            if len(array(forecast_dates)[month_idx]) < 6: 
                climate_dates_month = (climate_dates_month[0:-1]
                            if array(forecast_dates)[month_idx][0].day < 6
                            else climate_dates_month[1:]
                                       )
            if len(array(forecast_dates)[month_idx]) == 7: 
                climate_dates_month.append(climate_dates_month[-1])
            climate_dates_month_idx = [climate_dates.index(d)
                                       for d in climate_dates_month
                                       ]
            climate_mean_month = climate_mean[climate_dates_month_idx]
            #--------------------------------------------------------------
            #Replace hindcasts[idxrap, month_idx] by meanhindcast[month_idx]
            #for comparision with the mean of RAP forecasts, rather than with
            #the deepest of RAP forecasts
            modeloutput  = hindcasts[idxrap, month_idx]
            observations = smooth[array(q.start_idx) + array(month_idx)]
            #-- Compute the Bias --
            bias_array[m-1,lead-1, year - years[0]] = bias(modeloutput - observations)
            #-- Compute the RMSE --
            rmse_array[m-1,lead-1, year - years[0]] = rmse(modeloutput - observations)
            #-- Compute ACC --
            acc_array[m-1,lead-1, year - years[0]]  = \
                pearsonr(modeloutput  - climate_mean_month ,
                         observations - climate_mean_month)[0] 


#--------------------------------------------------------------
stats = {
    'bias_array'    : bias_array,
    'rmse_array'    : rmse_array,
    'acc_array'     : acc_array,
    'idxrap'        : idxrap,
    'hindcasts'     : hindcasts,
    'meanhindcast'  : meanhindcast,
    'forecast_dates': forecast_dates,
    'month_idx'     : month_idx,
    'rap_parameters': q.get_parameters(),
    'years'         : years,
    'months'        : months,
    'leads'         : leads,
    'Apert'         : Apert
}
fname = "Statistics_"+str(datetime.today().strftime('%Y-%m-%d-%H-%M'))+".pkl"
with open(fname, "wb") as statsfile:
    pickle.dump(stats, statsfile)
#--------------------------------------------------------------
