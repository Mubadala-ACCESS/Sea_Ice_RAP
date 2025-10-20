     --- Sea_Ice_RAP ---
Random Analog Prediction of Arctic Sea Ice Extent using NSIDC SIE data.

The core of the software is the "sea_ice" class in sea_ice_rap_lambda.py.

The best way to learn how to use the code is to interactively execute line by line the content of the file sea_ice_example_of_use.py in ipython or a similar interactive python shell, while reading the comments in that file. 

The following files use sea_ice_rap_lambda.py and were used to prepare the paper "Random Analog Prediction of Arctic Sea Ice Extent: a Benchmark for Seasonal Models." (by Faiq Raees & Francesco Paparella, submitted, 2025):

- sea_ice_explore_parameters.py explores the RAP parameter space and
                                saves a .pkl file containing the band-depth
				of the results.
				
- sea_ice_analyze_parameters.py uses the rusults of
                                sea_ice_explore_parameters.py to
                                identify the best parameters according
                                to the band-depth criterion discussed
                                in the paper.

- sea_ice_plot_forecasts.py     plots 5-months long hindcasts, starting
                                on Jun 1st, for a sample of 6 years in the
				satellite era.

- sea_ice_plot_september_averages.py  plots the september SIE average of
                                      the NSIDC data vs the RAP forecast.


- sea_ice_calc_statistics.py    computes the Bias, RMSE, and ACC of RAP
                                hindcasts for a large number of
				(target month, lead time) pairs. Saves the
				results in a .pkl file.

- sea_ice_plot_statistics.py    plots the content of the file saved by
                                sea_ice_calc_statistics.py
