This is our test site data files. Some are raw data, some are gap filled and/or disaggregated to hourly time step. 

Zhao et al. sites:
XDT: 35.7 N 94.13 E -- Also called "QT09" for AL temperature dataset
TGL: 33.07 N 91.94 E -- Also called "QT04" for AL temperature dataset
LDH: 31.82 N 91.74 E -- Also called "Ch04" for AL temperature dataset, and "QTB18" for borehole dataset
These three sites have coexisting meteorological stations, active layer temperature observations, and deep borehole temperature observations. These are very hard to find, and make them great study sites. Raw data is in Xcel files: "Active layer of ground temperature.xlsx", "Ground temperature dataset.xlsx", and "Meteorological dataset.xlsx". My gap filled data is located in folders with site names, under the "plots" and "time series" folders. Be careful, though... some of my gap filling processes aren't the most rigorous, so the data should be viewed for quality checks before using. We can meet and go over this if you'd like...

Zhao et al. only provides daily mean values in the raw output, so I disaggregated to hourly using the sub-daily variations in meteorological observations at the QOMS site from Ma et al. 

Ma et al. sites:
QOMS: 28.36 N 86.95 E
NAMORS: 30.77 N 90.98 E
These sites have meteorological stations with shallow ground temperature observations. They lack deep, borehole temperature observations. Raw data comes at an hourly timestep. Raw data is in the Xcel files titled "FLUX_<site>.xlsx", "RADM_<site_name>.xlsx", and "GRAD_<site_name>.xlsx". My gap filled datasets are in the "plots" and "time_series_Cloud_filled" folders. Again, be careful before using any of my gap filled datasets. 

CAMP sites:
Pyramid: 27.959 N 86.813 E
This site is just south of Mount Everest, but doesn't have any ground temperature observations. We may replace it with Wani et al., but I don't have the correct data for that site yet. 
