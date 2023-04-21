import rasterio as rio
import numpy as np
import matplotlib.pyplot as plt
from rasterio.plot import show
import argparse
import os



#########################
### Point of interest ###
#########################

# Read in the latitude and longitude of the point of interest
parser = argparse.ArgumentParser()
parser.add_argument('-lat',action="store",dest='lat',default=0.0)
parser.add_argument('-lon',action="store",dest='lon',default=0.0)
parser.add_argument('-model_depth',action="store",dest='model_depth',default=0.0)
args = parser.parse_args()

poi_lat = float(args.lat)
poi_lon = float(args.lon)
model_depth_m = float(args.model_depth)

###################
### File naming ###
###################

# Set directory path to the soil layer data
dir_path = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/POI_' + str(poi_lat).replace('.','_') + '_' + str(poi_lon).replace('.','_') + '/'

# Create directory for POI
if os.path.isdir(dir_path):
	os.system("rm -rf "+dir_path)
os.mkdir(dir_path)

# Find the SoilGrids tile name
lats = np.arange(21,62,2)
lat_i = np.argmin(np.abs(lats-poi_lat))
lower_lat = lats[lat_i] - 1
lons = np.arange(61,122,2)
lon_i = np.argmin(np.abs(lons-poi_lon))
lower_lon = lons[lon_i] - 1
lat_lon = str(lower_lat) + '_' + str(lower_lat+2) + '_' + str(lower_lon) + '_' + str(lower_lon+2)

# Check if SoilGrids tiles exist
if not os.path.isdir('/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/' + lat_lon):
	print "The SoilGrids tiles have not been generated for this location."
	print "Please download tiles from SoilGrids website."
	print "Exiting"
	exit()

#################
#### OPTIONS ####
################

# Plot the time series for visual inspection
show_plots = 0

# Save Plots
save_plots = 1

#######################
### Soil Properties ###
#######################

# These values are mostly from Charbeneau textbook, besides bedrock values. 
# Bedrock values are from Parajuli et al. (2017) -- Agricultural and Forest Management

# Hydraulic Conductivity [m/s]
K_silt = 6.9e-7
K_silty_loam = 1.3e-6
K_silty_clay_loam = 2.0e-7
K_silty_clay = 5.6e-8
K_clay = 5.6e-7
K_clay_loam = 7.2e-7
K_loam = 2.9e-6
K_sandy_clay_loam = 3.6e-6
K_sandy_clay = 4.1e-5
K_sandy_loam = 4.1e-5
K_clay_sand = 1.7e-7
K_silty_sand = 3.0e-7
K_sand = 8.2e-5
K_bedrock = 1.0e-10

# Saturated water content [-] (assumed equal to the porosity)
ts_silt = 0.46
ts_silty_loam = 0.45
ts_silty_clay_loam = 0.43
ts_silty_clay = 0.36
ts_clay = 0.38
ts_clay_loam = 0.41
ts_loam = 0.43
ts_sandy_clay_loam = 0.39
ts_sandy_clay = 0.38
ts_sandy_loam = 0.41
ts_clay_sand = 0.41
ts_silty_sand = 0.41
ts_sand = 0.43
ts_bedrock = 0.05

# Residual water content [-]
tr_silt = 0.034
tr_silty_loam = 0.067
tr_silty_clay_loam = 0.089
tr_silty_clay = 0.070
tr_clay = 0.048
tr_clay_loam = 0.062
tr_loam = 0.078
tr_sandy_clay_loam = 0.100
tr_sandy_clay = 0.100
tr_sandy_loam = 0.065
tr_clay_sand = 0.057
tr_silty_sand = 0.057
tr_sand = 0.045
tr_bedrock = 0.001

# van Genuchten N [-]
N_silt = 1.37
N_silty_loam = 1.41
N_silty_clay_loam = 1.23
N_silty_clay = 1.09
N_clay = 1.09
N_clay_loam = 1.31
N_loam = 1.56
N_sandy_clay_loam = 1.48
N_sandy_clay = 1.23
N_sandy_loam = 1.89
N_clay_sand = 2.28
N_silty_sand = 2.28
N_sand = 2.68
N_bedrock = 1.68

# van Genuchten alpha [1/m]
alpha_silt = 1.61
alpha_silty_loam = 2.00
alpha_silty_clay_loam = 1.00
alpha_silty_clay = 0.50
alpha_clay = 0.80
alpha_clay_loam = 1.89
alpha_loam = 3.57
alpha_sandy_clay_loam = 5.88
alpha_sandy_clay = 2.70
alpha_sandy_loam = 7.69
alpha_clay_sand = 12.34
alpha_silty_sand = 12.34
alpha_sand = 14.49
alpha_bedrock = 0.001

# Thermal properties of solids -- based on values reported by Bejan and Kraus (2003) and Dong et al. (2015)
kappa_solids = 2.8 # [W/(mK)]
C_solids = 2.0e6 # [J/m3K]

####################
### Soil Profile ###
####################

BD_profile = []
OC_profile = []
CF_profile = []
Clay_profile = []
Sand_profile = []
Silt_profile = []
K_profile = []
ts_profile = []
tr_profile = []
N_profile = []
alpha_profile = []
kappa_profile = []
C_profile = []

depths = [0,5,15,30,60,100,200] 	# [cm]

model_depth = model_depth_m * 100.0	# [cm]

# Find the lon and lat vectors
init_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/'+lat_lon+'/BD_0_5_'+lat_lon+'.tiff'
init_dataset = rio.open(init_file_name)
init = init_dataset.read(1, masked=True)
lat = np.zeros(np.shape(init)[0])
lon = np.zeros(np.shape(init)[1])
for i in range(len(lat)):
	point = init_dataset.xy(i,0)
	lat[i] = point[1]
for i in range(len(lon)):
	point = init_dataset.xy(0,i)
	lon[i] = point[0]

# Find the index of the POI
lon_dist = lon - poi_lon
lon_i = np.argmin(np.abs(lon_dist))
lat_dist = lat - poi_lat
lat_i = np.argmin(np.abs(lat_dist))

# Loop through depth profile
for i in range(len(depths)-1):

	# File name at depth
	d = str(depths[i]) + '_' + str(depths[i+1])

	# Read in all files at this depth
	BD_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/'+lat_lon+'/BD_'+d+'_'+lat_lon+'.tiff'
	BD_dataset = rio.open(BD_file_name)
	OC_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/'+lat_lon+'/OC_'+d+'_'+lat_lon+'.tiff'
	OC_dataset = rio.open(OC_file_name)
	CF_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/'+lat_lon+'/CF_'+d+'_'+lat_lon+'.tiff'
	CF_dataset = rio.open(CF_file_name)
	Clay_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/'+lat_lon+'/Clay_'+d+'_'+lat_lon+'.tiff'
	Clay_dataset = rio.open(Clay_file_name)
	Sand_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/'+lat_lon+'/Sand_'+d+'_'+lat_lon+'.tiff'
	Sand_dataset = rio.open(Sand_file_name)
	Silt_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/'+lat_lon+'/Silt_'+d+'_'+lat_lon+'.tiff'
	Silt_dataset = rio.open(Silt_file_name)

	# Convert file values according ISRIC website
	BD = 10.0*BD_dataset.read(1, masked=True)		# [kg/m^3]
	OC = OC_dataset.read(1, masked=True)/10.0 		# [kg/m^3]
	CF = CF_dataset.read(1, masked=True)/10.0 		# [%Vol]
	Clay = Clay_dataset.read(1, masked=True)/10.0 	# [%Vol] assuming uniform paritcle density
	Sand = Sand_dataset.read(1, masked=True)/10.0	# [%Vol] assuming uniform paritcle density
	Silt = Silt_dataset.read(1, masked=True)/10.0	# [%Vol] assuming uniform paritcle density

	# Find the variable values at the POI for this depth
	BD_val = BD[lat_i,lon_i]
	OC_val = OC[lat_i,lon_i]
	CF_val = CF[lat_i,lon_i]
	Clay_val = Clay[lat_i,lon_i]
	Sand_val = Sand[lat_i,lon_i]
	Silt_val = Silt[lat_i,lon_i]

	# print "BD = ",BD_val
	# print "OC = ", OC_val
	# print "CF = ", CF_val
	# print "Clay = ", Clay_val
	# print "Sand = ", Sand_val
	# print "Silt = ", Silt_val
	# print Clay_val+Sand_val+Silt_val

	# Append values at this depth
	BD_profile.append(BD_val)
	OC_profile.append(OC_val)
	CF_profile.append(CF_val)
	Clay_profile.append(Clay_val)
	Sand_profile.append(Sand_val)
	Silt_profile.append(Silt_val)

	# Classify the soil
	if np.round(Sand_val,-1) <= 10:
		if np.round(Clay_val,-1) == 0:
			# Silt
			K_profile.append(K_silt)
			ts_profile.append(ts_silt)
			tr_profile.append(tr_silt)
			N_profile.append(N_silt)
			alpha_profile.append(alpha_silt)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 10:
			# Silt
			K_profile.append(K_silt)
			ts_profile.append(ts_silt)
			tr_profile.append(tr_silt)
			N_profile.append(N_silt)
			alpha_profile.append(alpha_silt)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 20:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 30:
			# Silty Clay Loam
			K_profile.append(K_silty_clay_loam)
			ts_profile.append(ts_silty_clay_loam)
			tr_profile.append(tr_silty_clay_loam)
			N_profile.append(N_silty_clay_loam)
			alpha_profile.append(alpha_silty_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 40:
			# Silty Clay Loam
			K_profile.append(K_silty_clay_loam)
			ts_profile.append(ts_silty_clay_loam)
			tr_profile.append(tr_silty_clay_loam)
			N_profile.append(N_silty_clay_loam)
			alpha_profile.append(alpha_silty_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 50:
			# Silty Clay
			K_profile.append(K_silty_clay)
			ts_profile.append(ts_silty_clay)
			tr_profile.append(tr_silty_clay)
			N_profile.append(N_silty_clay)
			alpha_profile.append(alpha_silty_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) >= 60:
			# Clay
			K_profile.append(K_clay)
			ts_profile.append(ts_clay)
			tr_profile.append(tr_clay)
			N_profile.append(N_clay)
			alpha_profile.append(alpha_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 20:
		if np.round(Clay_val,-1) == 0:
			# Silt
			K_profile.append(K_silt)
			ts_profile.append(ts_silt)
			tr_profile.append(tr_silt)
			N_profile.append(N_silt)
			alpha_profile.append(alpha_silt)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 10:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 20:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 30:
			# Silty Clay Loam
			K_profile.append(K_silty_clay_loam)
			ts_profile.append(ts_silty_clay_loam)
			tr_profile.append(tr_silty_clay_loam)
			N_profile.append(N_silty_clay_loam)
			alpha_profile.append(alpha_silty_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 40:
			# Silty Clay Loam
			K_profile.append(K_silty_clay_loam)
			ts_profile.append(ts_silty_clay_loam)
			tr_profile.append(tr_silty_clay_loam)
			N_profile.append(N_silty_clay_loam)
			alpha_profile.append(alpha_silty_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) >= 50:
			# Clay
			K_profile.append(K_clay)
			ts_profile.append(ts_clay)
			tr_profile.append(tr_clay)
			N_profile.append(N_clay)
			alpha_profile.append(alpha_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 30:
		if np.round(Clay_val,-1) == 0:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 10:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 20:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 30:
			# Clay Loam
			K_profile.append(K_clay_loam)
			ts_profile.append(ts_clay_loam)
			tr_profile.append(tr_clay_loam)
			N_profile.append(N_clay_loam)
			alpha_profile.append(alpha_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 40:
			# Clay Loam
			K_profile.append(K_clay_loam)
			ts_profile.append(ts_clay_loam)
			tr_profile.append(tr_clay_loam)
			N_profile.append(N_clay_loam)
			alpha_profile.append(alpha_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) >= 50:
			# Clay
			K_profile.append(K_clay)
			ts_profile.append(ts_clay)
			tr_profile.append(tr_clay)
			N_profile.append(N_clay)
			alpha_profile.append(alpha_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 40:
		if np.round(Clay_val,-1) == 0:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 10:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 20:
			# Loam
			K_profile.append(K_loam)
			ts_profile.append(ts_loam)
			tr_profile.append(tr_loam)
			N_profile.append(N_loam)
			alpha_profile.append(alpha_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 30:
			# Clay Loam
			K_profile.append(K_clay_loam)
			ts_profile.append(ts_clay_loam)
			tr_profile.append(tr_clay_loam)
			N_profile.append(N_clay_loam)
			alpha_profile.append(alpha_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 40:
			# Clay Loam
			K_profile.append(K_clay_loam)
			ts_profile.append(ts_clay_loam)
			tr_profile.append(tr_clay_loam)
			N_profile.append(N_clay_loam)
			alpha_profile.append(alpha_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) >= 50:
			# Clay
			K_profile.append(K_clay)
			ts_profile.append(ts_clay)
			tr_profile.append(tr_clay)
			N_profile.append(N_clay)
			alpha_profile.append(alpha_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 50:
		if np.round(Clay_val,-1) == 0:
			# Silty Loam
			K_profile.append(K_silty_loam)
			ts_profile.append(ts_silty_loam)
			tr_profile.append(tr_silty_loam)
			N_profile.append(N_silty_loam)
			alpha_profile.append(alpha_silty_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 10:
			# Loam
			K_profile.append(K_loam)
			ts_profile.append(ts_loam)
			tr_profile.append(tr_loam)
			N_profile.append(N_loam)
			alpha_profile.append(alpha_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 20:
			# Loam
			K_profile.append(K_loam)
			ts_profile.append(ts_loam)
			tr_profile.append(tr_loam)
			N_profile.append(N_loam)
			alpha_profile.append(alpha_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 30:
			# Sandy Clay Loam
			K_profile.append(K_sandy_clay_loam)
			ts_profile.append(ts_sandy_clay_loam)
			tr_profile.append(tr_sandy_clay_loam)
			N_profile.append(N_sandy_clay_loam)
			alpha_profile.append(alpha_sandy_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 40:
			# Sandy Clay
			K_profile.append(K_sandy_clay)
			ts_profile.append(ts_sandy_clay)
			tr_profile.append(tr_sandy_clay)
			N_profile.append(N_sandy_clay)
			alpha_profile.append(alpha_sandy_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) >= 50:
			# Sandy Clay
			K_profile.append(K_sandy_clay)
			ts_profile.append(ts_sandy_clay)
			tr_profile.append(tr_sandy_clay)
			N_profile.append(N_sandy_clay)
			alpha_profile.append(alpha_sandy_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 60:
		if np.round(Clay_val,-1) <= 20:
			# Sandy Loam
			K_profile.append(K_sandy_loam)
			ts_profile.append(ts_sandy_loam)
			tr_profile.append(tr_sandy_loam)
			N_profile.append(N_sandy_loam)
			alpha_profile.append(alpha_sandy_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 30:
			# Sandy Clay Loam
			K_profile.append(K_sandy_clay_loam)
			ts_profile.append(ts_sandy_clay_loam)
			tr_profile.append(tr_sandy_clay_loam)
			N_profile.append(N_sandy_clay_loam)
			alpha_profile.append(alpha_sandy_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 40:
			# Sandy Clay
			K_profile.append(K_sandy_clay)
			ts_profile.append(ts_sandy_clay)
			tr_profile.append(tr_sandy_clay)
			N_profile.append(N_sandy_clay)
			alpha_profile.append(alpha_sandy_clay)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 70:
		if np.round(Clay_val,-1) <= 20:
			# Sandy Loam
			K_profile.append(K_sandy_loam)
			ts_profile.append(ts_sandy_loam)
			tr_profile.append(tr_sandy_loam)
			N_profile.append(N_sandy_loam)
			alpha_profile.append(alpha_sandy_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 30:
			# Sandy Clay Loam
			K_profile.append(K_sandy_clay_loam)
			ts_profile.append(ts_sandy_clay_loam)
			tr_profile.append(tr_sandy_clay_loam)
			N_profile.append(N_sandy_clay_loam)
			alpha_profile.append(alpha_sandy_clay_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 80:
		if np.round(Clay_val,-1) == 0:
			# Silty Sand
			K_profile.append(K_silty_sand)
			ts_profile.append(ts_silty_sand)
			tr_profile.append(tr_silty_sand)
			N_profile.append(N_silty_sand)
			alpha_profile.append(alpha_silty_sand)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 10:
			# Clay Sand
			K_profile.append(K_clay_sand)
			ts_profile.append(ts_clay_sand)
			tr_profile.append(tr_clay_sand)
			N_profile.append(N_clay_sand)
			alpha_profile.append(alpha_clay_sand)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
		if np.round(Clay_val,-1) == 20:
			# Sandy Loam
			K_profile.append(K_sandy_loam)
			ts_profile.append(ts_sandy_loam)
			tr_profile.append(tr_sandy_loam)
			N_profile.append(N_sandy_loam)
			alpha_profile.append(alpha_sandy_loam)
			kappa_profile.append(kappa_solids)
			C_profile.append(C_solids)
	if np.round(Sand_val,-1) == 90:
		# Sand
		K_profile.append(K_sand)
		ts_profile.append(ts_sand)
		tr_profile.append(tr_sand)
		N_profile.append(N_sand)
		alpha_profile.append(alpha_sand)
		kappa_profile.append(kappa_solids)
		C_profile.append(C_solids)

# Update parameter values in top soil for coarse fragements
for i in range(len(K_profile)):

	K_profile[i] = K_profile[i] * (1.0 - CF_profile[i]/100.0)
	ts_profile[i] = ts_profile[i] * (1.0 - CF_profile[i]/100.0)
	tr_profile[i] = tr_profile[i] * (1.0 - CF_profile[i]/100.0)

# Update parameter values in top soil for organic content
for i in range(len(K_profile)):

	f_oc = OC_profile[i]/130.0
	if f_oc > 1.0:
		print "The organic content fraction is greater than one!"
		exit()
	K_profile[i] = (1.0 - f_oc)*K_profile[i] + f_oc*2.8e-4
	ts_profile[i] = (1.0 - f_oc)*ts_profile[i] + f_oc*0.9
	kappa_profile[i] = (1.0 - f_oc)*kappa_profile[i] + f_oc*0.25
	C_profile[i] = (1.0 - f_oc)*C_profile[i] + f_oc*2.5e6

# Trim off 0 from depths array
depths = depths[1:len(depths)+1]

# Read in the depth to bedrock
bedrock_file_name = '/Volumes/SeagateBUP/Research/HMA_Permafrost/Sat_data/SoilGrids/Pelletier_bedrock/Global_Soil_Regolith_Sediment_1304/data/average_soil_and_sedimentary-deposit_thickness.tif'
bedrock_dataset = rio.open(bedrock_file_name)
bedrock = bedrock_dataset.read(1, masked=True)

# Find the lon and lat vectors
lat = np.zeros(np.shape(bedrock)[0])
lon = np.zeros(np.shape(bedrock)[1])
for i in range(len(lat)):
	point = bedrock_dataset.xy(i,0)
	lat[i] = point[1]
for i in range(len(lon)):
	point = bedrock_dataset.xy(0,i)
	lon[i] = point[0]

# Find the index of the POI
lon_dist = lon - poi_lon
lon_i = np.argmin(np.abs(lon_dist))
lat_dist = lat - poi_lat
lat_i = np.argmin(np.abs(lat_dist))

# Find the depth to bedrock at the POI
bedrock_val = bedrock[lat_i,lon_i]

if bedrock_val == -1:
	print "No value for depth to bedrock."
	print "Assuming depth to bedrock is 3 meters."
	bedrock_val = 3.0

# Give a buffer if depth to bedrock is shallow
if bedrock_val < 3.0:
	bedrock_val = 3.0

# Add point to profile from 2 meter depth down to bedrock
num_points = (np.round(100.0*bedrock_val,0) - depths[-1])/100.0

depths = [float(x) for x in depths]

for i in range(int(num_points)):

	depths.append((i+1)*100.0+200.0)

# Solve for the exponential parameters
K_beta = np.log(K_bedrock/K_profile[-1])/(np.round(100.0*bedrock_val,0)-200.0)
K_alpha = K_profile[-1]/(np.exp(200.0*K_beta))
ts_beta = np.log(ts_bedrock/ts_profile[-1])/(np.round(100.0*bedrock_val,0)-200.0)
ts_alpha = ts_profile[-1]/(np.exp(200.0*ts_beta))
tr_beta = np.log(tr_bedrock/tr_profile[-1])/(np.round(100.0*bedrock_val,0)-200.0)
tr_alpha = tr_profile[-1]/(np.exp(200.0*tr_beta))
N_beta = np.log(N_bedrock/N_profile[-1])/(np.round(100.0*bedrock_val,0)-200.0)
N_alpha = N_profile[-1]/(np.exp(200.0*N_beta))
alpha_beta = np.log(alpha_bedrock/alpha_profile[-1])/(np.round(100.0*bedrock_val,0)-200.0)
alpha_alpha = alpha_profile[-1]/(np.exp(200.0*alpha_beta))

# Append values at each new depth
for i in range(int(num_points)):

	d = depths[i+6]
	K_profile.append(K_alpha*np.exp(K_beta*d))
	ts_profile.append(ts_alpha*np.exp(ts_beta*d))
	tr_profile.append(tr_alpha*np.exp(tr_beta*d))
	N_profile.append(N_alpha*np.exp(N_beta*d))
	alpha_profile.append(alpha_alpha*np.exp(alpha_beta*d))
	kappa_profile.append(kappa_solids)
	C_profile.append(C_solids)

# Append bedroock values down to the model floor
num_points = (model_depth - np.round(100.0*bedrock_val,0))/100.0

for i in range(int(num_points)):

	depths.append((i+1)*100.0 + np.round(100.0*bedrock_val,0))
	K_profile.append(K_bedrock)
	ts_profile.append(ts_bedrock)
	tr_profile.append(tr_bedrock)
	N_profile.append(N_bedrock)
	alpha_profile.append(alpha_bedrock)
	kappa_profile.append(kappa_solids)
	C_profile.append(C_solids)

################
### Plotting ###
################

fig = plt.figure()
ax = fig.add_subplot(111)
plt.semilogx(K_profile,depths)
plt.gca().invert_yaxis()
plt.xlabel('K [m/s]')
plt.ylabel('Depth [cm]')
plt.grid()
a = plt.axes([0.6,0.2,0.2,0.45])
plt.plot(K_profile[0:6],depths[0:6])
plt.gca().invert_yaxis()
plt.yticks([5,15,30,60,100,200],fontsize=8)
plt.xticks([np.min(K_profile[0:6]),np.max(K_profile[0:6])],fontsize=8)
plt.grid()
if show_plots == 1:
	plt.show()
if save_plots == 1:
	fig.savefig(dir_path+'K_'+str(poi_lat).replace('.','_')+'_'+str(poi_lon).replace('.','_')+'.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.semilogx(ts_profile,depths)
plt.gca().invert_yaxis()
plt.xlabel('Saturated Water Content [-]')
plt.ylabel('Depth [cm]')
plt.grid()
a = plt.axes([0.6,0.2,0.2,0.45])
plt.plot(ts_profile[0:6],depths[0:6])
plt.gca().invert_yaxis()
plt.yticks([5,15,30,60,100,200],fontsize=8)
plt.xticks([np.min(ts_profile[0:6]),np.max(ts_profile[0:6])],fontsize=8)
plt.grid()
if show_plots == 1:
	plt.show()
if save_plots == 1:
	fig.savefig(dir_path+'ts_'+str(poi_lat).replace('.','_')+'_'+str(poi_lon).replace('.','_')+'.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.semilogx(tr_profile,depths)
plt.gca().invert_yaxis()
plt.xlabel('Residual Water Content [-]')
plt.ylabel('Depth [cm]')
plt.grid()
a = plt.axes([0.6,0.2,0.2,0.45])
plt.plot(tr_profile[0:6],depths[0:6])
plt.gca().invert_yaxis()
plt.yticks([5,15,30,60,100,200],fontsize=8)
plt.xticks([np.min(tr_profile[0:6]),np.max(tr_profile[0:6])],fontsize=8)
plt.grid()
if show_plots == 1:
	plt.show()
if save_plots == 1:
	fig.savefig(dir_path+'tr_'+str(poi_lat).replace('.','_')+'_'+str(poi_lon).replace('.','_')+'.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(N_profile,depths)
plt.gca().invert_yaxis()
plt.xlabel(r'van Genuchten $N$ [-]')
plt.ylabel('Depth [cm]')
plt.grid()
a = plt.axes([0.6,0.2,0.2,0.45])
plt.plot(N_profile[0:6],depths[0:6])
plt.gca().invert_yaxis()
plt.yticks([5,15,30,60,100,200],fontsize=8)
plt.xticks([np.min(N_profile[0:6]),np.max(N_profile[0:6])],fontsize=8)
plt.grid()
if show_plots == 1:
	plt.show()
if save_plots == 1:
	fig.savefig(dir_path+'N_'+str(poi_lat).replace('.','_')+'_'+str(poi_lon).replace('.','_')+'.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.semilogx(alpha_profile,depths)
plt.gca().invert_yaxis()
plt.xlabel(r'van Genuchten $\alpha$ [m/s]')
plt.ylabel('Depth [cm]')
plt.grid()
a = plt.axes([0.6,0.2,0.2,0.45])
plt.plot(alpha_profile[0:6],depths[0:6])
plt.gca().invert_yaxis()
plt.yticks([5,15,30,60,100,200],fontsize=8)
plt.xticks([np.min(alpha_profile[0:6]),np.max(alpha_profile[0:6])],fontsize=8)
plt.grid()
if show_plots == 1:
	plt.show()
if save_plots == 1:
	fig.savefig(dir_path+'alpha_'+str(poi_lat).replace('.','_')+'_'+str(poi_lon).replace('.','_')+'.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(kappa_profile,depths)
plt.gca().invert_yaxis()
plt.xlabel('Thermal Conductivity [W/(m K)]')
plt.ylabel('Depth [cm]')
plt.grid()
a = plt.axes([0.6,0.2,0.2,0.45])
plt.plot(kappa_profile[0:6],depths[0:6])
plt.gca().invert_yaxis()
plt.yticks([5,15,30,60,100,200],fontsize=8)
plt.xticks([np.min(kappa_profile[0:6]),np.max(kappa_profile[0:6])],fontsize=8)
plt.grid()
if show_plots == 1:
	plt.show()
if save_plots == 1:
	fig.savefig(dir_path+'kappa_'+str(poi_lat).replace('.','_')+'_'+str(poi_lon).replace('.','_')+'.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(C_profile,depths)
plt.gca().invert_yaxis()
plt.xlabel('Heat Capacity [J/(m3 K)]')
plt.ylabel('Depth [cm]')
plt.grid()
a = plt.axes([0.6,0.2,0.2,0.45])
plt.plot(C_profile[0:6],depths[0:6])
plt.gca().invert_yaxis()
plt.yticks([5,15,30,60,100,200],fontsize=8)
plt.xticks([np.min(C_profile[0:6]),np.max(C_profile[0:6])],fontsize=8)
plt.grid()
if show_plots == 1:
	plt.show()
if save_plots == 1:
	fig.savefig(dir_path+'C_'+str(poi_lat).replace('.','_')+'_'+str(poi_lon).replace('.','_')+'.png')

##############
### Saving ###
##############

# Convert depths from [cm] to [mm] before saving
depths = [10.0 * x for x in depths]
# Convert K from [m/s] to [mm/s] before saving
K_profile = [1000.0 * x for x in K_profile]
# Convert alpha from [m^-1] to [mm^-1] befor saving
alpha_profile = [x / 1000.0 for x in alpha_profile]

# Save profiles
np.savetxt(dir_path+'depths.txt',depths)
np.savetxt(dir_path+'K_profile.txt',K_profile)
np.savetxt(dir_path+'ts_profile.txt',ts_profile)
np.savetxt(dir_path+'tr_profile.txt',tr_profile)
np.savetxt(dir_path+'alpha_profile.txt',alpha_profile)
np.savetxt(dir_path+'N_profile.txt',N_profile)
np.savetxt(dir_path+'kappa_profile.txt',kappa_profile)
np.savetxt(dir_path+'C_profile.txt',C_profile)



