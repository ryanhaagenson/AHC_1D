# Enthalpy based formulation for 1D heat condution, vertical with phase change. 
# Semi-implicit formulation
# Apparent Heat Capacity approach
# Finite Difference method with Thomas Algorithm
# Water Content and Temperature Relationship Derived from Lovell (1956)
# Water Content Profile from Van Genuchten with Stipulated Water Table Depth

# Governing Equation (linear with lag of rho*H)
# d(rho*H)/dt - d/dz(k*dT/dz) = 0

import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
import os
# from numba import jit

###############################################################################
###############################################################################

print('1D Heat Conduction with Phase Change')
print('Apparent Heat Capcity Approach')
print('Water Content and Temperature Relationship Derived from Lovell (1956)') 
print('Water Content Profile from Van Genuchten with Stipulated Water Table Depth')
print(' ')

log_file = open('./run.log','w')

lines = []
lines.append('1D Heat Conduction with Phase Change\n')
lines.append('Apparent Heat Capcity Approach\n')
lines.append('Water Content and Temperature Relationship Derived from Lovell (1956)\n') 
lines.append('Water Content Profile from Van Genuchten with Stipulated Water Table Depth\n')
lines.append('\n')
log_file.writelines(lines)

###############################################################################
###############################################################################

### THOMAS ALGORITHM ###

# Define the Thomas Algorithm to solve system with tri-diagonal matrix
def thomas_algorithm(a,b,c,r):

	J = len(a)

	e = np.zeros(J)
	f = np.zeros(J)
	g = np.zeros(J)

	# Row 1
	e[0] = a[0]
	f[0] = b[0]
	g[0] = c[0]/f[0]

	for i in range(1,J):
		e[i] = a[i]
		f[i] = b[i] - e[i]*g[i-1]
		g[i] = c[i]/f[i]

	# Forward solve for y
	y = np.zeros(J)
	y[0] = r[0]/f[0]
	for i in range(1,J):
		y[i] = (r[i] - e[i]*y[i-1])/f[i]

	# Backward solve for x
	x = np.zeros(J)
	x[J - 1] = y[J - 1]
	for i in range(J-2,-1,-1):
		x[i] = y[i] - g[i]*x[i+1]

	return x

# jit_thomas = jit(nopython=True)(thomas_algorithm)

###############################################################################
###############################################################################

### MODEL PARAMETERS ###

# Read in the inputs file
if not(os.path.isfile('./AHC_params.inpts')):
	print('No input file found.')
	print('Exiting...')
	exit()

inpts = open(r'./AHC_params.inpts',"r")
file_lines = inpts.readlines()

# Model parameters are stored here
params = {}

for file_line in file_lines:

	if not(file_line[0] == '#' or file_line[0] == '\n'):

		key = file_line.split('=')[0]
		while key[-1] == ' ':
			key = key[0:len(key)-1]

		value = file_line.split('=')[1]
		while value[0] == ' ':
			value = value[1:len(value)]
		if ',' in value:
			value = value.split(',')
			value = [float(x) for x in value]
		else:
			try:
				value = float(value)
			except:
				value = str(value[0:len(value)-1])

		params[key] = value

# Now assign values to code objects
print('---MODEL PARAMETERS---')
lines = []
lines.append('---MODEL PARAMETERS---\n')

# Spatial Grid
print('Domain:')
lines.append('Domain:\n')
if 'Model_Depth' in params.keys(): 	# Depth of model domain [m]
	D   = params.get('Model_Depth')
	print('    Depth         = ' + str(D) + ' [m]')
	lines.append('    Depth         = ' + str(D) + ' [m]\n')

else:
	D   = 50.0
	print('    Depth         = ' + str(D) + ' [m] (DEFAULT)')
	lines.append('    Depth         = ' + str(D) + ' [m] (DEFAULT)\n')

if 'Node_Spacing' in params.keys(): 	# Spacing of Model Nodes [m]
	dz   = params.get('Node_Spacing')
	print('    Node Spacing  = ' + str(dz) + ' [m]')
	lines.append('    Node Spacing  = ' + str(dz) + ' [m]\n')

else:
	dz   = 0.1
	print('    Node Spacing  = ' + str(dz) + ' [m] (DEFUALT)')
	lines.append('    Node Spacing  = ' + str(dz) + ' [m] (DEFAULT)\n')
	
print(' ')
lines.append('\n')

z = np.linspace(0.0,D,int(D/dz+1))

# Time Parameters
print('Time:')
lines.append('Time:\n')

if 'Sim_Years' in params.keys(): 	# Years in Simulation [y]
	t_end_y   = params.get('Sim_Years')
	print('    Sim Years     = ' + str(t_end_y) + ' [y]')
	lines.append('    Sim Years     = ' + str(t_end_y) + ' [y]\n')

else:
	t_end_y   = 100
	print('    Sim Years     = ' + str(t_end_y) + ' [y] (DEFUALT)')
	lines.append('    Sim Years     = ' + str(t_end_y) + ' [y] (DEFUALT)\n')

dt = 86400 			# Time Step [seconds]
t_end_d = int(t_end_y*365) 	# End of simulation [days]
t_end   = t_end_d*86400 	# End of simulation [seconds]

print('')
lines.append('\n')

# Soil Grains/Rock Matrix Properties
print('Soil Grains/Rock Matrix:')
lines.append('Soil Grains/Rock Matrix:\n')

if 'Specific_Heat_Rock' in params.keys(): 	# Specific Heat of Rock [J/(kg*degC)]
	c_r   = params.get('Specific_Heat_Rock')
	print('    Specific Heat of Rock        = ' + str(c_r) + ' [J/(kg*degC)]')
	lines.append('    Specific Heat of Rock        = ' + str(c_r) + ' [J/(kg*degC)]\n')

else:
	c_r   = 790.0
	print('    Specific Heat of Rock        = ' + str(c_r) + ' [J/(kg*degC)] (DEFUALT)')
	lines.append('    Specific Heat of Rock        = ' + str(c_r) + ' [J/(kg*degC)] (DEFAULT)\n')

if 'Density_Rock' in params.keys(): 	# Density of Rock [kg/m^3]
	rho_r   = params.get('Density_Rock')
	print('    Density of Rock              = ' + str(rho_r) + ' [kg/m^3]')
	lines.append('    Density of Rock              = ' + str(rho_r) + ' [kg/m^3]\n')

else:
	rho_r   = 2690.0
	print('    Density of Rock              = ' + str(rho_r) + ' [kg/m^3] (DEFUALT)')
	lines.append('    Density of Rock              = ' + str(rho_r) + ' [kg/m^3] (DEFAULT)\n')
	
if 'Thermal_Conductivity_Rock' in params.keys(): 	# Thermal Conductivity of Rock [J/(m*s*degC)]
	k_r   = params.get('Thermal_Conductivity_Rock')
	print('    Thermal Conductivity of Rock = ' + str(k_r) + ' [J/(m*s*degC)]')
	lines.append('    Thermal Conductivity of Rock = ' + str(k_r) + ' [J/(m*s*degC)]\n')

else:
	k_r   = 4.5
	print('    Thermal Conductivity of Rock = ' + str(k_r) + ' [J/(m*s*degC)] (DEFUALT)')
	lines.append('    Thermal Conductivity of Rock = ' + str(k_r) + ' [J/(m*s*degC)] (DEFAULT)\n')

print(' ')
lines.append('\n')

# Water, Ice and Air Properties
print('Water:')
lines.append('Water:\n')

if 'Specific_Heat_Water' in params.keys(): 	# Specific Heat of Water [J/(kg*degC)]
	c_w   = params.get('Specific_Heat_Water')
	print('    Specific Heat of Water        = ' + str(c_w) + ' [J/(kg*degC)]')
	lines.append('    Specific Heat of Water        = ' + str(c_w) + ' [J/(kg*degC)]\n')

else:
	c_w   = 4.187e3
	print('    Specific Heat of Water        = ' + str(c_w) + ' [J/(kg*degC)] (DEFUALT)')
	lines.append('    Specific Heat of Water        = ' + str(c_w) + ' [J/(kg*degC)] (DEFAULT)\n')

if 'Density_Water' in params.keys(): 	# Density of Water [kg/m^3]
	rho_w   = params.get('Density_Water')
	print('    Density of Water              = ' + str(rho_w) + ' [kg/m^3]')
	lines.append('    Density of Water              = ' + str(rho_w) + ' [kg/m^3]\n')

else:
	rho_w   = 999.8
	print('    Density of Water              = ' + str(rho_w) + ' [kg/m^3] (DEFUALT)')
	lines.append('    Density of Water              = ' + str(rho_w) + ' [kg/m^3] (DEFAULT)\n')
	
if 'Thermal_Conductivity_Water' in params.keys(): 	# Thermal Conductivity of Water [J/(m*s*degC)]
	k_w   = params.get('Thermal_Conductivity_Water')
	print('    Thermal Conductivity of Water = ' + str(k_w) + ' [J/(m*s*degC)]')
	lines.append('    Thermal Conductivity of Water = ' + str(k_w) + ' [J/(m*s*degC)]\n')

else:
	k_w   = 0.58
	print('    Thermal Conductivity of Water = ' + str(k_w) + ' [J/(m*s*degC)] (DEFUALT)')
	lines.append('    Thermal Conductivity of Water = ' + str(k_w) + ' [J/(m*s*degC)] (DEFAULT)\n')

if 'Latent_Heat_Fusion' in params.keys(): 	# Latent Heat of Fusion [J/kg]
	L   = params.get('Latent_Heat_Fusion')
	print('    Latent Heat of Fusion         = ' + str(L) + ' [J/kg]')
	lines.append('    Latent Heat of Fusion         = ' + str(L) + ' [J/kg]\n')

else:
	L   = 334e3
	print('    Latent Heat of Fusion         = ' + str(L) + ' [J/kg] (DEFUALT)')
	lines.append('    Latent Heat of Fusion         = ' + str(L) + ' [J/kg] (DEFAULT)\n')

print(' ')
lines.append('\n')

print('Ice:')
lines.append('Ice:\n')

if 'Specific_Heat_Ice' in params.keys(): 	# Specific Heat of Ice [J/(kg*degC)]
	c_i   = params.get('Specific_Heat_Ice')
	print('    Specific Heat of Ice        = ' + str(c_i) + ' [J/(kg*degC)]')
	lines.append('    Specific Heat of Ice        = ' + str(c_i) + ' [J/(kg*degC)]\n')

else:
	c_i   = 2.108e3
	print('    Specific Heat of Ice        = ' + str(c_i) + ' [J/(kg*degC)] (DEFUALT)')
	lines.append('    Specific Heat of Ice        = ' + str(c_i) + ' [J/(kg*degC)] (DEFAULT)\n')

if 'Density_Ice' in params.keys(): 	# Density of Ice [kg/m^3]
	rho_i   = params.get('Density_Ice')
	print('    Density of Ice              = ' + str(rho_i) + ' [kg/m^3]')
	lines.append('    Density of Ice              = ' + str(rho_i) + ' [kg/m^3]\n')

else:
	rho_i   = 916.8
	print('    Density of Ice              = ' + str(rho_i) + ' [kg/m^3] (DEFUALT)')
	lines.append('    Density of Ice              = ' + str(rho_i) + ' [kg/m^3] (DEFAULT)\n')
	
if 'Thermal_Conductivity_Ice' in params.keys(): 	# Thermal Conductivity of Ice [J/(m*s*degC)]
	k_i   = params.get('Thermal_Conductivity_Ice')
	print('    Thermal Conductivity of Ice = ' + str(k_i) + ' [J/(m*s*degC)]')
	lines.append('    Thermal Conductivity of Ice = ' + str(k_i) + ' [J/(m*s*degC)]\n')

else:
	k_i   = 2.18
	print('    Thermal Conductivity of Ice = ' + str(k_i) + ' [J/(m*s*degC)] (DEFUALT)')
	lines.append('    Thermal Conductivity of Ice = ' + str(k_i) + ' [J/(m*s*degC)] (DEFAULT)\n')

print(' ')
lines.append('\n')

print('Air:')
lines.append('Air:\n')

if 'Specific_Heat_Air' in params.keys(): 	# Specific Heat of Air [J/(kg*degC)]
	c_a   = params.get('Specific_Heat_Air')
	print('    Specific Heat of Air        = ' + str(c_a) + ' [J/(kg*degC)]')
	lines.append('    Specific Heat of Air        = ' + str(c_a) + ' [J/(kg*degC)]\n')

else:
	c_a   = 1.0
	print('    Specific Heat of Air        = ' + str(c_a) + ' [J/(kg*degC)] (DEFUALT)')
	lines.append('    Specific Heat of Air        = ' + str(c_a) + ' [J/(kg*degC)] (DEFAULT)\n')

if 'Density_Air' in params.keys(): 	# Density of Air [kg/m^3]
	rho_a   = params.get('Density_Air')
	print('    Density of Air              = ' + str(rho_a) + ' [kg/m^3]')
	lines.append('    Density of Air              = ' + str(rho_a) + ' [kg/m^3]\n')

else:
	rho_a   = 1.3
	print('    Density of Air              = ' + str(rho_a) + ' [kg/m^3] (DEFUALT)')
	lines.append('    Density of Air              = ' + str(rho_a) + ' [kg/m^3] (DEFAULT)\n')
	
if 'Thermal_Conductivity_Air' in params.keys(): 	# Thermal Conductivity of Air [J/(m*s*degC)]
	k_a   = params.get('Thermal_Conductivity_Air')
	print('    Thermal Conductivity of Air = ' + str(k_a) + ' [J/(m*s*degC)]')
	lines.append('    Thermal Conductivity of Air = ' + str(k_a) + ' [J/(m*s*degC)]\n')

else:
	k_a   = 0.035
	print('    Thermal Conductivity of Air = ' + str(k_a) + ' [J/(m*s*degC)] (DEFUALT)')
	lines.append('    Thermal Conductivity of Air = ' + str(k_a) + ' [J/(m*s*degC)] (DEFAULT)\n')

print(' ')
lines.append('\n')

# Plotting Output
print('Plotting Output:')
lines.append('Plotting Output:\n')
if 'Plotting_Times' in params.keys(): 	# Plotting Times [y]
	plot_years   = params.get('Plotting_Times')
	print('    Plotting Times = ' + str(plot_years) + ' [y]')
	lines.append('    Plotting Times = ' + str(plot_years) + ' [y]\n')

else:
	plot_years = []
	print('No plotting times indicated.')
	lines.append('No plotting times indicated.')

if not(plot_years == []):
	plot_times = [int(x)*365*86400 for x in plot_years]

print(' ')
lines.append('\n')

# Soil/Rock Column Properties
print('---SOIL PARAMETERS---')
lines.append('---SOIL PARAMETERS---\n')

if 'Soil_File' in params.keys(): 	# Soil Parameters File [y]
	soil_file   = params.get('Soil_File')
	print('    Soil Parameters File = ' + str(soil_file))
	lines.append('    Soil Parameters File = ' + str(soil_file) + '\n')

else:
	print('No Soil Parameters File indicated.')
	lines.append('No Soil Parameters File indicated.')
	print('Exiting...')
	exit()

Depth_data = [] 		# Depth of soil input
phi_data = []			# Porosity [-]
theta_s_data = []		# Saturated Water Content [-]
theta_r_data = []		# Residual Water Content [-]
Alpha_data = []			# Van Genuchten Alpha [1/m]
N_data = []				# Van Genuchten N [-]
a_eta_data = []			# Product of parameter a and eta (Lovell, 1956)
b_lovell_data = [] 		# b parameter (Lovell, 1956)

with open('./'+soil_file,mode = 'r') as infile:
	csvFile = csv.DictReader(infile)

	for file_line in csvFile:
		if 'Depth' in file_line.keys():
			Depth_data.append(float(file_line["Depth"]))
		if 'Porosity' in file_line.keys():
			phi_data.append(float(file_line["Porosity"]))
		if 'theta_s' in file_line.keys():
			theta_s_data.append(float(file_line["theta_s"]))
		if 'theta_r' in file_line.keys():
			theta_r_data.append(float(file_line["theta_r"]))
		if 'Alpha' in file_line.keys():
			Alpha_data.append(float(file_line["Alpha"]))
		if 'N' in file_line.keys():
			N_data.append(float(file_line["N"]))
		if 'a_eta' in file_line.keys():
			a_eta_data.append(float(file_line["a_eta"]))
		if 'b' in file_line.keys():
			b_lovell_data.append(float(file_line["b"]))

# Assign default values to soil parameters not given in soil_file
if phi_data == []:
	phi_data = 0.25*np.ones(np.shape(Depth_data))
	print('    Porosity = 0.25 [-] (DEFUALT)')
	lines.append('    Porosity = 0.25 [-] (DEFUALT)\n')
else:
	print('    Porosity from soil file.')
	lines.append('    Porosity from soil file.\n')
if theta_s_data == []:
	theta_s_data = 0.25*np.ones(np.shape(Depth_data))
	print('    Saturated Water Content = 0.25 [-] (DEFUALT)')
	lines.append('    Saturated Water Content = 0.25 [-] (DEFUALT)\n')
else:
	print('    Saturated Water Content from soil file.')
	lines.append('    Saturated Water Content from soil file.\n')
if theta_r_data == []:
	theta_r_data = 0.05*np.ones(np.shape(Depth_data))
	print('    Residual Water Content = 0.05 [-] (DEFUALT)')
	lines.append('    Residual Water Content = 0.05 [-] (DEFUALT)\n')
else:
	print('    Residual Water Content from soil file.')
	lines.append('    Residual Water Content from soil file.\n')
if Alpha_data == []:
	Alpha_data = 0.5*np.ones(np.shape(Depth_data))
	print('    Van Genuchten Alpha = 0.5 [-] (DEFUALT)')
	lines.append('    Van Genuchten Alpha = 0.5 [-] (DEFUALT)\n')
else:
	print('    Van Genuchten Alpha from soil file.')
	lines.append('    Van Genuchten Alpha from soil file.\n')
if N_data == []:
	N_data = 2.0*np.ones(np.shape(Depth_data))
	print('    Van Genuchten N = 2.0 [-] (DEFUALT)')
	lines.append('    Van Genuchten N = 2.0 [-] (DEFUALT)\n')
else:
	print('    Van Genuchten N from soil file.')
	lines.append('    Van Genuchten N from soil file.\n')
if a_eta_data == []:
	a_eta_data = 0.1*np.ones(np.shape(Depth_data))
	print('    Lovell a*eta = 0.1 [degC^-1] (DEFUALT)')
	lines.append('    Lovell a*eta = 0.1 [degC^-1] (DEFUALT)\n')
else:
	print('    Lovell a*eta from soil file.')
	lines.append('    Lovell a*eta from soil file.\n')
if b_lovell_data == []:
	b_lovell_data = 0.2*np.ones(np.shape(Depth_data))
	print('    Lovell b = 0.1 [-] (DEFUALT)')
	lines.append('    Lovell b = 0.1 [-] (DEFUALT)\n')
else:
	print('    Lovell b from soil file.')
	lines.append('    Lovell b from soil file.\n')

# Interpolate values to z grid
phi = np.interp(z,Depth_data,phi_data)
theta_s = np.interp(z,Depth_data,theta_s_data)
theta_r = np.interp(z,Depth_data,theta_r_data)
alpha = np.interp(z,Depth_data,Alpha_data)
N = np.interp(z,Depth_data,N_data)
a_eta = np.interp(z,Depth_data,a_eta_data)
b_lovell = np.interp(z,Depth_data,b_lovell_data)

# Water content profile -- Van Genuchten
if 'Water_Table_Depth' in params.keys(): 	# Water Table Depth [m]
	WT   = params.get('Water_Table_Depth')
	print('    Water Table Depth = ' + str(WT) + ' [m]')
	lines.append('    Water Table Depth = ' + str(WT) + ' [m]\n')

else:
	WT   = 0.0
	print('    Water Table Depth = ' + str(WT) + ' [m] (DEFUALT)')
	lines.append('    Water Table Depth = ' + str(WT) + ' [m] (DEFAULT)\n')

psi 	= (z - WT) 						# Suction Head [m]
eta 	= np.array(theta_s) 			# Maximum unfrozen water content
eta[psi < 0] 	= theta_r[psi < 0] + (theta_s[psi < 0] - theta_r[psi < 0])/((1 + (alpha[psi < 0]*np.abs(psi[psi < 0]))**N[psi < 0])**(1-1/N[psi < 0]))

# Levell Model Parameters
a_lovell = a_eta/eta 					# Parameter (Lovell, 1956)
T_star = -(1./a_lovell)**(-1./b_lovell) # Implied Freezing Point Depression (Lovell, 1956)

print(' ')
lines.append('\n')

###############################################################################
###############################################################################

### BOUNDARY AND INITIAL CONDITOIONS ###

# Top BC -- read in values from file
print('---BOUNDARY AND INITIAL CONDITIONS---')
lines.append('---BOUNDARY AND INITIAL CONDITIONS---\n')
if 'GST_File' in params.keys(): 	# GST File Name [y]
	GST_file   = params.get('GST_File')
	print('    GST is transient, given in ' + str(GST_file))
	lines.append('    GST is transient, given in ' + str(GST_file) + '\n')

else:
	print('Error: No GST File Name indicated.')
	lines.append('Error: No GST File Name indicated.')
	print('Exiting...')
	exit()

GST = []
Day = []
Month = []

with open('./'+GST_file,mode = 'r') as infile:
	csvFile = csv.DictReader(infile)

	for file_line in csvFile:
		GST.append(file_line["GST"])
		Day.append(file_line["Day"])
		Month.append(file_line["Month"])

def get_GST(j):
	value = float(GST[j])
	return value

# Bottom BC
if 'Geo_Heat_Grad' in params.keys(): 	# Geothermal Heat Gradient [degC/m]
	geo   = params.get('Geo_Heat_Grad')
	print('    Geothermal Heat Gradient = ' + str(geo) + ' [degC/m]')
	lines.append('    Geothermal Heat Gradient = ' + str(geo) + ' [degC/m]\n')

else:
	geo   = 0.0
	print('    Geothermal Heat Gradient = ' + str(geo) + ' [degC/m] (DEFUALT)')
	lines.append('    Geothermal Heat Gradient = ' + str(geo) + ' [degC/m] (DEFAULT)\n')

# Initial Condition
if 'Init_Temp' in params.keys(): 	# Initial Ground Temp [degC/m]
	T_init   = params.get('Init_Temp')
	print('    Initial Ground Temp      = ' + str(T_init) + ' [degC]')
	lines.append('    Initial Ground Temp      = ' + str(T_init) + ' [degC]\n')

else:
	T_init   = 0.0
	print('    Initial Ground Temp      = ' + str(T_init) + ' [degC] (DEFUALT)')
	lines.append('    Initial Ground Temp      = ' + str(T_init) + ' [degC] (DEFAULT)\n')

T_i = T_init * np.ones(np.shape(z))

print('')
lines.append('\n')

log_file.writelines(lines)

###############################################################################
###############################################################################

### USER DEFINED FUNCTIONS ###

# Apparent Heat Capacity
def ARC(T_arc):

	# Calculate the unfrozen ARC for each grid point
	term1 = eta*rho_w*c_w
	term2 = (theta_s-eta)*rho_a*c_a
	term3 = (1.0-theta_s)*rho_r*c_r
	ARC_unfrozen = term1 + term2 + term3
	
	# Calculate the frozen ARC for each grid point
	# Supress "invalid value encountered in power" error by making positive T values negative
	# All ARC_frozen values of positive T will be overwritten anyway.
	T_calc = -1.0*np.abs(T_arc)
	term1 = (theta_s-eta)*c_a*rho_a
	term2 = -a_lovell*b_lovell*(-T_calc)**(-1.0-b_lovell)*T_calc*eta*c_i*rho_i
	term3 = (eta-a_lovell*(-T_calc)**(-b_lovell)*eta)*c_i*rho_i
	term4 = a_lovell*(-T_calc)**(-b_lovell)*eta*c_w*rho_w
	term5 = a_lovell*b_lovell*(-T_calc)**(-1.0-b_lovell)*eta*(L+T_calc*c_w)*rho_w
	term6 = (1.0-theta_s)*c_r*rho_r
	ARC_frozen = term1 + term2 + term3 + term4 + term5 + term6

	# Apply the correct ARC value at each point based on T_star
	ARC = ARC_unfrozen
	ARC[T_arc < T_star] = ARC_frozen[T_arc < T_star]

	return ARC

# Water content
def theta_w(T):
	theta_w = np.array(eta)
	theta_w[T < T_star] = eta[T < T_star]*a_lovell[T < T_star]*\
	(-T[T < T_star])**(-b_lovell[T < T_star])
	return theta_w

# Thermal conductivity
def k(T):
	k = (theta_w(T)*(k_w)**(0.5) + (eta - theta_w(T))*(k_i)**(0.5) \
		+ (phi - eta)*(k_a)**(0.5) + (1 - phi)*(k_r)**(0.5))**(2.0)
	return k

def k_half(T):
	# Arithmetic Average k at half node locations
	k_T = k(T)
	k_half = (k_T[0:len(k_T)-1] + k_T[1:len(k_T)])/2.0
	return k_half

def k_above(T):
	k_half_T = k_half(T)
	# Zero is a dummy value. No j+1/2 node for j = J
	k_above = np.append(k_half_T[0:len(k_half_T)],[0])
	return k_above

def k_below(T):
	k_half_T = k_half(T)
	# Zero is dummy value. No j-1/2 node for j = 0
	k_below = np.append([0],k_half_T[0:len(k_half_T)])
	return k_below

###############################################################################
###############################################################################

# Initialize
times = np.linspace(dt,t_end,t_end_d)
T_n = T_i
plot_count = 0
year_count = 0
annual_temps = np.zeros((len(T_n),365+1))
annual_temps[:,0] = z
annual_tw = np.zeros((len(T_n),365+1))
annual_tw[:,0] = z
annual_ti = np.zeros((len(T_n),365+1))
annual_ti[:,0] = z
annual_temp_mean = np.zeros((len(T_n),int(t_end_y+1)))
annual_temp_mean[:,0] = z
annual_tw_mean = np.zeros((len(T_n),int(t_end_y+1)))
annual_tw_mean[:,0] = z
annual_ti_mean = np.zeros((len(T_n),int(t_end_y+1)))
annual_ti_mean[:,0] = z
fig, (ax1,ax2,ax3) = plt.subplots(1,3)

# Time loop
for t in times:

	julian_day = int(round(t/(86400.0*365.0) % 1 *365,0))

	# Evaluate nonlinearities at previous time step
	ARC_T = ARC(T_n)
	k_above_T = k_above(T_n)
	k_below_T = k_below(T_n)

	# Get the ground surface temperature
	T_s = get_GST(julian_day)	

	# Build tridiagonal vectors a,b and c
	a = np.hstack((np.array([0.0]),
		-k_below_T[1:len(k_below_T)-1]/dz**2.0,
		np.array([-1.0/dz])))
	b = np.hstack((np.array([1.0]),
		ARC_T[1:len(ARC_T)-1]/dt+\
		k_below_T[1:len(k_below_T)-1]/dz**2.0+\
		k_above_T[1:len(k_above_T)-1]/dz**2.0,
		np.array([1.0/dz])))
	c = np.hstack((np.array([0.0]),
		-k_above_T[1:len(k_above_T)-1]/dz**2.0,
		np.array([0.0])))
	r = np.hstack((np.array([T_s]),
		np.multiply(ARC_T[1:len(ARC_T)-1],T_n[1:len(T_n)-1])/dt,
		np.array([geo])))

	# Solve system
	T_solve = thomas_algorithm(a,b,c,r)
	# T_solve = jit_thomas(a,b,c,r)
	
	# Update last solution
	T_n = T_solve

	# Update annual temps matrix	
	annual_temps[:,julian_day+1] = T_solve
	annual_tw[:,julian_day+1] = theta_w(T_solve)
	annual_ti[:,julian_day+1] = eta - theta_w(T_solve)

	# Calculate the annual mean temperature profile
	if julian_day == 364:
		annual_temp_mean[:,year_count+1] = np.mean(annual_temps,1)
		annual_tw_mean[:,year_count+1] = np.mean(annual_tw,1)
		annual_ti_mean[:,year_count+1] = np.mean(annual_ti,1)
		year_count += 1

	if t in plot_times:
		ax1.plot(T_solve,z,label='Year = ' + str(plot_years[plot_count]))
		ax2.plot(theta_w(T_solve),z,label='Year = ' + str(plot_years[plot_count]))
		ax3.plot(eta-theta_w(T_solve),z,label='Year = ' + str(plot_years[plot_count]))
		plot_count += 1

###############################################################################
###############################################################################

# Plotting
if not(plot_years == []):

	ax1.plot(T_star,z,'k--',label='Freezing Point')

	ax1.set(ylabel=str(r"Depth [m]"))
	ax1.set(xlabel=str(r"Temperature [$^{\circ}C$]"))
	ax1.grid()
	ax1.invert_yaxis()
	ax1.set_xlim([np.min(T_solve)-5.0,np.max(T_solve)+5.0])

	ax2.set(ylabel=str(r"Depth [m]"))
	ax2.set(xlabel=str(r"$\theta_{w}$ [$^{\circ}C$]"))
	ax2.grid()
	ax2.invert_yaxis()

	ax3.set(ylabel=str(r"Depth [m]"))
	ax3.set(xlabel=str(r"$\theta_{i}$ [$^{\circ}C$]"))
	ax3.grid()
	ax3.invert_yaxis()
	ax3.legend()
	
	fig.tight_layout()

	if not(os.path.isdir('./plots')):
		os.mkdir('./plots')
	plt.savefig('./plots/output_plot.png',dpi=600)

# Output
if not(os.path.isdir('./output')):
	os.mkdir('./output')
np.savetxt('./output/annual_mean_temps.txt',annual_temp_mean,delimiter=',',header='Depth, '+', '.join(['Year ' + str(x+1) for x in range(int(t_end_y))]))
np.savetxt('./output/final_year_temps.txt',annual_temps,delimiter=',',header='Depth, '+', '.join(['Day ' + str(x+1) for x in range(365)]))
np.savetxt('./output/annual_mean_theta_w.txt',annual_tw_mean,delimiter=',',header='Depth, '+', '.join(['Year ' + str(x+1) for x in range(int(t_end_y))]))
np.savetxt('./output/final_year_theta_w.txt',annual_tw,delimiter=',',header='Depth, '+', '.join(['Day ' + str(x+1) for x in range(365)]))
np.savetxt('./output/annual_mean_theta_i.txt',annual_ti_mean,delimiter=',',header='Depth, '+', '.join(['Year ' + str(x+1) for x in range(int(t_end_y))]))
np.savetxt('./output/final_year_theta_i.txt',annual_ti,delimiter=',',header='Depth, '+', '.join(['Day ' + str(x+1) for x in range(365)]))
