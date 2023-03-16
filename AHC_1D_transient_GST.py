# Enthalpy based formulation for 1D heat condution, vertical with phase change. 
# Semi-implicit formulation
# Essentially the Apparent Heat Capacity approach
# Finite Difference method with Thomas Algorithm

# Governing Equation (linear with lag of rho*H)
# d(rho*H)/dt - d/dz(k*dT/dz) = 0

# Numerical approach outlined in Hari's notes
# Assumptions: 
#     1. Full saturation throughout soil
#     2. Full transition between frozen/unfrozen water at 0 degC
#     3. Soil has a single porosity value

import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
import os
# from numba import jit

###############################################################################
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-q','--quiet',action="store",dest='quiet',default='n',help = 'Run in quiet mode [y/n]')
parser.add_argument('-GST','--GST_file',action="store",dest='GST_file',default='',help = 'GST file name')
parser.add_argument('-y','--sim_years',action="store",dest='sim_years',default='0',help = 'Number of simulation years')
parser.add_argument('-p','--plot_years',action="store",dest='plot_years',nargs='+',default=[],help='Years to plot solution.')
args = parser.parse_args()

quiet = str(args.quiet)
GST_file = str(args.GST_file)
sim_years = int(args.sim_years)
plot_years = list(args.plot_years)

if GST_file == '':
	print('Error: Define the GST_file.')
	print('Type python AHC_1D_transient.py -h for help.')
	print('Exiting...')
	exit()

if sim_years == 0:
	print('Error: Define number of simulation years.')
	print('Type python AHC_1D_transient.py -h for help.')
	print('Exiting...')
	exit()

###############################################################################
###############################################################################

if quiet == 'n':
	print('1D Heat Conduction with Phase Change')
	print('Apparent Heat Capcity Approach')
	print('Methodology from Hari Notes')
	print('')

log_file = open('./run.log','w')

lines = []
lines.append('1D Heat Conduction with Phase Change\n')
lines.append('Apparent Heat Capcity Approach\n')
lines.append('Methodology from Hari Notes\n')
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

# Rock Properties
phi		= 0.05		# Porosity [-]
c_r		= 790		# Specific Heat of Rock [J/(kg*degC)]
rho_r 	= 2690		# Rock Density [kg/m^3]
k_r 	= 4.5 		# Thermal Conductivity of Rock [J/(m*s*degC)]

# Water and Ice Properties
c_w 	= 4.187e3	# Specific Heat of Water [J/(kg*degC)]
c_i 	= 2.108e3	# Specific Heat of Ice [J/(kg*degC)]
rho_w 	= 999.8 	# Water density [kg/m^3]
rho_i 	= 916.8 	# Ice density [kg/m^3]
k_w 	= 0.58 		# Thermal Conductivity of Water [J/(m*s*degC)]
k_i 	= 2.18 		# Thermal Conductivity of Ice [J/(m*s*degC)]
L 		= 334e3 	# Latent Heat with Phase Change [J/kg]
B 		= 100		# Dummy model parameter

# Spatial Grid
D 		= 50 		# Depth of model domain [m]
dz 		= 0.1		# Node spacing [m]
z = np.linspace(0.0,D,int(D/dz+1))

# Time parameters
t_end_y = sim_years 	# End of simulatino [years]
t_end_d = t_end_y*365 	# End of simulation [days]
dt 		= 86400			# Time step [seconds]

t_end   = t_end_d*86400
times   = np.linspace(dt,t_end,t_end_d)
plot_times = [int(x)*365*86400 for x in plot_years]

if quiet == 'n':
	print('---MODEL PARAMETERS---')
	print('Soil Grains/Rock Matrix:')
	print('    porosity      = ' + str(phi) + ' [-]')
	print('    specific heat = ' + str(c_r) + ' [J/(kg*degC)]')
	print('    density       = ' + str(rho_r) + ' [kg/m^3]')
	print('    thermal cond. = ' + str(k_r) + ' [J/(m*s*degC)]')
	print('')

	print('Water:')
	print('    specific heat = ' + str(c_w) + ' [J/(kg*degC)]')
	print('    density       = ' + str(rho_w) + ' [kg/m^3]')
	print('    thermal cond. = ' + str(k_w) + ' [J/(m*s*degC)]')
	print('')

	print('Ice:')
	print('    specific heat = ' + str(c_i) + ' [J/(kg*degC)]')
	print('    density       = ' + str(rho_i) + ' [kg/m^3]')
	print('    thermal cond. = ' + str(k_i) + ' [J/(m*s*degC)]')
	print('')

	print('Domain:')
	print('    depth         = ' + str(D) + ' [m]')
	print('    node spacing  = ' + str(dz) + ' [m]')
	print('')

	print('Time:')
	print('    sim years     = ' + str(sim_years) + ' [y]')
	print('    time step     = ' + str(dt/86400.0) + ' [d]')
	print('')

lines = []
lines.append('---MODEL PARAMETERS---\n')
lines.append('Soil Grains/Rock Matrix:\n')
lines.append('    porosity      = ' + str(phi) + ' [-]\n')
lines.append('    specific heat = ' + str(c_r) + ' [J/(kg*degC)]\n')
lines.append('    density       = ' + str(rho_r) + ' [kg/m^3]\n')
lines.append('    thermal cond. = ' + str(k_r) + ' [J/(m*s*degC)]\n')
lines.append('\n')
lines.append('Water:\n')
lines.append('    specific heat = ' + str(c_w) + ' [J/(kg*degC)]\n')
lines.append('    density       = ' + str(rho_w) + ' [kg/m^3]\n')
lines.append('    thermal cond. = ' + str(k_w) + ' [J/(m*s*degC)]\n')
lines.append('\n')
lines.append('Ice:\n')
lines.append('    specific heat = ' + str(c_i) + ' [J/(kg*degC)]\n')
lines.append('    density       = ' + str(rho_i) + ' [kg/m^3]\n')
lines.append('    thermal cond. = ' + str(k_i) + ' [J/(m*s*degC)]\n')
lines.append('\n')
lines.append('Domain:\n')
lines.append('    depth         = ' + str(D) + ' [m]\n')
lines.append('    node spacing  = ' + str(dz) + ' [m]\n')
lines.append('\n')
lines.append('Time:\n')
lines.append('    sim years     = ' + str(sim_years) + ' [y]\n')
lines.append('    time step     = ' + str(dt) + ' [s]\n')
lines.append('\n')
log_file.writelines(lines)

###############################################################################
###############################################################################

### BOUNDARY AND INITIAL CONDITOIONS ###

# Top BC -- read in values from file
GST = []
Day = []
Month = []
with open("./GST_test.txt",mode = 'r') as infile:
	csvFile = csv.DictReader(infile)

	for lines in csvFile:
		GST.append(lines["GST"])
		Day.append(lines["Day"])
		Month.append(lines["Month"])

def get_GST(j):
	value = float(GST[j])
	return value

# Bottom BC
geo = 0.025		# geothermal gradient [degC/m]

# Initial Condition

T_init = -5.0
T_i = T_init * np.ones(np.shape(z))

if quiet == 'n':
	print('---BOUNDARY AND INITIAL CONDITIONS---')
	print('GST is transient, given in ' + str(GST_file))
	print('Geothermal Heat flux = ' + str(geo) + ' [degC/m]')
	print('Initial Temp = ' + str(T_init) + ' [degC]')
	print('')

lines = []
lines.append('---BOUNDARY AND INITIAL CONDITIONS---\n')
lines.append('GST is transient, given in ' + str(GST_file) + '\n')
lines.append('Geothermal Heat flux = ' + str(geo) + ' [degC/m]\n')
lines.append('Initial Temp = ' + str(T_init) + ' [degC]\n')
lines.append('\n')
log_file.writelines(lines)

###############################################################################
###############################################################################

### USER DEFINED FUNCTIONS ###

# Apparent Heat Capacity
def ARC(T):
	# d(pho*H)/dT
	term1=np.ones(len(T))*(1.0-phi)*rho_r*c_r
	term2=phi*rho_i*c_i*((-B/2.0)*T*(1.0/np.cosh(B*T))**2.0 \
		+ (0.5-0.5*np.tanh(B*T)))
	term3=phi*rho_w*((L+c_w*T)*(B/2.0)*(1.0/np.cosh(B*T))**2.0 \
		+ c_w*(0.5+0.5*np.tanh(B*T)))
	ARC=term1+term2+term3
	return ARC

# Thermal conductivity
def k(T):
	k = (1-phi)*k_r + phi*((0.5 - 0.5*np.tanh(B*T))*k_i + \
		(0.5 + 0.5*np.tanh(B*T))*k_w)
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
T_n = T_i
plot_count = 0
year_count = 0
annual_temps = np.zeros((len(T_n),365+1))
annual_temps[:,0] = z
annual_means = np.zeros((len(T_n),sim_years+1))
annual_means[:,0] = z

# Time loop
for t in times:

	julian_day = int(round(t/(86400.0*365.0) % 1 *365,0))

	# print("###################")
	# print("Current time =", t)
	# print("")

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

	# Calculate the annual mean temperature profile
	if julian_day == 364:
		annual_means[:,year_count+1] = np.mean(annual_temps,1)
		year_count += 1

	if t in plot_times:
		plt.plot(T_solve,z,label='Year = ' + str(plot_years[plot_count]))
		plot_count += 1

###############################################################################
###############################################################################

# Plotting
if not(plot_years == []):

	plt.ylabel(r"Depth [m]")
	plt.xlabel(r"Temperature [$^{\circ}C$]")
	plt.grid()
	plt.gca().invert_yaxis()
	plt.legend()
	if not(os.path.isdir('./plots')):
		os.mkdir('./plots')
	plt.savefig('./plots/output_plot.png',dpi=600)

# Output
if not(os.path.isdir('./output')):
	os.mkdir('./output')
np.savetxt('./output/annual_mean_temps.txt',annual_means,delimiter=',',header='Depth, '+', '.join(['Year ' + str(x+1) for x in range(sim_years)]))
np.savetxt('./output/final_year_temps.txt',annual_temps,delimiter=',',header='Depth, '+', '.join(['Day ' + str(x+1) for x in range(365)]))
