# Enthalpy based formulation for 1D heat condution, vertical with phase change. 
# Semi-implicit formulation
# Essentially the Apparent Heat Capacity approach
# Finite Difference method with Thomas Algorithm

# Governing Equation (linear with lag of rho*H)
# d(rho*H)/dt - d/dz(k*dT/dz) = 0

# Example problem from Hari's notes
# Initially geohermal (linear) gradient with T_top = -1 degC
# Surface temp suddenly increases to 0.5 degC

import numpy as np
import matplotlib.pyplot as plt
# from numba import jit
# from datetime import datetime

# start_time = datetime.now()

###############################################################################
###############################################################################

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

# Model Parameters

# Rock Properties
phi		= 0.02		# Porosity [-]
c_r		= 790		# Specific Heat of Rock [J/(kg*degC)]
rho_r 	= 2690		# Rock Density [kg/m^3]
k_r 	= 2.5 		# Thermal Conductivity of Rock [J/(m*s*degC)]

# Water and Ice Properties
c_w 	= 4.187e3	# Specific Heat of Water [J/(kg*degC)]
c_i 	= 2.108e3	# Specific Heat of Ice [J/(kg*degC)]
k_w 	= 0.58 		# Thermal Conductivity of Water [J/(m*s*degC)]
k_i 	= 2.18 		# Thermal Conductivity of Ice [J/(m*s*degC)]
L 		= 334e3 	# Latent Heat with Phase Change [J/kg]
rho_w 	= 999.8 	# Water density [kg/m^3]
rho_i 	= 916.8 	# Ice density [kg/m^3]
B 		= 100		# Dummy model parameter

# Spatial Grid
D 		= 50 		# Depth of model domain [m]
dz 		= 0.1		# Node spacing [m]
z = np.linspace(0.0,D,int(D/dz+1))

# Time parameters
t_end_d = 24000 	# End of simulation [days]
dt 		= 86400		# Time step [seconds]
t_end = t_end_d*86400
times = np.linspace(dt,t_end,t_end_d)

###############################################################################
###############################################################################

# Boundary Conditions

# Top BC
T_s = 1.0 		# Surface temperature [degC]
# Bottom BC
geo = 0.025		# geothermal gradient [degC/m]

# Initial Condition

T_init = -1.0
T_i = T_init * np.ones(np.shape(z))

###############################################################################
###############################################################################

# Define the apparent heat capacity and thermal conductivity functions

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

plot_times = np.array([150,300,600,1500,3000,6000,12000,24000,48000])*86400

# Time loop
for t in times:

	# print("###################")
	# print("Current time =", t)
	# print("")

	# Evaluate nonlinearities at previous time step
	ARC_T = ARC(T_n)
	k_above_T = k_above(T_n)
	k_below_T = k_below(T_n)

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

	if t in plot_times:
		print(t/86400.)
		plt.plot(T_solve,z)
		# print(np.min(T_solve))

plt.ylabel(r"Depth [m]")
plt.xlabel(r"Temperature [$^{\circ}C$]")
plt.grid()
plt.gca().invert_yaxis()
plt.show()

