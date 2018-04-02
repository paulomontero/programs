import numpy as np
import scipy as sp
import sys
from scipy.integrate import quad
from scipy.interpolate import interp2d

#cosmoparameter (adjust as needed for diferrent realizations)
omega_m = 0.3088
omega_r = 0.0 # this what it would be 9.134000000e-5, however at this redshift is irrelevant
omega_l = 1 - omega_m 
h = 0.6774

# factor of Hubble units -> km/Mpc/s
def hubble(z):
	return (h*100)*np.sqrt(omega_m*(1 + z)**3 + omega_r*(1 + z)**4 + omega_l)	

print hubble(0)
print hubble(1)
print hubble(2)
print hubble(10)

def integrand(x):
	return (1+x)/(hubble(x))**3

def Dratio(z_obs,z):
	D_obs = quad(integrand, z_obs, np.inf)
	print D_obs
	D = quad(integrand, z, np.inf)
	print D
	return D_obs[0]/D[0]

print Dratio(1, 21)

table = np.loadtxt('power_default1')

p_mxh   = table[:,2]
print p_mxh

#this does not work because I am doing it wrong, need for every z the k, intervals of 16 fixed
k_array = [2.094395e-02, 3.183814e-02, 4.770840e-02, 7.025257e-02, 9.773653e-02, 1.368335e-01, 1.929536e-01, 2.695639e-01, 3.777252e-01, 5.288205e-01, 7.401737e-01, 1.036455e+00, 1.451071e+00, 2.031429e+00, 2.784068e+00, 3.580407e+00]
print k_array[0]



z_vec = np.loadtxt('zarray.dat')
z_array = z_vec[::-1]
print z_array[0]
print z_array[83]

crosspower = interp2d(z_array, k_array, p_mxh, kind='linear')

crosspower_mxh = crosspower(7.9, 0.0997)
print crosspower_mxh


