import numpy as np
import scipy as sp
import sys
from scipy.integrate import quad
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
from scipy import optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import interpolate

#cosmoparameter (adjust as needed for diferrent realizations)
omega_m = 0.3088
omega_r = 0.0 # this what it would be 9.134000000e-5, however at the redshift involved this is irrelevant
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
	#print D_obs
	D = quad(integrand, z, np.inf)
	#print D
	# multiplying by prefactor
	return (D_obs[0]*hubble(z_obs))/(D[0]*hubble(z))

print Dratio(1, 21)

#load 21cmfast realization 
table = np.loadtxt('power_pop3_1')
#table = np.loadtxt('power_pop3_2')

p_mxh   = table[:,2]
print p_mxh

#this does not work because I am doing it wrong, need for every z the k, intervals of 16 fixed
k_array = [2.094395e-02, 3.183814e-02, 4.770840e-02, 7.025257e-02, 9.773653e-02, 1.368335e-01, 1.929536e-01, 2.695639e-01, 3.777252e-01, 5.288205e-01, 7.401737e-01, 1.036455e+00, 1.451071e+00, 2.031429e+00, 2.784068e+00, 3.580407e+00]
print k_array[0]
# remember that there is no h anywhere in the ouputs of 21cmfast



z_vec = np.loadtxt('zarray.dat')
#inverting array
z_array = z_vec[::-1]
print z_array[0]
print z_array[83]

crosspower = interp2d(z_array, k_array, p_mxh, kind='cubic')

crosspower_mxh = crosspower(7.9, 0.0997)
print crosspower_mxh

#PLOTTING STARTED
#important: z1 = 2.0, z2 =2.5, z3=3.0, z4=3.5, z5=4.0
#time to deal with the bias representative, do a piecewise linear guy.
#start with data from zobs = 2.5
x = np.array([6, 7, 8, 9, 10, 11, 12], dtype = float)
y_z2 = np.array([0.05903, 0.01993, 0, -0.002998, 0.001291, 0.006988, 0.009932], dtype = float)

tck_z2 = interpolate.splrep(x, y_z2, k=1, s=0)
xnew = np.linspace(5.9,34.7)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z2, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z2, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z2, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

#plt.show()

#next the others dpsi  guys, zobs = 2.0
y_z1 = np.array([0.05621, 0.01945, 0, -0.008235, -0.01222, -0.01311, -0.01468], dtype = float)

tck_z1 = interpolate.splrep(x, y_z1, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z1, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z1, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z1, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

#plt.show()

#next the others dpsi  guys zobs = 3.0
y_z3 = np.array([0.08588, 0.03052, 0, -0.006657, -0.001757, 0.005913, 0.01350], dtype = float)

tck_z3 = interpolate.splrep(x, y_z3, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z3, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z3, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z3, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

#plt.show()


#next the others dpsi  guys zobs = 3.5
y_z4 = np.array([0.13705, 0.05318, 0, -0.02108, -0.02326, -0.01762, -0.01140], dtype = float)

tck_z4 = interpolate.splrep(x, y_z4, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z4, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z4, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z4, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

#plt.show()


#next the others dpsi  guys zobs = 4.0
y_z5 = np.array([0.21727, 0.08893, 0, -0.04266, -0.05708, -0.05848, -0.05478], dtype = float)

tck_z5 = interpolate.splrep(x, y_z5, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z5, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z5, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z5, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

#plt.show()

# PLOTTING ENDED

#integral implementation with constant  psi derivative with  respect to z
#dpsi = 0.0715296
delta_z = 0.14
redshift_test = 5.90
k_test = 0.0977

#k_test = 0.0977/2


#37 in range for zmax = 11
#206 for zmax = 34.7
def crosspower_z2(redshift, k):
	result = 0 	
	zobs_z2 = 2.5
	for j in range(0,206):
	#for j in range(0,37):
		result += interpolate.splev(redshift, tck_z2, der = 1)*crosspower(redshift, k)*Dratio(zobs_z2, redshift)*delta_z
		#print 'did %d or %f \n' % (j,redshift) 	
		redshift = redshift + delta_z
	return result*(-1)

#print crosspower_final(redshift_test, k_test, zobs_test)
print 'Result for desired cross power spectrum at zobs = 2.5 is %f\n' % crosspower_z2(redshift_test, k_test) 


def crosspower_z1(redshift, k):
	result = 0 	
	zobs_z1 = 2.0
	for j in range(0,206):
	#for j in range(0,37):
		result += interpolate.splev(redshift, tck_z1, der = 1)*crosspower(redshift, k)*Dratio(zobs_z1, redshift)*delta_z
		#print 'did %d or %f \n' % (j,redshift) 	
		redshift = redshift + delta_z
	return result*(-1)

print 'Result for desired cross power spectrum at zobs = 2.0 is %f\n' % crosspower_z1(redshift_test, k_test)
 
def crosspower_z3(redshift, k):
	result = 0 	
	zobs_z3 = 3.0
	for j in range(0,206):
	#for j in range(0,37):
		result += interpolate.splev(redshift, tck_z3, der = 1)*crosspower(redshift, k)*Dratio(zobs_z3, redshift)*delta_z
		#print 'did %d or %f \n' % (j,redshift) 	
		redshift = redshift + delta_z
	return result*(-1)

print 'Result for desired cross power spectrum at zobs = 3.0 is %f\n' % crosspower_z3(redshift_test, k_test)

def crosspower_z4(redshift, k):
	result = 0 	
	zobs_z4 = 3.5
	for j in range(0,206):
	#for j in range(0,37):
		result += interpolate.splev(redshift, tck_z4, der = 1)*crosspower(redshift, k)*Dratio(zobs_z4, redshift)*delta_z
		#print 'did %d or %f \n' % (j,redshift) 	
		redshift = redshift + delta_z
	return result*(-1)

print 'Result for desired cross power spectrum at zobs = 3.5 is %f\n' % crosspower_z4(redshift_test, k_test)

def crosspower_z5(redshift, k):
	result = 0 	
	zobs_z5 = 4.0
	for j in range(0,206):
	#for j in range(0,37):
		result += interpolate.splev(redshift, tck_z5, der = 1)*crosspower(redshift, k)*Dratio(zobs_z5, redshift)*delta_z
		#print 'did %d or %f \n' % (j,redshift) 	
		redshift = redshift + delta_z
	return result*(-1)

print 'Result for desired cross power spectrum at zobs = 4.0 is %f\n' % crosspower_z5(redshift_test, k_test)

#here  I will interpolate the P_{m,\psi} (k, z_obs)
zobs_array = [2.0, 2.5, 3.0, 3.5, 4.0]
#crosspmxh_array =[]
crosspmxh_array = np.zeros(16*5)
a = 0
b = 0
for a in range(0,4):
	for b in range(0,16):
		if a==0:
			crosspmxh_array[b] = crosspower_z1(redshift_test,k_array[b]) 
			b = b + 1
		elif a==1:
			crosspmxh_array[b+16] = crosspower_z2(redshift_test,k_array[b]) 
			b = b + 1
		elif a==2:
			crosspmxh_array[b+16*2] = crosspower_z3(redshift_test,k_array[b]) 
			b = b + 1
		elif a==3:
			crosspmxh_array[b+16*3] = crosspower_z4(redshift_test,k_array[b]) 
			b = b + 1
		else:
			crosspmxh_array[b+16*4] = crosspower_z5(redshift_test,k_array[b]) 
			b = b + 1
	a = a + 1

print len(crosspmxh_array)
print len(k_array)
print len(zobs_array)

crosspower_final = interp2d(zobs_array, k_array, crosspmxh_array, kind='cubic')

print crosspower_final(2.2, 0.9999)
#this is where  matter matter starts
#comparison with matter power spectrum at observation redshift
table_z2 = np.loadtxt('lymanproject_z3_pk_nl.dat')

#convert from class h convertion
k_marrayz2   = table_z2[:,0]*h
p_mmz2 = table_z2[:,1]/(h**3)

#print len(k_marrayz2)
#print len(p_mmz2)
#print p_mmz2
#print k_marrayz2

autopower_z2 = interp1d(k_marrayz2, p_mmz2, kind='quadratic')

#get the delta instead of the power
def matterpower_z2(k):
	return ((k**3)*(autopower_z2(k)))/(2.0*(np.pi)**2)

print matterpower_z2(k_test)
print matterpower_z2(0.0970)

table_z1 = np.loadtxt('lymanproject_z2_pk_nl.dat')

#convert from class h convertion
k_marrayz1   = table_z1[:,0]*h
p_mmz1 = table_z1[:,1]/(h**3)

autopower_z1 = interp1d(k_marrayz1, p_mmz1, kind='quadratic')

#get the delta instead of the power
def matterpower_z1(k):
	return ((k**3)*(autopower_z1(k)))/(2.0*(np.pi)**2)


table_z3 = np.loadtxt('lymanproject_z4_pk_nl.dat')

#convert from class h convertion
k_marrayz3   = table_z3[:,0]*h
p_mmz3 = table_z3[:,1]/(h**3)

autopower_z3 = interp1d(k_marrayz3, p_mmz3, kind='quadratic')

#get the delta instead of the power
def matterpower_z3(k):
	return ((k**3)*(autopower_z3(k)))/(2.0*(np.pi)**2)


table_z4 = np.loadtxt('lymanproject_z5_pk_nl.dat')

#convert from class h convertion
k_marrayz4   = table_z4[:,0]*h
p_mmz4 = table_z4[:,1]/(h**3)

autopower_z4 = interp1d(k_marrayz4, p_mmz4, kind='quadratic')

#get the delta instead of the power
def matterpower_z4(k):
	return ((k**3)*(autopower_z4(k)))/(2.0*(np.pi)**2)


table_z5 = np.loadtxt('lymanproject_z6_pk_nl.dat')

#convert from class h convertion
k_marrayz5   = table_z5[:,0]*h
p_mmz5 = table_z5[:,1]/(h**3)

autopower_z5 = interp1d(k_marrayz5, p_mmz5, kind='quadratic')

#get the delta instead of the power
def matterpower_z5(k):
	return ((k**3)*(autopower_z5(k)))/(2.0*(np.pi)**2)


#comparison @ mu = 0
def comparison(cross,auto):
	return (-1*cross/auto)*200

final_z2 = comparison(crosspower_z2(redshift_test, k_test),matterpower_z2(k_test))

final_z1 = comparison(crosspower_z1(redshift_test, k_test),matterpower_z1(k_test))

final_z3 = comparison(crosspower_z3(redshift_test, k_test),matterpower_z3(k_test))

final_z4 = comparison(crosspower_z4(redshift_test, k_test),matterpower_z4(k_test))

final_z5 = comparison(crosspower_z5(redshift_test, k_test),matterpower_z5(k_test))
print 'Comparison of effect at zobs = 2.0 and k = %f is %e \n' % (k_test, final_z1)
print 'Comparison of effect at zobs = 2.5 and k = %f is %e \n' % (k_test, final_z2)
print 'Comparison of effect at zobs = 3.0 and k = %f is %e \n' % (k_test, final_z3)
print 'Comparison of effect at zobs = 3.5 and k = %f is %e \n' % (k_test, final_z4)
print 'Comparison of effect at zobs = 4.0 and k = %f is %e \n' % (k_test, final_z5)

#comparison @mu=0 but with the bias for redshifts that we have them. This is from Chris' paper
r2 = final_z2*(0.122/0.125)
r3 = final_z3*(0.156/0.201)

print 'Comparison of effect with the bias coefficients at zobs = 2.5 and k = %f is %e \n' % (k_test, r2)
print 'Comparison of effect with the bias coefficients at zobs = 3.0 and k = %f is %e \n' % (k_test, r3)

#here we put the bias coefficients back to the comparison
#first the flux bias I ignored the negative sign
bias_F_z1 = 0.12
bias_F_z2 = 0.18
bias_F_z3 = 0.27
bias_F_z4 = 0.37
bias_F_z5 = 0.55

#now the radiation bias
bias_G_z1 = (0.087 + 0.087 + 0.085 + 0.084 + 0.084 + 0.083)/6.0
bias_G_z2 = (0.148 + 0.148 + 0.146 + 0.145 + 0.145 + 0.144)/6.0
bias_G_z3 = (0.231 + 0.236 + 0.237 + 0.236 + 0.236 + 0.236)/6.0
bias_G_z4 = (0.332 + 0.346 + 0.354 + 0.356 + 0.356 + 0.357)/6.0
bias_G_z5 = (0.448 + 0.476 + 0.499 + 0.503 + 0.508 + 0.509)/6.0

print bias_G_z1
print bias_G_z2
print bias_G_z3
print bias_G_z4
print bias_G_z5

#ratios of the bias
bias_ratio_z1 = bias_G_z1/bias_F_z1
bias_ratio_z2 = bias_G_z2/bias_F_z2
bias_ratio_z3 = bias_G_z3/bias_F_z3
bias_ratio_z4 = bias_G_z4/bias_F_z4
bias_ratio_z5 = bias_G_z5/bias_F_z5

print bias_ratio_z1
print bias_ratio_z2
print bias_ratio_z3
print bias_ratio_z4
print bias_ratio_z5

#final answers with bias now
comp_z1 = bias_ratio_z1*final_z1
comp_z2 = bias_ratio_z2*final_z2
comp_z3 = bias_ratio_z3*final_z3
comp_z4 = bias_ratio_z4*final_z4
comp_z5 = bias_ratio_z5*final_z5

print 'Comparison effect with bias ratio at zobs = 2.0 and k = %f is %e \n' % (k_test, comp_z1)
print 'Comparison effect with bias ratio at zobs = 2.5 and k = %f is %e \n' % (k_test, comp_z2)
print 'Comparison effect with bias ratio at zobs = 3.0 and k = %f is %e \n' % (k_test, comp_z3)
print 'Comparison effect with bias ratio at zobs = 3.5 and k = %f is %e \n' % (k_test, comp_z4)
print 'Comparison effect with bias ratio at zobs = 4.0 and k = %f is %e \n' % (k_test, comp_z5)
