import numpy as np
import scipy as sp
import sys
from scipy.integrate import quad
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy import optimize
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
#table = np.loadtxt('power_default1')
#table = np.loadtxt('power_default2')
#table = np.loadtxt('power_default3')
#table = np.loadtxt('power_default4')
#table = np.loadtxt('power_default5')
#table = np.loadtxt('power_default6')
table = np.loadtxt('power_default7')
#table = np.loadtxt('power_default8')

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

#I want a table of crosspower_mxh as a function of redshift , for given k. Let's say k = 0.1
o=0
deltaplot_array = np.zeros(84)
for o in range(0,84):
	deltaplot_array[o] = crosspower(z_array[o], 0.1)
	o = o + 1	

np.savetxt('plotdeltamxh_default7.dat', np.transpose([z_array, deltaplot_array]), fmt='%e', delimiter=' ')









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

plt.show()

#next the others dpsi  guys, zobs = 2.0
y_z1 = np.array([0.05621, 0.01945, 0, -0.008235, -0.01222, -0.01311, -0.01468], dtype = float)

tck_z1 = interpolate.splrep(x, y_z1, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z1, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z1, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z1, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

plt.show()

#next the others dpsi  guys zobs = 3.0
y_z3 = np.array([0.08588, 0.03052, 0, -0.006657, -0.001757, 0.005913, 0.01350], dtype = float)

tck_z3 = interpolate.splrep(x, y_z3, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z3, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z3, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z3, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

plt.show()


#next the others dpsi  guys zobs = 3.5
y_z4 = np.array([0.13705, 0.05318, 0, -0.02108, -0.02326, -0.01762, -0.01140], dtype = float)

tck_z4 = interpolate.splrep(x, y_z4, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z4, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z4, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z4, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

plt.show()


#next the others dpsi  guys zobs = 4.0
y_z5 = np.array([0.21727, 0.08893, 0, -0.04266, -0.05708, -0.05848, -0.05478], dtype = float)

tck_z5 = interpolate.splrep(x, y_z5, k=1, s=0)

fig, axes = plt.subplots(2)
axes[0].plot(x, y_z5, 'x', label = 'data')
axes[0].plot(xnew, interpolate.splev(xnew, tck_z5, der = 0), label = 'Fit')
axes[1].plot(xnew, interpolate.splev(xnew, tck_z5, der = 1), label = '1st dev')

for ax in axes:
	ax.legend(loc = 'best')

plt.show()

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
#there was a typo here, need to fix this -> re run some realizations	
	zobs_z5 = 4.0
	for j in range(0,206):
	#for j in range(0,37):
		result += interpolate.splev(redshift, tck_z5, der = 1)*crosspower(redshift, k)*Dratio(zobs_z5, redshift)*delta_z
		#print 'did %d or %f \n' % (j,redshift) 	
		redshift = redshift + delta_z
	return result*(-1)

print 'Result for desired cross power spectrum at zobs = 4.0 is %f\n' % crosspower_z5(redshift_test, k_test)

#*************************************************************************************************************************************
#here  I will interpolate the P_{m,\psi} (k, z_obs)
zobs_array = [2.0, 2.5, 3.0, 3.5, 4.0]
#crosspmxh_array =[]
crosspmxh_array = np.zeros(16*5)
a = 0
b = 0
for a in range(0,5):
	b = 0
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
		elif a==4:
			crosspmxh_array[b+16*4] = crosspower_z5(redshift_test,k_array[b]) 
			b = b + 1
	a = a + 1

print len(crosspmxh_array)
print len(k_array)
print len(zobs_array)

crosspower_final = interp2d(zobs_array, k_array, crosspmxh_array, kind='cubic')

print 'The crosspower spectrum at observer redshift 2.2 and wavenumber 0.9999 is %f\n' % crosspower_final(2.2, 0.9999)


#************************************************************************************************************************
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


#old interpolator with cubic problems
#autopower_z2 = interp1d(k_marrayz2, p_mmz2, kind='cubic')
autopower_z2 = CubicSpline(k_marrayz2, p_mmz2)

#get the delta instead of the power
def matterpower_z2(k):
	return ((k**3)*(autopower_z2(k)))/(2.0*(np.pi)**2)

print matterpower_z2(k_test)
print 'Matter power spectrum P (not delta) for zobs =2.5 k= 0.09997 is %f' % autopower_z2(k_test) 

table_z1 = np.loadtxt('lymanproject_z2_pk_nl.dat')

#convert from class h convertion
k_marrayz1   = table_z1[:,0]*h
p_mmz1 = table_z1[:,1]/(h**3)


#autopower_z1 = interp1d(k_marrayz1, p_mmz1, kind='quadratic')
autopower_z1 = CubicSpline(k_marrayz1, p_mmz1)
#get the delta instead of the power
def matterpower_z1(k):
	return ((k**3)*(autopower_z1(k)))/(2.0*(np.pi)**2)


table_z3 = np.loadtxt('lymanproject_z4_pk_nl.dat')

#convert from class h convertion
k_marrayz3   = table_z3[:,0]*h
p_mmz3 = table_z3[:,1]/(h**3)

#autopower_z3 = interp1d(k_marrayz3, p_mmz3, kind='quadratic')
autopower_z3 = CubicSpline(k_marrayz3, p_mmz3)

#get the delta instead of the power
def matterpower_z3(k):
	return ((k**3)*(autopower_z3(k)))/(2.0*(np.pi)**2)


table_z4 = np.loadtxt('lymanproject_z5_pk_nl.dat')

#convert from class h convertion
k_marrayz4   = table_z4[:,0]*h
p_mmz4 = table_z4[:,1]/(h**3)

#autopower_z4 = interp1d(k_marrayz4, p_mmz4, kind='quadratic')
autopower_z4 = CubicSpline(k_marrayz4, p_mmz4)

#get the delta instead of the power
def matterpower_z4(k):
	return ((k**3)*(autopower_z4(k)))/(2.0*(np.pi)**2)


table_z5 = np.loadtxt('lymanproject_z6_pk_nl.dat')

#convert from class h convertion
k_marrayz5   = table_z5[:,0]*h
p_mmz5 = table_z5[:,1]/(h**3)

#autopower_z5 = interp1d(k_marrayz5, p_mmz5, kind='quadratic')
autopower_z5 = CubicSpline(k_marrayz5, p_mmz5)

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

#comparison @mu=0 but with the bias for redshifts that we have them
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

#here I want to plot the final cross-power as a function of k, for a given zobs, in this case zobs = 2.5, then z1 = zobs= 2.0
m=0
deltamm_array = np.zeros(16)
crossplot_array = np.zeros(16)
for m in range(0,16):
	#crossplot_array[m] = crosspower_z2(redshift_test, k_array[m])
	#crossplot_array[m] = crosspower_z1(redshift_test, k_array[m])
	#crossplot_array[m] = crosspower_z3(redshift_test, k_array[m])
	#crossplot_array[m] = crosspower_z4(redshift_test, k_array[m])	
	crossplot_array[m] = crosspower_z5(redshift_test, k_array[m])
	#deltamm_array[m] = matterpower_z2(k_array[m]) 
	#deltamm_array[m] = matterpower_z1(k_array[m])
	#deltamm_array[m] = matterpower_z3(k_array[m])
	#deltamm_array[m] = matterpower_z4(k_array[m])
	deltamm_array[m] = matterpower_z5(k_array[m])
	m = m + 1	

np.savetxt('plotprueba_default7z5.dat', np.transpose([k_array, crossplot_array, deltamm_array]), fmt='%e', delimiter=' ')

#I got what is needed to plot the ration Delta_{m,psi}/Delta_{m,m}, by modifying previous line

#****************************************************************************************************************************************
#Need to interpolate P_m(k,zobs)
automm_array = np.zeros(16*5)

aa = 0
bb = 0
for aa in range(0,5):
	bb = 0
	for bb in range(0,16): 
		if aa==0:
			automm_array[bb] = autopower_z1(k_array[bb]) 
			bb = bb + 1
		elif aa==1:
			automm_array[bb + 16] = autopower_z2(k_array[bb]) 
			bb = bb + 1
		elif aa==2:
			automm_array[bb + 16*2] = autopower_z3(k_array[bb])
			bb = bb + 1
		elif aa==3:
			automm_array[bb + 16*3] = autopower_z4(k_array[bb])
			bb = bb + 1
		elif aa==4:
			automm_array[bb + 16*4] = autopower_z5(k_array[bb])
			bb = bb + 1
	aa = aa + 1

autopower_final = interp2d(zobs_array, k_array, automm_array, kind='cubic')
#****************************************************************************************************************************************
#interpolation tests
print autopower_final(2.0,k_test), autopower_z1(k_test)
print autopower_final(2.5,k_test), autopower_z2(k_test)
print autopower_final(3.0,k_test), autopower_z3(k_test)
print autopower_final(3.5,k_test), autopower_z4(k_test)
print autopower_final(4.0,k_test), autopower_z5(k_test)

print crosspower_final(2.0,k_test), crosspower_z1(redshift_test,k_test)
print crosspower_final(2.5,k_test), crosspower_z2(redshift_test,k_test)
print crosspower_final(3.0,k_test), crosspower_z3(redshift_test,k_test)
print crosspower_final(3.5,k_test), crosspower_z4(redshift_test,k_test)
print crosspower_final(4.0,k_test), crosspower_z5(redshift_test,k_test)
#****************************************************************************************************************************************
#Ok let's try doing the 1D with the interpolations
def linear_integrand(kpar, kper, kmag, zobs):
	factor = kper*((1 + (kpar**2)/(kmag**2))**2)
	return factor*autopower_final(zobs, kmag)/(2.0*np.pi)

def second_integrand(kpar, kper, kmag, zobs):
	factor2 = kper*(1 + (kpar/kmag)**2)
	return factor2*crosspower_final(zobs, kmag)*(np.pi/kmag**3)

#********************************************************************************************************************************
#it worked, now make it a function
def integralisimo(k_par):
	deltak = 0.05
	ii = 0
	k_per = 0
	lzobs1 = 0
	szobs1 = 0
	lzobs2 = 0
	szobs2 = 0
	lzobs3 = 0
	szobs3 = 0
	lzobs4 = 0
	szobs4 = 0
	lzobs5 = 0
	szobs5 = 0
	k_max = 3.59
	for ii in range(0,71):
		k_per += 0.05 + deltak*ii
		k_mag = np.sqrt(k_par**2 + k_per**2)
		if (k_mag <= k_max): 
			lzobs1 += linear_integrand(k_par, k_per, k_mag, 2.0)*deltak
			szobs1 += second_integrand(k_par, k_per, k_mag, 2.0)*deltak
			lzobs2 += linear_integrand(k_par, k_per, k_mag, 2.5)*deltak
			szobs2 += second_integrand(k_par, k_per, k_mag, 2.5)*deltak
			lzobs3 += linear_integrand(k_par, k_per, k_mag, 3.0)*deltak
			szobs3 += second_integrand(k_par, k_per, k_mag, 3.0)*deltak
			lzobs4 += linear_integrand(k_par, k_per, k_mag, 3.5)*deltak
			szobs4 += second_integrand(k_par, k_per, k_mag, 3.5)*deltak
			lzobs5 += linear_integrand(k_par, k_per, k_mag, 4.0)*deltak
			szobs5 += second_integrand(k_par, k_per, k_mag, 4.0)*deltak
		ii = ii + 1

	return lzobs1, lzobs2, lzobs3, lzobs4, lzobs5, szobs1, szobs2, szobs3, szobs4, szobs5

#***************************************************************************************************************************************
onedcomp_array = np.zeros(10)
onedcomp_array = integralisimo(k_test)


onedcomp_z1 = bias_ratio_z1*onedcomp_array[5]/(onedcomp_array[0]*1.0)*(-200)
onedcomp_z2 = bias_ratio_z2*onedcomp_array[6]/(onedcomp_array[1]*1.0)*(-200)
onedcomp_z3 = bias_ratio_z3*onedcomp_array[7]/(onedcomp_array[2]*1.0)*(-200)
onedcomp_z4 = bias_ratio_z4*onedcomp_array[8]/(onedcomp_array[3]*1.0)*(-200)
onedcomp_z5 = bias_ratio_z5*onedcomp_array[9]/(onedcomp_array[4]*1.0)*(-200)

print 'Comparison effect for 1D with bias ratio at zobs = 2.0 and k = %f is %e \n' % (k_test, onedcomp_z1)
print 'Comparison effect for 1D with bias ratio at zobs = 2.5 and k = %f is %e \n' % (k_test, onedcomp_z2)
print 'Comparison effect for 1D with bias ratio at zobs = 3.0 and k = %f is %e \n' % (k_test, onedcomp_z3)
print 'Comparison effect for 1D with bias ratio at zobs = 3.5 and k = %f is %e \n' % (k_test, onedcomp_z4)
print 'Comparison effect for 1D with bias ratio at zobs = 4.0 and k = %f is %e \n' % (k_test, onedcomp_z5)
#****************************************************************************************************************************************
#time to get the data for plots, here we try to get the P1D_m(k,zobs)/P1D_{m,psi}(k,zobs)
mm = 0
pmm_arrayz1 = np.zeros(16)
pmm_arrayz2 = np.zeros(16)
pmm_arrayz3 = np.zeros(16)
pmm_arrayz4 = np.zeros(16)
pmm_arrayz5 = np.zeros(16)
mid_array = np.zeros(10)
pcrossplotz1_array = np.zeros(16)
pcrossplotz2_array = np.zeros(16)
pcrossplotz3_array = np.zeros(16)
pcrossplotz4_array = np.zeros(16)
pcrossplotz5_array = np.zeros(16)
for mm in range(0,16):
	mid_array = integralisimo(k_array[mm])
	pmm_arrayz1[mm] = mid_array[0]
	pmm_arrayz2[mm] = mid_array[1]
	pmm_arrayz3[mm] = mid_array[2]
	pmm_arrayz4[mm] = mid_array[3]
	pmm_arrayz5[mm] = mid_array[4]
	pcrossplotz1_array[mm] = mid_array[5]
	pcrossplotz2_array[mm] = mid_array[6]
	pcrossplotz3_array[mm] = mid_array[7]
	pcrossplotz4_array[mm] = mid_array[8]
	pcrossplotz5_array[mm] = mid_array[9]
	mm = mm + 1

np.savetxt('plot1Dprueba_default7z1.dat', np.transpose([k_array, pcrossplotz1_array, pmm_arrayz1]), fmt='%e', delimiter=' ')
np.savetxt('plot1Dprueba_default7z2.dat', np.transpose([k_array, pcrossplotz2_array, pmm_arrayz2]), fmt='%e', delimiter=' ')
np.savetxt('plot1Dprueba_default7z3.dat', np.transpose([k_array, pcrossplotz3_array, pmm_arrayz3]), fmt='%e', delimiter=' ')
np.savetxt('plot1Dprueba_default7z4.dat', np.transpose([k_array, pcrossplotz4_array, pmm_arrayz4]), fmt='%e', delimiter=' ')
np.savetxt('plot1Dprueba_default7z5.dat', np.transpose([k_array, pcrossplotz5_array, pmm_arrayz5]), fmt='%e', delimiter=' ')



#****************************************************************************************************************************************
#Here we start the process for finding out what happened at the 1D power spectrum, first the linear part. Let's ignore the bias ratio till the comparison is due
#def linear_integrand_z1(kpar, kper, kmag):
#	factor = kper*((1 + (kpar**2)/(kmag**2))**2)  
#	return factor*autopower_z1(kmag)/(2*np.pi)

#def linear_integrand_z2(kpar, kper, kmag):
#	factor = kper*((1 + (kpar**2)/(kmag**2))**2)  
#	return factor*autopower_z2(kmag)/(2*np.pi)

#def linear_integrand_z3(kpar, kper, kmag):
#	factor = kper*((1 + (kpar**2)/(kmag**2))**2)  
#	return factor*autopower_z3(kmag)/(2*np.pi)

#def linear_integrand_z4(kpar, kper, kmag):
#	factor = kper*((1 + (kpar**2)/(kmag**2))**2)  
#	return factor*autopower_z4(kmag)/(2*np.pi)

#def linear_integrand_z5(kpar, kper, kmag):
#	factor = kper*((1 + (kpar**2)/(kmag**2))**2)  
#	return factor*autopower_z5(kmag)/(2*np.pi)

#now the new term, again ignoring the bias and factor of 2
#def new_integrand_z1(kpar, kper, kmag):
#	factor =  factor = kper*(1 + kpar/kmag)
#	return factor*crosspower_z1(redshift_test,kmag)*(np.pi/kmag**3)

#def new_integrand_z2(kpar, kper, kmag):
#	factor =  factor = kper*(1 + kpar/kmag)
#	return factor*crosspower_z2(redshift_test,kmag)*(np.pi/kmag**3)

#def new_integrand_z3(kpar, kper, kmag):
#	factor =  factor = kper*(1 + kpar/kmag)
#	return factor*crosspower_z3(redshift_test,kmag)*(np.pi/kmag**3)

#def new_integrand_z4(kpar, kper, kmag):
#	factor =  factor = kper*(1 + kpar/kmag)
#	return factor*crosspower_z4(redshift_test,kmag)*(np.pi/kmag**3)

#def new_integrand_z5(kpar, kper, kmag):
#	factor =  factor = kper*(1 + kpar/kmag)
#	return factor*crosspower_z5(redshift_test,kmag)*(np.pi/kmag**3)

#delta_k = 0.01
#ii = 0
#lz1 = 0 
#lz2 = 0
#lz3 = 0
#lz4 = 0
#lz5 = 0
#sz1 = 0
#sz2 = 0 
#sz3 = 0
#sz4 = 0
#sz5 = 0
#k_per = 0
#let's do this for an specific k parallel
#k_par = k_test

#for ii in range(0,999):
#	k_per += 0.01 + delta_k*ii
#	k_mag = np.sqrt(k_par**2 + k_per**2)
#	lz1 += linear_integrand_z1(k_par, k_per, k_mag)
#	lz2 += linear_integrand_z2(k_par, k_per, k_mag)
#	lz3 += linear_integrand_z3(k_par, k_per, k_mag)
#	lz4 += linear_integrand_z4(k_par, k_per, k_mag)
#	lz5 += linear_integrand_z5(k_par, k_per, k_mag)
#	sz1 += new_integrand_z1(k_par, k_per, k_mag)
#	sz2 += new_integrand_z2(k_par, k_per, k_mag)
#	sz3 += new_integrand_z3(k_par, k_per, k_mag)
#	sz4 += new_integrand_z4(k_par, k_per, k_mag)
#	sz5 += new_integrand_z5(k_par, k_per, k_mag)
#	ii = ii + 1

#*******************************************************************************************************************************************

#final answers with bias now
#onedcomp_z1 = bias_ratio_z1*sz1/(lz1*1.0)*(-200)
#onedcomp_z2 = 2.0*bias_ratio_z2*sz2/(lz2*1.0)*(-1.0*100)
#onedcomp_z3 = 2.0*bias_ratio_z3*sz3/(lz3*1.0)*(-1.0*100)
#onedcomp_z4 = 2.0*bias_ratio_z4*sz4/(lz4*1.0)*(-1.0*100)
#onedcomp_z5 = 2.0*bias_ratio_z5*sz5/(lz5*1.0)*(-1.0*100)

#print 'Comparison effect for 1D with bias ratio at zobs = 2.0 and k = %f is %e \n' % (k_test, onedcomp_z1)
#print 'Comparison effect for 1D with bias ratio at zobs = 2.5 and k = %f is %e \n' % (k_test, onedcomp_z2)
#print 'Comparison effect for 1D with bias ratio at zobs = 3.0 and k = %f is %e \n' % (k_test, onedcomp_z3)
#print 'Comparison effect for 1D with bias ratio at zobs = 3.5 and k = %f is %e \n' % (k_test, onedcomp_z4)
#print 'Comparison effect for 1D with bias ratio at zobs = 4.0 and k = %f is %e \n' % (k_test, onedcomp_z5)
