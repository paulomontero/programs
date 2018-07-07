import numpy as np
import scipy as sp

#this script deals with the different time scales involved with the STAGE 3 part of our calculation


#****************************************************************************************************
#run dependant stuff
#temp = 2.67e+9 #temperature of the relevant patch of gas in K
total_rho = 2.8e+6 #total density in g/cm^3 (decided for central density value for now)
mach_number = 3.0 #dimensionless
mass_pbh = 2.0e+19 #in g 
c_sound = 1.2e+8 # in cm/s
x_relativity = 0.6 #dimensionless


#********************************************************************************************************
#useful definitions
R = 8.314e+7 #in erg/K/g
Kb = 1.381e-16 # in erg/K
#Ef = 4.8065e-7 # in erg
a_rad = 7.5657e-15 # in erg/cm^3/K^4
hbar = 1.054e-27 # erg*s
zbar = 7.0
abar = 14.0
e_charge = 1.602e-19 # in C
mass_e = 9.11e-28 # in g 
mass_nuc = 1.67e-24 # in g
alpha = 1.0/137.0 
c_light = 2.998e+10 # in cm/s
avogadro_number = 6.022e+23 # in 1/mol
epsilon = 8.854e-21 # in SI e-12 this has c^2*s^2/(kg*m^3)*(1 kg/ 1000 g)*(1 m/100 cm)^3 -> e-21


#****************************************************************************
#useful parameters
n_e = total_rho*zbar/(1.0*mass_nuc*abar)
R_c = (2.0*6.67e-8*mass_pbh)/(3.0*c_sound**2) # is for getting it in cm
t_code = R_c/(1.0*c_sound) # is for getting cm/s the speed, t code is in seconds 

print 'The distance unit in the code R_c is %e \n' % R_c 
print 'The time unit in the code is t_code % e \n' % t_code

radius_shell_133 = 0.215031*R_c

print 'The radius of shell 133 is %e \n' % radius_shell_133

Ef = (total_rho*3.0*zbar*np.pi**2*(hbar*c_light)**3/(1.0*mass_nuc*abar))**(1.0/3.0)
Ef_mev = Ef*624151
print 'Fermi energy is %e in erg or %e in MeV \n' % (Ef, Ef_mev) 

temp = 0.4*Ef/(1.0*Kb)
print 'Temperature is %e in K \n' % temp
#*************************************************************************************************************
#epsilon nuc part
#right now only for 1.75<T9<3.3 : i=1, T9<1.75 : i = 0, 3.3<T9<6.0 : i =2
qbar_array = [5.723e-6,5.11e-6, 5.069e-6]

x_array = np.linspace(1e6, 1e10, 1000)
i = 0
def nuc_energy(rho, T):
	T9 = T/1.0e+9	
	T9a = T9/(1.0 + 0.0396*T9)
	global i
	if (T9 < 1.75):
		i = 0
	elif (1.75 < T9 < 3.3):
		i = 1
	elif (3.3 < T9 < 6.0):
		i = 2
	lamb_nuc = 4.27e+26*(T9a**(5./6.)/T9**(3./2.))*np.exp(-84.165/(1.0*T9a**(1./3.)) - 2.12e-3*T9**3)
	#print 'The lambda nuclear is %e \n' % lamb_nuc
	return avogadro_number*rho*qbar_array[i]*(1.0/(14.0*2))**2*lamb_nuc/2.0

print 'The specific nuclear energy generation rate is %e in erg/g/s \n' % nuc_energy(total_rho, temp)



#*******************************************************************************************************************
#Screening stuff from the code T suggested, only considering weak screening
def screening_enhancement(rho, T):
	z2bar = abar*(6.0**2/28.0 + 8.0**2/28.0)
	print 'z2bar check %d \n' % z2bar
	pp = np.sqrt((rho/28.0)*(z2bar + zbar)/(1.0*T))
	print 'pp factor is %e \n' % pp
	qlamboz = 1.88e+8*pp/(1.0*T)
	print 'qlamboz factor %e \n' % qlamboz
	h12 = 6.0**2*qlamboz
	print 'h12 result is %e \n' % h12
	return np.exp(h12)

print 'The screening enhancement factor is %e \n' % screening_enhancement(total_rho, temp)

maho = nuc_energy(total_rho, temp)*screening_enhancement(total_rho, temp)
print 'The specific nuclear energy generation rate with screening enhancement is %e in erg/g/s \n' % maho

#Well leaving this for comparison, but it seems strong screening should be ready
def strong_screening_enhancement(rho, T):
	zhat = 12**(5.0/3.0) - 2*6**(5.0/3.0)
	zhat2 = 12**(5.0/12.0) - 2*6**(5.0/12.0)
	gamp = (2.27495e5/(1.0*T))*(rho*zbar/(1.0*abar))**(1.0/3.0)
#	print 'gamp is %e \n' % gamp
	tau12 = (4.248719e3/3.0)*(36**2*144/(24*T))**(1.0/3.0)
#	print 'Tau12 is %e \n' % tau12
	alpha12 = (36.0*gamp*(2.0/12.0)**(1.0/3.0))/(1.0*tau12)
#	print 'alpha12 is %e \n' % alpha12
	ss = tau12*(5.0/32.0 - alpha12*(0.014 + 0.0128*alpha12))
	vv = 36*gamp*(2/12)**(1.0/3.0)*alpha12*(0.0055 + alpha12*(-0.0098 + 0.0048*alpha12))
	cc = 0.896434*gamp*zhat - 3.4474*gamp**(1.0/4.0)*zhat2 - 0.5551*(np.log(gamp) + (5.0/3.0)*np.log(3.0))
	h12 = cc -alpha12**3*(ss + vv)
	rr = 1 - 0.0562*alpha12**3
#	print 'rr for strong screening is %e \n' % rr 
	if (rr > 0.77):
		xlgfac = rr
	elif (rr < 0.77):
		xlgfac = 0.77
	h12 = np.log(xlgfac) + h12
	return h12


print 'Degenerate screening factor is %e \n' % strong_screening_enhancement(total_rho, temp)


opa =  nuc_energy(total_rho, temp)*np.exp(strong_screening_enhancement(total_rho, temp))
print 'DEGENERATE The specific nuclear energy generation rate with DEGENERATE screening enhancement is %e in erg/g/s \n' % opa

#*******************************************************************************************************************
#conduction losses part
def cond_energy_loss(T, rho, x, r):
	gamma_cou = alpha*hbar*c_light*(4*np.pi*n_e/3.0)**(1./3.)/(Kb*T)
	#print 'Coupling parameter %e \n' % gamma_cou
	lamb_cou = np.log((2.0*np.pi*zbar/3.0)**(1./3.)*(3.0/2.0 + 3.0/1.0*gamma_cou)**(1./2.) - x**2/(2*(1 + x**2)))
	#print 'Lambda term %e \n' % lamb_cou	
	freq_ei = 4*e_charge**4*mass_e*zbar*lamb_cou*np.sqrt(1.0 + x**2)/(3.0*np.pi*hbar**3*16.0*np.pi**2*epsilon**2)
	#print 'The frequency of ei collisions %e in Hz \n' % freq_ei
	sigma_cond =np.pi**2*Kb**2*T*n_e/(3.0*mass_e*np.sqrt(1 + x**2)*freq_ei)
	#print 'The sigma_cond is %e in erg/T/t/L \n' % sigma_cond
	return 2*sigma_cond*T/(rho*r**2)
	#D_cond = (np.pi**3*hbar**3*n_e*Kb**2*T)/(4*e_charge**4*mass_e**2*zbar*(1 + x**2)*lamb_cou)
	#print 'The D_cond %e \n' % D_cond
	#kappa_cond = (4.0*a_rad*T**3)/(3.0*rho*D_cond)
	#print 'The kappa cond is %e \n' % kappa_cond
	#return (2.0*kappa_cond*T)/(1.0*rho*r**2)
#this is not working I suspect a big problem with units maybe try again with SI units
#possible solutions D_cond has the right units! Kappa cond is not actually kappa cond because units do not match. Also, use the correction for the frequency

print 'The energy loss by conduction in units of erg/g/s is %e \n' % cond_energy_loss(temp, total_rho, x_relativity, radius_shell_133)


#for plotting in the future
h12_T7_array = np.zeros(1000)
h12_T5_array = np.zeros(1000)
h12_T3_array = np.zeros(1000)
h12_T1_array = np.zeros(1000)
nuc_array = np.zeros(1000)
nuc_screen_array = np.zeros(1000)
cond_array = np.zeros(1000)
ii = 0
for ii in range (0,1000):
	h12_T7_array[ii] = strong_screening_enhancement(x_array[ii], 7e9)
	h12_T5_array[ii] = strong_screening_enhancement(x_array[ii], 5e9)
	h12_T3_array[ii] = strong_screening_enhancement(x_array[ii], 3e9)
	h12_T1_array[ii] = strong_screening_enhancement(x_array[ii], 1e9)	
	nuc_array[ii] = nuc_energy(x_array[ii], 7e9)
	nuc_screen_array[ii] = nuc_array[ii]*np.exp(h12_T7_array[ii])
	cond_array[ii] = cond_energy_loss(7e9, x_array[ii], x_relativity, radius_shell_133)
	ii = ii + 1



np.savetxt('plot_h12_rho.dat', np.transpose([x_array, h12_T7_array, h12_T5_array, h12_T3_array, h12_T1_array]), fmt='%e', delimiter=' ')
np.savetxt('plot_epsilon_rho.dat', np.transpose([x_array, cond_array, nuc_array, nuc_screen_array]), fmt='%e', delimiter=' ')

#**************************************************************************************************************************
#time burn part
#specific heat of ions in units of R
cp_ion = (5.0/28.0)*(R/(1.0*abar))

#specific heat of electrons in units of R
def cp_ele(T):
	return ((np.pi**2)/2.0)*(R/(1.0*abar))*(Kb*T/Ef)*(1 + ((np.pi**2)/3.0)*(Kb*T/Ef)**2)

dif_cp = np.abs(cp_ion - cp_ele(temp))
print 'The cp of ions is %e, the cp of electrons is %e, and the difference between both is %e \n' % (cp_ion, cp_ele(temp), dif_cp)


#the burning time

hola = (cp_ele(temp) + cp_ion)*temp
print 'CpT is %e \n' % hola

t_burn = (cp_ele(temp) + cp_ion)*temp/(1.0*nuc_energy(total_rho, temp))
print 'The burning time is equal to %e s \n' % t_burn

#************************************************************************************************************************
#kevin helmholtz instability time scale
#for this just take a look at the grad_vz.dat file, which is generated with grad_vz_withz.py and look at the highest grad vz, highest grad vz implies smaller time scale for the instability

#for shell 133 we have -9.841690
grad_vz = np.abs(-9.841690)
t_insta = 1.0/(1.0*grad_vz)*t_code

print 'The instability time scale is %e s \n' % t_insta



