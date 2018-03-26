import numpy as np
import sys
import pandas as pd
from classy import Class


#This is the script for the transfer function of linear pol, running class several times#
#dr = pd.DataFrame(columns=['wavenumber','redshift','E_2(z)'])
g = open('newtabletransfer.dat', 'w')


data = []   
i = 0
h = 0.6726
Omega_m = 0.308
# Got  the radiation density from equality of matter and radiation.
Omega_r = 9.134e-5
for i in range(0,2000):
	params = {
        	'output': 'vTk, pCl, tCl',
        	'A_s': 2.206e-9,
        	'n_s': 0.9652, 
        	'h': 0.6726,
        	'omega_b': 0.02222,
        	'omega_cdm': 0.1199,
        	'Omega_k' : 0,
        	'T_cmb' : 2.725,
        	'format' : 'class',
        	'background_verbose' : 1,
        	'thermodynamics_verbose' : 1,
        	'perturbations_verbose' : 1,
        	'primordial_verbose' : 1,
        	'spectra_verbose': 1,
        	'nonlinear_verbose' : 1,
        	'output_verbose' : 1,
		'transfer_verbose' : 1,
        	'tau_reio' : 0.079,
		'headers' : 'yes',
		'modes' : 's',
		'k_output_values' : 0.0001 + i*0.0005, 
        	'lensing_verbose' : 1,
		'radiation_streaming_trigger_tau_over_tau_k' : 10000
        	}
	cosmo = Class()
	cosmo.set(params)
	cosmo.compute()
#	kval = 0.0001 + i/(10000*1.0)
	kval = 0.0001 + i*0.0005
	data = cosmo.get_perturbations()['scalar']
	
	print kval

	print data[0]['a']	
	print len(data[0]['a'])

	out = np.zeros(shape=(len(data[0]['a']),4))
		
	out[:,0] = kval
	out[:,1]= 1.0/data[0]['a'] - 1
	# Here comes conformal time
	#out[:,2] = 3.086e19*np.sqrt(Omega_r + Omega_m*data[0]['a'])/(50*h*Omega_m)
	out[:,2] = 2*3.086e+19*np.sqrt(Omega_r + Omega_m/(1+out[:,1]))/(100*h*Omega_m)	
	out[:,3] = -5.0*(data[0]['pol0_g']+data[0]['pol2_g'])/(4.0*np.sqrt(6.0))
	
	print out[2,1]
	# reducing the redshift range. range from redshift 100 to 1000
	j = 0
	for j in range(0,len(data[0]['a'])):
		if (out[j,1] >= 100) and (out[j,1] <= 1089):
				# this doesn't work, the reason for this is the dataframe doesn't work like this. Trying again tomorrow.
#				df = pd.DataFrame({'wavenumber': out[j,0],'redshift': out[j,1],'conformal time': out[j,2],'E_2 (z)': out[j,3]})
#				dr = dr.append(df,ignore_index=True)
				g.write("%11.5E %11.5E %11.5E %11.5E\n" % (out[j,0], out[j,1], out[j,2], out[j,3]))
					
	j = j + 1
	# instead of redshift want conformal time. Do that soon. 

	# need to cut the output want to have redshifts from 100 - 1000
	i = i + 1
	cosmo.struct_cleanup()

	cosmo.empty()

# I oversampled the K's a more appropiate way is to + 0.0005 instead of +0.0001
#dr.to_csv('kpydata100.txt', sep='\t')
g.close()
cosmo.struct_cleanup()
cosmo.empty()
