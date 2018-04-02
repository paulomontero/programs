import sys
import os
import numpy as np

# this is the actual script to get power spectrums, put maximum redshift as input, you get it from the boxes
z_max = float(sys.argv[1])

i = 0
z_new = z_max
z_array = []

for i in range (0,84): 

	if (i == 2 or i == 5 or i == 21 or i == 26 or i == 45 or i == 69 or i == 75):
		z_new = z_new + 0.01		
		z_var = '%06.2f' % z_new 
		cmd= "./delta_ps updated_smoothed_deltax_z{z_var}_256_300Mpc xH_nohalos_z{z_var} power_default1_z{z_var}.dat {z_var}".format(z_var=z_var)
		z_array.append(z_new)
	
		os.system(cmd)
		z_new = z_new - 0.01
	else:
		z_var = '%06.2f' % z_new
	
        	cmd= "./delta_ps updated_smoothed_deltax_z{z_var}_256_300Mpc xH_nohalos_z{z_var} power_default1_z{z_var}.dat {z_var}".format(z_var=z_var)
		z_array.append(z_new)
	
		os.system(cmd)

	
	print 'Vamos por %d o z = %06.2f\n' % (i,z_new)
	z_new = (1.000 + z_new)/(1.020) - 1.000 # 1.02 is the zprimefactor in heat_params
	i = i + 1

np.savetxt('zarray.dat', z_array, fmt='%06.2f')
