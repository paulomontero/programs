import numpy as np
import sys
import os

#no need for mass of the PBH. Need for mach number only
mach_number = float(sys.argv[1])

table = np.loadtxt('flowdata.txt')

#solution for the conserved equation, z speed
sol_vz_array = np.zeros(1201)
r_final_array = np.zeros(1201)
r_initial_array = np.zeros(1201)


#grabbing everyone
i = 0
for i in range(0,1201):
	sol_vz_array[i] = (table[i + 1201*5000,6] + table[i + 1201*5000,7]/(table[i +1201*5000,5]*1.0)- 3.0)/(mach_number*1.0)
	r_final_array[i] = table[i + 1201*5000, 2]
 	r_initial_array[i] = table[i,2]
	i = i + 1			 
	#print j, E_init_array[j-133]

shell_array = np.zeros(1201)

a = 0
for a in range(0,1201):
	shell_array[a] = a
	a = a + 1 



np.savetxt('new_plot_vz_tstep5001.dat', np.transpose([shell_array, r_final_array, sol_vz_array, r_initial_array]), fmt='%f', delimiter=' ')


#differential part
last_table = np.loadtxt('new_plot_vz_tstep5001.dat')

delta_vz_array  = np.diff(last_table[:,2])
delta_r_array  = np.diff(last_table[:,1])
result_array = delta_vz_array/delta_r_array

ishell_array = np.zeros(1200)
b = 0
for b in range(0,1200):
	ishell_array[b] = b
	b = b + 1 
#print result_array 

np.savetxt('grad_vz.dat', np.transpose([ishell_array, result_array]), fmt='%f')

