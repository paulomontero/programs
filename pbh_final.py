import numpy as np

table = np.loadtxt('flowdata.txt')
#table = np.loadtxt('flowdata_test.txt')


radius_array = np.zeros(1201)
#radius_array = np.zeros(1201*2)

final_array = np.zeros(1201)
#final_array = np.zeros(1201*2)

shell_array = np.zeros(1201)
#shell_array = np.zeros(1201*2)
j = 0 
for j in range(0,5001):
	i = 0 
	for i  in range(0,1201):
	#for i  in range(0,1201*2):
		if (j == 0):
			radius_array[i] = table[i,2]

		if (j == 5000):
			final_array[i] = table[i + 1201*j,9]
			shell_array[i] = table[i + 1201*j,1]
			#final_array[i] = table[i + 1201*j*2,9]
			#shell_array[i] = table[i + 1201*j*2,1]
		i = i + 1

	j = j + 1

print final_array, shell_array, radius_array



np.savetxt('plot_temp_final.dat', np.transpose([shell_array, final_array, radius_array]), fmt='%f', delimiter=' ')
#np.savetxt('plot_temp_final_test.dat', np.transpose([shell_array, final_array, radius_array]), fmt='%f', delimiter=' ')
