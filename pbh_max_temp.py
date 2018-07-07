import numpy as np

table = np.loadtxt('flowdata.txt')
#table = np.loadtxt('flowdata_test.txt')

#**********************************************************************************
#define the row we want, basically which shell matters. Max should be close to 133, i.e. Rc/mach number 
#irow = 133

max_array = np.zeros(1201)
#max_array = np.zeros(1201*2)
#new_table = np.zeros((501,10))
temporary_array = np.zeros(5001)
radius_array = np.zeros(1201)
#radius_array = np.zeros(1201*2)


b = 0
irow = 0
i = 0
for irow in range (0,1201):
#for irow in range (0,1201*2):
	a = 0
	j = 0
	for j in range(0,5001):
		i = 0	
		for i in range(0,1201):
		#for i in range(0,1201*2):
			if (j == 0 and b == 0):
				radius_array[i] = table[i,2]
				
			#careful with the number here,  should be irow -1 
			if (i == irow):
				
				temporary_array[a] = table[i + 1201*j, 9]
			#	temporary_array[a] = table[i + 1201*j*2, 9]

#				print 'Vamos por a %d y i %d \n' %(a, i)

				a = a + 1
			i = i + 1
		j = j + 1
		b = 1
		max_array[irow] = np.amax(temporary_array)
	irow = irow + 1

m = 0
 
shell_array = np.zeros(1201)
#shell_array = np.zeros(1201*2)

for m in range (0,1201):
#for m in range (0,1201*2):
	shell_array[m] = m
	m = m + 1
print shell_array, radius_array


np.savetxt('plot_temp_max.dat', np.transpose([shell_array, max_array, radius_array]), fmt='%f', delimiter=' ')
#let's try this stuff
#np.savetxt('plot_temp_t.dat', new_table, fmt='%f', delimiter=' ')
