import numpy as np

table = np.loadtxt('flowdata.txt')

#**********************************************************************************
#define the row we want, basically which shell matters. Max should be close to 133, i.e. Rc/mach number 
irow = 133

new_table = np.zeros((5001,10))
i = 0
a = 0
j = 0
b = 0
for j in range(0,5001):
	i = 0	
	for i in range(0,1201):
		#careful with the number here,  should be irow -1 
		if (i == 133):
			new_table[a,:] = table[i + 1201*j,:]
			a = a + 1
			print 'Vamos por a %d y i %d \n' %(a, i)
		i = i + 1
	j = j + 1
	b = 1

print new_table

np.savetxt('plot_temp_t.dat', new_table, fmt='%f', delimiter=' ')
