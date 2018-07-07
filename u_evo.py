import numpy as np
import sys
import os

table = np.loadtxt('flowdata.txt')

t_array = np.zeros(5001)
u_shell_133_array = np.zeros(5001)
u_shell_1201_array = np.zeros(5001)
u_shell_600_array = np.zeros(5001)
u_shell_300_array = np.zeros(5001)
u_shell_900_array = np.zeros(5001)
u_shell_1000_array = np.zeros(5001)
i = 0
for i in range(0,5001):
	u_shell_133_array[i] = table[133 +i*1201,6]
	u_shell_300_array[i] = table[300 + i*1201,6]
	u_shell_600_array[i] = table[600 + i*1201,6]
	u_shell_900_array[i] = table[900 + i*1201,6]
	u_shell_1000_array[i] = table[1000 + i*1201,6]
	u_shell_1201_array[i] = table[1200 +i*1201,6]
	t_array[i] = i*4e-3
	i = i + 1



np.savetxt('plot_u_evo.dat', np.transpose([t_array, u_shell_133_array, u_shell_300_array, u_shell_600_array, u_shell_900_array, u_shell_1000_array, u_shell_1201_array]), fmt='%f', delimiter=' ')
