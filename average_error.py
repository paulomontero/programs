import numpy as np

realization_array = np.zeros(8)
realization_array[0] = 6.17
realization_array[1] = 6.28
realization_array[2] = 6.24
realization_array[3] = 6.25
realization_array[4] = 6.17
realization_array[5] = 6.46
realization_array[6] = 6.32
realization_array[7] = 6.43

average = (realization_array[0] + realization_array[1] + realization_array[2] + realization_array[3] + realization_array[4] + realization_array[5] + realization_array[6] + realization_array[7])/8.0

error = np.sqrt(((realization_array[0] - average)**2 + (realization_array[1] - average)**2 + (realization_array[2] - average)**2 + (realization_array[3] - average)**2 + (realization_array[4] - average)**2 + (realization_array[5] - average)**2 + (realization_array[6] - average)**2 + (realization_array[7] - average)**2)/(8.0))

print "The average for 1D, with bias ratio, of zobs = 4.0 is %f, the error bar is %f \n" % (average, error)
  
