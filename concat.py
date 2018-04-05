import numpy as np
import scipy as sp
import pandas as pd
import glob
import sys

# this guy combines all tables into one for one realization
path = r'/home/yoshi/21cmFAST-master/Log_files_default1/crosspower/'
allFiles = sorted(glob.glob(path + "/*.dat"))
frame =pd.DataFrame()
list_ = []

for file_ in allFiles:
	table_ = np.loadtxt(file_)
	df = pd.DataFrame(table_)
	list_.append(df)
	#df = pd.read_table(file_)
	#list_.append(df)

frame = pd.concat(list_)
#table1 = np.loadtxt('power_default1_z005.90.dat')
#df1 = pd.DataFrame(table1)

#table2 = np.loadtxt('power_default1_z006.04.dat')
#df2 = pd.DataFrame(table2)
#frames = [df1, df2]
#result = pd.concat(frames)

#MAKE SURE  TO CHANGE THE OUTPUT
np.savetxt('power_default1', frame, fmt='%e')


