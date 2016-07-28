#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

import numpy as np
import sys
import os
from plotting_functions import *

data_file = sys.argv[1]
dist_cutoff = float(sys.argv[2])

data = np.loadtxt(data_file)

nRows = len(data)
nCols = len(data[0])

if nRows != nCols:
	print 'This is not a residue-residue pair contact map dataset; number of rows and columns are not equal'
	sys.exit()

for i in range(nRows-1):
	for j in range(i,nCols):
		if data[i][j] < dist_cutoff:
			data[i][j] = float(1)
			data[j][i] = float(1)
		else: 
			data[i][j] = float(0)
			data[j][i] = float(0)

with open('functionalized_dist_matrix.dat','w') as f:
	np.savetxt(f,data)

matrix2d(data,'Residue Number','Residue Number','Physical Contact','contact','test_system')

