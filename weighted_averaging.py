#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
import os
from plotting_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

window_file = sys.argv[1]
nRes = int(sys.argv[2])

zeros = np.zeros
change_dir = os.chdir
square = np.square
sqrt = np.sqrt

# ----------------------------------------
# SUBROUTINES:

# ----------------------------------------
# MAIN PROGRAM:

window_data = np.loadtxt(window_file)

nWindows = len(window_data)

avg_array = zeros((nRes,nRes))
std_array = zeros((nRes,nRes))

nSteps = 0
for i in range(nWindows):
	change_dir('%03d.%03d.distance_matrix' %(window_data[i][0],window_data[i][1]))
	avg_file = '%03d.%03d.avg_distance_matrix.dat' %(window_data[i][0],window_data[i][1])
	std_file = '%03d.%03d.std_distance_matrix.dat' %(window_data[i][0],window_data[i][1])
	avg_data = np.loadtxt(avg_file)
	std_data = np.loadtxt(std_file)

	avg_data *= window_data[i][2]
	std_data = square(std_data)
	std_data *= window_data[i][2]

	avg_array += avg_data
	std_array += std_data
	nSteps += window_data[i][2]
	
	change_dir('..')

avg_array /= nSteps
std_array /= nSteps
std_array = sqrt(std_array)

for i in range(nRes-1):
	for j in range(i+1,nRes):
		avg_array[j][i] = avg_array[i][j]
		std_array[j][i] = std_array[i][j]

out1 = open('%03d.%03d.avg_distance_matrix.dat' %(window_data[0][0],window_data[-1][1]),'w')
out2 = open('%03d.%03d.std_distance_matrix.dat' %(window_data[0][0],window_data[-1][1]),'w')
for i in range(nRes):
	for j in range(nRes):
		out1.write('%f   ' %(avg_array[i][j]))
		out2.write('%f   ' %(std_array[i][j]))
	out1.write('\n')
	out2.write('\n')
out1.close()
out2.close()

matrix2d(avg_array,'Residue Number','Residue Number','Distance','avg','%03d.%03d' %(window_data[0][0],window_data[-1][1]),vmin=0.001,vmax=45,plt_title='Average COM-COM Residue Distance (Traj %03d - %03d)' %(window_data[0][0],window_data[-1][1]),cb_units='$\AA$')
matrix2d(std_array,'Residue Number','Residue Number','Distance','std','%03d.%03d' %(window_data[0][0],window_data[-1][1]),vmin=0.001,vmax=6,plt_title='Standard Deviation (Traj %03d - %03d)' %(window_data[0][0],window_data[-1][1]),cb_units='$\AA$')

