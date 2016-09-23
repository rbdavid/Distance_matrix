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
selection_file = sys.argv[2]
nProt_res = int(sys.argv[3])
system = sys.argv[4]

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

sels = open(selection_file,'r')
count = 0
for line in sels:
	if not line.startswith('#'): 
		count += 1

avg_array = zeros((nProt_res,count))
std_array = zeros((nProt_res,count))

nSteps = 0
for i in range(nWindows):
	change_dir('%03d.%03d.distance_matrix' %(window_data[i][0],window_data[i][1]))

	avg_data = np.loadtxt('%03d.%03d.%s.avg_dist_mtx.dat' %(window_data[i][0],window_data[i][1],system))
	std_data = np.loadtxt('%03d.%03d.%s.stdv_dist_mtx.dat' %(window_data[i][0],window_data[i][1],system))

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

with open('%03d.%03d.avg_distance_matrix.dat' %(window_data[0][0],window_data[-1][1]),'w') as f:
	np.savetxt(f,avg_array)

with open('%03d.%03d.std_distance_matrix.dat' %(window_data[0][0],window_data[-1][1]),'w') as f:
	np.savetxt(f,std_array)

matrix2d(avg_array,'Nucleic Selection','Protein Residue Number','Distance','%03d.%03d.pro-ssRNA' %(window_data[0][0],window_data[-1][1]),'avg',vmin=0.001,vmax=45.0,plt_title='Avg COM-COM Distance (Traj %03d - %03d)' %(window_data[0][0],window_data[-1][1]),cb_units='$\AA$',ylim=(0,nProt_res),xlim=(0,count))
matrix2d(std_array,'Nucleic Selection','Protein Residue Number','Distance','%03d.%03d.pro-ssRNA' %(window_data[0][0],window_data[-1][1]),'std',vmin=0.001,vmax=6.0,plt_title='St. Dev. COM-COM Distance (Traj %03d - %03d)' %(window_data[0][0],window_data[-1][1]),cb_units='$\AA$',ylim=(0,nProt_res),xlim=(0,count))

