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

change_dir = os.chdir

# ----------------------------------------
# SUBROUTINES:

# ----------------------------------------
# MAIN PROGRAM:

window_data = np.loadtxt(window_file)

nWindows = len(window_data)

for i in range(nWindows-1):
	change_dir('%03d.%03d.distance_matrix' %(window_data[i][0],window_data[i][1]))
	avg_file = '%03d.%03d.avg_distance_matrix.dat' %(window_data[i][0],window_data[i][1])
	std_file = '%03d.%03d.std_distance_matrix.dat' %(window_data[i][0],window_data[i][1])
	avg_data1 = np.loadtxt(avg_file).flatten()
	std_data1 = np.loadtxt(std_file).flatten()
	change_dir('..')

	for j in range(i+1,nWindows):
		change_dir('%03d.%03d.distance_matrix' %(window_data[j][0],window_data[j][1]))
		avg_file = '%03d.%03d.avg_distance_matrix.dat' %(window_data[j][0],window_data[j][1])
		std_file = '%03d.%03d.std_distance_matrix.dat' %(window_data[j][0],window_data[j][1])
		avg_data2 = np.loadtxt(avg_file).flatten()
		std_data2 = np.loadtxt(std_file).flatten()
		change_dir('..')
		
		plot_1d(avg_data1[:],avg_data2[:],'k.','Distance (Traj %03d - %03d)' %(window_data[i][0],window_data[i][1]),'Distance (Traj %03d - %03d)' %(window_data[j][0],window_data[j][1]),'avg.windows_%d_%d' %(i,j),'q_q_plot',draw_line=1)
		plot_1d(std_data1[:],std_data2[:],'k.','Distance (Traj %03d - %03d)' %(window_data[i][0],window_data[i][1]),'Distance (Traj %03d - %03d)' %(window_data[j][0],window_data[j][1]),'std.windows_%d_%d' %(i,j),'q_q_plot',draw_line=1)

