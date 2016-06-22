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

descriptor = sys.argv[1]	# descriptor to decide which trajectory group to analyze for all systems

sys_list = []
#sys_list.append(['AMBER_apo'])
#sys_list.append(['AMBER_atp'])
sys_list.append(['AMBER_ssrna'])
sys_list.append(['AMBER_ssrna_atp'])
sys_list.append(['AMBER_ssrna_adp_pi'])
sys_list.append(['AMBER_ssrna_adp'])
#sys_list.append(['AMBER_ssrna_pi'])

nSys = len(sys_list)

# ----------------------------------------
# SUBROUTINES:

# ----------------------------------------
# MAIN PROGRAM:

for i in range(nSys-1):
	avg1 = np.loadtxt('../../%s/Distance_matrix/%s.avg_distance_matrix.dat' %(sys_list[i][0],descriptor))
	std1 = np.loadtxt('../../%s/Distance_matrix/%s.std_distance_matrix.dat' %(sys_list[i][0],descriptor))

	for j in range(i+1,nSys):
		avg2 = np.loadtxt('../../%s/Distance_matrix/%s.avg_distance_matrix.dat' %(sys_list[j][0],descriptor))
		std2 = np.loadtxt('../../%s/Distance_matrix/%s.std_distance_matrix.dat' %(sys_list[j][0],descriptor))

		avg_data = avg1 - avg2
		std_data = std1 - std2

		out1 = open('AVG_dif.%s.%s.output' %(sys_list[i][0],sys_list[j][0]),'w')
		out2 = open('AVG_dif.%s.%s.hist.dat' %(sys_list[i][0],sys_list[j][0]),'w')
		count_array = np.array(len(avg1))
		for x in range(20,len(avg1)-1):
			for y in range(x+1,len(avg1)):
				if abs(avg_data[x][y]) > 5.0:
					out1.write('%s   %s   %03d (%d)   %03d (%d)   %f\n' %(sys_list[i][0],sys_list[j][0],x+1,x+168,y+1,y+168,avg_data[x][y]))
					count_array[x] += 1
					count_array[y] += 1
		for x in range(len(avg1)):
			out2.write('%03d   %03d   %d\n' %(x+1,x+168,count_array[x]))

		out1.close()
		out2.close()

		plt.bar(range(len(avg1)+168),count_array[:])

		matrix2d(abs(avg1-avg2),'Residue Number','Residue Number','Distance','%s.%s.%s' %(descriptor,sys_list[i][0],sys_list[j][0]),'difference_avg',cb_units='$\AA$',vmax=12.)
		matrix2d(abs(std1-std2),'Residue Number','Residue Number','Distance','%s.%s.%s' %(descriptor,sys_list[i][0],sys_list[j][0]),'difference_std',cb_units='$\AA$',vmax=2.0)

