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
sys_list.append(['AMBER_apo'])
sys_list.append(['AMBER_atp'])
sys_list.append(['AMBER_ssrna'])
sys_list.append(['AMBER_ssrna_atp'])
sys_list.append(['AMBER_ssrna_adp_pi'])
sys_list.append(['AMBER_ssrna_adp'])
sys_list.append(['AMBER_ssrna_pi'])

nSys = len(sys_list)

# ----------------------------------------
# SUBROUTINES:

# ----------------------------------------
# MAIN PROGRAM:

for i in range(nSys-1):
	avg1 = np.loadtxt('../../%s/Distance_matrix/%s.avg_distance_matrix.dat' %(sys_list[i],descriptor))
	std1 = np.loadtxt('../../%s/Distance_matrix/%s.std_distance_matrix.dat' %(sys_list[i],descriptor))

	for j in range(i+1,nSys):
		avg2 = np.loadtxt('../../%s/Distance_matrix/%s.avg_distance_matrix.dat' %(sys_list[j],descriptor))
		std2 = np.loadtxt('../../%s/Distance_matrix/%s.std_distance_matrix.dat' %(sys_list[j],descriptor))

		matrix2d(avg1-avg2,'Residue Number','Residue Number','Distance','%s.%s.%s' %(descriptor,sys_list[i],sys_list[j]),'difference_avg',cb_units='$\AA$')
		matrix2d(std1-std2,'Residue Number','Residue Number','Distance','%s.%s.%s' %(descriptor,sys_list[i],sys_list[j]),'difference_std',cb_units='$\AA$')




