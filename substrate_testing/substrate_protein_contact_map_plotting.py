#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

from plotting_functions import *
import sys
import os

start = int(sys.argv[1])
end = int(sys.argv[2])
step = int(sys.argv[3])
system = sys.argv[4]

change_dir = os.chdir

my_cmap = plt.cm.get_cmap('jet')

j = 0
for i in range(start,end,step):
	j += step
	change_dir('%03d.%03d.distance_matrix' %(i,j))
	avg_file = '%03d.%03d.%s.avg_dist_mtx.dat' %(i,j,system)
	std_file = '%03d.%03d.%s.stdv_dist_mtx.dat' %(i,j,system)
	
	avg_data = np.loadtxt(avg_file)
	std_data = np.loadtxt(std_file)

	nProt_res = len(avg_data)
	nSub_res = len(avg_data[0])

	matrix2d(avg_data,'Nucleic Selection','Protein Residue Number','Distance','avg.%s' %(system),'%03d.%03d.pro-ssRNA' %(1,5),cb_units='$\AA$',plt_title='Avg COM-COM Distance (Traj %03d - %03d)' %(1,5),xlim=(0,nSub_res),ylim=(0,nProt_res),vmax=45.0)
	matrix2d(std_data,'Nucleic Selection','Protein Residue Number','St. Dev. of Distance','stdev.%s' %(system),'%03d.%03d.pro-ssRNA' %(1,5),cb_units='$\AA$',plt_title='St. Dev. COM-COM Distance (Traj %03d - %03d)' %(1,5),xlim=(0,nSub_res),ylim=(0,nProt_res),vmax=6.0)

	print i, j

	change_dir('..')

