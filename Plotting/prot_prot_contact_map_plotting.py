#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

from plotting_functions import *
import sys
import os

start = int(sys.argv[1])
end = int(sys.argv[2])
step = int(sys.argv[3])

change_dir = os.chdir

j = 0
for i in range(start,end,step):
	j += step
	change_dir('%03d.%03d.distance_matrix' %(i,j))
	avg_file = '%03d.%03d.avg_distance_matrix.dat' %(i,j)
	std_file = '%03d.%03d.std_distance_matrix.dat' %(i,j)
	
	avg_data = np.loadtxt(avg_file)
	std_data = np.loadtxt(std_file)
	
	nRes = len(avg_data)
	
	if nRes != len(avg_data[0]):
		print 'length of matrix data is not nRes X nRes; something is fucked up'
		sys.exit()
	
	for x in range(nRes-1):
		for y in range(x+1,nRes):
			avg_data[y][x] = avg_data[x][y]
			std_data[y][x] = std_data[x][y]
	
	matrix2d(avg_data,'Prot Residue Number','Prot Residue Number','Distance','avg','test_system',vmin=0.001,vmax=45,plt_title='Avg COM-COM Distance (Traj %03d - %03d)' %(i,j),cb_units='$\AA$',xlim=(0,451),ylim=(0,451))
	matrix2d(std_data,'Prot Residue Number','Prot Residue Number','Distance','std','test_system',vmin=0.001,vmax=6,plt_title='St. Dev. of COM-COM Distance (Traj %03d - %03d)' %(i,j),cb_units='$\AA$',xlim=(0,451),ylim=(0,451))

	print i, j

	change_dir('..')

