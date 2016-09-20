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

#plt.pcolor(y,x,avg_data,cmap=my_cmap,vmin=0.001,vmax=30)
	plt.pcolor(avg_data,cmap=my_cmap,vmin=0.001,vmax=15)
	cb1 = plt.colorbar(extend='max',cmap=my_cmap)
	cb1.set_label(r'Average Distance ($\AA$)' )
	plt.title('Average COM-COM Distance (Traj %03d - %03d)' %(i,j))
	plt.ylabel(r'Protein Residue Number', size=14)
	plt.ylim((0,451))
	plt.xlabel(r'Nucleic Selection', size=14)
	plt.savefig('avg.%s.%03d.%03d.pcolor.png'%(system,i,j),dpi=300)
	plt.close()

	plt.pcolor(std_data,cmap=my_cmap,vmin=0.001,vmax=6.0)
	cb1 = plt.colorbar(extend='max',cmap=my_cmap)
	cb1.set_label(r'Standard Deviation of Distance ($\AA$)')
	plt.title('Standard Dev. COM-COM Distance (Traj %03d - %03d)' %(i,j))
	plt.ylabel(r'Protein Residue Number', size=14)
	plt.ylim((0,451))
	plt.xlabel(r'Nucleic Selection', size=14)
	plt.savefig('stdev.%s.%03d.%03d.pcolor.png'%(system,i,j),dpi=300)
	plt.close()

	#matrix2d(avg_data,'Nucleic Residue Number','Protein Residue Number','Distance','avg','AMBER_ssrna',vmin=0.001,vmax=25,plt_title='Average COM-COM Residue Distance (Traj %03d - %03d)' %(i,j),cb_units='$\AA$')
	#matrix2d(std_data,'Nucleic Residue Number','Protein Residue Number','Distance','std','AMBER_ssrna',vmin=0.001,vmax=4.0,plt_title='Standard Deviation (Traj %03d - %03d)' %(i,j),cb_units='$\AA$')
	print i, j

	change_dir('..')

