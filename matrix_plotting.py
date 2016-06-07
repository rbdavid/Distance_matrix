#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

from plotting_functions import *

avg_file = sys.argv[1]
std_file = sys.argv[2]

avg_data = np.loadtxt(avg_file)
std_data = np.loadtxt(std_file)

nRes = len(avg_data)

if nRes != len(avg_data[0]):
	print 'length of matrix data is not nRes X nRes; something is fucked up'
	sys.exit()

for i in range(nRes-1):
	for j in range(i+1,nRes):
		avg_data[j][i] = avg_data[i][j]
		std_data[j][i] = std_data[i][j]

matrix2d(avg_data,'Residue Number','Residue Number','avg','test_system')
matrix2d(std_data,'Residue Number','Residue Number','std','test_system')


#my_cmap = plt.cm.get_cmap('jet')
#my_cmap.set_under('w')
#
#plt.imshow(avg_data,cmap=my_cmap,vmin=0.001,interpolation='None',origin='lower')
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.xlabel('Residue Number')
#plt.ylabel('Residue Number')
#plt.savefig('avg_data.png')
#
#plt.imshow(std_data,cmap=my_cmap,vmin=0.001,interpolation='None',origin='lower')
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.xlabel('Residue Number')
#plt.ylabel('Residue Number')
#
#plt.savefig('std_data.png')

