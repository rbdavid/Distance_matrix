#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from distance_functions import *
import scipy.stats

# ----------------------------------------
# VARIABLE DECLARATION

pdb_file = sys.argv[1]
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

shapiro_test = True
nFrames = 2500

zeros = np.zeros
square = np.square
sqrt = np.sqrt
flush = sys.stdout.flush

important = 'protein'

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

# ----------------------------------------
# MAIN PROGRAM:

u = MDAnalysis.Universe(pdb_file)
u_important = u.select_atoms(important)

nRes = len(u_important.residues)
ffprint(nRes)

avg_matrix = zeros((nRes,nRes))
std_matrix = zeros((nRes,nRes))
temp_prot_com = zeros((nRes,3))
if shapiro_test:
	shapiro_array = zeros((nFrames,nRes,3))

nSteps = 0
count = 0
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
	nSteps += len(u.trajectory)

	for ts in u.trajectory:
		if ts.frame%1000 == 0:
			ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))

		for i in range(nRes):
			temp_prot_com[i] = u_important.residues[i].center_of_mass()
		
		if shapiro_test and ts.frame%10 == 0:
			for i in range(nRes):
				shapiro_array[count,i,:] = temp_prot_com[i]
			count += 1

		for i in range(nRes-1):
			for j in range(i+1,nRes):
				dist, dist2 = euclid_dist(temp_prot_com[i],temp_prot_com[j])
				avg_matrix[i,j] += dist
				std_matrix[i,j] += dist2
	start +=1

ffprint(nSteps)

avg_matrix /= nSteps
std_matrix /= nSteps
std_matrix = sqrt(std_matrix - square(avg_matrix))

with open('%03d.%03d.avg_distance_matrix.dat' %(int(sys.argv[3]),end),'w') as f:
	np.savetxt(f,avg_matrix)

with open('%03d.%03d.std_distance_matrix.dat' %(int(sys.argv[3]),end),'w') as f:
	np.savetxt(f,std_matrix)

if shapiro_test:
	with open('shapiro_wilks_test.dat','w') as f:
		for i in range(nRes):
			for j in range(3):
				W, P = scipy.stats.shapiro(shapiro_array[:,i,j])
				f.write('%d   %f   %f   ' %(j,W,P))
			f.write('\n')

