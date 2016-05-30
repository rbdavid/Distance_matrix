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
from MDAnalysis.analysis.align import *
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

pdb_file = sys.argv[1]
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
#out_file = sys.argv[5]

zeros = np.zeros


important = 'protein'

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

# ----------------------------------------
# MAIN PROGRAM:

u = MDAnalysis.Universe(pdb_file)
u_all = u.select_atoms('all')

u_important = u.select_atoms(important)

nRes = len(u_important.residues)

ffprint(nRes)

avg_matrix = zeros((nRes,nRes))
std_matrix = zeros((nRes,nRes))

nSteps = 0
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))

	for ts in u.trajectory:

		for i in range(nRes-1):
			res0 = u_important.residues[i]
			com0 = res0.center_of_mass()
			
			for j in range(i+1,nRes):
				res1 = u_important.residues[j]
				com1 = res1.center_of_mass()

				avg_matrix[i,j], std_matrix[i,j] += Euclid_distance(com0,com1,dist2=True) 

		nSteps += 1

avg_matrix /= nSteps
std_matrix /= nSteps
std_matrix -= avg_matrix
std_matrix = sqrt(std_matrix)

out1 = open('average_distance_matrix.dat','w')
out2 = open('std_distance_matrix.dat','w')

for i in range(nRes):
	for j in range(nRes):
		out1.write('%f   ' %(avg_matrix[i,j]))
		out2.write('%f   ' %(std_matrix[i,j]))
	out1.write('\n')
	out2.write('\n')
