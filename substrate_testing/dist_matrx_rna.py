#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# -----------------------------------
import MDAnalysis
import numpy as np
import sys
from numpy.linalg import *
from distance_functions import *

# -----------------------------------

pdb = sys.argv[1]
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
system = sys.argv[5]

zeros = np.zeros
square = np.square
sqrt = np.sqrt
flush = sys.stdout.flush

def ffprint(string):
        print '%s' %(string)
        flush()

pro = 'protein'
nucleic = 'nucleic or resname A5 A3 U5 U3 G5 G3 C5 C3'

# SIMPLE SELECTION STRINGS TO USE FOR THE DESIRED SELECTIONS
sugar = " name C5' H5' H5'' C4' H4' O4' C1' H1' C3' H3' C2' O2' HO2' "	# DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...
sugar_5= " name HO5' O5' or " + sugar
sugar_3= sugar + " O3' HO3' "
base = ' name N9 C8 H8 N7 C5 C6 N6 H61 H62 N1 C2 H2 N3 C4 O6 H1 H21 H22 H6 H5 N4 H41 H42 C2 O2 O4 H3'	# selection string that will select all appropriate atoms for any of the nucleic residues...

u = MDAnalysis.Universe(pdb)
protein_sel = u.select_atoms(pro)
nucleic_sel = u.select_atoms(nucleic)

nProt_res = protein_sel.n_residues
nNucl_res = nucleic_sel.n_residues

#sites = ['base','sugar','phosphate']
#selection_output = open('selection_order.output','w')
#rna_sel_list = []
#for i in range(len(sites)):
#	if sites[i] == 'base':
#		selection_output.write('### BASE ###')
#		sel_string = base 
#		for j in range(nNucl_res):
#			rna_sel_list.append(nucleic_sel.residues[j].select_atoms(base))
#
#	elif sites[i] == 'sugar':
#		selection_output.write('### SUGAR ###')
#		for j in range(nNucl_res):
#			if j == 0:
#				rna_sel_list.append(nucleic_sel.residues[j].select_atoms(sugar_5))
#			elif j == range(nNucl_res)[-1]:
#				rna_sel_list.append(nucleic_sel.residues[j].select_atoms(sugar_3))
#			else:
#				rna_sel_list.append(nucleic_sel.residues[j].select_atoms(sugar))
#
#	elif sites[i] == 'phosphate':
#		selection_output.write('### PHOSPHATE ###')
#		for j in range(nNucl_res):
#			if j == 0:
#				continue
#			else:
#				res_num = nucleic_sel.residues[j].resid
#				temp_O3 = nucleic_sel.residues[j-1][-1].index
#				if temp_O3.name != "O3'":
#					print 'Something is fucky'
#					sys.exit()
#				rna_sel_list.append(nucleic_sel.select_atoms(" (resid %s and name P OP1 OP2 O5') or bynum %s" %(res_num,temp_O3+1)))

selection_output = open('selection_order.output','w')
rna_sel_list = []
selection_output.write('# Number  Resname  Selection_Type #')
for i in range(nNucl_res):
	temp_resname = nucleic_sel.residues[i].resname
	selection_output.write('%d   %s   Base' %(i+1,temp_resname))
	rna_sel_list.append(nucleic_sel.residues[i].select_atoms(base))		# GETTING A DEPRECIATION WARNING; need to figure out what to do about this
	if temp_resname in ['A5','U5','C5','G5']:
		rna_sel_list.append(nucleic_sel.residues[i].select_atoms(sugar_5))
		selection_output.write('%d   %s   Sugar' %(i+1,temp_resname))
	elif temp_resname in ['A3','U3','C3','G3']:
		rna_sel_list.append(nucleic_sel.residues[i].select_atoms(sugar_3))
		selection_output.write('%d   %s   Sugar' %(i+1,temp_resname))
	else:
		rna_sel_list.append(nucleic_sel.residues[i].select_atoms(sugar))
		selection_output.write('%d   %s   Sugar' %(i+1,temp_resname))
	if i == 0:
		continue
	else:
		res_num = nucleic_sel.residues[i].resid
		temp_O3 = nucleic_sel.residues[i-1][-1]
		if temp_O3.name != "O3'":
			print "The last atom in residue %d is not named O3'"
			sys.exit()
		rna_sel_list.append(nucleic_sel.select_atoms("(resid %s and name P OP1 OP2 O5') or bynum %s" %(res_num,temp_O3.index+1)))
		selection_output.write('%d   %s   Phosphate' %(i+1,temp_resname))

selection_output.close()

nSubstrate_sels = len(rna_sel_list)
avg_matrix = zeros((nProt_res,nSubstrate_sels))
std_matrix = zeros((nProt_res,nSubstrate_sels))
temp_nucl_com = zeros((nSubstrate_sels,3))

nSteps = 0
while start <= end:
        u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
        nSteps += len(u.trajectory)
        for ts in u.trajectory:
                if ts.frame%1000 == 0:
                        ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))
		
                for i in range(nSubstrate_sels):
			temp_nucl_com[i] = rna_sel_list[i].center_of_mass()
		
		for i in range(nProt_res):
			for j in range(nSubstrate_sels):
				dist, dist2 = euclid_dist(protein_sel.residues[i].center_of_mass(), temp_nucl_com[j])
				avg_matrix[i,j] += dist
				std_matrix[i,j] += dist2
	start += 1

ffprint(nSteps)
avg_matrix /= nSteps
std_matrix /= nSteps
std_matrix = sqrt(std_matrix - square(avg_matrix))

with open('%03d.%03d.%s.avg_dist_mtx.dat' %(start,end,system),'w') as f:
	np.savetxt(f,avg_matrix)

with open('%03d.%03d.%s.stdv_dist_mtx.dat' %(start,end,system),'w') as f:
	np.savetxt(f,std_matrix)

