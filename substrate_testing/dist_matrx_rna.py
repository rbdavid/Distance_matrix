#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# -----------------------------------
import MDAnalysis
import numpy as np
import sys
from numpy.linalg import *
from distance_functions import *

# -----------------------------------

pdb = sys.argv[1]
prmtop = sys.argv[2]
traj_loc = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])
system = sys.argv[6]

zeros = np.zeros
square = np.square
sqrt = np.sqrt
flush = sys.stdout.flush

def ffprint(string):
        print '%s' %(string)
        flush()

pro = 'protein'
nucleic = 'nucleic or resname A5 A3 U5'

sites = ['base','sugar','phosphate']

# SIMPLE SELECTION STRINGS TO USE FOR THE DESIRED SELECTIONS
sugar = " name C5' H5' H5'' C4' H4' O4' C1' H1' C3' H3' C2' O2' HO2' "	# DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...
sugar_5= " name HO5' O5' or " + sugar
sugar_3= sugar + " O3' HO3' "
base = ' name N9 C8 H8 N7 C5 C6 N6 H61 H62 N1 C2 H2 N3 C4 O6 H1 H21 H22 H6 H5 N4 H41 H42 C2 O2 O4 H3'	# selection string that will select all appropriate atoms for any of the nucleic residues...

#phos = 'nucleic and name (P or OP1 or OP2)'
#phosphate = " name O5' or name P or name OP1 or name OP2 or name O3' "		# NOT WORKING SINCE THE O3' and O5' atoms are part of different residues... will need to figure this one out still...

u = MDAnalysis.Universe(prmtop,pdb)
protein_sel = u.select_atoms(pro)
nucleic_sel = u.select_atoms(nucleic)

nProt_res = protein_sel.n_residues
nNucl_res = nucleic_sel.n_residues

selection_output = open('selection_order.output','w')
rna_sel_list = []
for i in range(len(sites)):
	if sites[i] == 'base':
		selection_output.write('### BASE ###')
		sel_string = base 
		for j in range(nNucl_res):
			rna_sel_list.append(nucleic_sel.residues[j].select_atoms(base))

	elif sites[i] == 'sugar':
		selection_output.write('### SUGAR ###')
		for j in range(nNucl_res):
			if j == 0:
				rna_sel_list.append(nucleic_sel.residues[j].select_atoms(sugar_5))
			elif j == range(nNucl_res)[-1]:
				rna_sel_list.append(nucleic_sel.residues[j].select_atoms(sugar_3))
			else:
				rna_sel_list.append(nucleic_sel.residues[j].select_atoms(sugar))

	elif sites[i] == 'phosphate':
		selection_output.write('### PHOSPHATE ###')
		for j in range(nNucl_res):
			if j == 0:
				continue
			else:
				res_num = nucleic_sel.residues[j].resid
				temp_O3 = nucleic_sel.residues[j-1][-1].index
				if temp_O3.name != "O3'":
					print 'Something is fucky'
					sys.exit()
				rna_sel_list.append(nucleic_sel.select_atoms(" (resid %s and name P OP1 OP2 O5') or bynum %s" %(res_num,temp_O3+1)))

rna_sel_list = []
for i in range(nNucl_res):
	temp_resname = nucleic_sel.residues[i].resname
	print 'Resname:', temp_resname
	rna_sel_list.append(nucleic_sel.residues[i].select_atoms(base))		# GETTING A DEPRECIATION WARNING; need to figure out what to do about this
	print nucleic_sel.residues[i].select_atoms(base).atoms
	if temp_resname in ['A5','U5','C5','G5']:
		rna_sel_list.append(nucleic_sel.residues[i].select_atoms(sugar_5))
		print nucleic_sel.residues[i].select_atoms(sugar_5).atoms
	elif temp_resname in ['A3','U3','C3','G3']:
		rna_sel_list.append(nucleic_sel.residues[i].select_atoms(sugar_3))
		print nucleic_sel.residues[i].select_atoms(sugar_3).atoms
	else:
		rna_sel_list.append(nucleic_sel.residues[i].select_atoms(sugar))
		print nucleic_sel.residues[i].select_atoms(sugar).atoms

###	rna_sel_list.append(phosphate_sel) 		# NEED TO FIGURE OUT THE PHOSPHATE SELECTION STRING...

nSubstrate_sels = len(rna_sel_list)
avg_matrix = zeros((nProt_res,nSubstrate_sels))
std_matrix = zeros((nProt_res,nSubstrate_sels))

nSteps = 0
while start <= end:
        u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
        nSteps += len(u.trajectory)
        for ts in u.trajectory:
                if ts.frame%1000 == 0:
                        ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))
		temp_prot_com = zeros((nProt_res,3))
                for i in range(nProt_res):
			temp_prot_com[i] = protein_sel.residues[i].center_of_mass()
		
		for i in range(nProt_res):
			for j in range(nSubstrate_sels):
				dist, dist2 = euclid_dist(temp_prot_com[i], rna_sel_list[j].center_of_mass())
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

