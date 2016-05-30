#!/bin/bash

PDB_LOC='~/Projects/Molecular_Machines/Helicase_DNS3/Analysis/AMBER_apo/'
TRAJ_LOC='~/Projects/Molecular_Machines/Helicase_DNS3/Analysis/AMBER_apo/'
NPRODS=120

for ((prod=1;prod<=$NPRODS;prod+=5))
do
	((a=$prod+4))
	time ./matrix_calc.py $pdb_loc $traj_loc $prod $a 
done

