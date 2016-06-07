#!/bin/bash
#SBATCH --job-name=dist_cal.a_apo
#SBATCH --output=dist_cal.a_apo.output
#SBATCH --time=96:00:00 
#SBATCH --nodes=1
#SBATCH --exclusive

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/software/usr/gcc-4.9.2/lib64"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/software/usr/hpcx-v1.2.0-292-gcc-MLNX_OFED_LINUX-2.4-1.0.0-redhat6.6/ompi-mellanox-v1.8/lib"

export PYTHON_EGG_CACHE="./"

PDB_LOC='~/Projects/Molecular_Machines/Helicase_DNS3/Analysis/AMBER_apo/truncated.pdb'
TRAJ_LOC='~/Projects/Molecular_Machines/Helicase_DNS3/Analysis/AMBER_apo/Truncated'
NPRODS=150
NCPUS=20

prod=1
for ((i=1;i<=2;i++))
do
	j=1
	while ((j <= $NCPUS)) && ((prod <= $NPRODS))
	do
		echo $j $i $prod
		((a=$prod+4))
		printf -v x "%03d" $prod
		printf -v y "%03d" $a
		mkdir $x.$y.distance_matrix
		cd $x.$y.distance_matrix
		time ../matrix_calc.py $PDB_LOC $TRAJ_LOC $prod $a > dist_calc.output & 
		cd ../
		((j=$j+1))
		((prod=$prod+5))
	done
	wait
done

