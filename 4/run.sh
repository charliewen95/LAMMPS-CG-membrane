#!/bin/bash
#SBATCH -J LAMMPS
#SBATCH -o prot.o%j

##SBATCH -p skx-dev
##SBATCH -N 2
##SBATCH -n 96   # Total number of MPI tasks for  all nodes

#SBATCH -p skx
#SBATCH -N 4
#SBATCH -n 192   # Total number of MPI tasks for  all nodes

#SBATCH -t 24:00:00
#SBATCH -A TG-BIO240033
##SBATCH --mail-user=charliewen95@hotmail.com
##SBATCH --mail-type=ALL

LAMMPS="/home1/08207/tg875064/lammps/src/lmp_mpi"

DATA_FILE="2piezo.config"

BEAD_SIZE=7.5
STICKYNESS=7.0
TAIL_BEADS=3
BEND=1.0
#BEND1=2
#BEND2=20

N_STEPS=0

CTRL_FILE="standard_init.control"
OUTPUT_PREFIX="data.production.init"

mpiexec -np 192 ${LAMMPS} -in ${CTRL_FILE} -var DATA_FILE ${DATA_FILE} \
	-var DELTA_T       50     \
	-var T_DAMP        100000 \
	-var P_DAMP        100000 \
	-var P_EXT         0.0    \
	-var TRAJ_EVERY    10000   \
	-var N_STEPS       ${N_STEPS}       \
	-var OUTPUT_PREFIX ${OUTPUT_PREFIX} \
	-var BEAD_SIZE     ${BEAD_SIZE}     \
	-var STICKYNESS    ${STICKYNESS}    \
	-var TAIL_BEADS    ${TAIL_BEADS}    \
	-var BEND          ${BEND}         \
#	-var BEND1         ${BEND1}         \
#	-var BEND2         ${BEND2}
