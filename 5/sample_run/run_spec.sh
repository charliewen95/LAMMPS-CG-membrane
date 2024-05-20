#!/bin/bash
##SBATCH -D  /data/charlie/LAMMPS/post_process/jesper_code_eg
#SBATCH -J python-comp
#SBATCH --partition=defq
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#####SBATCH --share

exe=/home/jesper/bin/Spectrum

$exe traj=/data/charlie/LAMMPS/post_defense_final/BEND_1.25_kc/BEND_1.25/data.production1.lip.lammpstrj head_type=1 tail_type=3 max_k=10 max_l=10 gx=20 gy=20 scale=0.1 delta_q=0.05 start=0 stop=10000 | tee spectrum.log
mv spectrum.modified.txt spectrum.modified.txt
mv spectrum.standard.txt spectrum.standard.txt
