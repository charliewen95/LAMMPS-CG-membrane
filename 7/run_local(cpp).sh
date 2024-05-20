#!/bin/bash
##SBATCH -D  /home/charlie/
#SBATCH -J c++-comp
#SBATCH --partition=defq
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#####SBATCH --share

./1_SEARCH_EXTRACT > search.out <<-EOF
../data.production1.lip.lammpstrj
EOF

./2_rdf_height > height_rdf <<-EOF
extract.out
EOF

./3_binning_height > binned <<-EOF
height_rdf.out
EOF

echo "done"
