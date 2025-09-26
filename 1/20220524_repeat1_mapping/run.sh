#! /bin/bash

#
# EXAMPLE SCRIPT TO DO CG-MAPPING OF PIEZO2 REPEAT-1
# Jesper Madsen, May 24, 2022
#
# *** For illustrative purposes only ***

inpdb=piezo2_repeat1.pdb
scr_run=run.tcl
scr_parse=parse_csv.py
cg=cg.xyz

#use VMD, remember to cite!
vmd="/Applications/VMD_1.9.4a55-x86_64-Rev11.app/Contents/vmd/vmd_MACOSXX86_64"

# Parse selections, outputting the centres-of-mass of the corresponding CG beads:
python3 $scr_parse > $scr_run
$vmd $inpdb -dispdev text -e $scr_run

# Inspect, again using VMD
$vmd -m $cg $inpdb
