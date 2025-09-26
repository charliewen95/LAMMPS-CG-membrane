#!/usr/bin/env python2

import sys
sys.path.append('../src/')
import cgmap as cg
import mdtraj as md
import md_check as check

############################### config #####################################

input_traj = "6kg7-3_bb.pdb"
input_top  = "6kg7-3_bb.pdb"

output_traj = "6kg7_full-3_cg.pdb"
output_top  = "6kg7_full-3_cg.pdb"

reference_traj = "6kg7_bb.pdb"
reference_top  = "6kg7_bb.pdb"

output_dir ='./output/'
input_dir  ='./input/'
reference_dir ='./reference/'


############################### run ########################################

### pull in trajectories
trj = md.load(input_dir + input_traj,top=input_dir + input_top)

### define mapping based on knowledge of topology
### in this instance, map every residue into a single site
for a in trj.top.atoms: a.mass = a.element.mass
for a in trj.top.atoms: a.charge = 0

# first residue is SER148 (zero index'd)
name_lists = []
label_lists = []
molecule_types = []
resREF = 8
istart = 0
iend = 0
iname = "GLY"
molnum = 0


maxSize = len(list(trj.top.atoms))
stopFlag = False
tempMol = []
tempCGL = []
name_lists_key = []

cnt = 0
testNUM = resREF

for i, a in enumerate(trj.top.atoms) :
    resNAME = str(a.residue)[0:3]
    resNUM = int(str(a.residue)[3:7])
    aINDEX = a.index

    if resNAME not in name_lists_key :
        name_lists_key.append(resNAME)

####ADDED EVERY 3 RESIDUE    
    test = False
    if (testNUM == resNUM) :
        if (cnt % 3) == 0 :
            tmpstart = istart
    if (testNUM != resNUM) :
        if (cnt % 3) == 2 :
            test = True
        testNUM = resNUM 
        cnt += 1

#    if (resNUM != resREF) :
    if test :
        #first append name_lists and label
        iend = aINDEX - 1
        tempMol.append("index %d to %d" % (tmpstart, iend))
        tempCGL.append(iname)

        #then update things for next residue
        iname = resNAME
        istart = aINDEX
        if (resNUM < resREF) :
            #stopFlag = True
            molecule_types.append(int(molnum))
            name_lists.append(tempMol)
            label_lists.append(tempCGL)
            tempMol = []
            tempCGL = []
            molnum += 1
        resREF = resNUM

    # special case if last item
    if (i == (maxSize-1)) :
        iend = aINDEX
        tempMol.append("index %d to %d" % (istart, iend))
        tempCGL.append(iname)
        molecule_types.append(int(molnum))
        name_lists.append(tempMol)
        label_lists.append(tempCGL)

#name_lists=[['index 0 to 2', 'index 3 to 5', 'index 6 to 8', 'index 9 to 11', 'index 12 to 14', 'index 15 to 17', 'index 18 to 20', 'index 21 to 23', 'index 24 to 26', 'index 27 to 29', 'index 30 to 32', 'index 33 to 35', 'index 36 to 38', 'index 39 to 41', 'index 42 to 44', 'index 45 to 46', 'index 47 to 48', 'index 49 to 51', 'index 52 to 54', 'index 55 to 57', 'index 58 to 60', 'index 61 to 63', 'index 64 to 66', 'index 67 to 69', 'index 70 to 72', 'index 73 to 75', 'index 76 to 78', 'index 79 to 81', 'index 82 to 84', 'index 85 to 87', 'index 88 to 90', 'index 91 to 93', 'index 94 to 96', 'index 97 to 99', 'index 100 to 102', 'index 103 to 105', 'index 106 to 108', 'index 109 to 111', 'index 112 to 114', 'index 115 to 117', 'index 118 to 120', 'index 121 to 123', 'index 124 to 126', 'index 127 to 129', 'index 130 to 132', 'index 133 to 135', 'index 136 to 138']]

#label_lists=[['HHH', 'HHH', 'HHH', 'HHH', 'HHH', 'HHH', 'HHH', 'HHH', 'III', 'TTT', 'TTT', 'TTT', 'TTT', 'III', 'HHH', 'HHH', 'HHH', 'HHH', 'III', 'TTT', 'TTT', 'TTT', 'TTT', 'III', 'HHH', 'HHH', 'III', 'TTT', 'TTT', 'TTT', 'TTT', 'TTT', 'III', 'HHH', 'HHH', 'HHH', 'HHH', 'III', 'TTT', 'TTT', 'TTT', 'TTT', 'TTT', 'III', 'HHH', 'HHH', 'HHH']]

#label_lists=[['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'I', 'T', 'T', 'T', 'T', 'I', 'H', 'H', 'H', 'H', 'I', 'T', 'T', 'T', 'T', 'I', 'H', 'H', 'I', 'T', 'T', 'T', 'T', 'T', 'I', 'H', 'H', 'H', 'H', 'I', 'T', 'T', 'T', 'T', 'T', 'I', 'H', 'H', 'H']]

#molecule_types=[0]

#actual map command
print name_lists
print label_lists
print molecule_types

print "Lengths of all three lists should be equivalent: %d = %d = %d" % (len(name_lists), len(label_lists), len(molecule_types))

cg_trj = cg.map_unfiltered_molecules(            trj = trj,
                                      selection_list = name_lists, 
                                     bead_label_list = label_lists,
                                      molecule_types = molecule_types,
                                    mapping_function = "com")

cg_trj.save(output_dir + output_traj) 
cg_trj[0].save(output_dir + output_top)

############################### check results ###############################
# reloading results from disk.

cg_traj = cg_trj.load(output_dir + output_traj,top=output_dir + output_top) 
ref_cg_traj = cg_trj.load(reference_dir + reference_traj,
                          top=reference_dir + reference_top) 

result=check.md_content_equality(cg_traj,ref_cg_traj)

sys.exit(check.check_result_to_exitval(result))
