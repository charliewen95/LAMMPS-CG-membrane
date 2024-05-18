import sys, math
import numpy as np
import random

#
# Enforce "odd" N to allow zigzag stacking in case simple lattice bad initial guess.
#
# *     *    *    *    *    *    *    *    * 
#    *    *    *    *    *    *    *    *    *
# *     *    *    *    *    *    *    *    * 
#    *    *    *    *    *    *    *    *    *
#
# Separation between points on both rows and columns (via diagonal) should be the same.
# Therefore, points basically form equilateral triangles: offset for "indented" rows
# is therefore dx,dy = dr.[sqrt(2),sqrt(2)] where dr is the point separation.
#
# Note that passing dx = dy and x_adj = 0 gives simple lattice.
#

# Some modest modifications by Jesper Madsen, April 4th, 2017
# : Added feature to have lipid mixtures at variable ratio

# Some reducing of functionality by Jesper Madsen, April 4, 2022
# : ** DEBUG VERSION ONLY **
# : ** HARDCODED FOR 4-SITE **


def print_bilayer( N, atom_types, type_list, sigma, dx, dy, x_adj, delta, atom_index=1, lipid_index=1 ):
	atom_format = '%10d %10d %10d  %g %g %g'
	lipid_len = len(atom_types)


	# Upper monolayer
	#
	for j in range( 0, N ):
		y = (ylo + dy*j)

		for i in range( 0, N ):
			x = (xlo + dx*i)

			# Stagger overy other row
			if j%2 == 0: x += x_adj

                        if ( type_list[lipid_index-1] ):
        			for k in range( 0, lipid_len ):
        				z = +sigma * (float(lipid_len)-(0.5+k))
        				print atom_format % ( atom_index, lipid_index, atom_types[k], x,y,z )
        				atom_index += 1
        			lipid_index += 1
                        else:
        			for k in range( 0, lipid_len ):
        				z = +sigma * (float(lipid_len)-(0.5+k))
        				print atom_format % ( atom_index, lipid_index, atom_types[k], x,y,z )
        				atom_index += 1
        			lipid_index += 1
	#
	# Lower monolayer
	#
	for j in range( 0, N ):
		y = (ylo + dy*j)
		for i in range( 0, N ):
			x = (xlo + dx*i)

			# Stagger overy other row
			if j%2 == 0: x += x_adj

                        if ( type_list[lipid_index-1] ):
        			for k in range( 0, lipid_len ):
        				z = -sigma * (float(lipid_len)-(0.5+k))
        				print atom_format % ( atom_index, lipid_index, atom_types[k], x+x_adj,y,z )
        				atom_index += 1
        			lipid_index += 1
                        else:
        			for k in range( 0, lipid_len ):
        				z = -sigma * (float(lipid_len)-(0.5+k))
        				print atom_format % ( atom_index, lipid_index, atom_types[k], x+x_adj,y,z )
        				atom_index += 1
        			lipid_index += 1

	return atom_index, lipid_index

def print_bonds( N, atom_types, bonds, atom_index=1, bond_index=1 ):
	bond_format = '%10d  %10d  %10d %10d'
	lipid_len = len(atom_types)
	#
	# Upper monolayer
	#
	for i in range( 0, N ):
		for j in range( 0, N ):
			for bond in bonds:
				b_type, b_i, b_j = bond
				print bond_format % ( bond_index, b_type, atom_index+b_i, atom_index+b_j )
				bond_index += 1
			atom_index += lipid_len
	#
	# Lower monolayer
	#
	for i in range( 0, N ):
		for j in range( 0, N ):
			for bond in bonds:
				b_type, b_i, b_j = bond
				print bond_format % ( bond_index, b_type, atom_index+b_i, atom_index+b_j )
				bond_index += 1
			atom_index += lipid_len

	return atom_index, bond_index


def print_angles( N, atom_types, angles, type_list, atom_index=1, angle_index=1, lipid_index=1 ):
	angle_format = '%10d  %10d  %10d %10d %10d'
	lipid_len = len(atom_types)
	#
	# Upper monolayer
	#
	for i in range( 0, N ):
		for j in range( 0, N ):
                        if ( type_list[lipid_index-1] ):
        			for angle in angles:
	        			a_type, a_i, a_j, a_k = angle
				        print angle_format % ( angle_index, a_type, atom_index+a_i, atom_index+a_j, atom_index+a_k )
				        angle_index += 1
			        atom_index += lipid_len
                                lipid_index += 1

	#
	# Lower monolayer
	#
	for i in range( 0, N ):
		for j in range( 0, N ):
                        if ( type_list[lipid_index-1] ):
        			for angle in angles:
        				a_type, a_i, a_j, a_k = angle
        				print angle_format % ( angle_index, a_type, atom_index+a_i, atom_index+a_j, atom_index+a_k )
        				angle_index += 1
        			atom_index += lipid_len
                                lipid_index += 1



	return atom_index, angle_index

reference = 'Grime & Madsen (2019) arXiv:1910.05362'
sigma =  7.5           # default sigma = 0.75 nm, after B,P,&B
APL   = 70.0           # default APL = 0.7 nm^2, after max from B,P,&B

if len(sys.argv) < 4:
	print ''
	print 'Usage: %s N Lz APL [hex]' % ( sys.argv[0] )
	print ''
	print 'Where:'
	print ''
	print '  N   : make NxN lipid bilayer membrane on x,y plane'
	print '  Lz  : output simultion cell length on z axis'
	print '  APL : area per lipid in Angstrom^2'
	print '  hex : OPTIONAL switch to use hexagonal packing'
	print ''
	sys.exit( -1 )

N  = int( sys.argv[1] )
Lz = float( sys.argv[2] )
APL = float( sys.argv[3] )

n_lipids = N*N*2
#type_list = np.append( np.ones(n_i), np.zeros(n_ii) )
type_list=np.ones(n_lipids)
random.shuffle(type_list)
#print llist
#sys.exit(-1)

delta = math.sqrt(APL) # lipid x,y grid spacing in Angstrom via APL

# Standard lattice parameters
dx = delta
dy = delta
x_adj = 0.0

#
# Use hexagonal packing?
#
if( len(sys.argv) > 4 ):
	if (N%2 != 0): N += 1 # N should be even to allow staggering of rows
	theta = (30.0/180.0)*math.pi
	dy    = delta*math.cos(theta)
	x_adj = delta*math.sin(theta) # additional x stagger every second line!

Lx = dx*N
Ly = dy*N

total_atoms  = (2*N*N) * (4)
total_bonds  = (2*N*N) * (4-1) # total bonds in 2 monolayers with NxN 5-site lipids
total_angles = (2*N*N) * (4-2) # total angles in 2 monolayers with NxN 5-site lipids

xlo, xhi = -Lx/2, +Lx/2
ylo, yhi = -Ly/2, +Ly/2
zlo, zhi = -Lz/2, +Lz/2

atom_types  = [ 1 , 2 , 3 , 3 ]

# via Biophys J:
# tau = (mass*sig^2/eps)^(1/2) = ~3 ps
# so mass is ~47 g/mol with sig=7.5A and eps=kBT/0.85
#atom_mass = (3.0e-12)**2 * (300.0*1.3806e-23)/0.85 / (7.5e-10)**2 # mass in kg
#atom_mass = atom_mass * (1000.0*6.022e23) # mass in Daltons, or g/mol
atom_mass = 380

atom_masses = {}
for i in range(0,3): atom_masses[i+1] = atom_mass
#atom_masses[4] = 0.1
#atom_masses[5] = 0.1

bonds       = [ [1, 0,1],   [1, 1,2],   [1, 2,3] ] # [type, i,j] where i,j are ZERO-based atom indices in molecule
angles      = [ [1, 0,1,2], [1, 1,2,3],          ] # [type, i,j,k] where i,j,k are ZERO-based atom indices in molecule
lipid_len   = len( atom_types )

print '#LAMMPS data file, generated for lipid models in %s' % ( reference )
print '#Number of lipids: %d ' % ( n_lipids) 
print '%d atoms' % ( total_atoms )
print '3 atom types'
print '%d bonds' % ( total_bonds )
print '51 bond types'
print '%d angles' % ( total_angles )
print '1 angle types'
print ''
print '%g %g xlo xhi' % ( xlo, xhi )
print '%g %g ylo yhi' % ( ylo, yhi )
print '%g %g zlo zhi' % ( zlo, zhi )
print ''
print 'Masses'
print ''
for key in atom_masses:
	print '%10d  %g' % ( key, atom_masses[key] )
print ''
print 'Atoms'
print ''
print_bilayer( N, atom_types, type_list, sigma, dx, dy, x_adj, 1, 1 )
print ''
print 'Bonds'
print ''
print_bonds( N, atom_types, bonds, 1, 1 )
print ''
print 'Angles'
print ''
print_angles( N, atom_types, angles, type_list, 1, 1, 1 )
