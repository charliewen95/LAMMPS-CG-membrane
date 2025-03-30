# /************************************************************
# A simple script to analyze CG simulation of Piezo in vesicle
# Purpose of the script is to translate and rotate the system
# in order to place the vesicle at the origin and piezo at the
# "north pole". From averaged positions, it will be possible
# to extract the "dome" angle.
#
# NB: This approach is SLOW and STORAGE INEFFICIENT for large 
#     trajectories of large systems. Consider re-writing in a 
#     performance-language,  introducing a binning procedure
#     and/or density calculation.
#
#
#
# Jesper Madsen (jespermadsen@usf.edu)
# Sept., 2022
#
# /************************************************************

set out_sph_type3 [open "output_sphericalcoor_type3.dat" w]
set out_2dp_type3 [open "output_2dpolar_type3.dat" w]
set out_sph_pie [open "output_sphericalcoor_pie.dat" w]
set out_2dp_pie [open "output_2dpolar_pie.dat" w]
set out_sph_liptail [open "output_sphericalcoor_liptail.dat" w]
set out_2dp_liptail [open "output_2dpolar_liptail.dat" w]

# standard conversions, from wiki.tcl-lang.org
proc rect_to_spherical {x y z} {
	list [set rad [expr {hypot($x, hypot($y, $z))}]] [expr {atan2($y,$x)}] [expr {acos($z/($rad+1.0e-20))}]
}
proc spherical_to_rect {rad phi theta} {
	list [expr {$rad * cos($phi) * sin($theta)}] [expr {$rad * sin($phi) * sin($theta)}] [expr {$rad * cos($theta)}]
}

# we need this simplified version as well
proc spherical_to_2dpolar {rad theta} {
	list [expr {$rad*cos($theta)}] [expr {$rad*sin($theta)}]
}

# list-wise conversions, print to filehandles
proc lconvert_print { f1 f2 inlist } {
	foreach particle $inlist {
		set x [lindex $particle 0]
		set y [lindex $particle 1]
		set z [lindex $particle 2]
		set spherical [rect_to_spherical $x $y $z]
		puts $f1 "$spherical"
		set rad [lindex $spherical 0]
		set theta [lindex $spherical 2]
		set 2dpolar [spherical_to_2dpolar $rad $theta]
		puts $f2 "$2dpolar"
	}
	return 1
}

set num_steps [molinfo top get numframes]	 
#set num_steps 1 #test only first frame, for debugging

# ----- modify selections appropriately -----
set all [atomselect top "all"]
#set lipid_heads [atomselect top "name 3 and not resid 99999"]
set type3 [atomselect top "name 3"]
set piezo [atomselect top "resid 99999"]
set liptail [atomselect top "name 3 and not resid 99999"]
# -------------------------------------------

# loop over all frames in the trajectory
#for {set frame 0} {$frame < $num_steps} {incr frame} {
#for {set frame 0} {$frame < $num_steps} {set frame [expr {$frame+3}]} {
for {set frame 0} {$frame < $num_steps} {set frame [expr {$frame+1}]} {
	puts "#Frame: $frame"
	#  update selections
	$all frame $frame
	$type3 frame $frame
	$piezo frame $frame
	$liptail frame $frame
	$all update
	$type3 update
	$piezo update
	$liptail update
	
	# translate COM (of type3) to the origin, it's probably pretty close already
	set translate [measure center $type3]
	$all moveby [vecscale -1.0 $translate]

	# rotate so that piezo is at the (x) north pole {x,y,z} -> {x>0,0,0}
	set piezo_com [vecnorm [measure center $piezo]]
	$all move [transvecinv $piezo_com]

	# fetch and print spherical coordinates, notice we pass on {y z x} ON PURPOSE!
	# Recall that transvecinv(v) moves v to x. However, we'd like to change x -> z
	# so that we can use canonical notation in the coordinate conversions
	lconvert_print $out_sph_type3 $out_2dp_type3 [$type3 get {y z x}]
	lconvert_print $out_sph_pie $out_2dp_pie [$piezo get {y z x}]
	lconvert_print $out_sph_liptail $out_2dp_liptail [$liptail get {y z x}]
}

puts "# Done!"

quit
