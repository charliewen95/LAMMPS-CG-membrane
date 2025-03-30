vmd="/Applications/VMD_1.9.4a55-x86_64-Rev11.app/Contents/vmd/vmd_MACOSXX86_64"
# translate and rotate

$vmd r350B1.data.lip.lammpstrj -e shape_profile_redo.tcl -dispdev text

# prep csv
for tag in type3 pie liptail; do
	awk '{print $3,$1}' output_sphericalcoor_${tag}.dat | sed 's/ /,/' > theta_rad_${tag}.csv
done
