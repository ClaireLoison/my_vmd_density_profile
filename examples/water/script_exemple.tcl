mol new water_pure1.psf
mol addfile water_pure.1.skip1000.dcd waitfor all

proc pptable5file {file l1 l2 l3 l4 l5 } { foreach i1 $l1 i2 $l2  i3 $l3 i4 $l4  i5 $l5  { puts $file " [format %10.6f $i1]\t[format %10.6f $i2] \t[format %10.6f $i3] \t[format %10.6f $i4] \t[format %10.6f $i5]" }  }

source ~/Calculations/TOOLS/vmd_density_profile/orientation_profile.tcl

set wdir [density_dir_profile -rho direction -selection "water and resid 11 to 20" -average 1 -frame_from 18 -frame_to 19 -axisDef "water"]

set fileout [open "density_orientation.dat" w+]
# - Show output as a table
puts $fileout  "#**************************************************"
puts $fileout  "#output of density_dir_profile"
puts $fileout  "#column 1 | Bin breaks coordinates: (z, Angstroms) "
puts $fileout  "#column 2 | average number of oxygen "
puts $fileout  "#column 3 | average projection of water dipole along X"
puts $fileout  "#column 4 | average projection of water dipole along Y"
puts $fileout  "#column 5 | average projection of water dipole along Z"
puts $fileout  "#**************************************************"

pptable5file $fileout [lindex $wdir 0] [lindex $wdir 1]  [lindex $wdir 2] [lindex $wdir 3] [lindex $wdir 4] 
close $fileout

