# you can use this file by typing in a shell 
# vmd -dispdev text -e script_exemple.tcl
#######################
#=== preparations ! change the path to the path of the tcl package !

source ~/Calculations/TOOLS/my_vmd_density_profile/orientation_profile.tcl

proc pptable5file {file l1 l2 l3 l4 l5 } { foreach i1 $l1 i2 $l2  i3 $l3 i4 $l4  i5 $l5  { puts $file " [format %6.2f $i1]\t[format %6.2f $i2] \t[format %8.6f $i3] \t[format %8.6f $i4] \t[format %8.6f $i5]" }  }

#=== loading trajectory data from NAMD calculation

mol new water_pure1.psf
mol addfile water_pure.1.skip1000.dcd waitfor all

#=== analyse and output

set wdir [density_dir_profile -rho direction -selection "water and resid 11 to 20" -average 1 -frame_from 10 -frame_to 19 -axisDef "water"]

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


set wdir [density_dir_profile -rho direction -selection "water and resid 11 to 20" -average 1 -frame_from 10 -frame_to 19 -axisDef "generic" -refsel "name OH2" -com1sel "name OH2" -com2sel "name H1 H2" ]

set fileout [open "density_orientation_2.dat" w+]
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

exit
