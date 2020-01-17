# Core functions for computing density profiles. The Density Profile
# Tool is a VMD analysis plugin that computes 1-D projections of
# various atomic densities. The computation is performed in a single
# frame, a trajectory, or averaged over multiple frames.

# Toni Giorgino, ISIB-Consiglio Nazionale delle Ricerche, 2012
# http://multiscalelab.org/utilities/DensityProfileTool

# MODIFIED BY CLAIRE LOISON
# 16/01/2020

package provide density_orientation_profile 0.1
package require topotools
package require pbctools


# Declare the namespace for this dialog
namespace eval ::density_dir_profile:: {
    # Variables matching command line options
    variable dp_args
    variable dp_args_defaults {
	rho             number
	selection       all
	axis            z
	resolution      1
	Zsource         mass
	partial_charges 1
	frame_from      now
	frame_to        now
	frame_step      1
	average		0
	axisDir			2
	axisDef			water
    }
    array set dp_args $dp_args_defaults

    # List of args in "preferred" order
    variable dp_args_list {rho selection axis resolution Zsource partial_charges \
			       frame_from frame_to frame_step average axisDir axisDef}

}


# User-accessible proc
proc density_dir_profile { args } { return [eval ::density_dir_profile::density_dir_profile $args] }


# Help
proc ::density_dir_profile::density_dir_profile_usage { } {
    variable dp_args
    variable dp_args_list
    puts "VMD Density Profile tool.  Computes 1-D projections of various atomic densities."
    puts "The computation is performed in a single frame, a trajectory, or averaged"
    puts "over multiple frames. You probably want to wrap your trajectory."
    puts " "
    puts "Usage: density_profile <args>"
    puts "Args (with defaults):"
    foreach k $dp_args_list {
	puts "   -$k $dp_args($k)"
    }
    puts " "
    puts "Density  -rho  is one of  {number|mass|charge|electrons|direction}"
    puts "See documentation at http://multiscalelab.org/utilities/DensityProfileTool"
}


# Command line parsing (sets namespace variables). TODO: allow short
# substrings, e.g. -sel
proc ::density_dir_profile::parse_args {args} {
    variable dp_args
    foreach {a v} $args {
	if {![regexp {^-} $a]} {
	    error "Argument should start with -: $a"
	} 
	set a [string trimleft $a -]
	if {![info exists dp_args($a)]} {
	    error "Unknown argument: $a"
	} 
	set dp_args($a) $v
    }
}


# Main entry point. 
proc ::density_dir_profile::density_dir_profile {args} {
    variable dp_args
    variable dp_args_defaults
    array set dp_args $dp_args_defaults
    if {[llength $args]==0} {
	density_profile_usage
	return
    } 
    eval parse_args $args

    # parray dp_args

    # Compute the bare histogram
    set lhist [compute]

    # Reformat the histogram and return it
    set xbreaks [hist_to_xbreaks $lhist]
    set framelist [get_framelist]
    set values [hist_to_values $lhist]

    # If averaging or single-frame, flatten inner list
    if { $dp_args(average)==1 || [llength $framelist]==1 } {
	set values [average_sublists $values]
    }

    return [list $values $xbreaks]

}


# Convert histogram into lists of values.  If resolution=2 and
#  hist(2,-3)=0.23   ; hist(2,-1)=0.21;
#  hist(4,-2)=0.42   ; hist(4,0)=0.40 =>

# Convert histogram into list of breaks => {-6, -4, -2, 0}
proc ::density_dir_profile::hist_to_xbreaks {lhist} {
    variable dp_args
    array set hist $lhist
    lassign [get_x_range [array names hist]]  binmin binmax

    set xbreaks {}
    for {set bin $binmin} {$bin<=$binmax} {incr bin} {
	lappend xbreaks [expr $bin*$dp_args(resolution)]
    }
    return $xbreaks
}


# Values, return { {0.23 0} {0 0.42} {0.21 0} {0 0.40} }
proc ::density_dir_profile::hist_to_values {lhist} {
    variable dp_args
    array set hist $lhist
    lassign [get_x_range [array names hist]]  binmin binmax
    set framelist [get_framelist]

    # Outer cycle is on bins
    set v {}
    set nbins [expr $binmax-$binmin+1]

    for {set idx 0} {$idx<$nbins} {incr idx} {
	set bin [expr $idx+$binmin]
	set tmp {}
	# Inner cycle on frames
	foreach f $framelist {
	    lappend tmp $hist($f,$bin)
	}
	lappend v $tmp
    }
    return $v
}


# Average sublists: { {0.23 0} {0 0.42} {0.21 0} {0 0.40} } -> { 0.125
# 0.21 0.105 0.20 }. Also useful for flattening in case of one frame
proc ::density_dir_profile::average_sublists {vl} {
    set res {}
    foreach l $vl {
	lappend res [vecmean $l]
    }
    return $res
}



# Similar to average_sublists, but returns standard deviation of each bin.
proc ::density_dir_profile::stddev_sublists {vl} {
    set res {}
    foreach l $vl {
	lappend res [vecstddev $l]
    }
    return $res
}




# Compute and return the bare histogram. Relation between bins and
# coordinate is implicit; note that [0..$resolution) goes in bin 0,
# and so on. 
#
# This is the workhorse function, assumes namespace variables are
# properly set.
#
#  bin = [expr int(floor($x/$resolution))] 
#  xmin = $bin * $resolution
#
#  hist(frame,bin)
#
# The hist structure is contrieved because we don't know the bins
# until all of the frames have been scanned

proc ::density_dir_profile::compute { } {
    variable dp_args

    set resolution $dp_args(resolution)
    set axis $dp_args(axis)

    # Check if PBC box is set
    set area [transverse_area]
    set area_len [llength $area]
    if { $area_len == 1 } {
	if { $area == -1 } { 	
	    puts "Warning: No periodic cell information. Computing linear densities instead of volume densities."
	    set area [lrepeat [molinfo top get numframes] 1.0]
	} elseif { $area <= -2 } {
	    error "Cannot compute density, check arguments."
	}
    }
    
    # Values
    set rho [get_rho]
	puts "rho defined"
    
    # Build atomselection
    set as [atomselect top $dp_args(selection)]
    if { $dp_args(rho) == "direction" && $dp_args(axisDef) == "water" } {
			set newselection $dp_args(selection)
			set newselection [ append newselection " and name OH2" ]
		    set as [atomselect top $newselection]
    }
	puts "selection with [$as num] atoms"


    # Start loop over frames
    array unset hist
    set framelist [get_framelist]
    foreach f $framelist {
	puts "analyse frame $f"

	$as frame $f
	set xval [$as get $axis]
	puts $xval
	puts "end xval"
	
	if { $dp_args(rho) == "direction"  } {
		set dirval [ get_dir $as ]
		puts $dirval
		puts "end dirval"
	}
		
	# get area now and normalize density 
	# TODO profile
	set area_now [lindex $area $f]
	set rho_norm [vecscale [expr 1./$area_now/$resolution] $rho]

	# make histogram: hist(frame,bin)
	if { $dp_args(rho) == "direction"  } {
		puts "make histogram, [llength $xval ] [llength $rho_norm] [llength $dirval] "
	}
	
	#set atomnum 1
	if { $dp_args(rho) == "direction"  } {
	  foreach x $xval vn $rho_norm d $dirval {
	    # bin
	    #puts "analyse for atomnum $atomnum "
	    #set atomnum [ expr $atomnum + 1 ]
	    set bin [expr int(floor($x/$resolution))]
		#puts "frame $f : x = $x, vn = $vn , d = $d, bin = $bin"

	    if {! [info exists hist($f,$bin)] } { 
			set hist($f,$bin) 0.0
			set histd($f,$bin) 0.0
		}
		set hist($f,$bin) [expr $hist($f,$bin) + 1.0 ]
		set histd($f,$bin) [expr $histd($f,$bin) + $d]
		
	 }
    } else {
	  foreach x $xval vn $rho_norm  {
	    # bin
	    #puts "analyse for atomnum $atomnum "
	    #set atomnum [ expr $atomnum + 1 ]
	    set bin [expr int(floor($x/$resolution))]
		#puts "frame $f : x = $x, vn = $vn , d = $d, bin = $bin"

	    if {! [info exists hist($f,$bin)] } { set hist($f,$bin) 0.0 }
		set hist($f,$bin) [expr $hist($f,$bin) + $vn]
	 }
    
    }
    }
    $as delete
	
	puts "end histogram"
	
    # make bins for never-seen values
    fill_keys hist
	if { $dp_args(rho) == "direction"  } { 
	
		fill_keys histd 
		set framelist [get_framelist]
		
		lassign [get_x_range [array names hist]]  binmin binmax
		set nbins [expr $binmax-$binmin+1]
		puts "nbins = $nbins, binmin = $binmin, binmax = $binmax"
		
		foreach f $framelist {
			for {set idx 0} {$idx<$nbins} {incr idx} {
				set bin [expr $idx+$binmin]
				if { $hist($f,$bin) > 0.0 } {
					puts "renormalise $histd($f,$bin) by $hist($f,$bin) at frame $f and bin $bin"
					set hist($f,$bin) [expr $histd($f,$bin)/$hist($f,$bin)]
				}
			}
		}
		
		
		
	}
	
    # Return histogram 
	puts "Return histogram"
    return [array get hist]

}



# Sanity check on PBC. Return values
#   -1 if LINEAR densities can be still computed,
#   -2 if non-orthorhombic
#   -3 if wrong axis
#  else return a list of the transversal PBC area per frame
proc ::density_dir_profile::transverse_area { } {
    variable dp_args
    set axis $dp_args(axis)
    set pb_list [pbc get -all]
    set area_list {};		# to hold results

    foreach pb $pb_list {
	lassign $pb a b c alpha beta gamma
	# heuristic for unset box
	if {$a<2 || $b<2 || $c<2} {	
	    return -1;		# immediately
	} elseif {$alpha!= 90 || $beta!=90 || $gamma!=90} {
	    puts "Only orthorombic cells are supported"
	    return -2;		# immediately
	} else {
	    switch -- $axis {
		x { set area [expr $b*$c] }
		y { set area [expr $a*$c] }
		z { set area [expr $a*$b] }
		default {
		    puts "Wrong axis $axis, must be one of x,y or z"
		    set area -3
		}
	    }
	    lappend area_list $area
	}
    }
    return $area_list
}




# return the range over the 1st and 2nd dimension of a pseudo-2d array
# e.g. {2,3 5,4 2,4} -> {3 4}
proc ::density_dir_profile::get_x_range {kk} {    
    foreach k $kk {
	lappend xlist [lindex [split $k ,] 1]
    }
    set xlist [lsort -uniq -integer $xlist]
    set xmin [lindex $xlist 0]
    set xmax [lindex $xlist end]
    return [list $xmin $xmax]
}


# fill histogram keys so that there is one integer bin per each value
# between mi and max
# 
proc ::density_dir_profile::fill_keys arr {
    upvar $arr inp
    set keys [array names inp]
    lassign [get_x_range $keys]  xmin xmax
    set framelist [get_framelist]
    # puts "Filling frames $framelist, bins $xmin..$xmax"

    foreach f $framelist {
	for {set x $xmin} {$x<=$xmax} {incr x} {
	    if { ![info exists inp($f,$x)] } {
		set inp($f,$x) 0
	    }
	}
    }
}


# Auxiliary function, returns {from to step}, after possibly fixing
# "now". TODO first, last
proc ::density_dir_profile::get_framelist {} {
    variable dp_args
    set f $dp_args(frame_from)
    set t $dp_args(frame_to)
    set s $dp_args(frame_step)

    # Should be the atom selection frame?
    if { $f=="now" } { set f [molinfo top get frame] }
    if { $t=="now" } { set t [molinfo top get frame] }

    for {set i $f} {$i<=$t} {incr i $s} {
	lappend fl $i
    }
    return $fl
}


# Return the values of the selected property, as a list of one value
# per selected atom. 
#
#These will not change per-frame. => problem for direction ! 
# I shoud change that !
# 
proc ::density_dir_profile::get_rho {} {
    variable dp_args
    set as [atomselect top $dp_args(selection)]
	if { $dp_args(rho) == "direction" && $dp_args(axisDef) == "water" } {
			set newselection $dp_args(selection)
			set newselection [ append newselection " and name OH2" ]
		    set as [atomselect top $newselection]
    }
    
    if { [$as num]==0 } {
	$as delete
	error "Atom selection did not match any atom."
    }
    switch $dp_args(rho) {
	number { 
	    set tval [lrepeat [$as num] 1] 
	}
	mass { 
	    set tval [$as get mass] 
	}
	charge { 
	    set tval [$as get charge] 
	}
	electrons { 
	    set tval [getZ $as]
	    if {$dp_args(partial_charges)==1} {
		set pch [$as get charge]
		set tval [vecsub $tval $pch]
	    }
	}
	direction {
		#set tval [getDir $as $dp_args(axisDir) $dp_args(axisDef)]
		set tval [lrepeat [$as num] 1] 
	}
	default {
	    $as delete
	    error "Unknown rho, must be one of {number|mass|charge|electrons|direction}"
	}
    }
    $as delete
    return $tval
}


# Similar to [atomselect get ...] , but get Z number, based on the
# $Zsource option  . Requires an atomselection
proc ::density_dir_profile::getZ {as} {
    variable dp_args
    set attr $dp_args(Zsource)

    # If anything different than element was required, use
    # topotools. This will modify the "element" attribute, so make a
    # backup and restore when done.
    if { $attr == "element" } {
	set res [$as get atomicnumber]
    } else {
	set o_Z [$as get atomicnumber]
	topo guessatom element $attr
	set res [$as get atomicnumber]
	$as set atomicnumber $o_Z
    }

    # Z=0 or Z=-1 means unidentified
    if { [lsearch $res 0 ] != -1 || 
	 [lsearch $res -1] != -1} {
	error "Could not guess element for some atoms." 
    }

    return $res
}

# Similar to [atomselect get ...] , but get projection of  direction , based on the
# $Zsource option  . Requires an atomselection
proc ::density_dir_profile::get_dir {as} {
    variable dp_args
    #set dir $dp_args(axisDir)
    #set def $dp_args(axisDef)
    
    set def "water"
    set dir 2
    
    puts "def = $def"
    
#    set attr $dp_args(Zsource)

    if { $def != "water" } {
		error " the direction of water only is defined. \n Implement your molecule in orientation_profile.tcl to use it" 
    }
    
	set reslength [$as num]
	set residuesID [$as get resid]
	
	foreach myresid $residuesID {

			#puts "analyse residue $myresid"
			
			set myOH2 [atomselect top "name OH2 and resid $myresid"]
			set myH1 [atomselect top "name H1 and resid $myresid"]
			set myH2 [atomselect top "name H2 and resid $myresid"]
			set rest [atomselect top "not (name OH2 H1 H2 ) and resid $myresid"]
			
			#puts " atomnumbers [$myOH2 num]  [$myH1 num]  [$myH2 num] [$rest num]"
			
			if {([$myOH2 num] != 1) || ([$myH1 num] != 1) || ([$myH2 num] != 1) || ([$rest num] != 0)} {
				puts "number of OH2 in resid $myresid = [$myOH2 num] (should be 1)" 
				puts "number of H1 in resid $myresid = [$myH1 num] (should be 1)" 
				puts "number of H2 in resid $myresid = [$myH2 num] (should be 1)" 
				puts "number of rest atoms in resid $myresid = [$rest num] (should be 0)" 
				error "unexpected number of OH2, H1 and H2  and rest atoms. Check topology" 
			}
			#puts "test water ok"

			set coordO  [lindex [ $myOH2 get {x y z}] 0]
			#puts "coordO = $coordO"
			
			set coordH1 [lindex [ $myH1 get  {x y z}] 0]
			#puts "coordH1 = $coordH1"

			set coordH2  [lindex [ $myH2 get  {x y z}] 0]
			#puts "coordH2 = $coordH2"

			set dirvec  [ vecsub [vecadd $coordH1 $coordH2]  [vecadd $coordO $coordO] ]
			#puts "dirvec = $dirvec"

			set dirnorm [ vecnorm $dirvec ]
			#puts "dirnorm = $dirnorm"

			# attribute the direction on the oxygen only
			# the hydrogen have no weight here
			lappend res [ lindex $dirnorm  $dir]
			
			$myOH2 delete
			$myH1 delete
			$myH2 delete
			$rest delete
	}
	
    #puts "res = $res"
    return $res
}
puts "**************************************************"
puts "**************************************************"
puts "**************************************************"
puts "**************************************************"

proc pptable {l1 l2} { foreach i1 $l1 i2 $l2 { puts " [format %10.6f $i1]\t[format %10.6f $i2]" }  }

set wdens [density_dir_profile -selection "water and resid 11 to 20" ]
# - Show output as a table
puts "**************************************************"
puts "| Bin breaks coordinates: (z, Angstroms)"
puts "\t| Density * projection along z of water in each bin: (number * direction/A^3)"
pptable [lindex $wdens 1] [lindex $wdens 0]

set wdir [density_dir_profile -rho direction -selection "water and resid 11 to 20" ]
# - Show output as a table
puts "**************************************************"
puts "| Bin breaks coordinates: (z, Angstroms)"
puts "\t| Density * projection along z of water in each bin: (number * direction/A^3)"
pptable [lindex $wdir 1] [lindex $wdir 0]


