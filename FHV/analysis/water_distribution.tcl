set bw 2
set pi 3.1415
set r2 144
set vol [expr $pi*$r2*$bw]
set ifd 70
set ffd 180
set nbin [expr round(($ffd-$ifd)/$bw)]

set cid ccccc
mol addfile ../part6_${cid}.xtc waitfor all step 3
#mol addfile ../alligned${cid}.xtc waitfor all last 1
set ifrm 1
set ffrm [molinfo top get numframes]

##Setting directions
#set out1 [open dir_vec.dat "w"]
#for {set p 1} {$p < $np} {incr p} {
#set lat [expr asin(-1.0 + (2.0 * $p)/$np)]
#set lon [expr $ga*$p]
#set x [expr cos($lon)*cos($lat)]
#set y [expr sin($lon)*cos($lat)]
#set z [expr sin($lat)]
#set dir($p) "$x $y $z"
#puts $out1 "$x $y $z"
#}
#close $out1

##reading directions
#set in [open "directions.dat" "r"]
#set idata [read $in]
#set sidata [split $idata "\n"]
#set np [llength $sidata]
#for {set nd 1} {$nd < $np} {incr nd} {
#  set jj [expr $nd-1]
#  set dir($nd) [lindex $sidata $jj]
#}
#close $in


set npp 500
#initializing array
for {set p 1} {$p < $npp} {incr p} {
 for {set i 1} {$i<$nbin} {incr i} {
  set AccRes($p,$i) 0
  set ProRes($p,$i) 0
 }
}


set reference [atomselect top "backbone" frame 0]
for {set nf $ifrm} {$nf < $ffrm} {incr nf} {
 animate goto $nf
 set this_frame [atomselect top "backbone" ]
 set trans_mat [measure fit $this_frame $reference]
 #moving the COM of system to 0,0,0
 set protcen [atomselect top "protein"]
 set selall [atomselect top "all"]
 $selall move $trans_mat
 $selall moveby [vecscale -1 [measure center $protcen]]
 $protcen delete

 source beta.tcl
 
 #loop over all directions
 for {set p 1} {$p < $np} {incr p} {
  puts "$nf $p"
  $selall move [transvecinv $dir($p)]
  
  #loop along the cylinder to calculate the Number of water
  for {set i 1} {$i<$nbin} {incr i} {
   set ll [expr $i*$bw+$ifd]
   set ul [expr $ll+$bw]
#  set nw [atomselect top "name W and (y)^2+(z)^2 < $r2 and x < $j and x > $i"]
   set nw [atomselect top "water and (y)^2+(z)^2 < $r2 and x < $ul and x > $ll"]
   set AccRes($p,$i) [expr $AccRes($p,$i) + [expr [$nw num]/$vol]]
   $nw delete
   set pw [atomselect top "protein and (y)^2+(z)^2 < $r2 and x < $ul and x > $ll"]
   set ProRes($p,$i) [expr $ProRes($p,$i) + [expr [$pw num]/$vol]]
   $pw delete
  }
  $selall move [transvec $dir($p)]
 }
 $selall delete
 $this_frame delete
}
$reference delete


for {set p 1} {$p < $np} {incr p} {
 set out [open $cid-cyl_wat_$p.dat "w"]
 set out2 [open $cid-cyl_pro_$p.dat "w"]
 for {set i 1} {$i<$nbin} {incr i} {
  puts $out "[expr 0.5*$bw+$i*$bw] [expr $AccRes($p,$i)/($ffrm-$ifrm)]"
  puts $out2 "[expr 0.5*$bw+$i*$bw] [expr $ProRes($p,$i)/($ffrm-$ifrm)]"
 }
 close $out
}

set out1 [open dir_vec.dat "w"]
for {set p 1} {$p < $np} {incr p} {
  puts $out1 "$i $dir($p)"
}
close $out1
