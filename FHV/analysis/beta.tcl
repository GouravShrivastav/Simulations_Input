set sel [atomselect top all]
 $sel set beta -1
 $sel set occupancy -1
$sel delete

#selecting a subunits and assigning the same occupancy to a trimer and same beta to 5 consecutive trimers as they are forming a pentamer
set a 344
set b 311
set c 309
set d 21
set e 19
set f 20

set iid 0
set fid -1
set count 0

for {set i 0} {$i<60} {incr i} {
set bv [expr int($i/5) ]
set iid [expr $fid + 1]
set fid [expr $iid + $a - 1]
set sel [atomselect top "residue $iid to $fid"]
$sel set occupancy $count
$sel set beta $bv
$sel delete
set iid [expr $fid + 1]
set fid [expr $iid + $b - 1]
set sel [atomselect top "residue $iid to $fid"]
$sel set occupancy $count
$sel set beta $bv
$sel delete
set iid [expr $fid + 1]
set fid [expr $iid + $c - 1]
set sel [atomselect top "residue $iid to $fid"]
$sel set occupancy $count
$sel set beta $bv
incr count
$sel delete
}

set count 0

for {set i 0} {$i<60} {incr i} {
set bv [expr int($i/5) ]
set iid [expr $fid + 1]
set fid [expr $iid + $d - 1]
set sel [atomselect top "residue $iid to $fid"]
$sel set occupancy $count
$sel set beta $bv
$sel delete
set iid [expr $fid + 1]
set fid [expr $iid + $e - 1]
set sel [atomselect top "residue $iid to $fid"]
$sel set occupancy $count
$sel set beta $bv
$sel delete
set iid [expr $fid + 1]
set fid [expr $iid + $f - 1]
set sel [atomselect top "residue $iid to $fid"]
$sel set occupancy $count
$sel set beta $bv
incr count
$sel delete
}

#--------------------------------------------------------------------
#Now, based on occupancy and betavalues calculate the center and direction vector
#loop over frames
set C5 "0 1 2 3 4 5 6 7 8 9 10 11"
set C3 "{0  1  4} {0  1  6} {0  4  8}  {0  6  9} {0  8  9} {1  4 10} {1  6 11} {1 10 11} {2  3  5} {2  3  7} {2  5  8} {2  7  9} {2  8  9} {3  5 10} {3  7 11} {3 10 11} {4  5  8} {4  5 10} {6  7  9} {6  7 11}"
set C2 "{0 1} {0 4} {0 6} {0 8} {0 9} {1   4} {1   6} {1  10} {1  11} {2   3} {2   5} {2   7} {2   8} {2   9} {3   5} {3   7} {3  10} {3  11} {4   5} {4   8} {4  10} {5   8} {5  10} {6   7} {6   9} {6  11} {7   9} {7  11} {8   9} {10 11}"


 set np 1
 foreach i $C5 {
  set sel [atomselect top "protein and beta $i"]
  set dir($np) [measure center $sel]
  $sel delete
  incr np
 }
 foreach i $C3 {
  set sel [atomselect top "protein and beta $i"]
  set dir($np) [measure center $sel]
  $sel delete
  incr np
 }
 foreach i $C2 {
  set sel [atomselect top "protein and beta $i"]
  set dir($np) [measure center $sel]
  $sel delete
  incr np
 }

# quasi C3
 for {set i 0} {$i <60} {incr i} {
  set sel [atomselect top "occupancy $i"]
  set dir($np) [measure center $sel]
  incr np
  $sel delete
 }

# quasi C2
 for {set i 0} {$i <12} {incr i} {
 set sel [atomselect top "beta $i"]
 set occ [lsort -unique [$sel get occupancy]]
 $sel delete
  for {set j 0} {$j < 4} {incr j} {
   set oj [lindex $occ $j]
   set sel [atomselect top "occupancy $oj"]
   set jcm [measure center $sel]
   $sel delete
   for {set k [expr $j+1]} {$k < 5} {incr k} {
   set ok [lindex $occ $k]
    set sel [atomselect top "occupancy $ok"]
    set kcm [measure center $sel]
    $sel delete
    set v [vecsub $jcm $kcm]
    set vl [veclength $v]
    if {$vl > 48 && $vl < 58} {
     puts "$j $k $vl"
     set sel [atomselect top "occupancy $oj $ok"]
     set dir($np) [measure center $sel]
     $sel delete
     incr np
    }
   }
  }
 }

   
#center of alpha
 for {set i 0} {$i <60} {incr i} {
  set iid [expr $i*964]
  set fid [expr $iid+343]
  set sel [atomselect top "occupancy $i and residue $iid to $fid"]
  set dir($np) [measure center $sel]
  incr np
  $sel delete
 }
 puts "Total Points : $np"

