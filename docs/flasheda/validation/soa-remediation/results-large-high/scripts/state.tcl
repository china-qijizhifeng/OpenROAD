# Complete observable database state and geometric QoR, outside timed GPL.
set block [ord::get_db_block]
set dbu [$block getDbUnitsPerMicron]
set area 0
foreach inst [$block getInsts] {
 set master [$inst getMaster]
 set area [expr {$area+double([$master getWidth])*[$master getHeight]}]
 set box [$inst getBBox]
 puts "LARGE_CELL [$inst getName] [$box xMin] [$box yMin] [$inst getOrient] [$master getName] [$inst getPlacementStatus]"
}
foreach bterm [$block getBTerms] {
 foreach pin [$bterm getBPins] {
  foreach box [$pin getBoxes] {
   puts "LARGE_PORT [$bterm getName] [[$box getTechLayer] getName] [$box xMin] [$box yMin] [$box xMax] [$box yMax] [$pin getPlacementStatus]"
  }
 }
}
set hpwl 0
foreach net [$block getNets] {
 set terms {};set ports {}
 foreach iterm [$net getITerms] {lappend terms "[[$iterm getInst] getName]/[[$iterm getMTerm] getName]"}
 foreach bterm [$net getBTerms] {lappend ports [$bterm getName]}
 puts "LARGE_NET [$net getName] [$net getSigType] [lsort $terms] [lsort $ports]"
 if {[$net getSigType] in {POWER GROUND}} {continue}
 set xs {};set ys {}
 foreach iterm [$net getITerms] {
  lassign [$iterm getAvgXY] ok x y
  if {$ok} {lappend xs $x;lappend ys $y}
 }
 foreach bterm [$net getBTerms] {
  foreach bpin [$bterm getBPins] {
   foreach box [$bpin getBoxes] {
    lappend xs [expr {([$box xMin]+[$box xMax])/2}]
    lappend ys [expr {([$box yMin]+[$box yMax])/2}]
   }
  }
 }
 if {[llength $xs]>1} {
  set xs [lsort -integer $xs];set ys [lsort -integer $ys]
  set hpwl [expr {$hpwl+[lindex $xs end]-[lindex $xs 0]+[lindex $ys end]-[lindex $ys 0]}]
 }
}
puts "LARGE_AREA [expr {$area/($dbu*$dbu)}]"
puts "LARGE_HPWL [expr {double($hpwl)/$dbu}]"
puts "LARGE_DONE 1"
