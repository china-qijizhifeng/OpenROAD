# Incremental BonnRoute connects a new pin of a routed net and keeps its legal
# pin-to-pin paths (arXiv:2111.06169 Fig. 9; repair_antennas
# -keep_detailed_wires gives a routed net a diode this way).
source "helpers.tcl"
set test_name bonnroute_incremental_keep

read_lef Nangate45/Nangate45_tech.lef
read_lef Nangate45/Nangate45_stdcell.lef
read_def ta_ap_aligned.def
read_guides ta_ap_aligned.route_guide

set_bonnroute_dr_param profile 0
set_bonnroute_dr_param offtrack_access 0
set_bonnroute_dr_param postprocess 0
set_bonnroute_dr_param incremental 1
# Stops with DRT-9192 when a net with a new pin keeps it unconnected.
set_bonnroute_dr_param incremental_check 2
detailed_route -bonnroute -verbose 0

# A diode next to inv_M2_1 joins net_M2, whose routing stays.
set block [ord::get_db_block]
set net [$block findNet net_M2]
set master [[ord::get_db] findMaster ANTENNA_X1]
set diode [odb::dbInst_create $block $master keep_diode]
$diode setOrient MX
$diode setLocation 100700 19600
$diode setPlacementStatus PLACED
[$diode findITerm A] connect $net
detailed_route -bonnroute -verbose 0
if { [detailed_route_num_drvs] != 0 } {
  error "new pin routing has DRC violations"
}
# Two diodes on net_M3, one next to each inverter: each must be joined to its
# own gate (incremental_diode_node routes them one at a time; DRT-9193 stops
# the test when a diode lands on the other gate's node).
set net3 [$block findNet net_M3]
foreach {name x} {keep_diode_1 20900 keep_diode_2 180880} {
  set d [odb::dbInst_create $block $master $name]
  $d setOrient R0
  $d setLocation $x 100800
  $d setPlacementStatus PLACED
  [$d findITerm A] connect $net3
}
detailed_route -bonnroute -verbose 0
if { [detailed_route_num_drvs] != 0 } {
  error "multi-diode routing has DRC violations"
}
puts "pass"
