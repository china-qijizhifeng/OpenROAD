# A detailed_route call on a design whose pin access ran before reuses the
# access points of unchanged unique classes and the patterns of unchanged
# rows.  DRT_PA_INCREMENTAL_CHECK (set by the test) compares every such pin
# access with one from scratch and stops on a difference.
source "helpers.tcl"
set test_name pa_incremental

read_lef Nangate45/Nangate45_tech.lef
read_lef Nangate45/Nangate45_stdcell.lef
read_def ta_ap_aligned.def
read_guides ta_ap_aligned.route_guide

# Routes and returns {reused classes, classes, reused rows, rows} from the
# check's report.
proc route_with_check { label } {
  utl::redirectStringBegin
  detailed_route -verbose 0
  set log [utl::redirectStringEnd]
  set pattern {reused (\d+) of (\d+) unique classes and (\d+) of (\d+) rows;}
  append pattern { (\d+) of \d+ result lines differ}
  if { ![regexp $pattern $log -> rc nc rr nr diff] } {
    error "$label: the pin access was not checked"
  }
  puts "$label: classes $rc/$nc rows $rr/$nr differences $diff"
  return [list $rc $nc $rr $nr]
}

proc expect { label condition } {
  if { ![uplevel 1 [list expr $condition]] } {
    error "$label: expected $condition"
  }
}

detailed_route -verbose 0

# Nothing changed: every class and row is reused.
lassign [route_with_check unchanged] rc nc rr nr
expect unchanged {$rc == $nc && $rr == $nr && $nr > 0}

# A diode abutting inv_M2_1 joins net_M2: its class is new and the row of
# inv_M2_1 changes; the other rows are reused.
set block [ord::get_db_block]
set diode [odb::dbInst_create $block [[ord::get_db] findMaster ANTENNA_X1] \
  incremental_diode]
$diode setOrient MX
$diode setLocation 100700 19600
$diode setPlacementStatus PLACED
[$diode findITerm A] connect [$block findNet net_M2]
lassign [route_with_check diode] rc nc rr nr
expect diode {$rc > 0 && $rc < $nc && $rr > 0 && $rr < $nr}

# inv_M3_2 moves by one site and keeps its class: only its row changes.
[$block findInst inv_M3_2] setLocation 180500 100800
lassign [route_with_check moved] rc nc rr nr
expect moved {$rc == $nc && $rr == $nr - 1}

# inv_M2_2, the first instance of its class, moves: the class and the rows
# of its instances are computed again.
[$block findInst inv_M2_2] setLocation 100320 179200
lassign [route_with_check candidate] rc nc rr nr
expect candidate {$rc < $nc && $rr < $nr}

puts "pass"
