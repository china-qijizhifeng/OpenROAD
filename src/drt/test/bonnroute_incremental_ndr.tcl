# Incremental BonnRoute routes NDR nets with their wire type (route_ndr_nets 1).
source "helpers.tcl"
set test_name bonnroute_incremental_ndr

read_lef Nangate45/Nangate45_tech.lef
read_lef Nangate45/Nangate45_stdcell.lef
read_def ta_ap_aligned.def
read_guides ta_ap_aligned.route_guide

# This net spans several worker windows. The other net keeps initial routing
# so each ECO call takes the incremental BonnRoute path.
create_ndr -name WIDE_NDR -spacing {*2} -width {*2}
assign_ndr -ndr WIDE_NDR -net net_M2
set_bonnroute_dr_param profile 0
set_bonnroute_dr_param offtrack_access 0
set_bonnroute_dr_param postprocess 0
set_bonnroute_dr_param incremental 1
# Every incremental call checks that the routing space covers the design routing.
set_bonnroute_dr_param incremental_check 2
set_bonnroute_dr_param route_ndr_nets 1
detailed_route -bonnroute -verbose 0

foreach resident {1 0} {
  set_bonnroute_dr_param incremental_resident $resident
  set net [[ord::get_db_block] findNet net_M2]
  odb::dbWire_destroy [$net getWire]
  detailed_route -bonnroute -verbose 0
  if { [$net getWire] eq "NULL" } {
    error "NDR net was not routed (resident=$resident)"
  }
  if { [detailed_route_num_drvs] != 0 } {
    error "NDR net routing has DRC violations (resident=$resident)"
  }
}

puts "pass"
