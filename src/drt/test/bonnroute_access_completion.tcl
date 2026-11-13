# BonnRoute DR gives every pin of an instance a primary access path. The
# design is a cut-out of the placement of ibex_sky130hd around five
# sky130_fd_sc_hd__a222oi_1 cells (their six input pins share one met1
# track, five of them off the met2 tracks). With the preference for the
# access point pin access chose, the catalogue of pin B1 holds only paths
# that conflict with the primary paths of its neighbours B2 and A1, so the
# conflict-free selection leaves B1 without a primary path; its net then
# fails in BonnRoute and FlexDR repairs these cells last (ibex_sky130hd:
# 87-99% of the markers after iteration 4 on 0.58% of the instances).
# offtrack_complete_selection repeats the selection of such an instance
# with other catalogues (DRT-9240) until every pin has a primary path.
# bonnroute_access_completion.{def,guide} come from
# lead/ibexswing/root/kit/hk/tools/gen_access_completion.py.
source "helpers.tcl"
read_lef "sky130hd/sky130hd.tlef"
read_lef "sky130hd/sky130hd_std_cell.lef"
read_def "bonnroute_access_completion.def"
read_guides "bonnroute_access_completion.guide"

set_bonnroute_profile 2013-improved
set_bonnroute_dr_param offtrack_access_report 20
tee -variable route_log { detailed_route -bonnroute -verbose 1 }

# DRT-9240: instances with a pin without a primary path before the other
# catalogues, completed ones, and those that keep such a pin.
set found [regexp {selection completion: (\d+) instances with a pin without a primary path, (\d+) completed, (\d+) keep} \
  $route_log -> instances completed kept]
set without [regexp -inline {(\d+) pins with catalogues and no primary path} $route_log]
if {
  $found && $instances >= 1 && $completed == $instances && $kept == 0
  && [lindex $without 1] == 0 && [detailed_route_num_drvs] == 0
} {
  puts pass
} else {
  puts "fail: completion found=$found without=$without drvs=[detailed_route_num_drvs]"
}
