# The jumpered routes of repair_antennas_jumper_own_track_cugr.tcl read back
# from the database: no violation, the demand of the jumper layers met3 and
# met4 as before, no overflow.
set test_name repair_antennas_jumper_own_track_cugr
source "helpers.tcl"
read_liberty "sky130hs/sky130hs_tt.lib"
read_db [make_result_file "${test_name}.odb"]

set_placement_padding -global -left 2 -right 2
set_global_routing_layer_adjustment met2-met3 0.8
set_global_routing_layer_adjustment met4-met5 0.9
set_routing_layers -signal met1-met5

repair_antennas -jumper_only
report_congestion
