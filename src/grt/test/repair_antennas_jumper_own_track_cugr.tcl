# repair_antennas -jumper_only under CUGR on layers with one track per gcell
# edge (met4 and met5 with adjustment 0.9): a jumper takes the track with its
# wire; the pads of its own via stacks lie on that track and charge no edge
# the wire runs over. Every violation is repaired and the congestion report
# shows the jumper wires as the wire demand of met4. The routes read back
# from the database in a separate process keep the demand of the jumper
# layers met3 and met4, without overflow.
set test_name repair_antennas_jumper_own_track_cugr
source "helpers.tcl"
read_liberty "sky130hs/sky130hs_tt.lib"
read_lef "sky130hs/sky130hs.tlef"
read_lef "sky130hs/sky130hs_std_cell.lef"
read_def "gcd_sky130.def"

set_placement_padding -global -left 2 -right 2
set_global_routing_layer_adjustment met2-met3 0.8
set_global_routing_layer_adjustment met4-met5 0.9
set_routing_layers -signal met1-met5
global_route -use_cugr

check_antennas
repair_antennas -jumper_only
check_antennas
check_placement
report_congestion

set guide_file [make_result_file ${test_name}.guide]
write_guides $guide_file
diff_files ${test_name}.guideok $guide_file

write_db [make_result_file "${test_name}.odb"]
puts [exec [info nameofexecutable] -no_splash -no_init -exit \
  -threads [thread_count] \
  [file join [file dirname [info script]] "${test_name}_reload.tcl"] 2>@1]
