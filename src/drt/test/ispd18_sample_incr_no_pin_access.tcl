# Incremental detailed routing without pin access only re-runs the first
# iteration in the regions changed since the previous run.
source "helpers.tcl"

read_lef testcase/ispd18_sample/ispd18_sample.input.lef
read_def testcase/ispd18_sample/ispd18_sample.input.def
read_guides testcase/ispd18_sample/ispd18_sample.input.guide
set_routing_layers -signal Metal1-Metal9
detailed_route -verbose 0
set net [[ord::get_db_block] findNet "net1231"]
odb::dbWire_destroy [$net getWire]
set_debug_level DRT incremental 1
detailed_route -no_pin_access -verbose 0
set def_file [make_result_file ispd18_sample_incr_no_pin_access.def]
write_def $def_file
diff_files ispd18_sample.defok $def_file
