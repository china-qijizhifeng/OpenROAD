# Bundled LEF/DEF-only global placement: no invented timing constraints.
set_thread_count $::env(FLASHEDA_THREADS)
sta::set_thread_count 1
set design $::env(FLASHEDA_DESIGN)
set dir $::env(FLASHEDA_OPENROAD_SOURCE)/src/gpl/test
read_lef $dir/nangate45.lef
if {$design eq "macro03"} {read_lef $dir/swerv_wrapper_macro.lef}
read_def $dir/$design.def
if {$::env(FLASHEDA_ROUTABILITY) eq "1"} {
 set_routing_layers -signal metal2-metal6
 set_global_routing_layer_adjustment metal2-metal6 0.5
}
set begin [clock microseconds]
global_placement {*}$::env(FLASHEDA_GPL_ARGS)
puts "FLASHEDA_STAGE placement_ms [expr {([clock microseconds]-$begin)/1000.0}]"
source $::env(FLASHEDA_TASK)/state.tcl
