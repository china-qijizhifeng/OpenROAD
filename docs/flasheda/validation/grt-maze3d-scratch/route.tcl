set_thread_count $::env(GRT_THREADS)
read_lef $::env(GRT_LEF)
read_def $::env(GRT_DEF)
set_routing_layers -signal metal2-metal6
set_global_routing_layer_adjustment metal2-metal6 0.5
if { [info exists ::env(GRT_DEBUG)] } { set_debug_level GRT timer 1 }
proc cpu_ticks { } {
  set f [open /proc/self/stat r]
  set s [read $f]
  close $f
  set fields [split [string range $s [expr { [string last ")" $s] + 2 }] end]]
  return [list [lindex $fields 11] [lindex $fields 12]]
}
if { [info exists ::env(PERF_CTL)] } {
  set ctl [open $::env(PERF_CTL) w]
  set ack [open $::env(PERF_ACK) r]
  puts $ctl "enable"
  flush $ctl
  gets $ack response
}
lassign [cpu_ticks] user0 sys0
set start [clock microseconds]
global_route -verbose -congestion_iterations 50
set elapsed [expr { ([clock microseconds] - $start) / 1000000.0 }]
lassign [cpu_ticks] user1 sys1
if { [info exists ::env(PERF_CTL)] } {
  puts $ctl "disable"
  flush $ctl
  gets $ack response
  close $ctl
  close $ack
}
set user_ticks [expr { $user1 - $user0 }]
set sys_ticks [expr { $sys1 - $sys0 }]
puts "GRT_MEASURE wall_s $elapsed user_ticks $user_ticks sys_ticks $sys_ticks\
  threads $::env(GRT_THREADS)"
write_guides $::env(GRT_OUT)/route.guide
write_db $::env(GRT_OUT)/routed.odb
puts "GRT_DONE 1"
