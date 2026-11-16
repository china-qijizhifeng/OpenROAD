# The local jumper flag reaches the repair entry point and defaults to off.
source "helpers.tcl"
read_lef Nangate45/Nangate45.lef
read_def gcd.def

rename grt::have_routes grt::have_routes_real
proc grt::have_routes { } { return 1 }
rename grt::repair_antennas grt::repair_antennas_real
proc grt::repair_antennas { args } {
  set ::repair_args $args
  return 0
}

repair_antennas -jumper_only
if { [lindex $::repair_args 5] != 0 } {
  error "Local jumpers must default to off"
}
repair_antennas -jumper_only -local_jumpers -iterations 2 -ratio_margin 5
if { [lrange $::repair_args 1 end] ne {2 5 1 0 1} } {
  error "Incorrect repair argument translation: $::repair_args"
}
puts "pass"
