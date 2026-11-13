# FastRoute's slack updates after the first one of a run estimate only the
# parasitics of the nets whose planar route changed. The slacks and the net
# order of every update and the guides equal those of the full estimate in
# every update (the golden files are written with
# grt::set_incremental_partial_parasitics 0).
source "helpers.tcl"
read_lef "sky130hs/sky130hs.tlef"
read_lef "sky130hs/sky130hs_std_cell.lef"
read_liberty "sky130hs/sky130hs_tt.lib"
read_def "critical_nets_percentage.def"
read_sdc "critical_nets_percentage.sdc"

source "sky130hs/sky130hs.rc"
set_wire_rc -signal -layer "met2"
set_wire_rc -clock -layer "met5"

set_propagated_clock [all_clocks]
estimate_parasitics -placement

set_routing_layers -signal met1-met5 -clock met3-met5
# Congested enough for three slack updates.
set_global_routing_layer_adjustment met1-met5 0.84
set_thread_count 4

set_debug_level GRT partial_parasitics 1
tee -quiet -variable log {
  global_route -critical_nets_percentage 30
}

set guide_file [make_result_file critical_nets_partial_parasitics.guide]
write_guides $guide_file

# Slack and net order hashes of every slack update, and the estimate kind of
# every update.
set hash_file [make_result_file critical_nets_partial_parasitics.hashes]
set stream [open $hash_file w]
set estimates {}
foreach line [split $log "\n"] {
  if { [regexp {(partial slacks|net order) update=(\d+) nets=(\d+) hash=([0-9a-f]+)} \
          $line -> kind update nets hash] } {
    puts $stream "$kind $update $nets $hash"
  } elseif { [regexp {partial parasitics update=\d+ .* estimate=(\w+)} \
                $line -> estimate] } {
    lappend estimates $estimate
  }
}
close $stream
puts "estimates: $estimates"

check "Three slack updates" {
  llength $estimates
} 3

check "The first update estimates every net, the others the changed nets" {
  expr { $estimates eq {full incremental incremental} }
} 1

check "Slacks and net order of every update match the full estimate" {
  diff_files critical_nets_partial_parasitics.hashok $hash_file
} 0

check "Guides match the full estimate" {
  diff_files critical_nets_partial_parasitics.guideok $guide_file
} 0

exit_summary
