# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

# Source into a private OpenROAD process with its original DB, libraries and
# constraints already loaded. The caller must verify the prepared patch hash.
# Any error requires discarding the entire process. This is not a live-DB
# transaction and does not roll back earlier mutations. No checkpoint is saved.
namespace eval route_trials {}

proc route_trials::import_joint { patch_def } {
  set block [ord::get_db_block]
  set before_nets [$block getNets]
  foreach net $before_nets {
    set before_iterms($net) [$net getITerms]
  }

  # INCREMENTAL skips NETS in the tested importer; FLOORPLAN updates wires.
  read_def -floorplan_initialize $patch_def
  if { [$block getNets] ne $before_nets } {
    error "joint import changed net identities"
  }

  # definNet::end reverses existing ITerm lists even for unchanged connections.
  # Verify membership before restoring order, preserving hierarchical modnets.
  set restored_orders 0
  foreach net $before_nets {
    set prior $before_iterms($net)
    set current [$net getITerms]
    if { $current eq $prior } {
      continue
    }
    if { [lsort $current] ne [lsort $prior] } {
      error "joint import changed logical connections"
    }
    foreach terminal $current {
      $terminal disconnectDbNet
    }
    foreach terminal [lreverse $prior] {
      $terminal connect $net
    }
    if { [$net getITerms] ne $prior } {
      error "failed to restore original ITerm ordering"
    }
    incr restored_orders
  }
  if { ![drt::experimental_route_trial_cmd ""] } {
    error "joint connectivity rejected"
  }
  return $restored_orders
}
