// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute detailed routing on the drt database (Gester et al., TODAES
// 18(2), 2013, Sec. 3-4). Runs after pin access, guide processing and track
// assignment; writes complete net routings as initial routing, which
// FlexDR then cleans up (the paper's BR + external DRC clean-up flow).

#pragma once

#include <array>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "br/BrAccessFootprint.h"
#include "br/BrCore.h"
#include "br/BrDistanceRules.h"
#include "br/BrFastGrid.h"
#include "br/BrGeometry.h"
#include "br/BrPathCleanup.h"
#include "br/BrPinAccess.h"
#include "br/BrPotentialOracle.h"
#include "br/BrSameNetRules.h"
#include "br/BrShapeGrid.h"
#include "br/BrBusyThreads.h"
#include "br/BrOrderedBoxSchedule.h"
#include "br/BrTailPipeline.h"
#include "br/BrWireModel.h"
#include "frBaseTypes.h"
#include "frRTree.h"
#include "odb/geom.h"

namespace odb {
class dbBlock;
}

namespace utl {
class Logger;
}

namespace drt {

class frBlock;
class frBlockObject;
class frDesign;
class frNet;
class frInst;
class frInstTerm;
class frBTerm;
class frLayer;
class frViaDef;
struct RouterConfiguration;

struct BonnRouteDRParams
{
  // Via cost in pitches of the lower layer (gamma of Sec. 4.1).
  int via_cost = 4;
  // Jog cost factor (beta of Sec. 4.1).
  int jog_factor = 4;
  // Vertex search: jogs cost jog_factor_low (0: off) in GCells whose global
  // routing usage of the layer is at most jog_low_utilization percent of
  // its capacity; jog_factor elsewhere.
  int jog_factor_low = 0;
  int jog_low_utilization = 50;
  // Rip-up penalty per vertex, in pitches of the bottom layer.
  int ripup_cost = 40;
  int max_ripups = 3;
  // Sec. 3.3 / 4.4: rip-up level of the routed shapes of critical nets
  // (clock and NDR nets, routed first; critical_nets); other nets have
  // level 1. A rip-up sequence of the first pass removes shapes up to
  // first_pass_ripup_level, later passes up to 6 ("higher ripup effort").
  // Defaults 1 and 1: every routed shape is removable in every pass.
  int critical_ripup_level = 1;
  int first_pass_ripup_level = 1;
  // Sec. 4.4 / 5.1: a net whose connection fails inside its partition
  // region gets a rip-up sequence there (victims resident in the region, at
  // most first_pass_ripup_level); the sequences of a region settle at most
  // this percentage of the vertices its connections settled. 0: none (the
  // serial rip-up phase handles every failed net).
  int region_ripup_effort = 0;
  // D1 S2 (design/D1-connection-unit.md): a net whose connection fails in
  // a partition region keeps the connections it closed (committed, its
  // reservations stay active); later stages route only the terms not joined
  // yet. Nets still partial after the rip-up phase are removed and left to
  // FlexDR whole. 0: a failed net keeps nothing.
  bool keep_partial = true;
  // Parallel rip-up (Sec. 4.4 / 5.1): before each pass of the serial rip-up
  // queue, this many rounds of partition regions run the queued nets'
  // rip-up sequences in parallel, victims resident in the region; each
  // region gets the serial budget times its share of the queued nets. The
  // serial loop handles the rest. 0: serial only.
  int ripup_partitions = 0;
  // Rip-up sequences scheduled by net boxes (BrOrderedBoxSchedule.h): the
  // rip-up queue runs in chunks of ripup_box_chunk nets; a net's sequence
  // stays inside its routing area widened by ripup_box_margin tiles
  // (victims resident there) and starts when the earlier sequences of the
  // chunk whose boxes meet its box returned, so results equal the serial
  // order of the chunk at every thread count. A sequence that failed at
  // the limit of its box runs again with 4 x the margin, then over the
  // whole chip. Takes precedence over ripup_partitions. 0: serial queue.
  // Under an effort budget (final_ripup_effort) the first chunk of a pass
  // holds at most 16 sequences and a later chunk at most the remaining
  // budget divided by the mean effort of the pass's sequences so far, so a
  // chunk ends near the budget.
  bool ripup_boxes = false;
  int ripup_box_margin = 4;
  int ripup_box_chunk = 64;
  // ripup_box_chunk_retry: a sequence that failed at the limit of its box
  // with the margin runs again with 4 x the margin in the same chunk, after
  // the chunk's first sequences (in their order), scheduled by the larger
  // box; 0: in the next chunk. ripup_box_order: a chunk's sequences run in
  // layers of pairwise disjoint boxes (OrderedBoxSchedule::
  // disjointLayerOrder, layers formed in queue order), so a chain of waiting
  // sequences holds at most one sequence per layer; 0: queue order. Both
  // change the serial order the results equal, not their determinism.
  bool ripup_box_chunk_retry = false;
  bool ripup_box_order = false;
  // Off-track pin access: largest jog to a track, in pitches.
  int max_access_jog = 2;
  int max_access_options = 8;
  // Extra cost of a non-preferred access point, in pitches.
  int pin_cost = 2;
  // Routing area: guide tiles plus this many adjacent layers / tiles.
  int layer_dilation = 1;
  int retry_dilation = 2;
  // Width of the corridor: the guide area (path regions and the net's whole
  // area) includes this many tiles around the guide tiles. 0 keeps the
  // global routing corridor; the dilated area of a retry adds
  // retry_dilation on top.
  int guide_dilation = 0;
  // Vertex search: wires and jogs in the tiles around the guide tiles (the
  // guide_dilation and retry_dilation parts of the routing area) cost this
  // percentage more than in the guide tiles; 0: the same. With
  // use_intervals it is set to 0 (DRT-9131).
  int corridor_cost = 0;
  bool use_intervals = false;
  bool ripup = false;
  // Sec. 4.4: search area of a connection = GR tiles of a shortest path
  // along the guides from the source component to another component.
  bool component_regions = true;
  // Sec. 5.1 partition sequence; 0 selects batches of nets with disjoint
  // bounding boxes. partition_regions = 0 uses twice the thread count.
  bool partitions = true;
  int partition_regions = 0;
  // The speculative tail (tail_batch) takes over at the first partition
  // with fewer regions than this; 0 uses the thread count. The routing
  // depends on the thread count only through these two defaults: with both
  // parameters set, every thread count gives the same result.
  int tail_regions = 0;
  // Late partitions with fewer regions than threads are replaced by batches
  // of this many pending nets, searched in parallel on the unchanged
  // routing space and committed in net order; a net whose shapes conflict
  // with a commit earlier in its batch is searched again in the next batch.
  // Results depend on the net order and the batch size only. 0 keeps the
  // partition sequence to its end.
  int tail_batch = 32;
  // D6 (BrTailPipeline.h, TailBatchPolicy): batches of at least
  // tail_batch_min and at most tail_batch nets; the size halves after a
  // batch with more than tail_batch_high percent conflicts per attempt and
  // doubles below tail_batch_low percent. 0: every batch has tail_batch
  // nets. Barrier batches and the lookahead follow the same sizes.
  int tail_batch_min = 0;
  int tail_batch_high = 50;
  int tail_batch_low = 20;
  // D6 growth: above tail_batch, batches start at tail_batch nets and grow
  // up to tail_batch_max under the same rule (doubling below
  // tail_batch_low percent conflicts, halving above tail_batch_high, never
  // below tail_batch_min or, with tail_batch_min 0, tail_batch). The
  // deferred pass starts at tail_deferred_batch. 0, or at most the start
  // size: no growth.
  int tail_batch_max = 0;
  // Threads of the speculative batches; 0 uses all router threads. The
  // batch results do not depend on it (a check of that determinism).
  int tail_threads = 0;
  // Within a batch, start the searches with the largest estimated cost
  // first (longest processing time first); results do not depend on it.
  bool tail_order = true;
  // Pipeline instead of a barrier per batch (BrTailPipeline.h; off by
  // default: 27-49% shorter tails, but TNS 1-3% outside the measured noise
  // on sky130 designs): batches of
  // tail_batch nets keep their semantics (searched on the commits of the
  // earlier batches, decided in order), a conflicting net is searched again
  // tail_retry_delay batches later, and searches and decisions of nets with
  // disjoint tile boxes overlap across batches. Results depend on the net
  // order, tail_batch and tail_retry_delay only.
  bool tail_pipeline = false;
  int tail_retry_delay = 1;
  // Tail lookahead (BrTailPipeline.h, TailLookahead): the batches of the
  // barrier version (tail_batch nets, conflicting nets first in the next
  // batch), decided as soon as the searches and decisions of interacting
  // nets of the batch are done; idle threads search nets that certainly
  // belong to the next batch and whose boxes no undecided net of the
  // current batch can still write. Results equal the barrier version for
  // every thread count. Takes precedence over tail_pipeline.
  bool tail_lookahead = true;
  // Settled vertices a net may use in its speculative tail attempt; a net
  // that needs more (almost always one without a route, left to FlexDR) is
  // deferred to a final pass of the stage without a budget, so it does not
  // hold up the batches. 0 disables the budget.
  int tail_defer_settled = 300000;
  // Staged budgets: the first pass of a tail stage gives every net this
  // smaller budget, so a batch (whose wall time is its longest search)
  // holds short searches only; the nets over it are searched again in a
  // second pass with tail_defer_settled. 0, or a value at or above
  // tail_defer_settled: one pass with tail_defer_settled.
  int tail_first_settled = 100000;
  // The final pass of the stage for the deferred nets. 0 hands them to the
  // rip-up phase instead (final_ripup), whose rip-up mode search finds a
  // legal path at no penalty when one exists and does not exhaust the area
  // when none exists.
  bool tail_deferred_pass = true;
  // Settled vertices a net may use in the deferred pass; nets over it go to
  // the rip-up phase (most deferred nets have no route: jpeg 37 of 39).
  // 0: no budget.
  int tail_deferred_settled = 0;
  // Batch size of the deferred pass. Its searches are long and few of them
  // end with a route, so larger batches cost few conflicts and fewer
  // barriers. 0: tail_batch.
  int tail_deferred_batch = 128;
  // Sec. 5.1: a critical net routing step before the partitions; clock
  // nets are routed first, serially, on the whole chip.
  bool critical_nets = true;
  // Sec. 3.2 / 5.1: route nets with nondefault rules with their wire type
  // (NDR width, vias and spacing) in the critical-net step; 0 leaves them
  // to FlexDR as before.
  bool route_ndr_nets = true;
  // Sec. 3.6: wire types held by the fast grid ("a restricted set of
  // frequent wire types"): the default rule and the NDRs used by the most
  // nets, at most 5 (FastGrid::kMaxWireTypes), each by at least
  // fast_ndr_permille / 1000 of the nets; the others are answered by the
  // distance rule checker. Every fast wire type adds its fields to every
  // fast-grid build and update over the whole chip.
  int fast_wire_types = 5;
  int fast_ndr_permille = 10;
  // Interval property (ii) also for vias (Sec. 4.1); 0 relaxes vias per
  // vertex and lets intervals ignore via availability.
  bool interval_via_property = true;
  // Among interval events with equal key, process the one with the
  // smallest future cost first; Algorithm 4 (line 12) takes the largest.
  // Both are exact; this order reaches targets without labeling the whole
  // last flat layer of equal-length shortest paths.
  bool interval_nearest_first = false;
  // Sec. 4.4: bound on the depth of a rip-up sequence; a failed sequence
  // restores the original routing.
  int ripup_depth = 3;
  // Later passes over unrouted nets use more rip-up effort and a larger
  // routing area (Sec. 4.4, last sentence).
  int ripup_passes = 2;
  // Ahrens et al., arXiv:2111.06169, Sec. 5.2: legal track-assignment wires
  // become reservations of their net, usable by it at cost factor
  // discount_numerator / discount_denominator and blocked for other nets.
  bool reservations = false;
  int discount_numerator = 3;
  int discount_denominator = 4;
  // Lowest track-graph layer (0 = bottom routing layer) with reservations.
  int reserve_min_layer = 0;
  // Future cost: 0 = pi_H of the 2013 paper (Sec. 4.1); arXiv:2111.06169
  // potentials 1 = l1 (Table 3), 2 = simple model (Theorem 3), 3 = general
  // model (Section 4, guide area and discounts), 4 = Gester et al. 2013,
  // Sec. 4.1: the general (blockage-aware) model only for connections
  // whose guide path detours by at least potential_detour_tiles, pi_H
  // otherwise. Reservations or a potential select the vertex search.
  int potential = 0;
  int potential_detour_tiles = 2;
  // Debug: repeat every search that used a potential with pi_H and count
  // cost differences (a feasible potential has none).
  bool potential_check = false;
  // pi_H with the exact via-stack bound lb_via of Sec. 4.1 and the
  // smallest terminal cost per target box; 0 keeps dz * (cheapest via).
  bool exact_future_cost = true;
  // Sec. 4.2 wire spreading: extra cost (pitches of the layer) of a wire
  // vertex beside an occupied track, in GCells whose global routing
  // utilization (dbGCellGrid usage / capacity) is at most
  // spread_utilization percent. Selects the vertex search.
  bool spreading = true;
  // After routing, replace "via up, a wire on the layer above, via down"
  // (nothing else of the net on that wire) by a wrong-way wire on the lower
  // layer when it is legal against other nets: two vias fewer.
  bool collapse_bridges = true;
  // Also the reverse: a wire on the lower layer between two vias up (nothing
  // else of the net on it) becomes a wrong-way wire on the upper layer.
  bool collapse_lower_bridges = true;
  // Also spread on layers whose wires must stay on tracks in the preferred
  // direction (unidirectional, LEF58 right-way-on-grid-only).
  bool spread_unidirectional = false;
  int spread_cost = 1;
  int spread_utilization = 70;
  // Sec. 2.1 / 4.2: spread where global routing assigned the net extra
  // space (net property of `global_route -bonnroute` with extra_space 1);
  // without any such property the GCell congestion tiles are used.
  bool spread_from_extra_space = true;
  // Off-track pin access (Sec. 3.7, 3.8, 4.3, 4.4): per pin a catalogue of
  // shortest tau-feasible paths on the blockage grid to on-track points
  // within `offtrack_radius` pitches, primary paths per circuit class by
  // destructive bounding. 0 selects the previous access stubs from
  // OpenROAD access points and greedy reservations.
  bool offtrack_access = true;
  int offtrack_radius = 2;
  // Regenerate a pin's catalogue when fewer of its paths remain feasible
  // (Sec. 4.4 "If no or only a small number of access paths remains
  // feasible, we dynamically generate new access paths").
  int offtrack_regenerate = 2;
  // Minimum segment length tau in DBU; 0 derives it from the same-net
  // rules (Sec. 3.7). `offtrack_notch_tau` adds notch spacing + width.
  int offtrack_tau = 0;
  bool offtrack_notch_tau = false;
  // Blockage grid of a catalogue search: the Hanan lines alone, and one
  // search state per grid vertex, where tau is at most the manufacturing
  // grid (every coordinate is on that grid, so any two lines are at least
  // tau apart and every Hanan-grid path is tau-feasible). 0 keeps
  // Algorithm 3's lines at multiples of tau and the five direction states.
  bool offtrack_hanan_grid = true;
  // Among catalogue paths of equal cost, prefer the least preferred-
  // direction wire off the tracks of the search layers (lexicographic
  // tie-break; path costs unchanged): paths jog onto a track at the access
  // point instead of running beside the tracks. 0: ties resolve by the
  // search order.
  bool offtrack_on_track_ties = true;
  // D10: catalogues per pin class instead of per circuit class. The key of
  // a pin holds only its own search inputs that are equal between
  // instances of a master (local frame): master, pin, access points with
  // their vias and preferred marks, the instance's own shapes and the
  // tracks within the pin's search reach. Each class's catalogue is
  // generated at its first pin; every other pin maps it, keeps the paths
  // legal at the member, and generates its own catalogue when fewer than
  // min(offtrack_regenerate, catalogue size) remain. Selection per
  // instance. Needs access_classes.
  bool offtrack_pin_classes = false;
  // M3 statistics (verbose log, results unchanged): circuit classes with
  // one signature component left out, pin-key classes, and the blockage
  // grid lines of both grids.
  bool offtrack_access_stats = false;
  // D10 sampled comparison: every n-th member pin also generates its own
  // catalogue, compared with the mapped one (0: off; results unchanged).
  int offtrack_catalogue_sample = 0;
  // Access report (verbose log, results unchanged): wire sharing of the
  // catalogues and off-track wire and blocked tracks of the primary paths
  // after setup; pins without a primary path, pins of failed connects
  // (no feasible access option) and nets left to FlexDR, each listed up to
  // this many lines (0: off).
  int offtrack_access_report = 0;
  // With offtrack_hanan_grid: an instance whose selection leaves a pin with
  // a catalogue unassigned recomputes its catalogues on the Algorithm 3
  // grid (five states, no on-track ties) and solves again; the result that
  // assigns more pins is kept (ties: the Hanan one). 0: Hanan catalogues
  // only.
  bool offtrack_algorithm3_fallback = true;
  // An instance whose selection leaves a pin with a catalogue without a
  // primary path repeats the selection with other catalogues of its pins
  // (br::selectionCompletionVariants): every access point of a pin at equal
  // cost (the preference for the access point pin access chose fills a
  // catalogue with that point's paths), on the configured grid and on the
  // Algorithm 3 grid, then twice more with one more pitch of radius and
  // twice the paths. The first selection that gives every pin a primary
  // path is kept (DRT-9240). 2: the catalogues with all access points
  // alike at the catalogue size come before the Algorithm 3 fallback; 1:
  // all of them after it; 0: off. An instance whose selection is complete
  // keeps its catalogues and primary paths.
  int offtrack_complete_selection = 2;
  // Reuse a circuit class's catalogues and conflict-free solution for all
  // its instances; 0 solves every instance separately.
  bool access_classes = true;
  // Circuit signature area: the instance box expanded by this many pitches
  // of the lowest track layers (Sec. 4.3 "expanding the bounding area").
  // Paths mapped from the class representative are checked at each member.
  int access_class_margin = 0;
  // Off-track access also for pins of block (macro) instances; 0 keeps the
  // access stubs there (see OffTrack.md).
  bool offtrack_macro_pins = false;
  int access_node_limit = 20000;
  // Weights of the conflict-free solution (Sec. 4.3), in pitches of the
  // bottom track layer: per blocked track, per missing on-track
  // continuation direction; spreading per DBU two endpoints are closer than
  // two pitches.
  int access_track_weight = 1;
  int access_direction_weight = 1;
  int access_spread_weight = 1;
  // LEF58 end-of-line spacing (ENDTOEND), EOLKEEPOUT, CORNERSPACING and
  // cut class spacing tables in the routing space (Lef58Rules.md); 0 keeps
  // drt's EOL summary as the line-end rule and the plain cut spacing.
  bool lef58_rules = true;
  // Same-net cut spacing and via pad rules between the vias of a net: the
  // vertex search refuses a via whose cut or far pad violates them against
  // the via that entered the layer, or whose cut violates them against the
  // net's earlier cuts (Lef58Rules.md).
  bool same_net_vias = true;
  // Same-net post-processing of each connection (Sec. 4.4).
  bool postprocess = true;
  // Fill the gap between two separated pieces of one net's metal (wires,
  // via pads, patches, pins) that are closer than drt's spacing, with a
  // patch wire, when the fill is legal against other nets (same-net
  // spacing of via pads and pins, which the path post-processing does not
  // see).
  bool notch_fill = true;
  // Rip-up sequences (Sec. 4.4) for the nets the parallel phases leave,
  // also when `ripup` is off (2013-improved).
  bool final_ripup = true;
  // Effort budget of that phase in percent of the vertices settled by the
  // parallel phases (0: unlimited).
  int final_ripup_effort = 10;
  // The base of that budget leaves out the vertices settled by speculative
  // tail attempts decided as conflicts (searched again in a later batch), so
  // larger tail batches (tail_batch_max) do not lengthen the serial rip-up
  // phase by their re-searches. 0: every settled vertex counts.
  bool ripup_effort_conflicts_excluded = false;
  // A fixed budget of that phase in settled vertices, in place of the
  // percentage (for comparisons at an equal budget, e.g. the budget a
  // default run logged in DRT-9136). 0: final_ripup_effort.
  int64_t final_ripup_budget = 0;
  // Skip that phase when more than this percentage of the nets failed.
  int final_ripup_max_share = 2;
  // Catalogue paths start at the pin's OpenROAD access points (1) or
  // anywhere on the pin (0, as in the paper). FlexDR's clean-up needs every
  // pin contact at an access point, so 0 is incompatible with it (see
  // OffTrack.md; the final contact guard then drops affected nets).
  bool offtrack_from_access_points = true;
  // Prefer the access point pin access chose for each pin (its
  // conflict-free pattern): paths from the pin's other access points start
  // with a penalty of pin_cost pitches, as for the access stubs; 0 treats
  // all access points alike.
  bool offtrack_preferred_access_points = true;
  // Sec. 4.1 searches: per-worker label budget in MB (0: unlimited); a
  // search that would exceed it stops without a path (counted).
  int search_label_mb = 512;
  // Fast grid (Sec. 3.6) in the dense encoding: layers keep one word per
  // track-graph vertex (and via position) instead of runs, within
  // fast_grid_memory_mb MB of dense arrays, layers taken bottom-up; the
  // other layers keep the run encoding. Same words, same results.
  bool fast_grid_dense = false;
  int fast_grid_memory_mb = 8192;
  // FlexDR clean-up after BonnRoute. FlexDR's iterations 0-2 run in
  // incremental mode, which reroutes only nets without initial routing
  // (the nets BonnRoute left) and keeps every BonnRoute net fixed; iteration
  // 3 is the first that rips up nets at violations. After iteration 0 (which
  // routes the left nets and checks every tile), the clean-up runs the first
  // cleanup_incremental_iterations of the incremental iterations 1 and 2 and
  // continues with iteration 3; 2 is FlexDR's own sequence. With 1, the
  // left nets get one more reroute pass before iteration 3 and the result
  // stays within the spread of BonnRoute routings; 0 hands more markers to
  // iteration 3, which then rips up BonnRoute nets (more vias).
  int cleanup_incremental_iterations = 1;
  // In the incremental iterations, a FlexDR worker whose nets all carry
  // initial routing has nothing to reroute; it runs the design-rule check
  // (with FlexGC's patches) and writes back the same result without
  // building the maze grid graph. The routing result equals that of the
  // full worker (bit-identical result files).
  bool cleanup_check_only_workers = true;
  // Before FlexDR's iteration 0, run the connectivity check that otherwise
  // follows iteration 0 over BonnRoute's nets (all are marked modified): it
  // removes objects that are not on a connecting path (e.g. a pin-layer
  // segment next to a via whose pad already touches the pin) and shrinks
  // overshooting segments. Iteration 0 then checks the geometry FlexDR keeps;
  // the Recheck markers of this check are dropped, since iteration 0 checks
  // every tile.
  bool cleanup_normalize_initial_routing = false;
  // In FlexDR's marker-driven grid iterations (iteration 3 on), workers of
  // later checkerboard batches start ahead of their batch on idle threads;
  // a worker's result is kept when it read the design state its batch starts
  // from. The routing result equals that of the batch-by-batch run
  // (bit-identical result files); 0 runs batch by batch.
  bool cleanup_grid_workers_ahead = true;
  // FlexDR's clean-up after the first routing, iterations 3 to 16: how often
  // a worker reroutes one net (FlexDR's mazeEndIter, 8 in its own flow). A
  // worker that still has violations after that many reroutes of a net ends;
  // its result is written back when it has no more markers than before, and
  // the next iteration continues with shifted tiles. 0 keeps FlexDR's value.
  // From iteration 17 on FlexDR's value applies.
  // Default 0: with 4 the FlexDR time spread over ibex routings narrows on
  // the tuning routings but widens on held-out ones (round ga6).
  int cleanup_reroutes_per_net = 0;
  // Track assignment runs after BonnRoute instead of before it, for the nets
  // left to FlexDR only, and sees BonnRoute's routing as obstacles. FlexDR
  // reads track assignment only for those nets (their tile-boundary pins in
  // iteration 0). With reservations 1, BonnRoute reads track assignment
  // before routing, so track assignment keeps its place before BonnRoute.
  bool left_net_track_assignment = true;
  // Debug (Sec. 5.1): record the global routing tiles that the search and
  // the commit of every net read and write in the shape grid and the fast
  // grid, compare them with the net's tile box and its partition region, and
  // report how far beyond the net's core tile box each step reaches. The
  // routing is unchanged.
  bool partition_check = false;
  // Debug (dependency replay): with partition_check, write per region net
  // (partitions) and per rip-up box task the seconds, the schedule box and
  // the accessed tiles per structure, layer and read / write (tile boxes for
  // region nets, tiles for rip-up tasks) and the nets a rip-up task changed
  // to $RESULTS_DIR/bonnroute_dependency_trace.txt. Setting it sets
  // partition_check. The routing is unchanged.
  bool dependency_trace = false;
  // Ahrens et al., arXiv:2111.06169, Sec. 5.2 incremental routing: a
  // detailed_route -bonnroute call on a design that already carries routing
  // (the antenna-repair ECO) routes the nets that need routing against the
  // existing routing (rip-up and reroute included) and writes them; FlexDR
  // then only checks the changed regions and repairs markers
  // (BrIncremental.cpp). Default in every profile: the full-flow E2E (eco
  // line, round r1) gives final DRC 0 and no more antenna nets than the
  // FlexDR ECO on jpeg/ibex sky130hd, aes nangate45/asap7.
  bool incremental = true;
  // Keep the routing space alive across detailed_route calls and update it
  // with the design changes since the previous call; 0 rebuilds it from the
  // design and its routing in every incremental call.
  bool incremental_resident = true;
  // Input wires shorter than this many pitches of their layer that connect
  // only to illegal input wires get no reservation (Sec. 5.2).
  int reserve_short_pitches = 2;
  // In incremental calls (whatever the profile), legal input wires are
  // reservations of their net at the discount (Sec. 5.2: the input is a
  // legal detailed routing); the first call is unaffected. The search of
  // incremental calls is then the vertex search with discounted runs.
  bool incremental_reservations = true;
  // Future cost of incremental calls (values of `potential`); 3 is the
  // general model with discounts used in the paper's Tables 4 and 6.
  int incremental_potential = 3;
  // Rip-up sequences (ripup_depth, ripup_passes) for the nets that fail in
  // incremental calls, whatever `ripup` is for the first call.
  bool incremental_ripup = true;
  // Debug: at the start of every incremental call, count the routing shapes
  // of the design (wires, vias, patches of nets whose routing did not
  // change) that no shape of the routing space covers. 1 reports the count,
  // 2 also stops with an error when it is not 0.
  int incremental_check = 0;
  // C5 step 1 (design lead/cleanup/C5-revised.md): after routing, check
  // every via whose bottom pad lies on a pin shape of its own net with
  // FlexPA's via check (FlexGC, as FlexPA::isViaViolationFree), and with
  // every other via definition of the access point when it fails. Prints
  // statistics and one line per failing via; the routing is unchanged.
  int start_check = 0;
  // After FlexDR, extend routed metal on IO pins to the layer's minimum
  // area without counting pin-purpose shapes. Each patch is checked against
  // the final geometry with FlexGC before it is written. Summary DRT-9263.
  bool io_pin_area = true;
  // C5 step 2 (lead/cleanup/C5-revised.md, rule 3'): the searches place a
  // via whose bottom pad overlaps a pin shape of its own net only at an
  // access point whose FlexPA via definitions include the via, or where the
  // start check passes with same-net markers only; elsewhere on the pin the
  // path leaves the pin on the pin layer first. The check runs once per
  // instance class (master, orientation, offset in master coordinates, via
  // definition) before routing; the searches look the result up.
  bool planar_first = true;
  // Incremental calls: a new antenna-cell pin (diode of repair_antennas
  // -keep_detailed_wires) of a net whose kept routing joins all its other
  // pins connects on the lowest track layer to the kept metal of the nearest
  // other pin's node at that layer (the gate the diode serves), so that the
  // diode is on the gate's node for every antenna ratio of that layer and
  // above. Falls back to the unrestricted connection when that fails.
  bool incremental_diode_node = true;
  bool verbose = true;
  // Algorithm profile: 0 = 2013-improved (default), 1 = 2013-reimpl,
  // 2 = 2021 (2013-reimpl plus Ahrens et al., arXiv:2111.06169).
  // Setting `profile` assigns the bundle in applyProfile(); parameters set
  // afterwards override it.
  int profile = 0;

  bool set(const std::string& name, double value);
  void applyProfile(int which);
  const char* profileName() const;
};

class BonnRouteDR
{
 public:
  // block: source of the global routing congestion (dbGCellGrid) used by
  // wire spreading; may be nullptr.
  BonnRouteDR(frDesign* design,
              utl::Logger* logger,
              RouterConfiguration* router_cfg,
              const BonnRouteDRParams& params,
              odb::dbBlock* block = nullptr);
  ~BonnRouteDR();
  bool spreadLayer(int z) const;
  void buildJogTiles();

  void main();

  // Incremental mode (BrIncremental.cpp, arXiv:2111.06169 Sec. 5.2).
  // True if this routing space belongs to `design` and can be updated.
  bool residentFor(const frDesign* design) const;
  // Builds the routing space from the design and its existing routing
  // (first call without BonnRoute, or incremental_resident 0).
  void buildFromDesign();
  // After main(): turns the routes written by main() into input routing and
  // releases all reservations, so that later incremental calls start from
  // the routing space of the design as written.
  void keepRoutingSpace();
  // Updates the routing space with the design changes, routes the nets
  // that need routing against the existing routing and writes them;
  // `regions` receives the areas FlexDR checks afterwards.
  void routeIncremental(std::vector<odb::Rect>& regions);

 private:
  struct LayerInfo
  {
    frLayerNum fr = 0;
    frLayer* layer = nullptr;
    frCoord width = 0;
    frCoord pitch = 0;
    bool horizontal = true;
    // Via to the next track layer and its shapes relative to the origin.
    const frViaDef* up_via = nullptr;
    odb::Rect up_pad;    // bottom pad on this layer
    odb::Rect up_cut;    // cut on fr + 1
    odb::Rect down_pad;  // top pad of the via coming from below
    frCoord cut_spacing = 0;
    frCoord eol_width = 0;
    frCoord eol_space = 0;
  };
  // Routed geometry of one net, as written to frNet.
  // Routed geometry; wire_type selects the width (wires) and the NDR
  // spacing of the shapes (Sec. 3.2); access paths keep the default rule,
  // as drt tapers NDR wires at pins.
  struct OutWire
  {
    frLayerNum layer;
    odb::Point begin;
    odb::Point end;
    // A min-area extension beyond a via at `begin`: written as a patch wire
    // (FlexDR's form; its connectivity check removes dangling segments).
    bool patch = false;
    int wire_type = 0;
  };
  struct OutVia
  {
    const frViaDef* def;
    odb::Point point;
    int wire_type = 0;
  };
  // A shape on a drt layer; role selects its shape class (Sec. 3.1), ndr
  // the nondefault rule of a routed shape of an NDR net.
  struct RShape
  {
    frLayerNum layer;
    odb::Rect rect;
    br::ShapeRole role = br::ShapeRole::kFixed;
    int ndr = 0;
  };
  struct AccessOption
  {
    size_t vertex = 0;
    br::Cost cost = 0;
    int z = 0;
    std::vector<OutWire> wires;
    std::vector<OutVia> vias;
    // Planar access: the stub ends on the pin at `ap` (truncated end).
    bool planar = false;
    odb::Point ap;
    frLayerNum ap_layer = 0;
    // The path from the pin to the vertex as route nodes (drt layer
    // numbers); the post-processing combines it with on-track paths.
    std::vector<br::RouteNode> route;
    // Generated on the blockage grid (legality is checked for all shapes).
    bool offtrack = false;
    // Index in `vias` of the access point's own via, whose pin-layer pad
    // and cut pin access has checked; -1 if none.
    int pin_via = -1;
  };
  struct Term
  {
    std::vector<AccessOption> options;
    std::vector<RShape> reservation;
    bool reservation_active = false;
    int chosen = -1;
    // The pin, for dynamic generation of access paths.
    frInstTerm* inst_term = nullptr;
    frBTerm* bterm = nullptr;
    int primary = -1;
  };
  // Why a net is left to FlexDR (reportLeftNets): the step at which its
  // last connect() failed, the outcome of its last rip-up sequence, and the
  // output guard that dropped its routing.
  enum LeftStep : uint8_t
  {
    kLeftNone = 0,
    kLeftAccess,        // a pin without a feasible access option
    kLeftNoPath,        // a connection without a path in the net's area
    kLeftSearchBudget,  // a search stopped at its label or effort budget
    kNumLeftSteps
  };
  enum RipupOutcome : uint8_t
  {
    kRipupNotStarted = 0,
    kRipupBudget,        // the phase's effort budget was used up
    kRipupDepth,         // ripup_depth reached
    kRipupNoPath,        // no path in rip-up mode (LeftStep says why)
    kRipupVictimLimit,   // a victim was ripped up max_ripups times already
    kRipupVictimFailed,  // a victim could not be routed again
    kRipupVictimLevel,   // a victim above the pass's rip-up level
    kRipupVictimOutside,  // a victim not resident in the partition region
    kNumRipupOutcomes
  };
  // Rip-up sequences of one worker: where the running sequence failed first
  // (outcome and depth of that call) and counters.
  struct RipupState
  {
    uint8_t outcome = kRipupNotStarted;
    int outcome_depth = 0;
    int64_t sequences = 0;
    int64_t routed = 0;
    int64_t budget_skips = 0;
    int64_t events = 0;  // victims ripped up
    int64_t rollbacks = 0;
    int64_t plain_connects = 0;
    int64_t mode_connects = 0;
    double plain_seconds = 0;
    double mode_seconds = 0;
    // Victims ripped up, by their rip-up level (index 1..6).
    std::array<int64_t, br::kMaxRemovableLevel + 1> victims_by_level{};
    // Sequences that failed at a victim above the pass's level, at a victim
    // outside the partition region.
    int64_t victim_level_failures = 0;
    int64_t victim_outside_failures = 0;
    int64_t settled = 0;
  };
  enum LeftGuard : uint8_t
  {
    kGuardNone = 0,
    kGuardPinContact,      // a pin contact away from the pin's access points
    kGuardOnTrack,         // a wire off the tracks of an on-track-only layer
    kGuardMissingContact,  // a pin without a contact FlexDR recognises
  };
  // What blocks the access options of a pin at the end of routing, the
  // least severe kind over its options.
  enum AccessBlock : uint8_t
  {
    kAccessNoOption = 0,  // the pin has no access option
    kAccessFeasible,      // an option is legal now
    kAccessRouted,        // routed shapes of other nets only
    kAccessReservation,   // reservations of other nets
    kAccessFixed,         // blockages or pins of other nets
    kNumAccessBlocks
  };
  // D1 S2: the connections a net closed so far, from which connect
  // resumes: joined terms (term 0 always), terms whose access path is part
  // of a post-processed connection, chosen access options, track vertices
  // that carry metal (sources of later connections), and the geometry.
  struct PartialState
  {
    bool valid = false;
    int net = -1;  // the net whose connect produced it
    int joined = 0;
    std::vector<char> connected;
    std::vector<char> consumed;
    std::vector<int> chosen;
    std::vector<size_t> tree;
    std::vector<OutWire> wires;
    std::vector<OutVia> vias;
  };
  struct NetData
  {
    frNet* net = nullptr;
    int id = 0;
    uint8_t left_step = kLeftNone;
    int left_term = -1;
    uint8_t ripup_outcome = kRipupNotStarted;
    uint8_t left_guard = kGuardNone;
    std::vector<Term> terms;
    std::vector<OutWire> wires;
    std::vector<OutVia> vias;
    std::vector<RShape> shapes;  // committed shapes
    std::vector<RShape> pins;    // pin shapes in the shape grid
    std::vector<std::pair<frLayerNum, odb::Point>> pin_points;
    // planar_first: lower vertices of vias the searches may not place
    // (pad on an own pin, neither FlexPA-approved nor clean); sorted.
    std::vector<size_t> forbidden_pin_vias;
    // Legal input wires reserved for this net (arXiv:2111.06169 Sec. 5.2)
    // and the track-graph edges they cover.
    std::vector<RShape> input_reservation;
    bool input_reserved = false;
    std::vector<br::DiscountRun> discount_runs;
    // Guide tiles (z, tile x, tile y), deduplicated.
    std::vector<std::array<int, 3>> guide_tiles;
    // Tiles (z, tile x * tiles_y + tile y) with GR extra space for the net.
    std::vector<std::pair<int, int>> extra_space_tiles;
    // Routing read from the design (incremental mode, Sec. 5.2): routed
    // shapes of the net in the routing space while input_present.
    std::vector<RShape> input_shapes;
    bool input_present = false;
    // The input shapes were ripped up in this call: the frNet's routing is
    // replaced when the net is written.
    bool input_removed = false;
    // Terms built for this incremental call (access and guide tiles).
    bool prepared = false;
    // Kept pin-to-pin paths of the input routing are a terminal of the
    // net in this call (Fig. 9 of arXiv:2111.06169).
    bool keeps_input = false;
    // The frNet was deleted from the design.
    bool dead = false;
    // Fingerprints of the frNet's routing and pins at the last update.
    uint64_t routing_fp = 0;
    uint64_t terms_fp = 0;
    int ripups = 0;
    bool routed = false;
    bool failed = false;
    bool skip = false;
    // Sec. 4.4 / 5.1: the net's connections failed in the guide area; its
    // next attempt starts in the area dilated by retry_dilation.
    bool dilated_area = false;
    // Sec. 5.1: the net has a nondefault wire type, or a pin, reservation
    // or access path shape of a wide class; its updates and queries reach
    // the layer's whole rule distance (buildPartitionHalos, markWideNets).
    bool wide_shapes = false;
    // Highest track layer of the net's search area in this call (-1: no
    // limit); incremental_diode_node sets 0.
    int max_layer = -1;
    int64_t hpwl = 0;
    // Wire type (Sec. 3.2): 0 = default rule, else the net's NDR.
    int wire_type = 0;
    // keep_partial: shapes holds the committed shapes of partial_state.
    bool partial = false;
    PartialState partial_state;
  };

  // partition_check (Sec. 5.1): the steps of a net's unit of work, the
  // phases they run in, and per (phase, step, read / write) the number of
  // units by the tiles their accesses reach beyond the net's core tile box.
  enum FootprintStep
  {
    kStepAccess = 0,   // own-shape overlay, access option legality
    kStepRegenerate,   // dynamic generation of access paths
    kStepSearch,       // connections in the guide area
    kStepRetrySearch,  // the net again in the area dilated by retry_dilation
    kStepCommit,       // conflict check and insertion of the routed shapes
    kStepRelease,      // release of the reservations
    kNumFootprintSteps
  };
  enum FootprintPhase
  {
    kPhaseRegion = 0,  // a net of a partition region
    kPhaseTailSearch,  // speculative search of a tail attempt
    kPhaseTailCommit,  // decision and commit of a tail attempt
    kNumFootprintPhases
  };
  struct FootprintStats
  {
    static constexpr int kReachBins = 10;  // 0 .. 8 tiles, 9 or more
    using Histogram = std::array<int64_t, kReachBins>;
    std::array<int64_t, kNumFootprintPhases> units{};
    std::array<int64_t, kNumFootprintPhases> outside_box_units{};
    std::array<int64_t, kNumFootprintPhases> outside_region_units{};
    // [phase][step][0 reads, 1 writes]
    std::array<std::array<std::array<Histogram, 2>, kNumFootprintSteps>,
               kNumFootprintPhases>
        reach{};
    // The first units whose accesses left their partition region.
    std::vector<std::string> region_examples;
  };

  // Per-thread search state (Sec. 5.1: nets in disjoint areas are routed
  // in parallel).
  struct Worker
  {
    // partition_check: the accesses of the net being routed, its phase, and
    // the statistics of the units finished on this worker.
    br::AccessFootprint footprint;
    FootprintPhase footprint_phase = kPhaseRegion;
    FootprintStats footprint_stats;
    // Searching in a speculative batch: read-only on shared routing space.
    bool speculative = false;
    // Sec. 5.1: the partition region of the thread (nullptr: the whole
    // chip). A net continues in the dilated area only where that area's tile
    // box lies inside the region; otherwise connect() stops and sets
    // area_deferred, and the net continues in a later partition or the tail.
    const br::TileBox* limit = nullptr;
    bool area_deferred = false;
    // Rip-up mode: shapes up to this rip-up level may be removed.
    br::RipupLevel max_ripup_level = br::kMaxRemovableLevel;
    // region_ripup_effort: the rip-up sequences of this worker's regions,
    // the vertices they settled in the current region, and the rip-up
    // history increments of its regions (merged into the fast grid after
    // each partition; increments commute).
    RipupState ripup;
    int64_t region_ripup_settled = 0;
    std::vector<size_t> history_bumps;
    // keep_partial: progress of the last failed connectInRegion, partial
    // commits of the regions, connects that resumed a partial net.
    PartialState progress;
    int64_t partial_commits = 0;
    int64_t partial_resumes = 0;
    int64_t tail_partial_commits = 0;
    // Vertices settled by this worker's tail searches whose attempt was
    // decided as a conflict (ripup_effort_conflicts_excluded).
    int64_t conflict_settled = 0;
    int64_t area_deferrals = 0;
    // Settled-vertex budget of the current connect (0: none), the settled
    // count at its start, and whether a search was skipped for it.
    int64_t effort_budget = 0;
    int64_t effort_start = 0;
    bool effort_exceeded = false;
    // Spreading tiles of the net being routed (extra space mode).
    std::vector<std::vector<uint8_t>> spread_mask;
    br::Region region;
    std::unique_ptr<br::PathSearch> search;
    std::unique_ptr<br::IntervalSearch> interval_search;
    br::PotentialOracle potential;
    int64_t potential_fallbacks = 0;
    double potential_seconds = 0;
    int64_t potential_checks = 0;
    int64_t potential_mismatches = 0;
    // potential_check: pops and settled vertices of the checked searches
    // with the potential and of their pi_H repetitions (same state).
    int64_t potential_check_pops = 0;
    int64_t potential_check_settled = 0;
    int64_t potential_check_plain_pops = 0;
    int64_t potential_check_plain_settled = 0;
    // Guide path of the last buildPathRegion: planar tile steps minus the
    // tile distance of its ends.
    int path_detour = 0;
    // Scratch of the guide-path search for connection regions.
    std::vector<uint32_t> tile_guide;
    std::vector<uint32_t> tile_seen;
    std::vector<uint32_t> tile_target;
    std::vector<int64_t> tile_parent;
    uint32_t tile_epoch = 0;
    // Sec. 5.1: a thread changes the routing space only inside its region;
    // shape grid and fast grid are chunked by global routing tile, so the
    // worker updates them in place. Its configuration transitions are
    // memoised here.
    br::ShapeGridWriteCache cache;
    // The net being connected: its own shapes, the fast-grid overlay with
    // them removed and the reader of connect's own legality checks.
    std::vector<br::FastGridBuilder::OwnShape> own_shapes;
    br::FastGridOverlay overlay;
    br::FastGridReader reader;
    int64_t connects = 0;
    int64_t own_shape_count = 0;
    int64_t fast_grid_queries = 0;
    int64_t checker_queries = 0;
    int64_t aborted_searches = 0;
    br::CleanupStats cleanup;
    int64_t regenerations = 0;
    // connect() calls that failed because a pin had no feasible access
    // option even after dynamic generation.
    int64_t access_failures = 0;
    // (net, term) of each such failure (offtrack_access_report only).
    std::vector<std::pair<int, int>> access_failure_terms;
    int64_t regenerated_paths = 0;
    int64_t searches = 0;
    int64_t path_region_searches = 0;
    int64_t heap_pops = 0;
    int64_t held_pops = 0;
    int64_t settled = 0;
  };
  // Saved routing of a net, for restoring a failed rip-up sequence.
  struct SavedRouting
  {
    int net;
    bool was_routed;
    std::vector<OutWire> wires;
    std::vector<OutVia> vias;
    std::vector<std::pair<frLayerNum, odb::Point>> pin_points;
    std::vector<int> chosen;
    std::vector<char> reserved;
    bool input_reserved;
    bool was_partial = false;
    PartialState partial_state;
    bool input_present;
  };

  void buildTrackGraph();
  void buildFixed();
  void buildNets();
  void buildAccess(NetData& net);
  void reserveAccess(NetData& net);

  // Off-track pin access (Sec. 3.7-3.8, 4.3-4.4); see OffTrack.md.
  class Legality;
  struct AccessStats
  {
    int64_t pins = 0;
    int64_t searches = 0;
    int64_t paths = 0;
    int64_t empty = 0;
    int64_t fallback = 0;
    int64_t drc_rejected = 0;
    int64_t same_net_rejected = 0;
    int64_t offgrid_rejected = 0;
    int64_t grid_vertices = 0;
    int64_t settled = 0;
    double search_seconds = 0;
    int64_t instances = 0;
    int64_t classes = 0;
    int64_t class_reused = 0;
    int64_t reselected = 0;
    int64_t mapping_failed = 0;
    int64_t selections = 0;
    int64_t selection_nodes = 0;
    int64_t selection_optimal = 0;
    int64_t primaries = 0;
    double seconds = 0;
    // Why catalogues stayed empty, and the first such pins.
    int64_t empty_no_pin = 0;
    int64_t empty_no_targets = 0;
    int64_t empty_no_path = 0;
    int64_t empty_rejected = 0;
    int64_t macro_pins = 0;
    int64_t member_regenerated = 0;
    // offtrack_algorithm3_fallback: instances with a pin without a primary
    // path on the Hanan grid, and those whose Algorithm 3 catalogues assign
    // more pins.
    int64_t fallback_instances = 0;
    int64_t fallback_improved = 0;
    // offtrack_complete_selection: instances whose selection is complete
    // after the variants, instances that keep a pin without a primary path,
    // and the instances per kept variant (index into the variant list).
    int64_t completion_complete = 0;
    int64_t completion_incomplete = 0;
    std::array<int64_t, 8> completion_kept{};
    // D10 pin classes: classes, mapped paths of member pins and those
    // illegal at the member.
    int64_t pin_classes = 0;
    int64_t mapped_options = 0;
    int64_t mapped_dropped = 0;
    // Sampled comparison of a member pin's mapped catalogue with the one
    // generated at the member: sampled pins, equal (target, cost) sets,
    // pins whose whole search environment equals the representative's
    // (and of those, related by a translation, and their unequal sets).
    int64_t sample_pins = 0;
    int64_t sample_equal = 0;
    int64_t sample_mapped_only = 0;  // (target, cost) only in the mapped
    int64_t sample_own_only = 0;     // (target, cost) only in the own
    int64_t sample_same_environment = 0;
    int64_t sample_same_environment_translated = 0;
    int64_t sample_same_environment_translated_unequal = 0;
    // M3 (c): blockage-grid lines and vertices inside the window with the
    // grid in use and with the other grid (Hanan / Algorithm 3).
    int64_t grid_lines = 0;
    int64_t grid_lines_other = 0;
    int64_t grid_vertices_other = 0;
    int64_t search_states = 0;
    std::vector<std::string> empty_examples;
    // Thread seconds of the catalogue steps besides the search: on-track
    // targets, obstacle collection, blockage grid, exact path checks, and
    // the conflict-free selection.
    double target_seconds = 0;
    double blockage_seconds = 0;
    double grid_seconds = 0;
    double check_seconds = 0;
    double selection_seconds = 0;
    void add(const AccessStats& o);
  };
  // Wall seconds of the steps of buildOffTrackAccess (verbose log).
  struct AccessPhases
  {
    double pin_map = 0;
    double signatures = 0;
    double classification = 0;
    double class_catalogues = 0;
    double members = 0;
    double member_mapping = 0;
    double member_legality = 0;
    double member_catalogues = 0;
    double member_solve = 0;
    double member_reservations = 0;
    double top_level_pins = 0;
  };
  AccessPhases access_phases_;
  // Wall seconds of the setup steps of main() (verbose log).
  struct SetupPhases
  {
    double track_graph = 0;
    double routing_space = 0;
    double spreading = 0;
    double workers = 0;
    double nets = 0;
    double fixed_insert = 0;  // total: query, classes, shape-grid insert
    double fixed_query = 0;
    double fixed_classes = 0;
    double fast_grid_rebuild = 0;
    double compaction = 0;
    double same_net_rules = 0;
    double guide_tiles_and_stubs = 0;
    double offtrack_access = 0;
    double stub_reservations = 0;
    double input_wires = 0;
  };
  SetupPhases setup_phases_;
  // Fingerprints of the access catalogues, primaries and reservations after
  // setup, and of the routing before writeNets (verbose log).
  uint64_t accessDigest() const;
  uint64_t routingDigest() const;
  void buildSameNetRules();
  // pre_skip: net skip flags before buildAccess.
  void buildOffTrackAccess(const std::vector<char>& pre_skip);
  // offtrack_access_report: the catalogues and primary paths after setup,
  // and the access failures and left nets after routing.
  void reportOffTrackAccess() const;
  void reportAccessFailures() const;
  std::string termName(const Term& term) const;
  // Rectangles of the pin on its access layer (the lowest routing layer
  // that is a track layer or directly below the lowest one).
  // A start of catalogue paths: an access point with its up via, if any.
  struct AccessSource
  {
    frLayerNum layer = 0;
    odb::Point point;
    const frViaDef* via = nullptr;
    // The access point pin access chose for the pin (always for top-level
    // pins).
    bool preferred = true;
  };
  // The access points of the pin of `term` (FlexDR's: the instance's pin
  // access index, no-rotation transform; index 0 for top-level pins).
  // `preferred_only`: the access point pin access chose per pin (the
  // instance pin's access point), else all of the pin access index.
  void accessSources(const Term& term,
                     std::vector<AccessSource>& out,
                     bool preferred_only = false) const;
  // Wires of a layer that must stay on its tracks in the preferred
  // direction (multi-patterned, RECTONLY, RIGHTWAYONGRIDONLY: drt treats
  // them as unidirectional; FlexDR cannot route a wrong-way piece there).
  bool onTracksOnly(frLayerNum layer) const;
  bool wireOnTracks(const OutWire& wire) const;
  int64_t wrong_way_nets_ = 0;
  int64_t connected_vias_ = 0;
  int64_t notch_fills_ = 0;
  void accessPointsOf(const frBlockObject* pin,
                      std::vector<br::RouteNode>& out) const;
  // Pin contacts of the net's geometry that are not access points (wire
  // ends truncated on a pin and via pads at the via origin); see
  // br::pinContactViolations.
  // With `vias`, via pads are counted instead of truncated ends.
  int64_t pinContactViolations(const NetData& net,
                               const std::vector<OutWire>& wires,
                               const std::vector<OutVia>& vias,
                               bool count_vias = false) const;
  // True unless p lies on a pin of `net` on `layer` away from its access
  // points (no new nodes there in the post-processing).
  bool nodeAllowed(int net, frLayerNum layer, const odb::Point& p) const;
  // p is an access point (on `layer`) of a pin of the net there.
  bool atAccessPoint(const NetData& net,
                     frLayerNum layer,
                     const odb::Point& p) const;
  // Pins of the net without a contact FlexDR recognises (a wire end in
  // the net's pin points, or a via at an access point), with writeNets'
  // split `wires`; such a net is left to FlexDR.
  int pinsWithoutContact(const NetData& net,
                         const std::vector<OutWire>& wires) const;
  int64_t uncontacted_nets_ = 0;
  int64_t uncontacted_pins_ = 0;
  int64_t guard_nets_ = 0;
  int64_t guard_contacts_ = 0;
  // Via pads on pins away from access points (not a FlexDR term: retained
  // vias are never marked connected); reported only.
  int64_t guard_via_contacts_ = 0;
  bool pinShapes(const Term& term,
                 frLayerNum& layer,
                 std::vector<odb::Rect>& rects) const;
  // Catalogue of the pin's off-track access paths against the current
  // routing space (Sec. 4.3).
  // `extra_radius` pitches beyond offtrack_radius (dynamic generation).
  void offTrackCatalogue(int net,
                         const Term& term,
                         const Worker* worker,
                         std::vector<AccessOption>& out,
                         AccessStats& stats,
                         const br::CatalogueVariant& variant
                         = br::CatalogueVariant()) const;
  void candidateOf(const AccessOption& option,
                   int net,
                   br::AccessCandidate& candidate) const;
  // Chooses primary paths for `pins` (net ids per pin) by destructive
  // bounding; -1 for pins without a feasible path.
  std::vector<int> choosePrimaries(
      const std::vector<const std::vector<AccessOption>*>& catalogues,
      const std::vector<int>& nets,
      AccessStats& stats) const;
  bool regenerateAccess(NetData& net, int t, Worker& worker);
  // term.options := catalogue followed by the term's access stubs.
  static void keepStubs(Term& term, std::vector<AccessOption>& catalogue);
  bool optionLegal(const AccessOption& option,
                   int net,
                   const Worker* worker) const;
  const frViaDef* viaDefAbove(frLayerNum lower) const;
  // Route nodes -> wires and vias of the net; via defs of `hints` at their
  // points take precedence over the defaults.
  void routeToGeometry(const std::vector<br::RouteNode>& nodes,
                       const std::vector<OutVia>& hints,
                       NetData& net) const;
  // Sec. 4.4 post-processing of one connection.
  void postprocessConnection(std::vector<br::RouteNode>& nodes,
                             const std::vector<OutVia>& hints,
                             int net,
                             Worker& worker,
                             std::vector<OutWire>& patches) const;
  // `local`: the worker whose region contains the net (Sec. 5.1).
  // Releasing with `pending` collects the fast-grid recomputation there for
  // one later flush; without it the release flushes at once.
  void setReservations(NetData& net,
                       bool active,
                       Worker* local = nullptr,
                       br::FastGridBuilder::PendingRemovals* pending = nullptr);
  void reserveInputWires();
  void setInputReservation(NetData& net, bool active, Worker* local = nullptr);
  // Diff-net conflict with the shapes of an R-tree (input wires only).
  bool treeConflicts(const RTree<int>& tree,
                     const RShape& shape,
                     int net) const;
  // Cost unit of the search: discount_denominator when discounts are used.
  br::Cost unit() const { return graph_.wire_unit; }
  // Reservations need wire costs per DBU above 1 (discounted runs), which
  // only the vertex search supports; spreading and potentials work with
  // both searches.
  bool modernSearch() const { return params_.reservations; }
  // Tiles where global routing left room for spreading (Sec. 4.2).
  void buildSpreadTiles();
  void setReservation(Term& term,
                      int net,
                      bool active,
                      Worker* local = nullptr);
  void accessShapes(const AccessOption& option,
                    std::vector<RShape>& shapes) const;
  // extra_dilation enlarges the routing area of later rip-up passes.
  bool connect(NetData& net,
               bool ripup,
               Worker& worker,
               int extra_dilation = 0);
  bool connectWithoutOwnShapes(NetData& net,
                               bool ripup,
                               Worker& worker,
                               int extra_dilation);
  bool connectInRegion(NetData& net,
                       int dilation,
                       bool ripup,
                       const std::vector<std::vector<char>>& valid,
                       Worker& worker);
  // Measurement M5 (debug group bonnroute_search): one line per search of
  // connectInRegion; statistics only.
  void logSearchStatistics(const NetData& net,
                           const Worker& worker,
                           const br::SearchRequest& request,
                           const br::SearchResult& result,
                           bool path_region,
                           int dilation,
                           bool ripup) const;
  void buildRegion(const NetData& net, int dilation, br::Region& region);
  void buildGuideTiles(NetData& net);
  // Region of the guide path from the routed component to the nearest
  // unconnected terminal; false if the guides do not connect them.
  bool buildPathRegion(const NetData& net,
                       const std::vector<char>& connected,
                       const std::vector<size_t>& tree,
                       const std::vector<std::vector<char>>& valid,
                       Worker& worker);
  // Tile box of the net's guides, access options and pins.
  br::TileBox netCoreTileBox(const NetData& net) const;
  // Highest track layer of the net's routing area: guides and access
  // options plus layer_dilation.
  int netTopLayer(const NetData& net) const;
  // Highest track layer of the access paths of the net's pins: the track
  // layer above each pin's layer (offTrackCatalogue, buildAccess).
  int netAccessTopLayer(const NetData& net) const;
  // Halo of the net's box: the larger of area_halo_tiles_ at netTopLayer
  // and access_halo_tiles_ at netAccessTopLayer.
  int netHaloTiles(const NetData& net) const;
  // Tile box of the net's searches and commit in its guide area dilated by
  // `dilation` tiles, with halo, clipped to the tile grid.
  br::TileBox netAreaTileBox(const NetData& net, int dilation) const;
  // Tile box of everything the net's searches may touch, with halo.
  void netTileBox(const NetData& net,
                  int& xlo,
                  int& ylo,
                  int& xhi,
                  int& yhi) const;
  // netTileBox clipped to the tile grid.
  br::TileBox netClippedTileBox(const NetData& net) const;
  // partition_check: starts / ends recording the accesses of `net`'s unit
  // of work on `worker` (whose footprint the caller installed); `region`
  // is the partition region of the unit, if any.
  void beginFootprint(Worker& worker,
                      const NetData& net,
                      FootprintPhase phase,
                      const br::TileBox& box);
  void endFootprint(Worker& worker,
                    const NetData& net,
                    const br::TileBox* region);
  void reportFootprints();
  void routeParallel(const std::vector<int>& order, std::vector<int>& failed);
  // D6: the batch size policy of the speculative tail.
  br::TailBatchPolicy tailBatchPolicy() const;
  // DRT-9134 line of one tail pass; DRT-9135 when the counts are not
  // conserved.
  void reportTailOutcomes(const std::string& stage,
                          const br::TailOutcomes& outcomes) const;
  // settled_budget: stop starting sequences once the serial worker settled
  // this many vertices in the phase (0: unlimited).
  void routeRipupQueue(const std::vector<int>& failed,
                       int64_t settled_budget = 0);
  // ripup_partitions: one round of region-parallel rip-up sequences over
  // `queue` (pass `pass`); returns the settled vertices of the round.
  int64_t parallelRipupRound(const std::vector<int>& queue,
                             int pass,
                             int round,
                             int64_t settled_budget);
  // Parallel rip-up rounds: rounds, sequences, nets routed, rollbacks,
  // region budget stops, wall seconds.
  int64_t parallel_ripup_rounds_ = 0;
  int64_t parallel_ripup_sequences_ = 0;
  int64_t parallel_ripup_routed_ = 0;
  int64_t parallel_ripup_rollbacks_ = 0;
  int64_t parallel_ripup_budget_stops_ = 0;
  double parallel_ripup_seconds_ = 0;
  // ripup_boxes: one pass of the rip-up queue under OrderedBoxSchedule.
  // `queue` becomes the nets still unrouted; `settled` counts the vertices
  // the phase's sequences settled, against settled_budget.
  void routeRipupBoxes(std::vector<int>& queue,
                       int pass,
                       int64_t settled_budget,
                       int64_t& settled);
  // Box steps of a sequence: margin, 4 x margin, whole chip.
  static constexpr int kRipupBoxSteps = 3;
  struct RipupBoxStats
  {
    int64_t chunks = 0;
    std::array<int64_t, kRipupBoxSteps> sequences{};
    std::array<int64_t, kRipupBoxSteps> routed{};
    int64_t widened = 0;
    // ripup_box_chunk_retry: sequences run again in their own chunk.
    int64_t retried_in_chunk = 0;
    int64_t longest_chain = 0;
    // Sum of the chunks' longest chains (tasks) and critical paths (s).
    int64_t chain_sum = 0;
    double critical_path = 0;
    double wall = 0;
    double busy = 0;
    double capacity = 0;
    br::BusyThreads::Seconds busy_threads{};
  };
  RipupBoxStats ripup_box_stats_;
  // dependency_trace: the accesses of a footprint as trace text (tiles when
  // `tiles`, else tile boxes), and the trace file.
  static std::string traceAccesses(const br::AccessFootprint& footprint,
                                   bool tiles);
  void writeTrace(const std::vector<std::string>& lines);
  std::unique_ptr<std::ofstream> dependency_trace_out_;
  // try_plain: first try the connection without rip-up. The nets of the
  // rip-up queue have just failed that search, so their sequences start
  // with the rip-up mode search (which finds a legal path at no penalty
  // when one exists).
  // `worker` searches and `state` counts; `resident` (nullptr: the serial
  // phase) is the partition region the sequence stays in: victims must be
  // resident there and the worker's write cache is used.
  bool ripupSequence(int net,
                     int depth,
                     int dilation,
                     std::vector<SavedRouting>& journal,
                     bool try_plain,
                     Worker& worker,
                     RipupState& state,
                     const br::TileBox* resident);
  // region_ripup_effort: rip-up sequence of `net` inside `region` on the
  // region's thread; region_settled_start is the worker's settled count at
  // the region's start. True if the net is routed (and committed).
  bool regionRipup(NetData& net,
                   Worker& worker,
                   const br::TileBox& region,
                   int64_t region_settled_start);
  SavedRouting saveRouting(const NetData& net) const;
  // Rip-up sequences of the serial phase (outcome of the running sequence,
  // counters).
  RipupState serial_ripup_;
  // Kinds of shapes of other nets that `shape` violates: bit 0 blockages,
  // bit 1 pins, bit 2 reservations, bit 3 routed wires and vias.
  int blockerKinds(const RShape& shape, int net) const;
  AccessBlock accessBlock(const NetData& net, const Term& term) const;
  // Log: the nets left to FlexDR by cause; debug group bonnroute_left
  // prints one line per net.
  void reportLeftNets() const;
  int64_t input_wires_ = 0;
  int64_t input_reserved_ = 0;
  // `local`: the worker whose region contains the net (nullptr: serial).
  void restoreRouting(const SavedRouting& saved, Worker* local = nullptr);
  void pathToGeometry(const std::vector<size_t>& path, NetData& net) const;
  void shapesOf(const NetData& net, std::vector<RShape>& shapes) const;
  void viaShapes(const OutVia& via, std::vector<RShape>& shapes) const;
  RShape wireShape(const OutWire& wire) const;
  // Distance rules, shape classes and wire types from the drt tech (Sec.
  // 3.1, 3.2), the shape grid (Sec. 3.3), checker (3.4) and fast grid (3.6).
  void buildRoutingSpace();
  void rebuildFastGrid();
  // Shape-grid layer of a drt layer, or -1: track layers, their cuts, and
  // the pin layer and cut directly below the lowest track layer.
  int gridLayer(frLayerNum layer) const;
  // Thread-safe (per-thread class cache).
  br::ClassId shapeClass(const RShape& shape, int grid_layer) const;
  // Grid layer above a cut layer that holds its inter-layer projection, or
  // -1 (Sec. 3.2).
  int projectionLayer(int grid_layer) const;
  // Inserts / removes a shape in the shape grid (with the inter-layer
  // projection of a cut, Sec. 3.2) and updates the fast grid (Sec. 3.6).
  // `pending`: removals of a group, recomputed by FastGridBuilder::flush().
  void changeShape(const RShape& shape,
                   br::ShapeType type,
                   int owner,
                   bool insert,
                   Worker* local,
                   br::FastGridBuilder::PendingRemovals* pending = nullptr);
  // Sec. 4.4: the shapes of the components being connected (pins,
  // reservations), removed from routing space during the net's searches
  // through a fast-grid overlay.
  void ownShapes(const NetData& net,
                 std::vector<br::FastGridBuilder::OwnShape>& out) const;
  // Track-graph vertex kind usable by `net`: the fast grid, and where it
  // reports a blockage, the distance rule checker with the net's own
  // shapes treated as removed (Sec. 4.4).
  bool usable(size_t v, br::Kind kind, int net) const;
  void commitNet(NetData& net, Worker* local = nullptr);
  // Sec. 5.1: a sequence of balanced partitions into thread regions.
  // initial_regions: regions of the first partition, 0 = parameter;
  // stage names the step in the log.
  void routePartitions(const std::vector<int>& order,
                       std::vector<int>& failed,
                       int initial_regions = 0,
                       const char* stage = "partition");
  void routeSpeculative(const std::vector<int>& order,
                        std::vector<int>& failed,
                        const char* stage);
  // One pass of routeSpeculative: nets over the settled budget (0: none)
  // go to `deferred`, in net order. routeSpeculativePipeline schedules the
  // pass with TailLookahead (tail_lookahead) or TailPipeline.
  void routeSpeculativeBatches(const std::vector<int>& order,
                               std::vector<int>& failed,
                               std::vector<int>& deferred,
                               const std::string& stage,
                               int threads,
                               int64_t budget);
  void routeSpeculativePipeline(const std::vector<int>& order,
                                std::vector<int>& failed,
                                std::vector<int>& deferred,
                                const std::string& stage,
                                int threads,
                                int64_t budget);
  // Estimated search effort of a tail net (orders the searches of a batch).
  double tailSearchEstimate(const NetData& net) const;
  // Routing digest of all nets (FNV-1a over net id, routed flag, wires,
  // vias and pin points).
  uint64_t tailRoutingDigest() const;
  void partitionRegions(const std::vector<int>& nets,
                        const std::vector<std::array<int, 4>>& boxes,
                        int count,
                        bool shifted,
                        std::vector<std::array<int, 4>>& regions) const;
  void ripupNet(NetData& net, Worker* local = nullptr);
  // Removes the committed shapes of `net` from the routing space.
  void uncommitShapes(NetData& net, Worker* local);
  // keep_partial: whether `progress` joined more terms than the net's
  // committed partial routing (term 0 alone without one).
  static bool progressAdds(const NetData& net, const PartialState& progress)
  {
    return progress.valid && progress.net == net.id
           && progress.connected.size() == net.terms.size()
           && progress.joined > (net.partial ? net.partial_state.joined : 1);
  }
  // keep_partial: commits `state` as the net's partial routing.
  void commitPartial(NetData& net, const PartialState& state, Worker* local);
  // Rip-up victims: nets reported by the distance rule checker as blocking
  // the net's new shapes (Sec. 3.4, 4.2).
  void collectVictims(const NetData& net, std::vector<int>& victims) const;
  // Distance rule check of a shape against everything but `net`'s own
  // shapes (the adapter for off-track paths and pin access).
  bool shapeConflicts(const RShape& shape,
                      int net,
                      const Worker* worker = nullptr) const;
  // Conflicts with routed shapes and reservations of other nets only
  // (blockages and pins are not looked at): the pin-layer pad and cut of an
  // access point's via, whose fixed surroundings pin access checked.
  bool routedConflicts(const RShape& shape, int net) const;
  frCoord spacing(int z, frCoord width1, frCoord width2, frCoord prl) const;
  void writeNets();
  // start_check: FlexPA's via check on the written vias (statistics only).
  void checkPathStarts();
  // FlexPA's via check (FlexPA::isViaViolationFree) of via `def` at p on the
  // pin `term` (frInstTerm or frBTerm) of `net`; with own_net_only, only
  // markers whose sources are all `net` count.
  bool viaCheckClean(frBlockObject* term,
                     frNet* net,
                     const frViaDef* def,
                     const odb::Point& p,
                     bool own_net_only) const;
  // planar_first: NetData::forbidden_pin_vias of every net, or of the nets
  // `only` (incremental calls: the nets prepared in the call).
  void buildPinViaLegality(const std::vector<int>* only = nullptr);
  // The access point of `term` (frInstTerm or frBTerm) at p on layer, or
  // nullptr.
  const frAccessPoint* accessPointAt(const frBlockObject* term,
                                     frLayerNum layer,
                                     const odb::Point& p) const;
  // notch_fill for one written net; fills of earlier nets are in `fills`.
  void fillNotches(frNet* net,
                   int id,
                   std::vector<std::pair<frLayerNum, odb::Rect>>& fills);
  // collapse_bridges: replace via up, one wire, via down by a wrong-way
  // wire on the lower layer where that wire is legal and neither via lies
  // on a pin of the net (pin contacts stay as routed). Returns the count.
  int64_t collapseBridges();
  struct BridgeStats
  {
    int64_t candidates = 0;
    int64_t other_wire = 0;
    int64_t patch = 0;
    int64_t other_via_or_pin = 0;
    int64_t pin_contact = 0;  // a bridge via on a pin of the net
    int64_t illegal = 0;
    int64_t illegal_routed = 0;
    int64_t illegal_fixed = 0;
    int64_t collapsed = 0;
  } bridge_stats_, lower_bridge_stats_;

  // Incremental mode (BrIncremental.cpp).
  struct FixedEntry
  {
    RShape shape;
    br::ShapeType type;
    int owner;
  };
  // Fixed shapes of an instance in the routing space and the state they
  // were read in.
  struct InstRecord
  {
    std::string name;
    odb::Rect bbox;
    int orient = 0;
    const void* master = nullptr;
    std::vector<const frNet*> term_nets;
    std::vector<FixedEntry> shapes;
    bool seen = false;
  };
  struct SyncStats
  {
    int64_t instances_changed = 0;
    int64_t nets_reread = 0;
    int64_t nets_added = 0;
    int64_t nets_removed = 0;
    int64_t shapes_removed = 0;
    int64_t shapes_inserted = 0;
    std::vector<odb::Rect> inst_regions;
    std::unordered_set<const frInst*> changed_insts;
  };
  // buildFixed hook: records the fixed shapes of instances.
  void recordFixedShape(const frBlockObject* obj,
                        const RShape& shape,
                        br::ShapeType type,
                        int owner);
  void recordInstanceState(const frInst* inst, InstRecord& record) const;
  void readInstanceShapes(const frInst* inst, InstRecord& record);
  void removeInstanceShapes(InstRecord& record, SyncStats& stats);
  void syncInstances(SyncStats& stats);
  void syncNets(SyncStats& stats);
  // Routed shapes of the frNet's own geometry (wires, vias, patches).
  void routingShapes(const frNet* net, std::vector<RShape>& shapes) const;
  // The routed shape of a patch wire (or notch fill) as routingShapes reads
  // it from the design.
  RShape patchShape(frLayerNum layer, const odb::Rect& rect) const;
  // incremental_check: design routing shapes outside the routing space.
  void checkRoutingSpace(int64_t& nets, int64_t& uncovered) const;
  static uint64_t routingFingerprint(const frNet* net);
  static uint64_t termsFingerprint(const frNet* net);
  int addNetData(frNet* net);
  // Wire type of the net's nondefault rule, 0 for the default rule.
  int ndrType(const frNet* net) const;
  void setInputShapes(NetData& net,
                      bool present,
                      br::FastGridBuilder::PendingRemovals* pending = nullptr);
  void insertDesignRouting();
  bool routableType(const NetData& net) const;
  // Terms, access catalogues and guide tiles of a net routed in this call;
  // false if BonnRoute cannot route it (FlexDR does).
  bool prepareNet(NetData& net,
                  const std::unordered_set<const frInst*>& changed_insts);
  // Rip-up hook: a victim of an incremental call is prepared before it is
  // removed; false keeps it in place.
  bool prepareVictim(int net);
  // Sec. 5.2 / Fig. 9: pieces of a net's input routing (wires, vias,
  // patches) classified into kept pin-to-pin paths, reservations and
  // dropped pieces.
  struct InputPiece
  {
    std::vector<RShape> shapes;
    const frBlockObject* obj = nullptr;
    int kind = 0;  // 0 wire, 1 via, 2 patch
  };
  struct InputClassification
  {
    std::vector<InputPiece> pieces;
    std::vector<char> kept;
    std::vector<char> reserved;
    // Component label per term (net.terms order) and per piece.
    std::vector<int> term_label;
    std::vector<int> piece_label;
    bool complete = false;
    bool any_illegal = false;
  };
  struct ClassifyStats
  {
    int64_t nets = 0;
    int64_t pieces = 0;
    int64_t illegal = 0;
    int64_t kept = 0;
    int64_t reserved = 0;
    int64_t dropped = 0;
    int64_t discount_runs = 0;
  };
  void classifyInput(const NetData& net, InputClassification& out) const;
  // Applies a classification: dropped pieces leave the routing space,
  // reserved ones become reservations of the net (with discount runs on
  // their tracks), kept ones stay.
  void applyClassification(NetData& net,
                           const InputClassification& c,
                           ClassifyStats& stats);
  // Terms of a component joined by kept pieces become one terminal whose
  // options include the track vertices on the kept pieces.
  void mergeKeptComponents(NetData& net, const InputClassification& c);
  // incremental_diode_node: new diode pins of a net, the gates they serve
  // and the targets of one diode's connection.
  struct DiodeConnection
  {
    std::vector<std::pair<const frInstTerm*, const frInstTerm*>> pairs;
    // Access options of the paired gates before the merge.
    std::vector<AccessOption> gate_options;
    std::unordered_set<size_t> targets;
    // Per pair: the gate's own options before the merge, the first cap
    // (lowest violating layer from the antenna checker, else 0), the cap the
    // diode was routed under (-1: unrestricted), whether the checker paired it.
    std::vector<std::vector<AccessOption>> gate_options_by_pair;
    std::vector<int> start_cap;
    std::vector<int> achieved;
    std::vector<char> from_checker;
  };
  // A gate pin with an antenna violation (ant::AntennaChecker) and the
  // track layer of its lowest violating layer.
  struct GateViolation
  {
    std::string inst;
    std::string term;
    int z = 0;
  };
  // Antenna checker violations of the nets `ids` (their odb routing).
  void antennaViolationsOf(
      const std::vector<int>& ids,
      std::unordered_map<int, std::vector<GateViolation>>& out) const;
  // Routes the diodes of a net one at a time, each against its own gate's
  // node under a rising cap; commits the union of the connections.
  bool routeDiodeNet(int id, const InputClassification& c, DiodeConnection& dc);
  std::unordered_map<int, DiodeConnection> diode_connections_;
  // Before mergeKeptComponents: false if the net has no new diode pin or
  // its kept routing is not one component.
  bool planDiodeConnection(const NetData& net,
                           const InputClassification& c,
                           const std::vector<GateViolation>* violations,
                           DiodeConnection& out) const;
  // Targets of the diode connections with the search capped at track layer
  // `cap`: track vertices at or below it on kept metal joined to a paired
  // gate through layers at or below it, and the gates' options there.
  void diodeTargets(const InputClassification& c,
                    DiodeConnection& dc,
                    int cap) const;
  // The pins a and b of `net` are joined by its metal on the layers at or
  // below track layer `cap` (antenna node of that layer).
  bool joinedOnLowLayers(const frNet* net,
                         const frInstTerm* a,
                         const frInstTerm* b,
                         int cap) const;
  std::unordered_map<int, InputClassification> classified_;
  std::unordered_map<const frInst*, InstRecord> inst_records_;
  bool record_fixed_ = false;
  bool incremental_active_ = false;
  // writeNets() writes only the nets marked here (empty: all routed nets).
  std::vector<char> write_mask_;
  const frDesign* built_design_ = nullptr;
  const frBlock* built_block_ = nullptr;
  const std::unordered_set<const frInst*>* victim_changed_insts_ = nullptr;
  // Bulk updates of the incremental mode change the shape grid only; the
  // fast grid is rebuilt afterwards.
  bool grid_only_ = false;
  // Incremental calls: wire costs scaled to the discount denominator and
  // the vertex search set up (enableIncrementalDiscounts).
  bool incremental_discounts_ = false;
  void enableIncrementalDiscounts();
  // Components of the net's routing and pins in the design (touching metal
  // on a layer joins two objects); 1 for a connected net. Debug check of
  // the nets written by an incremental call.
  int routingComponents(const frNet* net) const;
  void gridOnlyChange(const RShape& shape,
                      br::ShapeType type,
                      int owner,
                      bool insert);

  int tileX(frCoord x) const;
  int tileY(frCoord y) const;

  frDesign* design_;
  utl::Logger* logger_;
  RouterConfiguration* router_cfg_;
  BonnRouteDRParams params_;

  br::TrackGraph graph_;
  // Interning is internally locked (classes appear during routing).
  mutable br::ShapeClassTable shape_classes_;
  uint64_t class_generation_ = 0;
  // Wire type 0: default widths and vias (the fast-grid wire type); then
  // one wire type per non-default rule (Sec. 3.2).
  std::vector<br::WireType> wire_types_;
  // Per wire type and track layer: wire width and via definition; nets per
  // wire type; wire type of each NDR. Wire types 1.. are the NDRs, ordered
  // by their number of nets (the first ones are fast, Sec. 3.6).
  std::vector<std::vector<frCoord>> type_widths_;
  std::vector<std::vector<const frViaDef*>> type_vias_;
  std::vector<int> type_nets_;
  std::unordered_map<const frNonDefaultRule*, int> ndr_type_;
  int fast_types_ = 1;
  br::ShapeGrid shape_grid_;
  br::DistanceRuleChecker checker_;
  br::FastGrid grid_;
  br::FastGridBuilder fast_builder_;
  br::ShapeGridWriteCache serial_cache_;
  // Shape-grid layers per drt layer (-1: not covered); number of track
  // grid layers (2 z + 1 for z track layers); the low layers follow.
  std::vector<int> grid_of_fr_;
  int track_grid_layers_ = 0;
  // Tiles between a net's routing area and its partition box border, so
  // that its routing-space updates and queries stay inside (Sec. 5.1): the
  // value for a routing area on all layers.
  int halo_tiles_ = 1;
  // Largest distance at which an update or query of a shape-grid layer
  // reads or writes (rule reach, model, one cell), per layer; and the same
  // with the rule reach of the narrow classes, which holds for nets whose
  // own shapes are all narrow on a layer without wide fast-grid models:
  // their block and tile accesses use narrowReach(), wide shapes come from
  // the locked wide-shape lists.
  std::vector<int64_t> need_by_grid_layer_;
  std::vector<int64_t> narrow_need_by_grid_layer_;
  // The halo per net (buildPartitionHalos). area_halo_tiles_[z]: routing
  // area up to track layer z; such a net reads and writes the grid layers
  // up to the cut layer above z and the low layers, so the interaction
  // distances of the layers above z, the largest ones, do not enter.
  // access_halo_tiles_[z]: pins with access paths up to track layer z
  // (reservation release and dynamic generation of access paths around the
  // pins).
  std::vector<int> area_halo_tiles_;
  std::vector<int> access_halo_tiles_;
  // The same tables for nets without wide shapes (NetData::wide_shapes).
  std::vector<int> narrow_area_halo_tiles_;
  std::vector<int> narrow_access_halo_tiles_;
  // Routed shapes of a wide class inserted or removed by a partition thread
  // for a net that was not marked wide_shapes (expected: none).
  std::atomic<int64_t> unexpected_wide_updates_{0};
  void buildPartitionHalos();
  // Sets NetData::wide_shapes from the net's wire type, pins, reservations
  // and the shapes of its access options (end of setup).
  void markWideNets();
  // NetData::wide_shapes of one net (as markWideNets).
  bool netHasWideShapes(const NetData& net) const;
  // Largest distance at which two shapes of the routing space interact
  // (any layer, incl. via-layer rules), plus one pitch.
  int64_t interaction_reach_ = 0;
  size_t label_limit_ = 0;
  std::atomic<int64_t> inconsistent_updates_{0};
  // Debug group bonnroute_connect: one line per search of a connection
  // (phase, area, settled vertices) and one per partition region.
  bool trace_connections_ = false;
  std::vector<Worker> workers_;
  std::vector<LayerInfo> layers_;  // per z
  std::vector<int> z_of_fr_;       // routing layer -> z, -1
  std::vector<int> cut_z_of_fr_;   // cut layer -> lower z, -1
  frCoord tile_x0_ = 0;
  frCoord tile_y0_ = 0;
  frCoord tile_dx_ = 1;
  frCoord tile_dy_ = 1;
  int tiles_x_ = 1;
  odb::Rect die_;
  odb::dbBlock* block_ = nullptr;
  std::vector<std::vector<uint8_t>> spread_tiles_;
  std::vector<std::vector<uint8_t>> jog_low_tiles_;
  int64_t jog_low_tile_count_ = 0;
  // Global routing extra space (Sec. 2.1) per net, read from odb.
  void loadExtraSpace();
  int64_t extra_space_nets_ = 0;
  int64_t extra_space_tiles_ = 0;
  int64_t spread_tile_count_ = 0;
  int tiles_y_ = 1;

  std::vector<NetData> nets_;
  std::unordered_map<const frNet*, int> net_index_;
  // Same-net rules and tau per track layer z; cleanup layers per drt layer.
  std::vector<br::SameNetLayerRules> same_net_rules_;
  std::vector<br::GeoCoord> tau_;
  std::vector<br::CleanupLayer> cleanup_layers_;
  AccessStats access_stats_;
  // Statistics.
  int64_t ripup_budget_skipped_ = 0;
  double partition_wall_ = 0;
  double partition_busy_ = 0;
  double partition_capacity_ = 0;
};

}  // namespace drt
