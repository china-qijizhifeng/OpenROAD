# BonnRoute on-track routing space: data structures

This directory implements the routing-space data structures of Gester et al.,
*BonnRoute: Algorithms and Data Structures for Fast and Good VLSI Routing*,
TODAES 18(2), 2013, Sections 3.1–3.6, and connects them to the on-track
searches (Sec. 4.1, 4.2) and to `BonnRouteDR`. Every module except
`BonnRouteDR` is free of drt dependencies and tested standalone against a
brute-force reference.

| Module | Paper | Contents |
| ----- | ----- | ----- |
| `BrGeometry.h` | – | `Rect`, `Point`, `AppendOnlyTable` (stable, append-only storage) |
| `BrWireModel.{h,cpp}` | 3.1, 3.2 | `SpacingTable`, `LineEndRule`, `LayerRules`, `ShapeClass`, `ShapeClassTable`, `Requirement`, `shapesViolate` (pairwise reference rule), `WireModel`, `ViaModel`, `WireType`, `makeWireType` |
| `BrShapeGrid.{h,cpp}` | 3.3 | `ShapeType`, ripup levels, `CellShape`, `CellConfigurationTable`, `ShapeGrid` (cells, AVL interval trees per row and tile chunk, owners, shared cells), `maximalCellSize` |
| `BrDistanceRules.{h,cpp}` | 3.4 | `DistanceRuleChecker`: `check`, `query(l, W)`, `shapeLevel`; maximal intervals; blocking-net report |
| `BrTrackGraph.h` | 3.5 | `TrackLayer`, `TrackGraph`, shape `Kind`s of the search |
| `BrFastGrid.{h,cpp}` | 3.6 | `FastGrid` (64-bit words per interval, edge query bits, via-layer words), `FastGridBuilder` (keeps it equal to the checker and answers the rest exactly) |
| `BrCore.h`, `BrCore.cpp`, `BrInterval.cpp` | 4.1, 4.2 | vertex and interval searches reading the fast grid; `PagedLabels` search scratch |
| `BonnRouteDR.{h,cpp}` | 3–5 | rules, classes and wire types from the drt tech; shape grid, checker and fast grid in routing |

Shape-grid layer numbering: wiring layer `z` of the track graph is layer
`2z`, the via layer between `z` and `z + 1` is `2z + 1`.

## A1 Wire models, via models, wire types, shape classes (Sec. 3.1, 3.2)

* `WireModel` = rectangle relative to the stick figure + `ClassId`; the metal
  of a wire is the Minkowski sum of its stick figure and the rectangle.
  `ViaModel` = bottom pad, cut, top pad and, when the via layer above has an
  inter-layer via rule, the projection of the cut onto that layer.
  `WireType` = per wiring layer (preferred, non-preferred) wire model and per
  via layer a via model. `BonnRouteDR::buildRoutingSpace` builds wire type 0
  from the default widths and default vias and one wire type per NDR (NDR
  width, NDR preferred via).
* `ShapeClass` = (layer, kind {metal, cut, cut projection}, width row, run
  column in x, run column in y, line-end flag). Two shapes of one class have
  the same requirements to every other class.
* `SpacingTable`: widths × widths × run-length, "a nondecreasing function of
  their widths and common run-length". The drt adapter samples
  `frLayer::getMinSpacingValue` at every threshold `t` and `t + 1` of the PRL
  table, the two-widths table and the range rules, so bucket borders match
  drt's lookup for both strict and non-strict comparisons
  (`BonnRouteWireModelTest` checks this against random rule functions).

### Rule treatment (Sec. 3.1 "conservative or optimistic modeling")

| drt rule | Treatment |
| ----- | ----- |
| Spacing PRL table, two-widths table, spacing range | Sampled exactly at their thresholds; a table that is not nondecreasing is replaced by its monotone closure (conservative). |
| Common run-length | Bounded by the smaller of the two shapes' own extents (class run columns): conservative, since the true overlap is at most that. |
| Run-length of a preferred-direction wire along its direction | Unbounded (`kUnboundedRun`): wires may be extended later ("exact run-lengths are often not known during the construction of a route"); conservative. |
| Corner-to-corner spacing | Euclidean distance with the run-length-0 column. |
| End-of-line (drt's `drEolSpacingConstraint`: width, space, within) | Line-end extension rectangles (Fig. 2) at both ends of every preferred-direction wire and via pad whose end edge is shorter than the EOL width: pessimistic as in the paper ("always assume wires to have line-ends in preferred direction"). Jogs carry no extension (optimistic, as in the paper). Pins and blockages carry none (optimistic). Parallel-edge and LEF58 EOL variants are not modelled. |
| Cut spacing | Largest cut spacing of the layer, edge to edge (a centre-to-centre rule is thereby modelled conservatively). Adjacent-cuts, parallel-overlap, cut-area and LEF58 cut spacing tables are not modelled. |
| Inter-layer cut spacing (`SPACING ... LAYER`) | Projection of the cut onto the next higher via layer (Sec. 3.2); spacing = larger of the rules in both directions between the two cut layers; projections do not interact with each other. |
| NDR spacing | Not modelled (NDR wire types carry widths and vias only); NDR nets remain skipped by BonnRouteDR as before. |
| Same-net rules (notch, short edge, min area, min step) | Not in the shape grid (diff-net only); Workstream B. |
| Corner spacing, spacing-table influence, width tables, LEF58 metal spacing tables | Not modelled; left to the FlexDR clean-up (the paper also cleans up with an external router, Sec. 5.3). |

Point models versus swept wires: the answer at a location is that of the
zero-dimensional stick figure (the model rectangle with its line-end
extensions). A wire on an edge is the swept shape (model ⊕ segment, the class
of the model, line-ends at its two ends), which is how a committed wire is
stored. Every violation of a swept wire is a violation of a point model on the
segment, so constant point levels on an edge imply its level.

## A2 Shape grid (Sec. 3.3)

* Cells: per layer a tensor grid of equal cells with one neighbour per side.
  **Cell size** (paper: "small enough such that shapes of different nets
  cannot be legally present in the same cell"): the largest `d` with
  `2 d² ≤ s²` for the smallest diff-net spacing `s` of the layer (for a via
  layer with projections also the projection rule and the cut spacing of the
  layer below) — two points of one half-open `d × d` cell are closer than
  `s`. BonnRouteDR reduces `d` to a divisor of the track pitch when one is at
  least `d / 2` and anchors the grid so that tracks lie at cell centres;
  equal geometry on different tracks then yields equal configurations.
* Cell configurations: the multiset of shape intersections relative to the
  cell's anchor (its lower-left corner) with class, type and the line-end ends
  of the original shape inside the piece. `CellConfigurationTable` hash-conses
  them (configuration 0 = empty, never stored in a tree). The table is
  append-only (configurations are never reclaimed); writers memoise
  transitions (configuration, piece, insert/remove) → configuration per thread.
* Intervals: maximal runs of cells with equal (configuration, owner) in
  preferred direction, in an **AVL tree per row** keyed by the first cell.
  *Choice*: rows are split at global-routing-tile borders into chunks with
  their own AVL tree, node pool and shared-cell table (see A5 concurrency);
  an interval never crosses a chunk. Insert / remove of a shape is
  O((k + 1) log n) per row chunk, k = intervals it meets.
* Owner: "for each nonempty interval we store the net that the shapes of this
  interval belong to". Blockage pieces have no owner; pins, reservations,
  wires and vias are owned. A cell that would hold owned shapes of two nets
  (only possible for inputs that already violate rules, e.g. overlapping pins
  of different nets) becomes a one-cell *shared cell* whose parts are stored
  per owner; the statistics report their number.
* **Ripup levels** (paper: "type and net identifier of a shape are used to
  assign a ripup level"): 0 legal, 7 non-removable (blockage, pin, other nets'
  reservation), wires and vias `1 + net priority` clamped to 1..6
  (`setNetRipupLevel`, default 1). The fast grid stores "the lowest level at
  which the shape becomes legal" = maximum level of the violating shapes.
* Removal is atomic (a dry run checks every piece first).

## A3 Distance rule checking (Sec. 3.4)

* `check(l, models, options)` → per answer group the level (or allowed /
  not allowed), the maximal horizontal and vertical interval of locations
  with the same answer, and on request the nets whose removal makes a blocked
  group legal with their ripup levels. `query(z, l, W)` is the paper's
  interface (per wire type: preferred, non-preferred, via up, via down; a via
  model's shapes are one group). `shapeLevel` checks one placed shape with the
  reference pairwise rule on the stored pieces.
* Only `ShapeGrid::collect` is used: per layer it returns the pieces of all
  intervals meeting the window (models bloated by the layer's interaction
  distance) and the area **known** from them — per row from the end of the
  preceding interval to the start of the following one, limited to the tile
  chunks the window meets and to 128 cells (`kKnownHorizon`) beyond it.
  *Choice of maximality*: an interval ends where the answer changes or where a
  location's window would leave the known area; result flags tell which. The
  answer as a function of the location on a line is computed exactly from the
  violation spans of all (model, piece) pairs (Minkowski differences of the
  spacing regions and line-end extensions with the line).
* Own shapes: `ignored_owner` treats a net's owned shapes as removed (Sec.
  4.4), used for off-track checks; the fast grid itself is net-agnostic.

## A4 Fast grid (Sec. 3.6)

* Wiring-layer word: 12 bits per fast wire type (4 kinds × 3-bit level),
  five wire types in bits 0–59; BonnRouteDR keeps wire type 0 (the default
  rule) in it. Via-layer word (intervals along the lower layer's tracks over
  its via positions): 6 bits per wire type (cut, projection).
* *Choice for the edge bit*: bit 60 marks the preferred-direction edge to the
  next position and bit 61 the jog edge to the next track, both at the lower
  endpoint (the paper places "an extra bit at one of the vertices"). A bit is
  set **exactly** when, for some fast wire type, the swept wire's level exceeds
  the maximum of its endpoints' levels; the search then asks the checker for
  that edge. An edge's level is max(endpoints, swept wire).
* Storage: per global routing tile the runs of its tracks (`uint16` start and
  64-bit word per run; a track without runs is all legal), allocated on the
  first non-zero word. Tiles are written by one thread at a time.
* `FastGridBuilder` recomputes, after every shape-grid change, the tracks,
  jog pairs and edges the changed shape can affect, from the checker's
  maximal along-track intervals (one query per interval, not per vertex).
  Wire types outside the fast grid and marked edges are answered exactly
  through the same checker (`LegalityOracle`). Searches count both kinds of
  answers (`SearchResult::fast_grid_queries` / `checker_queries`).
* Rip-up history (Sec. 4.2 cost increase) is kept sparsely beside the grid.

## A5 Integration (BonnRouteDR)

* `buildRoutingSpace`: layer rules from the drt tech, classes, wire types,
  shape-grid geometry, checker, fast grid. `buildFixed` inserts blockages and
  pins (bulk), then the fast grid is built in parallel (vertex data, then jog
  bits) and compacted.
* Reservations (pin access primary paths, input wires) are shape-grid shapes
  of type reservation owned by their net: non-removable for other nets,
  released on connect, restored by rip-up rollback (`setReservation`,
  `restoreRouting`).
* Own shapes (Sec. 4.4: "The shapes of the components that contributed to S
  and T, including reserved pin access paths, are removed temporarily from
  routing space"): `connect` removes the net's pins and active reservations
  before its searches and reinserts them afterwards; committed reservations
  are then released by the caller as before.
* Rip-up victims: `collectVictims` asks the checker for the nets blocking each
  new shape (and its cut projection) with `ignored_owner` = the net.
* Off-track checks keep the adapter `shapeConflicts(shape, net, worker)`
  (Workstream B calls it): checker `shapeLevel` of the classified shape and its
  projection, ignoring the net's own shapes.
* Concurrency (Sec. 5.1): threads route nets whose tile boxes lie in disjoint
  regions. Shape-grid chunks and fast-grid tiles are per global routing tile
  and lazily allocated per tile, so every tree, node pool, run list and tile
  is written by the one thread whose region contains it; updates are applied
  in place, no deferral or merge is needed, and results do not depend on the
  thread count (configuration ids may differ between runs; they are content
  addressed and nothing depends on their values). `netTileBox` adds
  `halo_tiles_` tiles beyond the largest routing area, computed from the
  interaction distance (`2 · (reach + model + cell) + pitch`, rounded up to
  tiles), so all reads and writes of a net stay in its region.
  `BonnRoutePartitionTest` checks equality with a serial run; run it with
  `-fsanitize=thread` as well.
* Search scratch: per-vertex labels live in 1024-label pages allocated when a
  search touches them (`PagedLabels`), so a die-spanning window costs a page
  table, not 29 B per window vertex per worker.

## Deviations and their reasons

* Row chunks per global routing tile (not one AVL tree per full row): needed
  for race-free in-place updates by partition threads; complexity per
  operation unchanged up to the chunking.
* Two edge bits instead of one: they fit into the four spare bits and name
  the edge, so a set bit costs one checker query instead of one per incident
  edge.
* Fast-grid wire types: only the default rule is fast; NDR nets are not routed
  by BonnRouteDR (unchanged from before).
* Interval knowledge horizon of 128 cells and chunk borders bound interval
  lengths; queries stay O(window) instead of following empty rows arbitrarily
  far.
* Shared cells exist only for inputs that already violate diff-net rules.

## Tests and measurements

`src/drt/test`: `bonnroute_wire_model_test`, `bonnroute_shape_grid_test`,
`bonnroute_distance_rules_test`, `bonnroute_fast_grid_test`,
`bonnroute_partition_test`, `bonnroute_core_test`; all randomized against
brute-force references that keep every shape in a vector and apply
`shapesViolate` pairwise with the same classes. `bonnroute_routing_space_bench
[seed] [extent] [searches]` prints shape-grid and fast-grid statistics,
memory, the fraction of path-search legality queries answered by the fast
grid and query / update times against an R-tree with the same rule.
