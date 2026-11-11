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
* **Interaction distances** (the window widths of Sec. 3.3/3.4).
  `distance(a, b)` is the largest distance at which shapes of two classes
  violate a rule (spacing, line-end extension, 1 for an overlap rule).
  Per layer `ShapeClassTable` keeps three bounds: `ruleReach` (= `reach`)
  = the table's largest value (partition halos, shape-grid border, callers
  that classify shapes after choosing a window), `classReach` = the largest
  `distance` between any two classes interned so far, and `narrowReach` =
  the same over narrow classes only. Both grow monotonically inside
  `intern()` (under its lock, atomics for the readers). A class is **wide**
  when its width row's largest spacing is at least `wideRatioPercent()`
  (default 200 %) of the narrowest row's; wide shapes go to a separate index
  (A2), so a macro blockage or wide strap no longer widens every window of
  its layer. `reachOf(c)` = `classReach` for a wide class, `narrowReach`
  otherwise; a window around models extends by the largest `reachOf` of the
  models. Nangate45 metal2-6 have tables up to 1.5 um (3000 DBU) for 1.5 um
  wide shapes; with only wires and pins on a layer the windows extend by
  the 0.07 um narrow spacing instead (`bonnroute_routing_space_bench` with
  `BR_RULES=n45`). The criterion only affects speed: exactness holds for
  every threshold, and the tests draw it at random.

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
* **Wide shapes** (*choice*, see A1): kept whole (clipped to the grid area)
  in per-tile-chunk lists instead of the cells; a shape is listed in every
  chunk it meets and reported once, from the first chunk of the queried
  area. `collect(window)` adds the wide shapes within
  `classReach - narrowReach`
  of the window and of the known rows it returns. A window that extends at
  least `narrowReach` beyond a model therefore holds every shape the model
  can violate, for every location whose window lies in the known area.
  The lists are written by the thread owning the chunk like the trees; a
  reader whose wider area reaches into another region's chunk takes the
  list's lock and only sees shapes outside its own models' reach.

## A3 Distance rule checking (Sec. 3.4)

* `check(l, models, options)` → per answer group the level (or allowed /
  not allowed), the maximal horizontal and vertical interval of locations
  with the same answer, and on request the nets whose removal makes a blocked
  group legal with their ripup levels. `query(z, l, W)` is the paper's
  interface (per wire type: preferred, non-preferred, via up, via down; a via
  model's shapes are one group). `shapeLevel` checks one placed shape with the
  reference pairwise rule on the stored pieces.
* Only `ShapeGrid::collect` is used: per layer it returns the pieces of all
  intervals meeting the window (models bloated by their interaction
  distance `reachOf`, A1), the wide shapes near it (A2) and the area
  **known** from them — per row from the end of the
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
* `FastGridBuilder` keeps the fast grid equal to the checker ("continuously
  updated data computed by the distance rule checking module") with exact
  incremental updates. All of them use the checker's own violation spans
  (`appendViolationSpans`, the Minkowski regions of model and shape on a
  track line, as the old `forEachBlocked` marking but with the paper's
  rules):
  * *Insert* (`shapeInserted`): level(v) = max(old level, level of the new
    shape where it violates v), from the new shape alone; no shape-grid
    query. Edge bits: with S' = max(S, new level if the swept wire is hit),
    M' >= M and S <= M where the bit was clear, a clear bit becomes set
    exactly if the swept wire is hit and the new level exceeds M'; a set bit
    stays set while M' = M and is rechecked exactly otherwise (rare).
  * *Remove* (`shapeRemoved`): only levels equal to the removed shape's level
    where it violates, and marked edges whose swept wire it hits, can change.
    Only these candidates (with the edge bits at both sides and the jog bits
    of the pairs below and above) are recomputed exactly
    (`recomputeCandidates`): one shape-grid window for all of them, pieces
    prefiltered exactly by their interaction distance, unchanged neighbours
    read from the stored words. A group of removals (a net's rip-up, its
    released reservations) is recomputed once (`PendingRemovals`, `flush`).
  * The initial build computes whole tracks the same way (one window per
    track block, no per-interval queries).
  Wire types outside the fast grid and marked edges are answered exactly
  through the same checker (`LegalityOracle`). Searches count both kinds of
  answers (`SearchResult::fast_grid_queries` / `checker_queries`).
* Lookups: the searches read the fast grid through a `FastGridReader`
  (per search object) with run cursors: per track the last run and the run
  array of its tile, direct mapped, so sequential access along a track and
  its neighbouring tracks costs a comparison per lookup and a step to the
  next run at a run border. The searches pass (z, i, j) instead of vertex
  ids. The interval search builds its intervals from the same cursors, so
  consecutive vertices of a fast-grid run are read in O(1) each.
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
  routing space"). *Choice*: an exact overlay instead of physical removal and
  reinsertion. Before a net's searches, `buildOverlay` computes the words of
  every fast-grid position the net's own shapes (pins, active reservations,
  their cut projections) can influence, with those shapes ignored (the
  removal candidates of all own shapes, recomputed once per layer), and the
  searches of the net read the fast grid through it; every checker query of
  the connection (stub and landing-pad checks, marked edges, other wire
  types, off-track legality, victims) uses `ignored_owner` = the net. The
  answers are exactly those with the shapes removed (`BonnRouteFastGridTest`
  compares the overlay with the brute-force reference ignoring the owner),
  but the routing space is not changed: no shape-grid writes, no fast-grid
  restoration afterwards. Committed reservations are then released by the
  caller as before. The log reports the own shapes per connect.
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
  rule interaction distance (`2 · (ruleReach + model + cell) + pitch`,
  rounded up to tiles), so all reads and writes of a net stay in its region;
  the halo does not depend on the classes interned later.
  `BonnRoutePartitionTest` checks equality with a serial run; run it with
  `-fsanitize=thread` as well.
* Search scratch: per-vertex labels live in 1024-label pages allocated when a
  search touches them (`PagedLabels`), so a die-spanning window costs a page
  table, not 29 B per window vertex per worker.

## Merge with the potentials, spreading, critical nets and off-track work

* Input-wire reservations (arXiv:2111.06169) and pin-access reservations are
  both net-owned `kReservation` shapes in the shape grid (`setReservation`,
  `setInputReservation`); both are removed while their net is connected
  and restored by rip-up rollback.
* Potentials and the `potential_check` fallback guard are unchanged; they
  only use the track graph and the search request.
* Wire spreading reads the fast grid: a neighbouring track is occupied where
  its wire vertex is not legal (`PathSearch::spreadCost`, level != 0), the
  core-test reference uses the same expression. Because the net's own shapes
  are removed during its searches, "occupied" means occupied by other nets,
  as before.
* Critical nets first and the empty-partition skip are kept; every partition
  (the critical one on the whole chip by one thread) updates the routing
  space in place, so the sequence semantics are those of the deferred
  version: a thread sees its own commits at once and nothing of other
  regions.
* Low layers: the pin layer and cut directly below the lowest track layer are
  further shape-grid layers (after the track layers, outside the fast grid).
  They replace the previous read-only R-tree of those layers, so pins there
  are owned shapes (removed while connecting) and access reservations there
  block other nets.
* Off-track legality (`BonnRouteDR::Legality`, `optionLegal`,
  `regenerateAccess`, `postprocessConnection`) goes through the shape grid
  and checker. Track-graph usability for a given net
  (`BonnRouteDR::usable`): the fast grid, and where it reports a blockage the
  exact checker with that net's shapes ignored
  (`FastGridBuilder::exactLevel`), since the fast grid is net-agnostic.
* Search scratch: `search_label_mb` (default 512 MB per worker) bounds the
  labels of one search; a search that exceeds it stops without a path
  (`SearchResult::aborted`, counted and reported with the peak scratch) and
  the net goes the usual way of a failed connection. Between searches at most
  64 MB of pages stay allocated per worker.

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
* Windows extend by the interaction distance of the classes that occur, and
  wide shapes live outside the cells (A1, A2); the paper sizes windows by
  the rules. Same answers, smaller windows.
* No inter-layer via rule between the cut below the lowest track layer and the
  first track via layer is modelled.

## Tests and measurements

`src/drt/test`: `bonnroute_wire_model_test`, `bonnroute_shape_grid_test`,
`bonnroute_distance_rules_test`, `bonnroute_fast_grid_test`,
`bonnroute_partition_test`, `bonnroute_core_test`; all randomized against
brute-force references that keep every shape in a vector and apply
`shapesViolate` pairwise with the same classes. `bonnroute_routing_space_bench
[seed] [extent] [searches]` prints shape-grid and fast-grid statistics,
memory, the fraction of path-search legality queries answered by the fast
grid and query / update times against an R-tree with the same rule.
