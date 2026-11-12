# Off-track wiring and pin access

This note describes how `BrSameNetRules.h`, `BrBlockageGrid.*`,
`BrOffTrackSearch.*`, `BrPinAccess.*`, `BrPathCleanup.*` and their use in
`BonnRouteDR` implement Gester et al., *BonnRoute: Algorithms and Data
Structures for Fast and Good VLSI Routing*, TODAES 18(2), 2013, Sections 3.7,
3.8, 4.3 and the post-processing part of 4.4. Every place where the paper
leaves a choice open is listed with the choice made and the reason.

Sources. The paper text was available. The cited references were not:
Nieberg, *Gridless pin access in detailed routing*, DAC 2011 (paywalled), and
Maßberg and Nieberg, *Rectilinear paths with minimum segment lengths*, Discrete
Appl. Math. 161(12), 2013 (open archive, but not retrievable from the build
hosts). From the latter only the abstract was verified: an extended Hanan grid
with O(n^2) storage, an adapted shortest path algorithm, O(n^4 log n) time and
O(n^4) space. Everything that the paper attributes to these references (the
mapping of same-net rules to minimum segment lengths, the proof of
Theorem 3.2, the details of destructive bounding) is reconstructed here and,
where possible, checked empirically by `//src/drt/test:bonnroute_offtrack_test`.

## Same-net rules to tau (Section 3.7)

`minimumSegmentLength` maps the rules of a wiring layer to one minimum
segment length tau.

- Short edge (LEF `MINSTEP`, also `LEF58_MINSTEP`): tau_step = min_step +
  width. A centre-line segment of length L between two bends has boundary
  edges of length L - w, L or L + w; with L >= min_step + w every boundary
  edge except the line-end edges is at least min_step long, and line-end
  edges are adjacent only to such edges. The test rasterizes random paths,
  traces the boundary polygon and finds no adjacent pair of short edges with
  this tau (and violations in about a fifth of the control paths with
  segments only >= w).
- Notch: tau_notch = notch spacing + width (the legs of a U with base L are
  L - w apart), with the layer's minimum spacing as notch spacing. It is only
  added to tau with `offtrack_notch_tau 1`, because it also forbids short Z
  jogs, which have no notch. By default, catalogue paths containing a U-notch
  are discarded (exact check) and the post-processing collapses notches.
- Minimum area is not mapped to tau. It constrains a connected polygon, not
  each segment, and access paths are attached to a pin or continued on-track.
  The post-processing extends isolated via-to-via pieces instead.
- Without a mapped rule there is no minimum segment length: tau is the
  manufacturing grid (`offtrack_tau` overrides it), and every tau is rounded
  up to the grid. An earlier version used the wire width, which has no rule
  behind it: on Nangate45 (no MINSTEP, no AREA) it forbade the short jog
  from an access point to the nearest metal2 track whenever that track was
  closer than one width, and forced detours over the neighbouring track.
- The Sec. 3.8 alternative ("short edge avoidance ... by introducing
  corresponding edges in G or removing edges") is available as
  `OffTrackRequest::arc_filter` but unused: tau_step already guarantees the
  short-edge rule.

## Algorithm 3 (Section 3.8)

`blockageGridCoordinates` is Algorithm 3 for one axis. The test compares it
with a line-by-line transcription on 4000 random inputs.

- Line 7 reads `x_i + lambda * u`; `u` is read as tau.
- `k in N` includes 0. Then [c_min, c_max] is the maximal run of Hanan lines
  around i with consecutive gaps below 4 tau. Excluding 0 would always extend
  the run across the neighbouring gap, however wide, and break the paper's
  O(n^4) vertex bound.
- The sentinels -inf and inf end the runs and generate no coordinates.
- s and t are Hanan lines as in line 1; for pin access, "s" is the pin, so
  all pin rectangle borders are Hanan lines, and the on-track target points
  and the search window border are added like t.
- Per axis at most 4 m^2 + m coordinates for m distinct Hanan lines; the test
  asserts this bound and reports the growth (`vertex growth` table).

## Hanan grid when tau is at most the manufacturing grid (D9)

Algorithm 3 adds lines at multiples of tau so that some shortest path has
every segment at least tau long (Theorem 3.2). With `offtrack_hanan_grid 1`
(default) `BlockageGridInput::grid_unit` is the manufacturing grid, and
`blockageGridCoordinates` returns the sorted, unique Hanan lines when tau is
at most the grid unit and every Hanan line is on the grid; otherwise it runs
Algorithm 3.

- Why it suffices: obstacles are rounded outward to the grid, and access
  points, targets and the window are on it, so two distinct lines are at
  least one grid unit, hence at least tau, apart. Every path of the Hanan
  grid is tau-feasible, and the Hanan grid contains a shortest
  obstacle-avoiding path for positive axis weights; that path is a shortest
  tau-feasible path.
- One state per vertex: `OffTrackRequest::merge_direction_states` lets the
  search use one state per grid vertex with arcs to the four neighbours
  when there is no arc filter and, on every layer with wires, consecutive
  grid lines inside the window are at least tau apart. A path of that graph
  is a tau-feasible grid path, and every path of G is one of its paths, so
  the distances are those of G. The incoming-direction states serve only
  the tau constraint (the cost model has no bend cost).
- Same costs, other paths: among paths of equal cost the returned one may
  differ; with a target limit, the target vertex settled at a tied limit
  cost may differ. `testHananGrid` compares the Hanan grid with one state
  against the brute force and against Algorithm 3 with five states (cost of
  every target), on lattices of unit 1-3 with tau up to the unit, and with
  tau above it (Algorithm 3 kept); `testSearchMatchesReference` compares the
  merged search with the five-state reference on its 6000 requests.
- `offtrack_hanan_grid 0` restores Algorithm 3 and five states.

## Multi-layer blockage grid

Section 3.8: "we have a blockage grid for each wiring and via layer.
Coordinates are also induced by blockages of neighboring layers. Vertices of
two adjacent wiring layers are connected if a via is allowed here."
`buildBlockageGrid`:

- H_z = borders of the blockages of wiring layers z - 1, z, z + 1 and of the
  via layers below and above z, plus s/t points, pin borders, window.
- X_z = Algorithm 3 (H_z, tau_z) for x and y.
- Via layer grid z|z+1 = X_z u X_{z+1}; wiring layer grid = X_{z-1} u X_z u
  X_{z+1}, so every via-layer vertex is a vertex of both wiring layers.
- A superset of the grid keeps Theorem 3.2 (the endpoints of a shortest
  tau-feasible path stay grid vertices, and G still represents each segment,
  see below). The two-layer brute-force test with via blockages confirms
  equality on 1000 instances.

Via blockages are the Minkowski differences of other nets' shapes (bloated by
spacing) with the bottom pad, cut and top pad of the via model; a via arc
additionally requires the exact legality check of its three shapes.

## Path-preserving digraph G and Theorem 3.2

`OffTrackSearch` runs Dijkstra on G without building it explicitly.

- States per grid vertex: the paper's four incoming directions plus one
  "no incoming direction" state for sources and via landings. The paper
  says "up to four vertices"; the fifth state realizes "each segment has
  length at least tau" for the first segment at s and the first segment after
  a via (a turn arc in every direction). With it every segment, first and
  last included, has length >= tau, and a path may also leave a layer without
  wire (via directly at the pin, stacked vias).
- Straight arcs go to the neighbouring grid vertex in the incoming
  direction; turn arcs go to the nearest grid vertex at distance >= tau in
  each perpendicular direction; via arcs go to the fresh state of the same
  point on the adjacent layer. An arc exists only if its segment avoids the
  open interiors of the obstacles (the paper's tau-feasibility; touching an
  obstacle border is allowed).
- Theorem 3.2 is checked on random instances: 3000 single-layer instances
  (0-8 random rectangles, tau in {1, 2, 3, 4, 5, 7}) and 1000 two-layer
  instances (planar-only or via-only layers, via blockages, via costs 0-5).
  In every instance the cost on G equals a brute-force Dijkstra over (point,
  direction, run length) on the integer grid with a margin of 4 tau + 4
  around the instance, and every returned path is re-verified as
  tau-feasible and obstacle-free with the reported cost.
- Sources are, by default (`offtrack_from_access_points 1`), the pin's
  OpenROAD access points: points, as in the theorem. With
  `offtrack_from_access_points 0` they are regions (all grid vertices inside
  the pin rectangles, pin borders as Hanan lines), the paper's "anywhere on
  the pin"; that extension is not covered by the brute-force test, and it is
  incompatible with FlexDR clean-up (see "Pin contacts" below).
- A pin layer without wires (below the track graph) whose vias are allowed
  only at the access points is a `points_only` layer of the blockage grid:
  its blockages and the via blockages next to it induce no coordinates, and
  it keeps only the given points. The two-layer brute-force test includes
  such instances (vias restricted to the given points) with equal results.
  The restriction matters: with vias allowed anywhere, a down-and-up via
  pair resets the incoming direction and can shorten paths, and the reduced
  grid would miss those.
- Search costs are DBU of wire weighted by the on-track jog factor for the
  non-preferred direction (1 for the preferred direction), plus one wire
  width per via (a tie-breaker). Theorem 3.2 is stated for length; the
  brute-force test also covers axis-weighted instances (about half of the
  layers, weights 1-4) with equal results. Before, both directions cost 1,
  and catalogue paths carried about twice the non-preferred metal2 of the
  access stubs (generated gcd-like rows: 206 vs 109 DBU per path). The on-track cost
  model (jog factor for non-preferred wires, `via_cost` pitches per via, cost
  unit) gives the access option's cost used by the on-track search and the
  selection. With a via cost of four pitches the search settled about three
  times as many states on the generated rows.
- Dijkstra uses an admissible, consistent future cost (l1 distance to the
  lattice spanned by the target coordinates of a layer plus via costs), the
  same idea as the on-track search's pi of Sec. 4.1. It changes no distance;
  targets are settled in increasing distance. The test compares it with plain
  Dijkstra on every Theorem 3.2 instance.
- `max_targets` (2 x `max_access_options`) stops the search once that many
  target vertices are settled: the catalogue consists of the nearest
  on-track points. Some of them may still fail the exact checks.
- Implementation (exact): the result of a run, including the order of the
  paths and every counter, is a function of the request. The heap holds
  (f, state) pairs ordered by f only and sees the pushes and pops of a
  `std::priority_queue` in the order of the formulation above; ties
  therefore resolve as they always did. Per state one 16-byte label
  (distance, 32-bit parent, stamp) and 16-byte heap items; a window of
  2^32 / 5 or more grid vertices is not searched. Per grid line: the nearest
  line at distance >= tau in both directions (turn arcs), the prefix
  counts of unit segments meeting an obstacle (a segment between two lines
  is free iff the counts at its ends are equal; built on first use of the
  line), the l1 distance to the nearest target coordinate of each target
  layer (future cost), and for on-track-only layers whether arcs may start
  there (`OffTrackLayer::tracks`). `BonnRouteOffTrackTest` keeps the
  previous search verbatim and compares results, order and counters on
  6000 random requests (one to three layers, tau 0 to 5, points-only
  layers, source costs, duplicate targets, target and cost limits, via
  rules, arc filters, on-track-only layers).

## Diff-net legality adapter

`br::DiffNetLegality` is the only interface through which the off-track code
asks diff-net questions: `allowed` (exact check of one shape), `blockages`
(Minkowski obstacles for a model) and `conflict` (two candidate shapes).
`BonnRouteDR::Legality` implements it with the shape grid and the distance
rule checker of Sec. 3.3-3.4 (see [DataStructures.md](DataStructures.md)),
which also cover the pin layer and the cut directly below the lowest track
layer (pins on M1 with routing from M2): `allowed` is `shapeConflicts`
(the checker's placed-shape level, the asking net's own shapes ignored),
`conflict` the pairwise distance rule with the same shape classes, and
`blockages` the stored pieces of the window bloated by the largest
requirement of their class pair plus their line-end extensions, minus the
model (conservative: an l-infinity box around the Euclidean region).
Off-track shapes are classified by their geometry without line-ends of
their own.

## Pin access catalogues and circuit classes (Section 4.3)

- Access layer: the lowest layer of the pin's access points among the track
  layers and the routing layer directly below the lowest track layer (with
  `offtrack_from_access_points 0`: the lowest layer with pin shapes). The
  search uses that layer (planar wires only if it is a track layer) and the
  next routing layer. A via at an access point is the access point's own via
  (OpenROAD's choice); pin access checked its pin-layer pad and cut, so only
  its top pad is checked against the routing space, as for the previous
  access stubs. Other vias use the default via with all three shapes.
- Targets: track-graph vertices within `offtrack_radius` pitches (Chebyshev,
  around the pin's bounding box) on these layers that the fast grid allows
  for a preferred wire. Search window: bounding box of pin and targets plus
  2 tau + pitch.
- Catalogue ("several DRC clean paths connecting it to on-track points within
  a small radius"): one shortest tau-feasible path per target, at most
  `max_access_options`, after the exact diff-net check of every shape and a
  same-net notch check.
- Circuit classes ("identify geometrically equal situations on the chip area
  (up to translation, mirroring, and rotation)"): the signature of an
  instance consists of, in the instance's local frame (inverse of its
  placement transform as a lattice isometry), all fixed shapes of the access
  layers and cuts within the instance box expanded by `access_class_margin`
  pitches (default 0: the cell box itself; clipped to it), its access
  points labelled with their relation to the instance's access pins
  (net group or foreign), the track coordinates of the lowest four track
  layers in that area, the via models' shapes, the access pins with their net
  groups, and the radius. Equal signatures map onto each other by
  M = place_member o place_representative^-1. Via models enter in local
  orientation, so an asymmetric via separates orientations that would map
  its pads incorrectly.
- Each class's catalogues and conflict-free solution are computed once on a
  representative against fixed shapes and pins (preprocessing), in parallel.
  Each member maps them; its mapped primary paths are checked against the
  routing space including the reservations of earlier instances (the paper
  is silent on conflicts between neighbouring circuits); if one fails, the
  member's own selection is re-solved on its legal mapped paths. A mapping
  whose endpoint is not a track vertex (possible only through track
  irregularities outside the signature area) falls back to per-instance
  generation. `access_classes 0` gives every instance its own class.
  Since the signature area is the cell box by default, neighbours can
  differ between members: mapped paths are re-checked at each member, and a
  pin without a legal mapped path gets its own catalogue ("member pins
  regenerated" in the log). On the generated rows the classes drop from 480
  to 8 for 480 instances (margin 0) with no mapped path illegal at a member;
  the earlier area (the search reach, about 2 radius pitches + 2 tau + 3
  pitches) gave one class per instance, as on tinyRocket and aes.
- Pins without any catalogue path keep the previous access stubs (counted as
  fallback); top-level pins get a catalogue each and their best legal path.
- Selection conflicts: every candidate shape is classified once (role
  fixed, no line-end ends, as the pairwise legality check classifies it);
  shape pairs farther apart in x or y than the layer's rule reach are not
  checked, since no rule applies beyond it (`BonnRouteOffTrackTest` checks
  that no such pair violates a rule).
- Setup log (verbose): `BonnRoute off-track access phases` gives the wall
  seconds of pin map, signatures, classification, class catalogues (with
  the thread seconds of targets, obstacle collection, blockage grid,
  search, path checks and selection), members and top-level pins;
  `BonnRoute DR access digest` and `routing digest` fingerprint the
  catalogues, primaries and reservations after setup and the routing
  before writeNets, for comparing runs.

## Catalogues per pin class (D10)

With `offtrack_pin_classes 1` the unit of reuse is one pin instead of one
instance with all its pins.

- Key of a pin (instance's local frame, as the circuit signature): master
  and pin index; the pin's access points with their vias and, with
  `offtrack_preferred_access_points 1`, their preferred marks; the shapes of
  the instance itself (pins and obstructions, marked same net as the pin or
  not) and the tracks of the two lowest track layers within the pin's search
  reach (access-point box plus `offtrack_radius` + 3 pitches + 2 tau); the
  via models and the catalogue parameters. Other pins' preferred marks and
  connections, neighbours and power shapes are not in the key.
- Catalogues: each class's catalogue is generated at its first pin (instance
  order) in that pin's own surroundings, in parallel. Then, in parallel per
  instance and before any reservation, every other pin maps the catalogue
  (`map_option`), keeps the paths `optionLegal` at the member, and generates
  its own catalogue when the mapping fails or fewer than
  min(`offtrack_regenerate`, catalogue size) paths remain; the instance's
  conflict-free selection follows. In instance order each instance keeps its
  selection where the primary paths are legal against earlier reservations,
  else re-solves on its legal paths, and reserves.
- Legality equals the circuit-class flow: every kept path is checked at the
  member. Completeness: paths blocked at the representative by its own
  neighbours are missing at members unless they regenerate.
- `offtrack_catalogue_sample n` checks every n-th member pin against a
  catalogue generated at the member: equal (target, cost) sets, sets with
  mapped-only or own-only entries, and pins whose whole surroundings within
  the reach equal the representative's (all shapes, labelled); for those
  related by a translation the sets must be equal.
- `offtrack_access_stats 1` (M3) logs the circuit classes with one signature
  component left out (net groups, preferred marks, shapes not of the
  instance, tracks), the pin-key classes, and the blockage-grid lines and
  vertices inside the windows with the grid in use and with the other grid.

## Conflict-free solution by destructive bounding

The paper names the criteria ("spreading of endpoints, number of blocked
tracks, directions of feasible on-track continuation, and length") and the
method ("a branch-and-bound based enumeration technique called destructive
bounding") without formulas. Implemented in `selectConflictFree`:

- Cost of a choice (one path or none per pin): sum of path costs +
  `access_direction_weight` pitches per missing on-track continuation
  (along the track both ways, via up, read from the fast grid) +
  `access_track_weight` pitches per track in the union of tracks the paths
  block for other nets' wires + `access_spread_weight` per DBU by which two
  endpoints on one layer are closer than two pitches + a dominating penalty
  per pin without a path. Two paths of different nets in conflict make the
  choice infeasible (pairwise check through the adapter).
- Destructive bounding as in destructive improvement (Klein and Scholl,
  EJOR 1999): fix a hypothetical bound F (initially the sum of the pins'
  individual lower bounds); enumerate depth-first, at every node removing
  from all open pins' domains the conflicting paths and those whose lower
  bound exceeds F; an empty domain refutes the node. If F is refuted, F
  becomes the smallest exceeding lower bound. The first complete choice found
  is optimal for the cost model; the test checks this against exhaustive
  enumeration on 3000 instances.
- `access_node_limit` bounds the enumeration; beyond it the greedy choice is
  taken and counted as not optimal.

## Connections (Section 4.4)

- Source and target sets: the endpoints of the catalogue paths of each pin
  (item 2 of Sec. 4.4; item 1 is the on-track search itself).
- Checking: before each connection every catalogue path is checked against
  earlier routed nets through the adapter (all shapes, worker-local commits
  included). With fewer than `offtrack_regenerate` feasible paths the pin's
  catalogue is generated again against the current routing space and the new
  paths are added.
- Reservations: as before, the primary paths are reserved as complete
  geometry and released when the net is connected.

## Same-net post-processing (Section 4.4)

After each on-track path, the connection is assembled as one chain of route
nodes: target catalogue path (pin to track vertex), on-track path and, for
the first connection, the source catalogue path reversed. `cleanupConnection`
processes each layer run between vias; the chain ends and via positions stay
fixed.

1. Repeated points and interior collinear points are removed (spikes
   included: the removed metal is dangling).
2. A segment shorter than tau is removed by shifting a perpendicular
   neighbour by it (the short segment merges into the next parallel segment,
   fix A/B), or, for a U with a short base, by collapsing the shorter leg
   (fix C). A change is kept only if its new metal passes the adapter and it
   reduces the run's short inner segments and notches without adding short
   end segments.
3. A via-to-via piece below the minimum area is extended beyond one of its
   vias along its end segment (or the preferred direction for a bare via
   stack), if legal; the extension is emitted as an extra wire.

The paper does not list the operations; these are the local repairs whose
legality can be checked with the same adapter. Track vertices of the
on-track path become sources of later connections only if they still carry
metal after the post-processing.

Blockage boxes and the search window are grown outward to the
manufacturing grid, so all blockage-grid coordinates are on the grid (the
round-8 runs had Off Grid markers on metal2 only with off-track access).

Pin-layer pad and cut of an access point's via: pin access checked them
against the fixed surroundings, so they are checked here against other
nets' routed shapes and reservations only (`routedConflicts`; neighbouring
access vias), the top pad against everything. Before, they were not
checked at all, and catalogue paths could start at access points whose
pads conflict with a neighbour's chosen access via (metal1 spacing).

Macro pins: fakeram45 pins (metal3, 0.07 x 0.07 at the macro edge) abut the
macro obstruction on their layer, which also covers metal1/metal2 under
them. Every path from such a pin starts inside the obstruction's spacing
halo, so the search finds none ("no tau-feasible path", 210 pins on
tinyRocket). OpenROAD's pin access treats the abutment as intended; these
pins keep the access stubs (`offtrack_macro_pins 0`).

Layers whose wires must stay on tracks in the preferred direction (drt's
`isUnidirectional`: multi-patterned or RECTONLY, plus RIGHTWAYONGRIDONLY;
asap7 M1-M3): the catalogue search removes every planar arc that is not a
preferred-direction arc along a track (Sec. 3.8, "removing edges that induce
a violation"), the post-processing makes no jog shifts or collapses there,
and min-area extensions are preferred-direction only. A wrong-way or
off-track piece there is unroutable for FlexDR: a worker boundary through it
becomes a boundary pin FlexDR cannot reach on that layer (aes_asap7,
`DRT-0255 ... cannot find path`, source "aps" on the box edge on M2). The
final check in `writeNets` leaves any net with such a wire to FlexDR
("wrong-way or off-track wires" in the log). The test compares the filtered
search with a brute force restricted the same way.

Access points: all access points of the pin's access index, those pin
access did not choose starting with a penalty of `pin_cost` pitches (an
initial distance of their sources; `offtrack_preferred_access_points 1`),
as the access stubs price them. An earlier version used the chosen access
point only: on asap7, where M2 is on-track-only, every catalogue path of a
pin then lay on one M2 track, neighbouring primaries on the same track
blocked them, dynamic generation found the same paths again ("5732
generations added 0 paths") and 6596 of 17338 nets failed in BonnRoute
(948 without off-track access). Access points off the manufacturing grid
are never used (they put the via and the grid lines off it).

Every pin keeps its access stubs after its catalogue paths (the catalogue
first, so primary indices are unchanged): a pin whose catalogue paths are
all blocked later still has the stubs' alternatives. Dynamic generation
searches two pitches beyond `offtrack_radius` against the current routing
space. The log counts connections that failed because a pin had no
feasible access option left.

Min-area extensions are written as FlexDR patch wires anchored at their via
(`frPatchWire`, origin at the via, offset box on the manufacturing grid).
Written as path segments, FlexDR's connectivity check removed them as
dangling and reported Recheck markers (asap7: 14,308 on M2 for 16,968
extensions), and their lengths were off the grid (Off Grid on sky130 met1).

Via connection marks: a via pad at an access point of a pin of the net is
marked bottom/top connected, as FlexDR marks its own vias. Unmarked, FlexDR's
connectivity check (`buildPin2epMap`) does not see the pin connected and may
remove the access via as a dangling object (Recheck on li1 on sky130).

## Pin contacts and FlexDR (DRT-1231)

FlexDR's clean-up (`initNets_searchRepair_pin2epMap`, `initNet_term_helper`
in `FlexDR_init.cpp`) makes a pin a term of a worker wherever a truncated
path-segment end or a connected via pad of the net lies on the pin, and then
requires one of the pin's access points (`getPinAccess(inst pin access
index)`) inside the worker's route box; otherwise it stops with
`DRT-1231 Pin ... does not have access point`. The previous access stubs
always touched pins at access points. The producers of pin-touching
geometry and how each keeps contacts at access points:

| Producer | Contact | Handling |
| ----- | ----- | ----- |
| Catalogue paths (preprocessing) | first node on the pin; access via | start at the access points; on a pin layer without wires vias only there (`points_only`) |
| Dynamic regeneration (Sec. 4.4) | same | same code path as the catalogue |
| Circuit-class mapping | mapped first node | mapped access points are the member's access points for equal signatures; the guard catches the rest |
| Fallback stubs | jog start / access via | start at an access point (unchanged) |
| Planar `pin_points` (truncated ends) | truncated wire end | only the first node of a planar option, an access point |
| Post-processing (jog shift, U collapse, min-area extension) | moved bends, extension tips | chain ends and vias never move; no new node or tip on a pin of the net away from its access points (`CleanupRequest::node_allowed`) |
| On-track paths | vias landing on pins of the net | unchanged; covered by the guard |

Final guard in `writeNets`, aligned with FlexDR: a pin becomes a worker
term only at truncated path-segment ends (`pin2epMap`, and
`initNets_segmentTerms` for border segments) and at via pads marked
bottom/top connected. FlexDR sets those marks only on vias its maze creates,
at access points (`FlexDR_maze.cpp`, `checkViaConnectivityToAP` requires
`hasAccessPoint`); vias of the initial routing keep the default false. So a
net is left to FlexDR only if a truncated end (a planar pin point) lies on a
pin of the net away from its access points (`br::pinContactViolations`); via
pads on pins away from access points (e.g. on-track vias on top-level or
macro pins) are only counted in the log. The earlier guard also dropped nets
for such via pads, even with off-track access off.
`offtrack_from_access_points 0` produces such contacts by construction and
is therefore incompatible with FlexDR clean-up; the guard then drops most
nets.

## Parameters and A/B

`offtrack_access 0 postprocess 0` restores the previous access (OpenROAD
access points with short stubs, greedy reservations) and geometry exactly.
`offtrack_access 1` (default) enables everything above;
`postprocess` can be switched independently. See `src/drt/README.md`.
