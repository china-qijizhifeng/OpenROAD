# BonnRoute implementation scope and paper fidelity

This adapter implements parts of Gester et al., *BonnRoute: Algorithms and Data
Structures for Fast and Good VLSI Routing*, TODAES 18(2), 2013. It is an
experimental GR and initial-DR implementation combined with CUGR and FlexDR.
Neither its default configuration nor its explicit interval/rip-up mode is a
complete reproduction of the paper. The measurements section refers to code at
`450f22fcc924cbc4084f45596fbf30762b1db207`; the tables describe the current
source.

## Models of connectivity, capacity and geometry

| Area | Paper | Current adapter and possible consequence |
| ----- | ----- | ----- |
| GR pin connectivity | Branches may reach different cells of the same pin and use the pin's existing conductivity to connect. | Multiple access cells are candidates, but the output chooses one cell per pin. It cannot represent arbitrary same-pin cell forests; this can require extra routing. |
| GR capacity estimation (Section 2.5) | Count usable track-graph vertices, account for prerouted single-tile nets and estimated intra-tile connections of longer nets, and statistically estimate intermediate-layer capacity loss from stacked vias. | Import CUGR edge capacity minus existing demand. Via-pad demand is accounted for, but this is not the paper's capacity estimator. The difference can affect congestion prediction and detours. |
| DR geometry (Sections 3.1–3.6) | Shape grid of cell configurations in interval trees, distance rule checker with maximal intervals and blocking-net report, fast grid with edge query bits. | Implemented as described in [`src/br/DataStructures.md`](src/br/DataStructures.md), including the per-rule modelling (some drt rules are modelled conservatively, optimistically or not at all and are left to FlexDR). |
| DR off-track geometry (Section 3.8) | Blockage grid and path-preserving digraph for off-track paths. | Absent; OpenROAD access points and short local access stubs, checked by the distance rule checker. |

Source entry points are [`BonnRouteGlobal::oracle` and `run`](../grt/src/cugr/src/BonnRoute.cpp),
and [`BonnRouteDR::changeShape`, `shapeConflicts` and `buildAccess`](src/br/BonnRouteDR.cpp).
Multiple pin access candidates and swept-edge legality checks are implemented;
their limitations should not be described as their complete absence.

These are verified model differences. Their individual effects on runtime,
wirelength, vias and final routing quality have not been isolated by controlled
experiments. They cannot yet be assigned shares of the measured regression.

## DR execution

| Area | Paper | Adapter |
| ----- | ----- | ----- |
| Interval search (Section 4.1, Algorithm 4) | Among intervals active at the same reduced distance, choose by maximum future cost; relax each neighbouring interval once from the member of J_I(delta) with maximum future cost. | Equal-key events are ordered by the maximum future cost of their flat run; label values come from prefix/suffix minima in O(log labels); jogs to a neighbouring interval are relaxed once from the representative. Property (ii) includes vias (`interval_via_property 1`), so via-heterogeneous runs split; with `0` vias are relaxed per settled vertex. Only the l1 future cost and unit wire costs are supported. |
| Pin access reservations (Sections 4.3–4.4) | Conflict-free primary access paths reserved before routing; removed temporarily for the connecting components; released once connected. Paths can be generated on the fly. | The first conflict-free option of each pin is reserved as complete geometry (stub wires and via) with exact per-owner counts and released when the net is committed. Options are chosen per pin greedily, not by the paper's per-circuit-class destructive bounding, and no paths are generated on the fly. |
| Connection regions (Section 4.4) | Route along a shortest global-routing path from the source component to another component; the area is the tiles of that path plus adjacent layers. | `component_regions 1`: breadth-first tile path over the guides from the current tree to the nearest unconnected terminal; the search runs in those tiles plus `layer_dilation` layers and falls back to all guide tiles if it fails. |
| Rip-up sequences (Section 4.4) | Bounded depth; restore the original situation if the sequence fails; later passes with more effort and larger areas. | Recursive sequence bounded by `ripup_depth` with a journal of every touched net; failure restores all of them. `ripup_passes` later passes enlarge the area by one tile each. Sequences are serial. |
| Parallel routing (Section 5.1) | Sequence of balanced partitions; threads change routing only inside their region; fewer regions later; critical nets first. | `partitions 1`: weighted-quantile column/row partitions of pending nets' boxes, alternating shifted cuts, region count halved every two partitions down to one. A net is routed in a partition only if its box (with halo) lies inside one region. There is no separate critical-net step. |

The default leaves unresolved nets to FlexDR (`ripup 0`) and uses the vertex
search (`use_intervals 0`). Section 5.2 of the 2013 paper uses an external
router for DRC cleanup, so the combination itself is compatible; the counts
of nets left to FlexDR are reported per run.

## arXiv:2111.06169 (reservations and potentials)

`reservations 1` implements Section 5.2 of Ahrens et al.: preferred-direction
track-assignment wires on track whose pieces (cut at GCell boundaries) conflict
neither with fixed shapes, pin reservations nor other nets' input wires become
reservations of their net. They block other nets and cost
`discount_numerator / discount_denominator` of normal wire cost for the owner;
costs use an integral unit of `discount_denominator`. The future cost hook
(`FutureCostOracle`) accepts the potentials of Sections 2–4; with `potential 0`
the l1 bound uses the minimum of normal and discounted rates and stays
feasible. Paper statements about multi-labeling, path post-processing and
backup rip-up searches are not modelled.

## Future cost and external cleanup

Section 4.1 explicitly describes DR's Dijkstra-based search with a future cost
similar to A*. A future cost in DR therefore does not by itself depart from the
paper. The adapter has a geometric future-cost lower bound, but not the full
blockage-aware future cost. A proposed A* change to the GR oracle is a separate
variant whose fidelity and guarantees require their own justification.

Section 5.2 explicitly combines BonnRoute with an external industrial router
for DRC cleanup. Using FlexDR is compatible with that integration pattern.
The amount of connection and rerouting work delegated to FlexDR still matters:
final DRC=0 establishes the combined flow's final result, not BonnRoute-alone
legality, completeness or paper fidelity.

## What the measurements establish

At the measured code revision, switching DR defaults to vertex search with
serial rip-up disabled reduces total PA+GR+DR from 231.842 to 108.356 seconds
on tinyRocket and from 188.863 to 151.465 seconds on AES Nangate45. More
unresolved nets go to FlexDR. These are configuration and combined-flow gains;
they do not demonstrate an optimization of a faithful reproduction.

The current flow remains 50.3% and 121.0% slower than original upstream
OpenROAD on those designs. The upstream comparison includes other fork
changes, so the entire difference cannot be attributed solely to BonnRoute.
Current AES CPU sampling identifies the GR oracle call stack as 92.59% of
user-mode samples for the profiled PA+GR process. This locates work in the
current implementation; it does not prove that missing paper features caused
the slowdown or that the paper algorithm itself is slower.

Implementation optimizations should preserve the selected algorithm's search
order, cost model and update semantics. GR A*, changed reuse rules, altered
price batches and adaptive stopping require separate algorithm-variant
comparisons. A fidelity experiment must also record the unresolved differences
above rather than identify a parameter switch as a complete paper mode.
