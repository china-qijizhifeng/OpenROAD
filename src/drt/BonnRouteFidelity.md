# BonnRoute implementation scope and paper fidelity

This adapter implements parts of Gester et al., *BonnRoute: Algorithms and Data
Structures for Fast and Good VLSI Routing*, TODAES 18(2), 2013. It is an
experimental GR and initial-DR implementation combined with CUGR and FlexDR.
Neither its default configuration nor its explicit interval/rip-up mode is a
complete reproduction of the paper. The comparisons below describe routing
code measured at `450f22fcc924cbc4084f45596fbf30762b1db207`.

## Models of connectivity, capacity and geometry

| Area | Paper | Current adapter and possible consequence |
| ----- | ----- | ----- |
| GR pin connectivity | Branches may reach different cells of the same pin and use the pin's existing conductivity to connect. | Multiple access cells are candidates, but the output chooses one cell per pin. It cannot represent arbitrary same-pin cell forests; this can require extra routing. |
| GR capacity estimation (Section 2.5) | Count usable track-graph vertices, account for prerouted single-tile nets and estimated intra-tile connections of longer nets, and statistically estimate intermediate-layer capacity loss from stacked vias. | Import CUGR edge capacity minus existing demand. Via-pad demand is accounted for, but this is not the paper's capacity estimator. The difference can affect congestion prediction and detours. |
| DR geometry (Sections 3.4 and 3.8) | Query permitted wire/via models and maximal usable intervals, handle rule exceptions, and use a blockage graph for off-track paths. | Check wire/jog edges and endpoints, via geometry and spacing; use OpenROAD access points and short local access stubs. The full shape-grid exception checker and off-track blockage-grid search are absent. Some feasible connections or rule interactions can fall outside this model. |

Source entry points are [`BonnRouteGlobal::oracle` and `run`](../grt/src/cugr/src/BonnRoute.cpp),
and [`BonnRouteDR::markShape`, `shapeConflicts` and `buildAccess`](src/br/BonnRouteDR.cpp).
Multiple pin access candidates and swept-edge legality checks are implemented;
their limitations should not be described as their complete absence.

These are verified model differences. Their individual effects on runtime,
wirelength, vias and final routing quality have not been isolated by controlled
experiments. They cannot yet be assigned shares of the measured regression.

## DR execution differs as well

| Area | Paper | Current adapter |
| ----- | ----- | ----- |
| Interval search (Section 4.1, Algorithm 4) | Among intervals active at the same reduced distance, choose by maximum future cost; propagate using the maximum-future-cost representative specified by the algorithm. | Label envelopes and dominance removal exist, but equal-key events use interval ID. Processing still visits settled vertices, including individual via relaxations; the representative shortcut applies only to restricted jog cases. Label values are obtained by scanning labels. |
| Pin access reservations (Sections 4.3–4.4) | Select a conflict-free set of complete primary paths from access catalogues and reserve them; reservations of participating components are temporarily removed, then released when connected. New paths can be generated as needed. | Reserve the top pad of the first option's access vias in the fixed index. This does not implement the complete path-reservation lifecycle or dynamic off-track access generation. |
| Connection regions (Section 4.4) | Select a path between components along the global routing and construct the search region from that path's tiles. | Grow a tree starting at terminal 0, searching the union of all net guides and access-point tiles, with dilation. |
| Rip-up sequences (Section 4.4) | Bound sequence depth and restore the original state if the sequence fails. | Use a serial victim queue with per-net rip-up limits. Commit a new route after removing victims; there is no sequence-wide rollback if victim rerouting fails. |
| Parallel routing (Section 5.1) | Use a sequence of balanced spatial partitions, change the partitions, and reduce region/thread counts later to permit longer connections. Route critical nets first. | Batch nets with disjoint full-net bounding boxes plus halos, search in parallel, then commit in net order. The rip-up queue uses a single worker. This does not implement the paper's partition sequence. |

Source entry points are [`IntervalSearch::Event`](src/br/BrCore.h),
[`IntervalSearch::value` and `process`](src/br/BrInterval.cpp), and
[`BonnRouteDR::reserveAccess`, `buildRegion`, `connectInRegion`, `routeParallel`
and `routeRipupQueue`](src/br/BonnRouteDR.cpp).

The default `use_intervals=0`, `ripup=0` bypasses two paper mechanisms. Setting
both parameters to `1` selects the existing implementations but does not close
the differences above. Differential tests against Dijkstra establish results
on the tested core graphs; they do not establish the paper's interval-search
complexity or full geometric completeness.

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
