# BonnRoute implementation scope and paper fidelity

BonnRoute GR (`global_route -bonnroute`, inside CUGR) and DR
(`detailed_route -bonnroute`, before FlexDR clean-up) implement Gester et al.,
*BonnRoute: Algorithms and Data Structures for Fast and Good VLSI Routing*,
TODAES 18(2), 2013. `set_bonnroute_profile` selects the algorithm variant:

| Profile | Purpose |
| ----- | ----- |
| `2013-reimpl` | The 2013 paper as written, as far as implemented (this document lists every remaining difference). |
| `2013-improved` | Default. The 2013 algorithms with measured optimizations that change search order or scheduling (vertex search, speculative tail batches, l1 GR bound, deterministic price batches, unresolved nets to FlexDR). |
| `2021` | Ahrens et al., arXiv:2111.06169, bulk routing (Table 3): `2013-reimpl` with the general-model potential for every search. See [2021 TODO](#2021-todo). |

Both 2013 profiles share every data structure and every exact speed-up
(cache-line alignment, per-class reach windows, tile-grouped fast-grid
updates, per-batch reservation flushes, exact pi_H). Exactness of those
speed-ups is checked by the differential tests in `src/drt/test` and
`src/grt/test`.

## 2013 paper, section by section

| Section | Status in `2013-reimpl` |
| ----- | ----- |
| 2.1 GR graph, pins as vertex sets | Implemented on OpenROAD GCells; pins connect through zero-cost cliques over the cells of their access points. |
| 2.1–2.2 resources, extra space | Edge capacities, the wirelength/via objective, extra space s(n,e) on a grid of 0.5 tracks with convex power (switched wire capacitance with coupling) and yield (critical area) resources; the oracle takes the optimum s per edge (`extra_space`, see `src/grt/README.md`). DR wire spreading uses the extra space per net (`spread_from_extra_space`). |
| 2.2 oracle (Algorithm 1) | Path composition with goal-oriented Dijkstra (l1 and landmark bounds). |
| 2.3 resource sharing (Algorithm 2) | t = 125, epsilon = 1, reuse rule, concurrent price updates. |
| 2.4 rounding, rip-up and reroute | Randomized rounding, re-choose, CUGR maze clean-up. |
| 2.5 capacity estimation | Usable track vertices, single-tile net prerouting (track-graph shortest paths), intra-tile Steiner estimates (rectilinear spanning trees), statistical stacked-via correction. |
| 3.1–3.2 rules, wire and via models, wire types | Spacing tables, line-end model as in the paper, one wire type per NDR; see [`src/br/DataStructures.md`](src/br/DataStructures.md) for the treatment of each drt rule. |
| 3.3 shape grid | Cells, configuration table, interval trees per row and GCell chunk, owners and rip-up levels. |
| 3.4 distance rule checking | Allowed models, maximal intervals, blocking-net report. |
| 3.5 track graph | OpenROAD tracks. Track optimization: [TODO](#2013-todo). |
| 3.6 fast grid | Intervals of 12-bit words per wire type, eight rip-up levels, edge query bits, via-layer words. |
| 3.7–3.8 same-net rules, blockage grid | Minimum segment length tau from the tech, Algorithm 3, path-preserving digraph ([`src/br/OffTrack.md`](src/br/OffTrack.md)). |
| 4.1–4.2 on-track search | Interval search (Algorithm 4 tie order, property (ii) with vias) with pi_H, and pi_P (general potential) for connections whose guide path detours; wire spreading costs from global routing congestion (vertices with a spreading cost are single-vertex intervals). |
| 4.2 rip-up mode, history | Rip-up penalties growing with history at a vertex. |
| 4.3 off-track pin access | Catalogues per circuit class, conflict-free primary paths by destructive bounding, reservations, dynamic regeneration. Paths start at OpenROAD access points, which FlexDR clean-up requires (`offtrack_from_access_points 0` starts anywhere on the pin); macro pins keep access stubs. |
| 4.4 connecting nets | Guide-path connection regions, temporary removal of own shapes, same-net post-processing, rip-up sequences with rollback, later passes with larger areas. |
| 5.1 parallelization | Balanced partition sequence with fewer regions later; critical nets first: clock nets and, with `route_ndr_nets 1`, NDR (wide-wire) nets. |
| 5.2 external clean-up | FlexDR cleans up; logs report the nets and violations it receives. |

## 2013 TODO

- **Wire types in DR (Sec. 3.2, 3.6):** NDR wire extensions and minimum cut counts of NDR vias are not modelled (`route_ndr_nets 1` routes NDR nets with one wire type per NDR).
- **Tile size (Sec. 2.1):** tiles for 50–100 wires per layer; `gcell_scale` is an experiment, and DR on coarse guides is slow.
- **Track optimization (Sec. 3.5, Müller 2009):** maximal usable track length.

## 2021 TODO

The `2021` profile sets `2013-reimpl` plus the general-model potential
(`potential 3`) for every search, as in the bulk routing experiment (Table 3).
Reservations (`reservations 1`, discount 3/4) are the tool of the paper's
incremental routing (Sec. 5.2, Tables 4-6), whose input is a legal detailed
routing. Fed with track assignment wires in bulk routing, they block other
nets and send more nets to the serial rip-up phase; the paper suggests
reserving only higher layers for such input (`reserve_min_layer`). Open items
of arXiv:2111.06169:

- Discount rates in the interval search (reservations run in the vertex search).
- Via reservations; dropping short wires that only connect to illegal input wires.
- The incremental scenario: keep pin-to-pin paths of legal detailed wires, replace the rest by global wires, re-route with reservations.
- Region-dependent costs (layer-range penalties).
- Theorem 7 and Theorem 14 query structures (the paper uses the Sec. 5.1 variant in practice, as implemented).

## Measurements

Results are generated from run logs; see the PR description and the E2E
reports. Results with BonnRoute GR or DR describe the combined flow with
FlexDR clean-up; the logs report how many nets and violations FlexDR
receives.
