# Detailed Routing

The Detailed Routing (`drt`) module in OpenROAD is based on the open-source
detailed router, TritonRoute. TritonRoute consists of several main 
building blocks, including pin access analysis, track assignment,
initial detailed routing,  search and repair, and a DRC engine.
The initial development of the
[router](https://vlsicad.ucsd.edu/Publications/Conferences/363/c363.pdf)
is inspired by the [ISPD-2018 initial detailed routing
contest](http://www.ispd.cc/contests/18/).  However, the current framework
differs and is built from scratch, aiming for an industrial-oriented scalable
and flexible flow.

TritonRoute provides industry-standard LEF/DEF interface with
support of [ISPD-2018](http://www.ispd.cc/contests/18/) and
[ISPD-2019](http://www.ispd.cc/contests/19/) contest-compatible route
guide format.

## Commands

```{note}
- Parameters in square brackets `[-param param]` are optional.
- Parameters without square brackets `-param2 param2` are required.
```

### Detailed Route

This command performs detailed routing.

Developer arguments
- `-or_seed`, `-or_k`

Distributed arguments
- `-distributed` , `-remote_host`, `-remote_port`, `-shared_volume`, `-cloud_size`

Pin Access arguments
- `-via_in_pin_bottom_layer`, `-via_in_pin_top_layer`, `-via_access_layer`, `-min_access_points`

```tcl
detailed_route 
    [-output_maze filename]
    [-output_drc filename]
    [-output_guide_coverage filename]
    [-drc_report_iter_step step]
    [-db_process_node name]
    [-disable_via_gen]
    [-droute_end_iter iter]
    [-via_in_pin_bottom_layer layer]
    [-via_in_pin_top_layer layer]
    [-via_access_layer layer]
    [-or_seed seed]
    [-or_k k]
    [-bottom_routing_layer layer]
    [-top_routing_layer layer]
    [-verbose level]
    [-distributed]
    [-remote_host rhost]
    [-remote_port rport]
    [-shared_volume vol]
    [-cloud_size sz]
    [-clean_patches]
    [-no_pin_access]
    [-bonnroute]
    [-min_access_points count]
    [-save_guide_updates]
    [-repair_pdn_vias layer]
    [-single_step_dr]
```

#### Options

| Switch Name | Description |
| ----- | ----- |
| `-output_maze` | Path to output maze log file (e.g. `output_maze.log`). |
| `-output_drc` | Path to output DRC report file (e.g. `output_drc.rpt`). |
| `-output_guide_coverage` | Path to output guide coverage file (e.g. `sample_coverage.csv`). |
| `-drc_report_iter_step` | Report DRC on each iteration which is a multiple of this step. The default value is `0`, and the allowed values are integers `[0, MAX_INT]`. |
| `-db_process_node` | Specify the process node. |
| `-disable_via_gen` | Option to diable via generation with bottom and top routing layer. The default value is disabled. | 
| `-droute_end_iter` | Number of detailed routing iterations. The default value is `-1`, and the allowed values are integers `[1, 64]`. |
| `-via_in_pin_bottom_layer` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-via_in_pin_top_layer` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-via_access_layer` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-min_access_points` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-or_seed` | Refer to developer arguments [here](#developer-arguments). |
| `-or_k` | Refer to developer arguments [here](#developer-arguments). |
| `-bottom_routing_layer` | Deprecated. Use the `set_routing_layers` command instead. |
| `-top_routing_layer` | Deprecated. Use the `set_routing_layers` command instead. |
| `-verbose` | Sets verbose mode if the value is greater than 1, else non-verbose mode (must be integer, or error will be triggered.) |
| `-distributed` | Refer to distributed arguments [here](#distributed-arguments). |
| `-remote_host` | Refer to distributed arguments [here](#distributed-arguments). |
| `-remote_port` | Refer to distributed arguments [here](#distributed-arguments). |
| `-shared_volume` | Refer to distributed arguments [here](#distributed-arguments). |
| `-cloud_size` | Refer to distributed arguments [here](#distributed-arguments). |
| `-clean_patches` | Clean unneeded patches during detailed routing. | 
| `-no_pin_access` | Disables pin access for routing. |
| `-bonnroute` | Run the experimental BonnRoute initial detailed router after track assignment, then run FlexDR routing and DRC repair. |
| `-save_guide_updates` | Flag to save guides updates. |
| `-repair_pdn_vias` | This option is used for PDKs where M1 and M2 power rails run in parallel. |

#### Developer arguments

Some arguments that are helpful for developers are listed here.  

| Switch Name | Description |
| ----- | ----- |
| `-or_seed` | Random seed for the order of nets to reroute. The default value is `-1`, and the allowed values are integers `[0, MAX_INT]`. | 
| `-or_k` | Number of swaps is given by $k * sizeof(rerouteNets)$. The default value is `0`, and the allowed values are integers `[0, MAX_INT]`. |

### BonnRoute Profile

Select one of two BonnRoute algorithm profiles for global
(`global_route -bonnroute`) and detailed (`detailed_route -bonnroute`)
routing. Parameters set afterwards with `set_bonnroute_param` or
`set_bonnroute_dr_param` override the profile.

```tcl
set_bonnroute_profile name
```

| Name | Meaning |
| ----- | ----- |
| `2013-improved` | Default. Gester et al. 2013 with measured optimizations: vertex search, speculative tail batches, unresolved nets to FlexDR; GR oracle with the l1 goal-oriented bound and deterministic frozen-price batches. |
| `2021` | `2013-reimpl` plus the bulk routing setting of Ahrens et al., arXiv:2111.06169 (Table 3): the general-model potential for every DR search, in the interval search. Reservations with the 3/4 discount (`set_bonnroute_dr_param reservations 1`) belong to that paper's incremental routing, which starts from a legal detailed routing; see the 2021 TODO list in [the fidelity notes](BonnRouteFidelity.md#2021-todo). |
| `2013-reimpl` | Gester et al. 2013 as written: interval search with property (ii) for vias and the Algorithm 4 tie order, rip-up sequences, the partition sequence to its end; GR oracle with landmark goal orientation and concurrent price updates. Results vary between runs. See [the fidelity notes](BonnRouteFidelity.md) for what remains different. |

### BonnRoute Detailed Routing Parameters

Set BonnRoute detailed routing parameters before `detailed_route -bonnroute`.
Unknown names are errors.

```tcl
set_bonnroute_dr_param name value
```

| Name | Default | Meaning |
| ----- | ----- | ----- |
| `via_cost` | `4` | Via cost in pitches of the lower routing layer. |
| `jog_factor` | `4` | Non-preferred direction wire cost multiplier. |
| `ripup_cost` | `40` | Base rip-up penalty in pitches of the bottom routing layer. |
| `max_ripups` | `3` | Maximum rip-ups per net. |
| `max_access_jog` | `2` | Maximum pin-to-track jog in pitches (access stubs, `offtrack_access 0`). |
| `max_access_options` | `8` | Maximum access options (catalogue paths) retained per pin. |
| `pin_cost` | `2` | Extra cost of a non-preferred access point, in pitches (access stubs). |
| `offtrack_access` | `1` | Off-track pin access (Sections 3.8, 4.3): per pin a catalogue of shortest minimum-segment-length paths on the blockage grid, primary paths per circuit class by destructive bounding; `0` uses OpenROAD access points with short stubs and greedy reservations. |
| `offtrack_radius` | `2` | Radius of the on-track endpoints of catalogue paths around the pin, in pitches. |
| `offtrack_regenerate` | `2` | Generate a pin's catalogue again against the current routing when fewer of its paths are feasible (Section 4.4). |
| `offtrack_tau` | `0` | Minimum segment length in DBU; `0` derives it per layer from the same-net rules (Section 3.7). |
| `offtrack_notch_tau` | `0` | `1` adds notch spacing plus width to the derived minimum segment length. |
| `access_classes` | `1` | Share catalogues and conflict-free solutions within circuit classes; `0` solves every instance separately. |
| `access_node_limit` | `20000` | Node limit of destructive bounding per circuit; beyond it the greedy choice is used. |
| `access_track_weight` / `access_direction_weight` / `access_spread_weight` | `1` / `1` / `1` | Weights of the conflict-free solution: pitches per blocked track, pitches per missing on-track continuation, cost per DBU two endpoints are closer than two pitches. |
| `postprocess` | `1` | Same-net post-processing of each connection (Section 4.4). `offtrack_access 0` with `postprocess 0` is the previous behaviour. |
| `access_class_margin` | `0` | Circuit-class signature area: the instance box plus this many pitches (Section 4.3). Paths mapped from the class representative are re-checked at every member; a pin left without a legal path gets its own catalogue. |
| `offtrack_macro_pins` | `0` | `1` builds off-track catalogues also for pins of block (macro) instances; by default they keep the access stubs (their pins abut the macro obstruction, so no minimum-segment-length path leaves them). |
| `offtrack_preferred_access_points` | `1` | Prefer the access point pin access chose for each pin: paths from its other access points start with a penalty of `pin_cost` pitches, as for the access stubs; `0` treats all access points alike. Access points off the manufacturing grid are never used. |
| `offtrack_from_access_points` | `1` | Start catalogue paths at the pin's OpenROAD access points. `0` starts anywhere on the pin, as in the paper; it is incompatible with the FlexDR clean-up, whose workers need every pin contact at an access point (the contact guard then leaves the affected nets to FlexDR). |
| `layer_dilation` | `1` | Extra routing layers around a guide. |
| `retry_dilation` | `2` | Extra guide tiles for a retry. |
| `use_intervals` | `0` | Use vertex search by default; `1` selects interval label sets. |
| `ripup` | `0` | Leave unresolved nets to FlexDR by default; `1` enables serial rip-up and rerouting. |
| `interval_via_property` | `1` | Interval property (ii) of Section 4.1 also for vias; `0` relaxes vias per vertex. |
| `interval_nearest_first` | `0` | Among interval events with equal key, take the smallest future cost first instead of the largest (Algorithm 4, line 12). Exact either way; `1` stops at targets without labeling the whole last layer of equal-length paths. |
| `component_regions` | `1` | Search each connection first inside the GR tiles of a shortest guide-tile path to another component (Section 4.4). |
| `ripup_depth` | `3` | Maximum depth of a rip-up sequence; a failed sequence restores the original routing (Section 4.4). |
| `ripup_passes` | `2` | Rip-up passes over unrouted nets; each later pass enlarges the routing area by one tile. |
| `partitions` | `1` | Route in a sequence of balanced thread partitions with decreasing region counts (Section 5.1); `0` uses batches of nets with disjoint boxes. |
| `critical_nets` | `1` | Route clock nets and, with `route_ndr_nets 1`, nets with nondefault rules first, serially on the whole chip (critical net routing step of Section 5.1). |
| `route_ndr_nets` | `1` | Route nets with nondefault rules with their wire type: NDR width, NDR preferred vias and NDR spacing (Sections 3.2, 3.6); `0` leaves them to FlexDR. |
| `fast_wire_types` | `5` | Wire types held by the fast grid: the default rule and the NDRs with the most nets, at most 5; other wire types are answered by the distance rule checker. |
| `fast_ndr_permille` | `10` | An NDR becomes a fast wire type only if at least this many per mille of the nets use it; every fast wire type adds its fields to every fast-grid build and update. CTS clock NDRs (a few nets) stay with the checker. |
| `tail_batch` | `32` | Once a partition would have fewer regions than threads, route the remaining nets in batches of this size: searched in parallel on the unchanged routing space, committed in net order, and searched again in the next batch on a conflict with an earlier commit of the same batch. `0` keeps the partition sequence to its end. |
| `tail_threads` | `0` | Threads of the speculative batches; `0` uses all router threads. The batch results are the same for every thread count. |
| `tail_order` | `1` | Start the searches of a batch with the largest estimated effort (guide tiles times terminals) first. Every search of a batch reads the same routing space, so results are the same either way; `0` starts them in net order. |
| `tail_pipeline` | `0` | `1` runs the batches as a pipeline: each batch keeps its semantics (searched on the commits of the earlier batches, decided in net order), and searches and commits of nets with disjoint tile boxes (including the interaction halo) overlap across batches instead of waiting at a barrier per batch. A conflicting net is searched again `tail_retry_delay` batches later. Results depend on the net order, `tail_batch` and `tail_retry_delay` only, not on the thread count. The default `0` uses a barrier per batch and searches a conflicting net again at the front of the next batch: in the E2E runs the pipeline shortened the tail by a further 12-33% (aes_nangate45 most), with TNS 1-3% worse than the noise of the barrier batches on the sky130 designs. |
| `tail_retry_delay` | `1` | With `tail_pipeline 1`, the batch distance at which a net that conflicted with an earlier commit of its batch is searched again. |
| `tail_defer_settled` | `300000` | Settled vertices a net may use in its speculative attempt (checked before each of its searches). A net that needs more, almost always one without a route that goes to FlexDR anyway, is deferred to a final pass of the stage without a budget, in net order, so its long search does not hold up the other nets. `0` disables the budget. |
| `partition_regions` | `0` | Regions of the first partition; `0` means twice the thread count. |
| `reservations` | `0` | Reserve legal track-assignment wires for their nets and discount them (Ahrens et al., arXiv:2111.06169, Section 5.2). Selects the vertex search. |
| `discount_numerator` / `discount_denominator` | `3` / `4` | Cost factor of reserved edges for their own net. |
| `reserve_min_layer` | `0` | Lowest track-graph layer whose input wires become reservations. |
| `potential` | `0` | Future cost: `0` l1 bound of the 2013 paper; arXiv:2111.06169 potentials `1` l1 (Table 3), `2` simple model (Theorem 3), `3` general model on the search area with discounts (Section 4); `4` uses the general model only for connections whose guide path detours (Gester et al. 2013, Section 4.1), the l1 bound otherwise. |
| `exact_future_cost` | `1` | The l1 future cost uses the exact via-stack cost to the nearest target layer (`lb_via` of Section 4.1) and the cheapest terminal cost of each target box; `0` uses the layer distance times the cheapest via. |
| `notch_fill` | `1` | Fill gaps between separated pieces of one net's metal on a layer (wires, via pads, patches, pins) that are closer than drt's spacing with patch wires, when the fill is legal against other nets: same-net spacing of via pads and pins, which the path post-processing does not model. |
| `final_ripup` | `1` | Rip-up sequences (`ripup_depth`, `ripup_passes`, Section 4.4) for the nets the parallel phases leave, also when `ripup` is `0`; nets left to FlexDR are routed in its iteration 0 against the fixed BonnRoute routing. |
| `final_ripup_effort` | `10` | Effort budget of that phase (without `ripup`) in percent of the vertices settled by the parallel phases; nets beyond it are left to FlexDR, smallest nets first. `0`: unlimited. |
| `potential_check` | `0` | Debug: repeat each search that used a potential with the l1 bound and report how many shortest path costs differ (expected 0). |
| `potential_detour_tiles` | `2` | Detour, in GCells beyond the tile distance of its ends, from which `potential 4` uses the general model. |
| `spreading` | `0` | Wire spreading (Section 4.2): in GCells whose global routing utilization is at most `spread_utilization` percent, a wire vertex beside an occupied track costs `spread_cost` pitches more. Reads the congestion that `global_route` stores in the database. Selects the vertex search. |
| `spread_from_extra_space` | `1` | With `spreading 1`, spread a net where BonnRoute global routing assigned it extra space (`set_bonnroute_param extra_space 1`, Section 2.1); designs without that information use the GCell congestion tiles. |
| `spread_cost` / `spread_utilization` | `1` / `70` | Spreading cost in pitches of the bottom layer; utilization threshold in percent. |
| `cleanup_incremental_iterations` | `2` | FlexDR clean-up: how many of FlexDR's incremental iterations 1 and 2 follow iteration 0 (which routes the nets BonnRoute left and checks every tile). Incremental iterations reroute only nets without BonnRoute routing and keep every BonnRoute net fixed; iteration 3, the first that rips up nets at violations, follows the last one. Range `0`–`2`. Fewer incremental iterations pay off when few nets are left to FlexDR; designs with many left nets clean up faster with both. |
| `cleanup_check_only_workers` | `0` | FlexDR clean-up, incremental iterations 0–2: a worker whose nets all carry BonnRoute's routing has no net to reroute and runs only the design-rule check (with FlexDR's patches), without the maze grid graph. The result equals the full worker's. |
| `verbose` | `1` | Print coverage and search statistics. |

The interval search maintains reduced-distance label envelopes on tracks,
removes dominated labels, and propagates jog labels between intervals
(Gester et al., TODAES 18(2), 2013, Section 4.1, Algorithm 4). Via edges are
relaxed per vertex because adjacent layers can have different track pitches.
The fast grid checks wire and jog edges as well as their endpoints, including
obstacles that lie strictly between track vertices.

With `offtrack_access 1` pin access follows Sections 3.7, 3.8 and 4.3:
same-net rules are mapped to a minimum segment length per layer, each pin gets
a catalogue of shortest minimum-segment-length paths on the blockage grid
(Algorithm 3) to nearby on-track points, and a conflict-free primary path per
pin is chosen per circuit class by destructive bounding and reserved. Each
connection is combined with its access paths and post-processed for same-net
rules (Section 4.4). The interpretation choices are listed in
[OffTrack.md](src/br/OffTrack.md). Diff-net legality of off-track geometry
still uses the adapter's spacing model; the paper's full shape-grid exception
checker and blockage-aware future cost are not implemented.
Unrouted nets and rule violations are handled by FlexDR. Final DRC, wirelength
and runtime therefore describe the combined flow. The log reports BonnRoute's
own routed/skipped/failed counts, search effort and elapsed time separately.

By default, nets unresolved after the partitioned phase go straight to
FlexDR. With `ripup 1`, each unresolved net starts a rip-up sequence: victims
are removed, the new route is committed and the victims are reconnected,
recursively up to `ripup_depth`; if any step fails, every net touched by the
sequence is restored. Primary pin access paths are reserved as complete
geometry before routing and released when their net is connected.

See the [paper-fidelity notes](BonnRouteFidelity.md) for what remains
different from the paper. Configuration comparisons in the combined flow are
tradeoffs of that flow, not evidence about a faithful reproduction.

The database-independent search core has a Bazel regression target:

```shell
bazelisk test //src/drt/test:bonnroute_core_test
```

It checks obstacle-edge geometry and compares interval and vertex searches
with a plain Dijkstra reference on generated graphs, including returned path
costs, pin ownership, restricted regions, via edges and rip-up penalties.

The off-track modules have their own randomized test and a measurement
binary on generated standard-cell rows:

```shell
bazelisk test //src/drt/test:bonnroute_offtrack_test
bazelisk run //src/drt/test:bonnroute_offtrack_bench -- 12 40 1
```

The test compares Algorithm 3 with a literal transcription, shortest paths on
the path-preserving digraph with brute-force shortest minimum-segment-length
paths (Theorem 3.2, one and two layers), destructive bounding with exhaustive
enumeration, circuit signatures under the eight lattice isometries, and the
post-processing invariants.

### Detailed Route Debugging

The following command and arguments are useful when debugging error
messages from `drt` and to understand its behavior.

```tcl
detailed_route_debug 
    [-pa]
    [-ta]
    [-dr]
    [-maze]
    [-net name]
    [-pin name]
    [-box x1 y1 x2 y2]
    [-iter iter]
    [-pa_markers]
    [-dump_dr]
    [-dump_dir dir]
    [-snapshot_dir dir]
    [-dump_last_worker]
    [-pa_edge]
    [-pa_commit]
    [-write_net_tracks]
```

#### Options

| Switch Name | Description |
| ----- | ----- |
| `-pa` | Enable debug for pin access. |
| `-ta` | Enable debug for track assignment. |
| `-dr` | Enable debug for detailed routing. |
| `-maze` | Enable debug for maze routing. | 
| `-net` | Enable debug for net name. |
| `-pin` | Enable debug for pin name. |
| `-box` | Set the box for debugging given by lower left/upper right coordinates. |
| `-worker` | Debugs routes that pass through the point `{x, y}`. |
| `-iter` | Specifies the number of debug iterations. The default value is `0`, and the accepted values are integers `[0, MAX_INT`. |
| `-pa_markers` | Enable pin access markers. |
| `-dump_dr` | Filename for detailed routing dump. |
| `-dump_dir` | Directory for detailed routing dump. |
| `-snapshot_dir` | Directory for snapshots produced if the debug level of `DRT snapshot` is set to 1. |
| `-pa_edge` | Enable visibility of pin access edges. |
| `-pa_commit` | Enable visibility of pin access commits. |
| `-write_net_tracks` | Enable writing of net track assigments. |

### Check Pin Access

This function checks pin access.

```tcl
pin_access
    [-db_process_node name]
    [-bottom_routing_layer layer]
    [-top_routing_layer layer]
    [-via_access_layer layer]
    [-via_in_pin_bottom_layer layer]
    [-via_in_pin_top_layer layer]
    [-min_access_points count]
    [-verbose level]
    [-distributed]
    [-remote_host rhost]
    [-remote_port rport]
    [-shared_volume vol]
    [-cloud_size sz]
```

#### Options

| Switch Name | Description |
| ----- | ----- |
| `-db_process_node` | Specify process node. |
| `-bottom_routing_layer` | Deprecated. Use the `set_routing_layers` command instead. |
| `-top_routing_layer` | Deprecated. Use the `set_routing_layers` command instead. |
| `-via_access_layer` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-via_in_pin_bottom_layer` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-via_in_pin_top_layer` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-min_access_points` | Refer to pin access arguments [here](#pin-access-arguments). |
| `-verbose` | Sets verbose mode if the value is greater than 1, else non-verbose mode (must be integer, or error will be triggered.) |
| `-distributed` | Refer to distributed arguments [here](#distributed-arguments). |
| `-remote_host` | Refer to distributed arguments [here](#distributed-arguments). |
| `-remote_port` | Refer to distributed arguments [here](#distributed-arguments). |
| `-shared_volume` | Refer to distributed arguments [here](#distributed-arguments). |
| `-cloud_size` | Refer to distributed arguments [here](#distributed-arguments). |

#### Pin Access Arguments

We have compiled all pin access related arguments in this section.

| Switch Name | Description |
| ----- | ----- |
| `-via_access_layer` | Top via access layer name. Any pin on this layer or below can only be accessed through a via, disallowing planar access. |
| `-via_in_pin_bottom_layer` | Via-in pin bottom layer name. Used in conjunction with via_in_pin_top_layer to specify the layer range for via-in-pin rules. |
| `-via_in_pin_top_layer` | Via-in pin top layer name. Used in conjunction with via_in_pin_bottom_layer to specify the layer range for via-in-pin rules. For pins within this layer range, any via enclosure must be fully contained within the pin boundary and cannot extend beyond the pin edges. |
| `-min_access_points` | Minimum number of access points per pin. |

#### Distributed Arguments

We have compiled all distributed arguments in this section.

```{note}
Additional setup is required. Please refer to this [guide](./doc/Distributed.md).
```

| Switch Name | Description |
| ----- | ----- |
| `-distributed` | Enable distributed mode with Kubernetes and Google Cloud. |
| `-remote_host` | The host IP. |
| `-remote_port` | The value of the port to access from. |
| `-shared_volume` | The mount path of the nfs shared folder. |
| `-cloud_size` | The number of workers. |

### Set Routing Watermark Strength

This command sets the cost multiplier the detailed router applies to wiring
against a layer's preferred direction on nets tagged by `set_routing_watermark`.
A larger value suppresses wrong-way wiring more strongly on those nets, which is
what carries the routing watermark. Preferred-direction wiring, pin access and
untagged nets are unaffected. This is router configuration and does not
persist in the database, so call it in the same process as `detailed_route`,
before it.

```tcl
set_routing_watermark_strength
    strength
```

#### Options

| Switch Name | Description |
| ----- | ----- |
| `strength` | Finite cost multiplier from `0` to `10000` inclusive. `1` tags nets without biasing them. Defaults to `100`. |

Scaled edge costs and the path costs built from them saturate at the router's
maximum cost rather than wrapping. Very large strengths make distinct expensive
paths indistinguishable, which is why the range is bounded.

### Get Routing Watermark Strength

This command returns the cost multiplier currently in effect.

```tcl
get_routing_watermark_strength
```

## Useful Developer Commands

If you are a developer, you might find these useful. More details can be found in the [source file](./src/TritonRoute.cpp) or the [swig file](./src/TritonRoute.i).

| Command Name | Description |
| ----- | ----- |
| `detailed_route_set_default_via` | Set default via. |
| `detailed_route_set_unidirectional_layer` | Set unidirectional layer. |
| `step_dr` | Refer to function `detailed_route_step_drt`. | 
| `check_drc` | Refer to function `check_drc_cmd`. |



## Example scripts

Example script demonstrating how to run TritonRoute on a sample design of `gcd`
in the Nangate45 technology node.

```shell
# single machine example 
./test/gcd_nangate45.tcl

# distributed example
./test/gcd_nangate45_distributed.tcl
```

## Regression tests

There are a set of regression tests in `./test`. For more information, refer to this [section](../../README.md#regression-tests). 

Simply run the following command:

```shell
bazelisk test //src/drt/...
```

## Limitations

## References

Please cite the following paper(s) for publication:

1.   A. B. Kahng, L. Wang and B. Xu, "TritonRoute: The Open Source Detailed
    Router", IEEE Transactions on Computer-Aided Design of Integrated Circuits
    and Systems (2020), doi:10.1109/TCAD.2020.3003234. [(.pdf)](https://ieeexplore.ieee.org/ielaam/43/9358030/9120211-aam.pdf)
1.   A. B. Kahng, L. Wang and B. Xu, "The Tao of PAO: Anatomy of a Pin Access
    Oracle for Detailed Routing", Proc. ACM/IEEE Design Automation Conf., 2020,
    pp. 1-6. [(.pdf)](https://vlsicad.ucsd.edu/Publications/Conferences/377/c377.pdf)

## Authors

TritonRoute was developed by graduate students Lutong Wang and
Bangqi Xu from UC San Diego, and serves as the detailed router in the
[OpenROAD](https://theopenroadproject.org/) project.

## License

BSD 3-Clause License. See [LICENSE](../../LICENSE) file.
