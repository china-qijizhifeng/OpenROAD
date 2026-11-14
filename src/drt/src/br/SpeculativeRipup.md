# Ordered speculative rip-up searches

`speculative_ripup N` precomputes the root rip-up-mode connect for at most N
queue entries against one committed routing snapshot. The default is zero.
Incremental routing, rip-up boxes, region-parallel rip-up and zero-depth
sequences retain their existing execution path.

Set it after selecting the routing profile, before detailed routing:

```tcl
set_bonnroute_profile 2013-improved
set_bonnroute_dr_param speculative_ripup 8
detailed_route -bonnroute
```

The value is the queue batch size, clamped to 0..256; routing threads are still
selected independently. It is an opt-in experiment and does not change the
default profile settings.

Each search owns a NetData copy and worker scratch. It reads the routing space
and collects its history increments without applying them. The main thread
waits for the whole batch before making queue decisions. It then accepts or
reexecutes each root search in the original order, applies its history, removes
victims, commits the root, and recursively reconnects victims as before.

A candidate is reusable only if its root revision is unchanged and no earlier
sequence wrote any tile read by its search. Read footprints include shape-grid,
wide-shape, fast-grid and history queries. Writes include insertions, removals,
reservation updates, history increments and rollback. Spatial boxes deliberately
combine all layers and all steps, accepting false invalidations. A removed
obstacle can improve the optimum without making the old path illegal, so path
legality alone does not establish reuse. Failed searches follow the same rule.
The serial decision phase records only writes; its reads do not invalidate
another search and are excluded from footprint instrumentation.
Any root or victim in a sequence journal increments its revision, including on
rollback, because metadata such as the rip-up counter can change.

The phase's settled-vertex budget counts the accepted root search and its
serial recursive work exactly once. Invalidated, skipped and budget-excluded
precomputations are reported as wasted work and do not advance that budget.
Worker statistics include physical search work. The speculative summary also
reports the logical settled count, reuse, invalidation and parallel search time.

Regression coverage uses both vertex and interval searches and both fast-grid
encodings. It covers a still-legal but suboptimal path after obstacle removal,
a failed search made feasible by removal, distant writes that allow reuse,
rollback and root revisions, and a history-only change of the optimal path.
Real route validation compares the same binary with the switch off/on, fixed
threads and input, including budget boundaries and recursive failures.
