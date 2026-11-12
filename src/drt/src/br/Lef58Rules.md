# LEF58 rules of asap7 in the BonnRoute routing space

The routing space (shape grid, distance rule checker, fast grid; Gester et
al. 2013, Sec. 3) is BonnRoute's DRC: whatever it does not model is left to
the FlexDR clean-up. This note lists the asap7 rules behind the FlexDR
iteration-0 markers after BonnRoute on aes_asap7 (attribution
`e2e/attrib/REPORT.md`, `p1/REPORT.md`, `p1/REPORT-p1b.md`), how drt checks
them, their parameters in `test/asap7/asap7_tech_1x_201209.lef` (DATABASE
MICRONS 4000: 1 nm = 4 DBU) and how the routing space treats them.

## Inventory (base 81526af603)

| marker class (aes_asap7) | drt check | asap7 rule (um / DBU) | nets | shape kinds | BonnRoute at 81526af603 |
| --- | --- | --- | --- | --- | --- |
| Lef58EolKeepOut M2 723, M3 364, M4 256, M1 47, M5 26 (iteration 0); 30 / 37 / 92 / 54 / 6 in the BonnRoute output | `FlexGC_eol.cpp` `checkMetalEOLkeepout_main` / `_helper`, `getEolKeepOutQueryBox` | M1 `EOLKEEPOUT 0.01825 EXTENSION 0.0 0.0 0.031` (EOL width 73, backward 0, side 0, forward 124), whole-rectangle form. M2, M3 `EOLKEEPOUT 0.025 EXTENSION 0.0 0.0125 0.031 CORNERONLY` (100; 0; 50; 124). M4, M5 `EOLKEEPOUT 0.025 EXTENSION 0.048 0.02425 0.048 CORNERONLY` (100; 192; 97; 192). M6, M7 `EOLKEEPOUT 0.05 EXTENSION 0.048 0.03225 0.048 CORNERONLY` (200; 192; 129; 192) | both: drt exempts only the same polygon (`getPin`) and fixed-fixed pairs, so separate polygons of one net count | EOL side (any direction, shorter than the EOL width, convex corners) of a wire end, via pad, pin or obstruction against: CORNERONLY: any polygon vertex of another polygon strictly inside the region; otherwise any rectangle of another polygon overlapping it | Folded into drt's EOL summary (`drEolSpacingConstraint`: space = max(EOL space, ENDTOEND, keep-out forward), within = max(EOL within, keep-out side), width = min EOL width - 1) and modelled as the line-end extension (Fig. 2) at the preferred-direction ends of preferred wires and via pads. Backward part, corner form, jogs, pins and obstructions not modelled; M1 width 72 = summary width 72, `edge < 72` false: M1 pads and pins carry no line-end (drt: 72 < 73). Same net not modelled |
| Lef58CutSpacingTable V3 280-450, V4 20 (343 / 5 in the BonnRoute output) | `FlexGC_cut.cpp` `checkLef58CutSpacingTbl` / `_main` / `_helper` | V3-V7 `LEF58_CUTCLASS` (V3 0.018 x 0.024 = 72 x 96; V4 Vx 96 x 96; V5 Vx 96 x 128; multi-cut classes up to 16 cuts) and `LEF58_SPACINGTABLE DEFAULT 0.034` (136) with every entry `-` (= DEFAULT): every cut pair needs 136 edge to edge (Euclidean). V3 has no plain `SPACING` (commented out) | both: no SAMENET / SAMEMETAL table, so the table applies to cuts of one net | cut against cut; in the attribution same-net VIA34 pairs 72 apart (adjacent M3 tracks: an M4 bridge) and 96 apart (adjacent M4 tracks: an M3 bridge) | Cut spacing = the layer's plain `SPACING`, else the default via's cut `minDXDY` (V3 72, V4 96); cut classes and the table not modelled. Same net not modelled |
| Lef58SpacingEndOfLine M2 88, M3 90, M4 8-15 | `FlexGC_eol.cpp` `checkMetalEndOfLine_eol*` (`_hasEol_getQueryBox`, `_endToEndHelper`, `_hasParallelEdge`) | M2 `SPACING 0.018 ENDOFLINE 0.025 WITHIN 0.02 ENDTOEND 0.031 PARALLELEDGE 0.025 WITHIN 0.02` (space 72, EOL width 100, within 80, end-to-end 124, parallel edge 100 / 80); M3 same with WITHIN 0.0125 (50); M4, M5 `SPACING 0.024 ENDOFLINE 0.025 WITHIN 0.04 ENDTOEND 0.04` (96, 100, 160, 160); M6, M7 `SPACING 0.032 ENDOFLINE 0.0375 WITHIN 0.04 ENDTOEND 0.04` (128, 150, 160, 160) | diff-net (and same net between polygons) | EOL side against an opposite-facing side of another polygon closer than the space and within `within` laterally; end-to-end: both sides EOL, overlapping laterally, closer than the end-to-end space | Same summary line-end as above (on M4: 192 deep, +-160 wide, at preferred wires and via pads); pins, obstructions and jogs carry none |
| Corner Spacing M1 319, M3 316 (in the BonnRoute output; FlexDR clears them in iteration 0) | `FlexGC_main.cpp` `checkMetalCornerSpacing_main` | M1 `CORNERSPACING CONVEXCORNER CORNERONLY 0.01 WIDTH 0.018 SPACING 0.018` (72); M2, M3 `CONVEXCORNER WIDTH 0.0 SPACING 0.02` (80); M4-M7 `CONVEXCORNER CORNERONLY 0.048 WIDTH 0.0 SPACING 0.04` (160). drt keeps CORNERONLY's within but does not use it; equal x/y spacing: violation if max(dx, dy) < spacing | both (no same-net exemption in drt; fixed-fixed exempt) | convex corner of a polygon against a rectangle diagonally beyond it whose nearest corner is a polygon corner; on M1 a BonnRoute VIA12 pad against a neighbouring cell's obstruction | Euclidean diagonal with the run-length-0 spacing (72 on M1-M3): optimistic for 72 <= l2 and max(dx, dy) < 72 / 80 |
| same-net keep-out and cut spacing between via pads / cuts of one net (VIA23 pads 32 apart on M2 / M3, VIA34 pads 96 apart on M4, VIA34 cuts 72 / 96 apart) | the checks above with both shapes on one net | as above | same net | via pad / via pad on the layer where the two vias are not connected; cut / cut | not modelled: the search removes the net's own shapes (Sec. 4.4) and the post-processing models notch, short edge and minimum area only |

## Model (`set_bonnroute_dr_param lef58_rules 1`, the default)

`BrDrtRules` builds the rules from the drt tech for BonnRouteDR and for the
differential test. A routing layer with LEF58 end-of-line spacing,
EOLKEEPOUT or CORNERSPACING rules gets *side rules* and a corner spacing in
place of drt's EOL summary (the summary line-end rule stays only for LEF58
EOL extension rules); other layers, and every layer with `lef58_rules 0`,
keep the previous model.

Real sides and corners. Every metal shape class knows which of its four
sides are polygon sides (`ShapeClass::sides_real`) and, per EOL width of the
layer, which of them are end-of-line sides (`eol_sides`: shorter than the
width; both corners of a rectangle side are convex). A routed
preferred-direction wire ends at a via whose pad covers the wire's end
square on asap7 (`LayerRules::wire_ends_covered`, checked per layer for
every wire type's pads), so the wire's ends are no polygon sides and its
flanks only are real: the pad carries the end-of-line side and the corners,
as in the merged polygon drt checks. Without such pads, wires keep their
ends. Via pads, jogs, patches (role via pad on such layers), pins and
obstructions are rectangles with real sides; a pin or obstruction narrower
than the EOL width has end-of-line sides in both directions. A stored piece
of a shape (shape-grid cell) marks its sides that are cell or clip borders
(`PlacedShape::line_end_ends` bits 2-5, kept on side-rule layers only), so
the rules see the sides and corners of the shape, not of its cells.

Rules (`appendSideTests`; each is a pair of rectangles that must not overlap
strictly, one from each shape, so the distance rule checker and the fast
grid compute their violation spans along a line exactly as for the line-end
extensions):

| rule | test | drt conditions not modelled (conservative) |
| --- | --- | --- |
| EOLKEEPOUT (`kRegionOverlap`) | region `forward` beyond an EOL side, `backward` behind it, `side` beyond its ends, overlapping the other shape | EXCEPTWITHIN, CLASS |
| EOLKEEPOUT CORNERONLY (`kRegionCorner`) | a real corner of the other shape strictly inside that region | EXCEPTWITHIN, CLASS |
| ENDOFLINE spacing (`kFacing`) | a real side of the other shape facing the EOL side, `0 < gap < space`, laterally overlapping the side bloated by `within` | PARALLELEDGE, MINLENGTH, ENCLOSECUT, "nothing in between", hasRoute; the larger of space, wrong-direction space and end-PRL spacing |
| ENDTOEND (`kFacingEol`) | as above with the facing side an EOL side, `gap < end-to-end space`, lateral overlap of both sides bloated by the extension | OTHERENDWIDTH |
| CORNERSPACING CONVEXCORNER | a real corner of each shape, the other shape in its quadrant, `max(dx, dy) < spacing` (drt with equal x/y spacing) | width table (largest value), EXCEPTEOL; CORNERTOCORNER taken as the larger Chebyshev region |
| LEF58_SPACINGTABLE (cut) | per cut-class pair the largest entry (SIDE / END, both PRL entries, both orders), DEFAULT without a class; edge to edge, max with the plain cut spacing | centre-to-centre entries taken edge to edge, PRL entries |

A piece whose side ends at a cell border extends its region by one unit
across that border, so the open regions of the pieces of one side cover the
open region of the whole side (a corner exactly on the border lies strictly
inside one of them). `LayerRules::sideReach` adds this unit.

Interaction distances, windows, halos and the shape-grid border include the
side rules (`LayerRules::reach`, `ShapeClassTable::distance`).

## Same net (`same_net_vias 1`, vertex search)

The searches ignore the net's own shapes (Sec. 4.4), so same-net rules
between the vias of a net were never seen. The markers are almost all
*bridges*: the path leaves a layer by a via, runs one or two positions on
the next layer and returns on the same via layer (a track change on a
unidirectional layer). On asap7 V3-V6 every two vias of one net at
neighbouring via positions violate the 136 DBU cut spacing, and the two
pads on the layer the path left (not connected there) violate the M4
keep-out (VIA34 pads 96 apart) or the M3 spacing / keep-out (VIA23 pads 32
apart).

Choice: the vertex search checks a via when it relaxes it. It walks back
from the settled vertex along its layer (parent pointers, at most
`SameNetVias::reach` away) to the via by which the path entered the layer;
if both vias are on the same via layer, the via is refused when the two
cuts, or the two pads on the other layer (unless they touch), violate the
routing-space rules (`pairAllowed`, the diff-net rule model, conservative
for same-net pairs). Every via is also checked against the net's cuts
outside the search (`allowed`: earlier connections and the access vias of
connected terms). The check is local (one walk of a few vertices per
settled vertex with via edges, only on layers with vias) and makes every
committed path keep these rules. It is not exact for the shortest path: a
label-setting search keeps one label per vertex, so a vertex reached
through a refused bridge can hide a costlier path through it that would
have been allowed; such a connection costs more or fails (the net goes to
FlexDR as before). Rejecting in post-processing or `pathToGeometry` would
need a new search per conflict (each retry finds the next bridge on the
next track), and a state per vertex for "entered by a via" would double the
label memory of every search.

Not covered: pads against the net's pads outside the search (their
connectivity is not known during the search), a via pair whose path leaves
the layer in between (down, a detour two layers away, up again near the
first via), the target term's access via against the path, and the
interval search (`use_intervals 1`, profile 2013-reimpl).
