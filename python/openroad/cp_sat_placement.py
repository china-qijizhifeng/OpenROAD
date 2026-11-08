# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Opt-in, pre-route local placement using finite CP-SAT joint choices.

The selected single-height cells stay on their current horizontal rows.
Position and horizontal mirroring are optimized jointly. HPWL is a geometric
proxy, not a timing, routability, or signoff guarantee. Re-estimate parasitics
after an accepted change before running timing analysis.
"""

from dataclasses import dataclass
import json
from pathlib import Path

from ._placement_model import Candidate, Cell, Problem, evaluate, overlaps, solve

_MIRROR = {"R0": "MY", "MY": "R0", "MX": "R180", "R180": "MX"}


@dataclass(frozen=True)
class PlacementResult:
    status: str
    hpwl_before: int
    hpwl_after: int
    moved_instances: tuple[str, ...] = ()
    mirrored_instances: tuple[str, ...] = ()
    displacement_dbu: int = 0
    candidate_count: int = 0
    solver_seconds: float = 0.0
    solver_version: str = ""
    reason: str = ""


def _box(rect):
    return rect.xMin(), rect.yMin(), rect.xMax(), rect.yMax()


def _signal(net):
    return net is not None and net.getSigType() not in ("POWER", "GROUND")


def _pin_location(terminal):
    valid, x, y = terminal.getAvgXY()
    if not valid:
        raise ValueError(f"Missing pin geometry: {terminal.getName()}")
    return x, y


def _nets(instances):
    nets = {}
    for inst in instances:
        for terminal in inst.getITerms():
            net = terminal.getNet()
            if _signal(net):
                nets[net.getName()] = net
    for net in nets.values():
        if (
            net.isSpecial()
            or net.getWire() is not None
            or len(net.getGuides())
            or len(net.getCapNodes())
        ):
            raise ValueError(
                f"CP-SAT placement requires unrouted, unextracted signal nets: {net.getName()}"
            )
    return dict(sorted(nets.items()))


def _anchors(nets, selected):
    anchors = {}
    for name, net in nets.items():
        points = []
        for terminal in net.getITerms():
            inst = terminal.getInst()
            if inst.getName() not in selected:
                if not inst.isPlaced():
                    raise ValueError(f"Unplaced boundary instance: {inst.getName()}")
                points.append(_pin_location(terminal))
        for terminal in net.getBTerms():
            valid, x, y = terminal.getFirstPinLocation()
            if not valid:
                raise ValueError(f"Unplaced boundary port: {terminal.getName()}")
            points.append((x, y))
        if points:
            anchors[name] = (
                (min(p[0] for p in points), min(p[1] for p in points)),
                (max(p[0] for p in points), max(p[1] for p in points)),
            )
        else:
            anchors[name] = ()
    return anchors


def _row_for(inst, rows):
    master = inst.getMaster()
    site = master.getSite()
    bbox = _box(inst.getBBox())
    for row in rows:
        if row.getDirection() != "HORIZONTAL" or row.getSite() != site:
            continue
        rowbox = _box(row.getBBox())
        spacing = row.getSpacing()
        if (
            spacing > 0
            and spacing == site.getWidth()
            and rowbox[1] == bbox[1]
            and rowbox[3] == bbox[3]
            and rowbox[0] <= bbox[0] < bbox[2] <= rowbox[2]
            and (bbox[0] - row.getOrigin()[0]) % spacing == 0
        ):
            return row
    raise ValueError(
        f"Expected a single-height cell on a contiguous horizontal row: {inst.getName()}"
    )


def _pin_sums(inst, terminal):
    """Keep unrounded centroids so reflected pins agree with dbITerm::getAvgXY."""
    ox, oy = inst.getOrigin()
    orient = inst.getOrient()
    sx = -1 if orient in ("MY", "R180") else 1
    sy = -1 if orient in ("MX", "R180") else 1
    xx = yy = count = 0
    for pin in terminal.getMTerm().getMPins():
        for shape in pin.getGeometry():
            xx += sx * (shape.xMin() + shape.xMax()) + 2 * ox
            yy += sy * (shape.yMin() + shape.yMax()) + 2 * oy
            count += 2
    if not count:
        raise ValueError(f"Missing pin geometry: {terminal.getName()}")
    return xx, yy, count


def _truncate(numerator, denominator):
    return (
        numerator // denominator if numerator >= 0 else -((-numerator) // denominator)
    )


def _candidates(inst, row, opendp, obstacles, displacement_sites, allow_mirroring):
    old_x, old_y, right, top = _box(inst.getBBox())
    width, height = right - old_x, top - old_y
    spacing = row.getSpacing()
    left_pad, right_pad = (
        opendp.padLeft(inst) * spacing,
        opendp.padRight(inst) * spacing,
    )
    terminals = sorted(
        (t for t in inst.getITerms() if _signal(t.getNet())),
        key=lambda t: t.getMTerm().getName(),
    )
    pin_sums = [_pin_sums(inst, t) for t in terminals]
    old_orient = inst.getOrient()
    orientations = [old_orient]
    if allow_mirroring and inst.getMaster().getSymmetryY():
        orientations.append(_MIRROR[old_orient])
    rowbox = _box(row.getBBox())
    candidates = []
    for step in range(-displacement_sites, displacement_sites + 1):
        x = old_x + step * spacing
        if x < rowbox[0] or x + width > rowbox[2]:
            continue
        occupied = (x - left_pad, old_y, x + width + right_pad, top)
        if any(overlaps(occupied, obstacle) for obstacle in obstacles):
            continue
        for orient in orientations:
            mirrored = orient != old_orient
            pins = tuple(
                (
                    _truncate(
                        (
                            n * (x + old_x + width) - px
                            if mirrored
                            else px + n * (x - old_x)
                        ),
                        n,
                    ),
                    _truncate(py, n),
                )
                for px, py, n in pin_sums
            )
            candidates.append(
                Candidate(x, old_y, orient, abs(x - old_x), mirrored, pins)
            )
    # The original pose is candidate zero whenever permitted by the model.
    candidates.sort(key=lambda c: (c.displacement, c.mirrored, c.x, c.orient))
    return Cell(
        inst.getName(),
        width,
        height,
        left_pad,
        right_pad,
        tuple(t.getNet().getName() for t in terminals),
        tuple(candidates),
    )


def _extract(design, names, displacement_sites, allow_mirroring, max_candidates):
    block = design.getBlock()
    if block is None:
        raise ValueError("Load a design before CP-SAT placement")
    instances = []
    for name in names:
        inst = block.findInst(name)
        if inst is None:
            raise ValueError(f"Unknown instance: {name}")
        if (
            not inst.isCore()
            or inst.isFixed()
            or inst.isDoNotTouch()
            or not inst.isPlaced()
        ):
            raise ValueError(f"Expected a movable, placed core instance: {name}")
        if inst.getRegion() is not None or inst.getGroup() is not None:
            raise ValueError(
                f"Region/group-constrained instances are not supported: {name}"
            )
        if inst.getOrient() not in _MIRROR:
            raise ValueError(f"Unsupported orientation: {name}")
        if any(len(t.getPrefAccessPoints()) for t in inst.getITerms()):
            raise ValueError(f"Run CP-SAT placement before pin access: {name}")
        instances.append(inst)
    opendp = design.getOpendp()
    rows = list(block.getRows())
    selected_rows = [_row_for(inst, rows) for inst in instances]
    nets = _nets(instances)
    anchors = _anchors(nets, set(names))
    obstacles = [_box(b.getBBox()) for b in block.getBlockages()]
    for inst in block.getInsts():
        if inst.getName() in names:
            continue
        if not inst.isPlaced():
            raise ValueError(f"Unplaced instance: {inst.getName()}")
        left, bottom, right, top = _box(inst.getBBox())
        site = inst.getMaster().getSite()
        if site is not None:
            left -= opendp.padLeft(inst) * site.getWidth()
            right += opendp.padRight(inst) * site.getWidth()
        obstacles.append((left, bottom, right, top))
    cells = []
    count = 0
    for inst, row in zip(instances, selected_rows):
        cell = _candidates(
            inst, row, opendp, obstacles, displacement_sites, allow_mirroring
        )
        count += len(cell.candidates)
        if count > max_candidates:
            raise ValueError(
                "Candidate budget exceeded; reduce selected instances or displacement"
            )
        cells.append(cell)
    return Problem(tuple(cells), anchors), instances, nets


def _actual_hpwl(nets):
    total = 0
    for net in nets.values():
        points = [_pin_location(t) for t in net.getITerms()]
        for terminal in net.getBTerms():
            valid, x, y = terminal.getFirstPinLocation()
            if not valid:
                raise ValueError(f"Unplaced boundary port: {terminal.getName()}")
            points.append((x, y))
        if points:
            total += max(p[0] for p in points) - min(p[0] for p in points)
            total += max(p[1] for p in points) - min(p[1] for p in points)
    return total


def _apply(design, instances, problem, choices, nets, expected_hpwl):
    saved = [(inst, tuple(inst.getOrigin()), inst.getOrient()) for inst in instances]
    try:
        for inst, cell, choice in zip(instances, problem.cells, choices):
            candidate = cell.candidates[choice]
            inst.setOrient(candidate.orient)
            inst.setLocation(candidate.x, candidate.y)
        # Includes rules not encoded in the small spatial model (e.g. edge spacing).
        design.getOpendp().checkPlacement(False)
        if _actual_hpwl(nets) != expected_hpwl:
            raise ValueError("OpenDB pin HPWL disagrees with the decoded model")
    except BaseException:
        for inst, origin, orient in saved:
            inst.setOrient(orient)
            inst.setOrigin(*origin)
        raise


def optimize_local_placement(
    design,
    instances,
    *,
    max_displacement_sites=4,
    allow_mirroring=True,
    timeout_ms=10_000,
    max_cells=16,
    max_candidates=1024,
    problem_path=None,
):
    """Jointly place named cells in a loaded OpenROAD Design, before routing.

    Select explicit instances to bound the search; all other objects remain
    fixed. Distances/results use integer DBU, displacement bounds use sites.
    Existing placement must pass checkPlacement. Routed/extracted affected nets,
    multi-height or region/group-constrained selected cells are rejected.

    The candidate with minimum affected-net pin-centroid HPWL is chosen, then
    minimum total displacement, then minimum mirror count. Only a proven optimum
    of this finite model with a strict HPWL improvement is applied. A solver
    timeout/unknown/unsat leaves the database unchanged. Failed physical checks
    restore the original origins/orientations and return rejected; unexpected
    exceptions also restore the placement before propagating.

    problem_path optionally exports the finite model as versioned JSON for
    reproducibility and future independently verified SAT backends.
    """
    if isinstance(instances, str):
        raise ValueError("Pass a sequence of instance names, not a single string")
    names = tuple(instances)
    if (
        not names
        or any(not isinstance(name, str) for name in names)
        or len(set(names)) != len(names)
    ):
        raise ValueError(
            "instances must be a nonempty sequence of distinct instance names"
        )
    if type(allow_mirroring) is not bool:
        raise ValueError("allow_mirroring must be a boolean")
    for name, value, minimum in (
        ("max_displacement_sites", max_displacement_sites, 0),
        ("timeout_ms", timeout_ms, 1),
        ("max_cells", max_cells, 1),
        ("max_candidates", max_candidates, 1),
    ):
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if len(names) > max_cells or (2 * max_displacement_sites + 1) > max_candidates:
        raise ValueError(
            "Search budget exceeded; reduce the window or raise the explicit budget"
        )
    problem, selected, nets = _extract(
        design, names, max_displacement_sites, allow_mirroring, max_candidates
    )
    design.getOpendp().checkPlacement(False)
    before = _actual_hpwl(nets)
    if problem_path is not None:
        Path(problem_path).write_text(json.dumps(problem.to_dict(), indent=2) + "\n")
    solution = solve(problem, timeout_ms)
    common = dict(
        hpwl_before=before,
        hpwl_after=before,
        candidate_count=sum(len(c.candidates) for c in problem.cells),
        solver_seconds=solution.seconds,
        solver_version=solution.solver_version,
    )
    if solution.status != "optimal":
        return PlacementResult(solution.status, reason=solution.reason, **common)
    hpwl, displacement, _ = evaluate(problem, solution.choices)
    if hpwl >= before:
        return PlacementResult("unchanged", **common)
    try:
        _apply(design, selected, problem, solution.choices, nets, hpwl)
    except RuntimeError as error:
        return PlacementResult("rejected", reason=str(error), **common)
    chosen = [
        (cell, cell.candidates[k]) for cell, k in zip(problem.cells, solution.choices)
    ]
    common["hpwl_after"] = hpwl
    return PlacementResult(
        "applied",
        **common,
        displacement_dbu=displacement,
        moved_instances=tuple(
            cell.name for cell, candidate in chosen if candidate.displacement
        ),
        mirrored_instances=tuple(
            cell.name for cell, candidate in chosen if candidate.mirrored
        ),
    )


def main():
    """Checkpoint-to-checkpoint entry point using the installed OpenROAD wheel."""
    import argparse
    from dataclasses import asdict
    from openroad import Design, Tech

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_odb", type=Path)
    parser.add_argument("output_odb", type=Path)
    parser.add_argument("--instances", nargs="+", required=True)
    parser.add_argument("--max-displacement-sites", type=int, default=4)
    parser.add_argument("--timeout-ms", type=int, default=10_000)
    parser.add_argument("--no-mirroring", action="store_true")
    parser.add_argument("--problem", type=Path)
    args = parser.parse_args()
    if args.input_odb.resolve() == args.output_odb.resolve():
        parser.error("Input and output checkpoints must differ")
    tech = Tech()
    design = Design(tech)
    design.readDb(str(args.input_odb))
    result = optimize_local_placement(
        design,
        args.instances,
        max_displacement_sites=args.max_displacement_sites,
        allow_mirroring=not args.no_mirroring,
        timeout_ms=args.timeout_ms,
        problem_path=args.problem,
    )
    print(json.dumps(asdict(result), indent=2))
    if result.status not in ("applied", "unchanged"):
        return 2
    design.writeDb(str(args.output_odb))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
