# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Finite placement choices, independent of OpenDB and the SAT/SMT backend."""

from dataclasses import asdict, dataclass
from itertools import combinations
import time


@dataclass(frozen=True)
class Candidate:
    x: int
    y: int
    orient: str
    displacement: int
    mirrored: bool
    # One integer pin centroid for each entry in Cell.pin_nets.
    pins: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class Cell:
    name: str
    width: int
    height: int
    pad_left: int
    pad_right: int
    pin_nets: tuple[str, ...]
    candidates: tuple[Candidate, ...]


@dataclass(frozen=True)
class Problem:
    cells: tuple[Cell, ...]
    # Fixed boundary pin centroids, compressed to bounding-box corners.
    anchors: dict[str, tuple[tuple[int, int], ...]]

    def to_dict(self):
        return {"schema": "openroad.local-placement.v1", "units": "DBU", **asdict(self)}


@dataclass(frozen=True)
class Solution:
    status: str
    choices: tuple[int, ...] = ()
    reason: str = ""
    seconds: float = 0.0
    solver_version: str = ""


def bounds(cell, candidate):
    return (
        candidate.x - cell.pad_left,
        candidate.y,
        candidate.x + cell.width + cell.pad_right,
        candidate.y + cell.height,
    )


def overlaps(a, b):
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def evaluate(problem, choices):
    """Check a decoded model and compute its objective without calling the solver."""
    if len(choices) != len(problem.cells):
        raise ValueError("Wrong number of placement choices")
    selected = []
    pins = {name: list(points) for name, points in problem.anchors.items()}
    for cell, choice in zip(problem.cells, choices):
        if not 0 <= choice < len(cell.candidates):
            raise ValueError(f"Invalid candidate for {cell.name}")
        candidate = cell.candidates[choice]
        selected.append((cell, candidate))
        for name, point in zip(cell.pin_nets, candidate.pins):
            pins.setdefault(name, []).append(point)
    for (cell_a, a), (cell_b, b) in combinations(selected, 2):
        if overlaps(bounds(cell_a, a), bounds(cell_b, b)):
            raise ValueError(f"Overlapping candidates: {cell_a.name}, {cell_b.name}")
    hpwl = sum(
        max(p[0] for p in points)
        - min(p[0] for p in points)
        + max(p[1] for p in points)
        - min(p[1] for p in points)
        for points in pins.values()
        if points
    )
    return (
        hpwl,
        sum(c.displacement for _, c in selected),
        sum(c.mirrored for _, c in selected),
    )


def _build_model(problem, cp_model):
    model = cp_model.CpModel()
    choices, x_intervals, y_intervals, displacement, mirrors = [], [], [], [], []
    pins = {name: list(points) for name, points in problem.anchors.items()}
    domains = {name: list(points) for name, points in problem.anchors.items()}
    for index, cell in enumerate(problem.cells):
        if not cell.candidates:
            model.add(False)
            continue
        choice = model.new_int_var(0, len(cell.candidates) - 1, f"choice_{index}")
        choices.append(choice)

        def choose(label, values):
            value = model.new_int_var(min(values), max(values), f"{label}_{index}")
            model.add_element(choice, values, value)
            return value

        x = choose("x", [c.x for c in cell.candidates])
        y = choose("y", [c.y for c in cell.candidates])
        x_intervals.append(
            model.new_fixed_size_interval_var(
                x - cell.pad_left,
                cell.width + cell.pad_left + cell.pad_right,
                f"x_interval_{index}",
            )
        )
        y_intervals.append(
            model.new_fixed_size_interval_var(y, cell.height, f"y_interval_{index}")
        )
        displacement.append(
            choose("displacement", [c.displacement for c in cell.candidates])
        )
        mirrors.append(choose("mirrored", [int(c.mirrored) for c in cell.candidates]))
        for pin, net in enumerate(cell.pin_nets):
            px = choose(f"pin_{pin}_x", [c.pins[pin][0] for c in cell.candidates])
            py = choose(f"pin_{pin}_y", [c.pins[pin][1] for c in cell.candidates])
            pins.setdefault(net, []).append((px, py))
            domains.setdefault(net, []).extend(c.pins[pin] for c in cell.candidates)
    if x_intervals:
        model.add_no_overlap_2d(x_intervals, y_intervals)
    spans = []
    for index, (net, points) in enumerate(pins.items()):
        if not points:
            continue
        for axis in (0, 1):
            low = min(p[axis] for p in domains[net])
            high = max(p[axis] for p in domains[net])
            minimum = model.new_int_var(low, high, f"n{index}_{axis}_min")
            maximum = model.new_int_var(low, high, f"n{index}_{axis}_max")
            model.add_min_equality(minimum, [p[axis] for p in points])
            model.add_max_equality(maximum, [p[axis] for p in points])
            spans.append(maximum - minimum)
    return model, choices, (sum(spans), sum(displacement), sum(mirrors))


def solve(problem, timeout_ms):
    """Minimize HPWL, displacement, then mirrors with a shared time budget.

    CP-SAT remains optional for the rest of OpenROAD. Each objective is fixed
    only after OPTIMAL; FEASIBLE/UNKNOWN never yields an applicable placement.
    """
    try:
        import ortools
        from ortools.sat.python import cp_model
    except ImportError as error:
        raise ImportError(
            "Local placement requires OR-Tools; install the pinned dependency "
            "from bazel/requirements_cp_sat_lock.txt"
        ) from error

    started = time.monotonic()
    model, choices, objectives = _build_model(problem, cp_model)
    validation_error = model.validate()
    if validation_error:
        raise ValueError(f"Invalid CP-SAT placement model: {validation_error}")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    values = []
    for objective in objectives:
        remaining = timeout_ms / 1000 - (time.monotonic() - started)
        if remaining <= 0:
            return Solution(
                "unknown",
                reason="Time budget exhausted",
                seconds=time.monotonic() - started,
                solver_version=ortools.__version__,
            )
        solver.parameters.max_time_in_seconds = remaining
        model.minimize(objective)
        status = solver.solve(model)
        elapsed = time.monotonic() - started
        if status == cp_model.INFEASIBLE:
            return Solution(
                "unsat", seconds=elapsed, solver_version=ortools.__version__
            )
        if status == cp_model.MODEL_INVALID:
            raise ValueError(
                f"Invalid CP-SAT placement model: {solver.response_stats()}"
            )
        if status != cp_model.OPTIMAL:
            return Solution(
                "unknown",
                reason=solver.status_name(status),
                seconds=elapsed,
                solver_version=ortools.__version__,
            )
        value = solver.value(objective)
        values.append(value)
        model.add(objective == value)
    selected = tuple(solver.value(choice) for choice in choices)
    if evaluate(problem, selected) != tuple(values):
        raise RuntimeError("Decoded model does not match the optimization objective")
    return Solution(
        "optimal",
        selected,
        seconds=time.monotonic() - started,
        solver_version=ortools.__version__,
    )
