# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Compare the CP-SAT optimizer against exhaustive, independent enumeration."""

from itertools import combinations, product
import random
import unittest

from _placement_model import Candidate, Cell, Problem, evaluate, solve


def enumerate_best(problem):
    best = None
    for choices in product(*(range(len(c.candidates)) for c in problem.cells)):
        rectangles = []
        pins = {name: list(points) for name, points in problem.anchors.items()}
        displacement = mirror_count = 0
        for cell, index in zip(problem.cells, choices):
            pose = cell.candidates[index]
            rectangles.append(
                (
                    pose.x - cell.pad_left,
                    pose.y,
                    pose.x + cell.width + cell.pad_right,
                    pose.y + cell.height,
                )
            )
            displacement += pose.displacement
            mirror_count += int(pose.mirrored)
            for net, point in zip(cell.pin_nets, pose.pins):
                pins.setdefault(net, []).append(point)
        if any(
            max(a[0], b[0]) < min(a[2], b[2]) and max(a[1], b[1]) < min(a[3], b[3])
            for a, b in combinations(rectangles, 2)
        ):
            continue
        length = 0
        for points in pins.values():
            if points:
                xs, ys = zip(*points)
                length += max(xs) - min(xs) + max(ys) - min(ys)
        value = length, displacement, mirror_count
        if best is None or value < best:
            best = value
    return best


def generated_problems():
    yield Problem((), {})
    yield Problem((), {"boundary": ((-5, -2), (7, 3)), "empty": ()})
    yield Problem((Cell("none", 2, 1, 0, 0, (), ()),), {})
    # Non-overlap must also prove infeasibility, not only empty domains.
    pose = Candidate(0, 0, "R0", 0, False, ())
    yield Problem(
        (Cell("a", 2, 1, 0, 0, (), (pose,)), Cell("b", 2, 1, 0, 0, (), (pose,))), {}
    )
    adjacent = Candidate(2, 0, "R0", 0, False, ())
    yield Problem(
        (Cell("a", 2, 1, 0, 1, (), (pose,)), Cell("b", 2, 1, 0, 0, (), (adjacent,))), {}
    )
    for seed in range(40):
        rng = random.Random(seed)
        cells = []
        for index in range(1 + seed % 3):
            nets = () if seed % 5 == 0 else ("shared", "shared", f"private{index}")
            poses = []
            for k in range(2 + seed % 3):
                x, y = rng.randrange(-3, 6), rng.randrange(-1, 2)
                mirrored = bool(k % 2)
                pins = ((x, y), (x + 2, y + 1), (x + 1, y)) if nets else ()
                poses.append(
                    Candidate(
                        x,
                        y,
                        "MY" if mirrored else "R0",
                        abs(x - index * 2),
                        mirrored,
                        pins,
                    )
                )
            cells.append(
                Cell(str(index), 2, 1, seed % 2, (seed // 2) % 2, nets, tuple(poses))
            )
        anchors = {} if seed % 5 == 0 else {"shared": ((-4, -2), (8, 3))}
        yield Problem(tuple(cells), anchors)


class PlacementModelTest(unittest.TestCase):
    def test_optimum_and_unsat_match_exhaustive_enumeration(self):
        for index, problem in enumerate(generated_problems()):
            with self.subTest(problem=index):
                expected = enumerate_best(problem)
                result = solve(problem, 10_000)
                if expected is None:
                    self.assertEqual(result.status, "unsat")
                else:
                    self.assertEqual(result.status, "optimal", result.reason)
                    self.assertEqual(evaluate(problem, result.choices), expected)


if __name__ == "__main__":
    unittest.main()
