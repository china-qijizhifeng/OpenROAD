# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Exercise local CP-SAT placement against real OpenDB and detailed-placement checks."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from openroad import Design, Tech
from openroad import cp_sat_placement
import odb
from openroad._placement_model import Solution

HERE = Path(__file__).resolve().parent


def snapshot(block):
    poses = {
        inst.getName(): (
            tuple(inst.getOrigin()),
            inst.getOrient(),
            inst.getPlacementStatus(),
        )
        for inst in block.getInsts()
    }
    connectivity = {
        net.getName(): (
            tuple(
                sorted(
                    (t.getInst().getName(), t.getMTerm().getName())
                    for t in net.getITerms()
                )
            ),
            tuple(sorted(t.getName() for t in net.getBTerms())),
        )
        for net in block.getNets()
    }
    ports = {
        term.getName(): tuple(term.getFirstPinLocation()) for term in block.getBTerms()
    }
    return poses, connectivity, ports


class CpSatPlacementTest(unittest.TestCase):
    def setUp(self):
        self.tech = Tech()
        self.tech.readLef(str(HERE / "cp_sat_test.lef"))
        self.design = Design(self.tech)
        self.design.readDef(str(HERE / "cp_sat_test.def"))
        self.block = self.design.getBlock()
        self.design.getOpendp().checkPlacement(False)

    def test_packed_row_requires_joint_swap_and_survives_database_roundtrip(self):
        before = snapshot(self.block)
        single = cp_sat_placement.optimize_local_placement(
            self.design, ["A"], allow_mirroring=False
        )
        self.assertEqual(single.status, "unchanged")
        self.assertEqual(snapshot(self.block), before)
        with tempfile.TemporaryDirectory() as directory:
            problem = Path(directory) / "problem.json"
            result = cp_sat_placement.optimize_local_placement(
                self.design, ["A", "B"], allow_mirroring=False, problem_path=problem
            )
            self.assertEqual(result.status, "applied")
            self.assertEqual((result.hpwl_before, result.hpwl_after), (10000, 6000))
            self.assertEqual(set(result.moved_instances), {"A", "B"})
            self.assertEqual(result.displacement_dbu, 4000)
            self.assertEqual(tuple(self.block.findInst("A").getLocation()), (4000, 0))
            self.assertEqual(tuple(self.block.findInst("B").getLocation()), (2000, 0))
            after = snapshot(self.block)
            self.assertEqual(before[1:], after[1:])
            for name in ("L", "R", "M", "N"):
                self.assertEqual(before[0][name], after[0][name])
            self.assertEqual(
                json.loads(problem.read_text())["schema"], "openroad.local-placement.v1"
            )
            path = Path(directory) / "placed.odb"
            self.design.writeDb(str(path))
            restored_tech = Tech()
            restored = Design(restored_tech)
            restored.readDb(str(path))
            restored.getOpendp().checkPlacement(False)
            self.assertEqual(snapshot(restored.getBlock()), after)

    def test_pin_geometry_drives_horizontal_mirroring_on_both_row_orientations(self):
        before = snapshot(self.block)
        for name, expected in (("M", "R180"), ("N", "MY")):
            with self.subTest(instance=name):
                result = cp_sat_placement.optimize_local_placement(
                    self.design, [name], max_displacement_sites=0
                )
                self.assertEqual(result.status, "applied")
                self.assertEqual(result.moved_instances, ())
                self.assertEqual(result.mirrored_instances, (name,))
                self.assertEqual(result.hpwl_before - result.hpwl_after, 1500)
                self.assertEqual(self.block.findInst(name).getOrient(), expected)
                self.design.getOpendp().checkPlacement(False)
        self.assertEqual(snapshot(self.block)[1:], before[1:])

    def test_external_instance_pin_is_a_fixed_anchor(self):
        net = odb.dbNet_create(self.block, "external_anchor")
        self.block.findInst("N").findITerm("A").connect(net)
        self.block.findInst("R").findITerm("A").connect(net)
        before = snapshot(self.block)
        result = cp_sat_placement.optimize_local_placement(
            self.design, ["N"], allow_mirroring=False
        )
        self.assertEqual(result.status, "applied")
        self.assertEqual((result.hpwl_before, result.hpwl_after), (8000, 4000))
        self.assertEqual(tuple(self.block.findInst("N").getLocation()), (6000, 4000))
        after = snapshot(self.block)
        self.assertEqual(after[0]["R"], before[0]["R"])
        self.assertEqual(after[1:], before[1:])
        self.design.getOpendp().checkPlacement(False)

    def test_padding_limits_motion_before_fixed_neighbor(self):
        inst = self.block.findInst("N")
        guard = odb.dbInst_create(self.block, inst.getMaster(), "upper_guard")
        guard.setLocation(6000, 4000)
        guard.setPlacementStatus("LOCKED")
        self.design.getOpendp().setPadding(inst, 0, 1)
        self.design.getOpendp().checkPlacement(False)
        before = snapshot(self.block)
        result = cp_sat_placement.optimize_local_placement(
            self.design, ["N"], allow_mirroring=False
        )
        # The finite model leaves one empty site. This library has no one-site
        # cell, so real DPL rejects that gap and the transaction must roll back.
        self.assertEqual(result.status, "rejected")
        self.assertIn("DPL-0033", result.reason)
        self.assertEqual(snapshot(self.block), before)
        self.design.getOpendp().checkPlacement(False)

        # Two-site padding permits a legal improvement without a one-site gap.
        inst.setLocation(0, 4000)
        self.design.getOpendp().setPadding(inst, 0, 2)
        self.design.getOpendp().checkPlacement(False)
        result = cp_sat_placement.optimize_local_placement(
            self.design, ["N"], allow_mirroring=False
        )
        self.assertEqual(result.status, "applied")
        self.assertEqual(tuple(inst.getLocation()), (2000, 4000))
        self.assertEqual(result.hpwl_before - result.hpwl_after, 2000)
        self.assertEqual(tuple(guard.getLocation()), (6000, 4000))
        self.design.getOpendp().checkPlacement(False)

    def test_solver_unknown_and_unsat_do_not_mutate_design(self):
        before = snapshot(self.block)
        for status in ("unknown", "unsat"):
            with self.subTest(status=status), patch.object(
                cp_sat_placement,
                "solve",
                return_value=Solution(status, reason="injected status"),
            ):
                result = cp_sat_placement.optimize_local_placement(
                    self.design, ["A", "B"]
                )
                self.assertEqual(result.status, status)
                self.assertEqual(snapshot(self.block), before)

    def test_failed_post_apply_checks_restore_exact_origins_and_orientations(self):
        before = snapshot(self.block)
        original_check = type(self.design.getOpendp()).checkPlacement

        def reject_changed_pose(opendp, *args):
            original_check(opendp, *args)
            if snapshot(self.block) != before:
                raise RuntimeError("injected downstream placement rule failure")

        with patch.object(
            type(self.design.getOpendp()), "checkPlacement", reject_changed_pose
        ):
            result = cp_sat_placement.optimize_local_placement(
                self.design, ["M"], max_displacement_sites=0
            )
        self.assertEqual(result.status, "rejected")
        self.assertEqual(snapshot(self.block), before)
        self.design.getOpendp().checkPlacement(False)

        original_hpwl = cp_sat_placement._actual_hpwl
        calls = 0

        def mismatching_hpwl(nets):
            nonlocal calls
            calls += 1
            return original_hpwl(nets) + int(calls == 2)

        with patch.object(cp_sat_placement, "_actual_hpwl", mismatching_hpwl):
            with self.assertRaisesRegex(ValueError, "HPWL disagrees"):
                cp_sat_placement.optimize_local_placement(
                    self.design, ["M"], max_displacement_sites=0
                )
        self.assertEqual(snapshot(self.block), before)
        self.design.getOpendp().checkPlacement(False)

    def test_unsupported_instance_states_and_routed_nets_fail_before_mutation(self):
        inst = self.block.findInst("A")
        cases = (
            (
                lambda: inst.setPlacementStatus("LOCKED"),
                lambda: inst.setPlacementStatus("PLACED"),
            ),
            (lambda: inst.setDoNotTouch(True), lambda: inst.setDoNotTouch(False)),
            (
                lambda: inst.setPlacementStatus("UNPLACED"),
                lambda: inst.setPlacementStatus("PLACED"),
            ),
        )
        for index, (mutate, restore) in enumerate(cases):
            with self.subTest(state=index):
                mutate()
                before = snapshot(self.block)
                with self.assertRaises(ValueError):
                    cp_sat_placement.optimize_local_placement(self.design, ["A"])
                self.assertEqual(snapshot(self.block), before)
                restore()
        wire = odb.dbWire_create(self.block.findNet("right"))
        before = snapshot(self.block)
        with self.assertRaisesRegex(ValueError, "unrouted"):
            cp_sat_placement.optimize_local_placement(self.design, ["A"])
        self.assertEqual(snapshot(self.block), before)
        odb.dbWire_destroy(wire)


if __name__ == "__main__":
    unittest.main()
