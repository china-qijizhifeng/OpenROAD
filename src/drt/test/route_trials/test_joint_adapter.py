# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Portable preparation tests; these do not substitute for native DRC/STA."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

import joint_adapter as adapter

TECH = """LAYER metal2
  TYPE ROUTING ;
  DIRECTION VERTICAL ;
END metal2
LAYER metal3
  TYPE ROUTING ;
  DIRECTION HORIZONTAL ;
END metal3
VIA via23 DEFAULT
  LAYER metal2 ;
  RECT -0.1 -0.1 0.1 0.1 ;
  LAYER cut2 ;
  RECT -0.05 -0.05 0.05 0.05 ;
  LAYER metal3 ;
  RECT -0.1 -0.1 0.1 0.1 ;
END via23
"""
SOURCE = """VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
DESIGN fixture ;
UNITS DISTANCE MICRONS 1000 ;
DIEAREA ( 0 0 ) ( 300 300 ) ;
TRACKS X 0 DO 31 STEP 10 LAYER metal2 ;
TRACKS Y 0 DO 31 STEP 10 LAYER metal2 ;
TRACKS X 0 DO 31 STEP 10 LAYER metal3 ;
TRACKS Y 0 DO 31 STEP 10 LAYER metal3 ;
COMPONENTS 0 ;
END COMPONENTS
NETS 3 ;
- first ( driver Z ) ( sink A ) + USE SIGNAL
  + ROUTED metal2 ( 10 20 0 ) ( 20 20 5 )
  NEW metal2 ( 20 20 ) via23
  NEW metal3 ( 20 20 ) ( 20 50 ) ( 80 50 ) ( 80 20 )
  NEW metal2 ( 80 20 ) via23
  NEW metal2 ( 80 20 ) ( 90 20 ) ;
- second ( another Z ) ( other A ) + USE SIGNAL
  + ROUTED metal2 ( 10 100 ) ( 20 100 )
  NEW metal2 ( 20 100 ) ( 80 100 )
  NEW metal2 ( 80 100 ) ( 90 100 ) ;
- untouched ( unchanged Z ) + USE SIGNAL
  + ROUTED metal2 ( 200 200 ) ( 200 240 ) ;
END NETS
END DESIGN
"""


def edits():
    return [
        {
            "net": "first",
            "ports": [[2, 20, 20], [2, 80, 20]],
            "remove_segments": [
                [3, [20, 20], [20, 50]],
                [3, [20, 50], [80, 50]],
                [3, [80, 50], [80, 20]],
            ],
            "remove_vias": [["via23", [20, 20]], ["via23", [80, 20]]],
            "add_atomic_segments": [[2, [20, 20], [80, 20]]],
            "add_vias": [],
        },
        {
            "net": "second",
            "ports": [[2, 20, 100], [2, 80, 100]],
            "remove_segments": [[2, [20, 100], [80, 100]]],
            "remove_vias": [],
            "add_atomic_segments": [
                [2, [20, 100], [20, 110]],
                [2, [20, 110], [80, 110]],
                [2, [80, 110], [80, 100]],
            ],
            "add_vias": [],
        },
    ]


class AdapterTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "original.def"
        self.source.write_text(SOURCE)
        self.tech = self.root / "tech.lef"
        self.tech.write_text(TECH)
        self.output = self.root / "prepared"

    def proposal(self, change=lambda proposal: None):
        data = {
            "schema": "joint_route_proposal_v1",
            "original_def_sha256": adapter.sha(self.source),
            "edits": edits(),
        }
        change(data)
        data["proposal_sha256"] = adapter.digest(data)
        path = self.root / "proposal.json"
        path.write_text(json.dumps(data))
        return path

    def prepare(self, path=None):
        return adapter.prepare(
            self.source, path or self.proposal(), self.tech, self.output
        )

    def test_joint_budget_allows_one_net_to_grow(self):
        before = adapter.sha(self.source)
        receipt = self.prepare()
        self.assertEqual([row["ports"] for row in receipt["metrics"]], [2, 2])
        self.assertEqual(
            [row["new_length"] - row["old_length"] for row in receipt["metrics"]],
            [-60, 20],
        )
        self.assertEqual(
            sum(row["added_vias"] - row["removed_vias"] for row in receipt["metrics"]),
            -2,
        )
        self.assertEqual(adapter.sha(self.source), before)
        self.assertEqual(
            set(dict(adapter.nets(self.output / "patch.def"))), {"first", "second"}
        )

    def test_later_net_failure_emits_nothing(self):
        path = self.proposal(
            lambda p: p["edits"][1]["ports"].__setitem__(1, [2, 70, 100])
        )
        before = adapter.sha(self.source)
        with self.assertRaisesRegex(ValueError, "boundary"):
            self.prepare(path)
        self.assertFalse(self.output.exists())
        self.assertEqual(adapter.sha(self.source), before)

    def test_later_duplicate_removal_emits_nothing(self):
        def change(p):
            row = p["edits"][1]
            row["remove_segments"].append(copy.deepcopy(row["remove_segments"][0]))

        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.prepare(self.proposal(change))
        self.assertFalse(self.output.exists())

    def test_real_internal_branch_is_required_boundary(self):
        self.source.write_text(
            SOURCE.replace(
                "NEW metal2 ( 80 100 ) ( 90 100 ) ;",
                "NEW metal2 ( 80 100 ) ( 90 100 ) NEW metal2 ( 50 100 ) ( 50 90 ) ;",
            )
        )
        with self.assertRaisesRegex(ValueError, "boundary"):
            self.prepare()
        self.assertFalse(self.output.exists())

        def change(p):
            row = p["edits"][1]
            row["ports"].append([2, 50, 100])
            row["add_atomic_segments"] = [[2, [20, 100], [80, 100]]]

        receipt = self.prepare(self.proposal(change))
        self.assertEqual(receipt["metrics"][1]["ports"], 3)

    def test_retained_extensions_and_logical_connections_preserved(self):
        self.prepare()
        tokens = dict(adapter.nets(self.output / "patch.def"))["first"]
        prefix, segments, _ = adapter.parse(
            tokens, adapter.tech_vias(self.tech), strict=True
        )
        self.assertEqual(prefix, "( driver Z ) ( sink A ) + USE SIGNAL".split())
        self.assertIn((2, (10, 20), (20, 20), 0, 5), segments)

    def test_reversed_extension_roundtrip(self):
        tokens = (
            "( i A ) + USE SIGNAL + ROUTED metal2 ( 10 80 40 ) ( 10 20 0 ) ;".split()
        )
        route = adapter.parse(tokens, {}, strict=True)
        emitted = adapter.emit_net("sample", route, {})
        self.assertEqual(adapter.parse(emitted.split()[2:], {}, strict=True), route)
        self.assertEqual(route[1], [(2, (10, 20), (10, 80), 0, 40)])

    def test_tampering_and_wrong_original_are_rejected(self):
        path = self.proposal()
        proposal = json.loads(path.read_text())
        proposal["edits"][0]["ports"][0][1] += 1
        path.write_text(json.dumps(proposal))
        with self.assertRaisesRegex(ValueError, "hash"):
            self.prepare(path)
        with self.assertRaisesRegex(ValueError, "original"):
            self.prepare(
                self.proposal(lambda p: p.update(original_def_sha256="0" * 64))
            )
        self.assertFalse(self.output.exists())

    def test_joint_wire_growth_is_rejected(self):
        def change(p):
            p["edits"][1]["add_atomic_segments"] = [
                [2, [20, 100], [20, 140]],
                [2, [20, 140], [80, 140]],
                [2, [80, 140], [80, 100]],
            ]

        with self.assertRaisesRegex(ValueError, "Joint wire length"):
            self.prepare(self.proposal(change))
        self.assertFalse(self.output.exists())

    def test_track_and_cycle_rejections(self):
        row = edits()[0]
        row["add_atomic_segments"] = [[2, [20, 21], [80, 21]]]
        with self.assertRaisesRegex(ValueError, "off track"):
            adapter.validate_new_tracks(
                self.source, self.tech, row, adapter.tech_vias(self.tech)
            )
        segments = [
            adapter.canonical_segment(2, a, b)
            for a, b in [
                ((0, 0), (10, 0)),
                ((10, 0), (10, 10)),
                ((0, 10), (10, 10)),
                ((0, 0), (0, 10)),
            ]
        ]
        with self.assertRaisesRegex(ValueError, "cycle"):
            adapter.connected_tree(segments, [], [(2, 0, 0), (2, 10, 10)], {})

    def test_audit_detects_untouched_net_changes(self):
        self.prepare()
        patched = dict(adapter.nets(self.output / "patch.def"))
        original = dict(adapter.nets(self.source))
        actual = self.root / "candidate.def"

        def emit():
            actual.write_text(
                "NETS 3 ;\n"
                + "\n".join(
                    "- " + name + " " + " ".join(patched.get(name, tokens))
                    for name, tokens in original.items()
                )
                + "\nEND NETS\n"
            )

        emit()
        receipt = self.output / "transaction.json"
        self.assertTrue(
            adapter.audit(self.source, actual, self.tech, receipt)["all_equal"]
        )
        original["untouched"] = [
            "250" if token == "240" else token for token in original["untouched"]
        ]
        emit()
        result = adapter.audit(self.source, actual, self.tech, receipt)
        self.assertFalse(result["all_equal"])
        self.assertEqual(result["mismatches"][0]["net"], "untouched")

    def test_existing_output_is_not_overwritten(self):
        self.output.mkdir()
        marker = self.output / "marker"
        marker.write_text("keep")
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.prepare()
        self.assertEqual(marker.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
