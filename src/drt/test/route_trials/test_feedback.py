# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

import math
import unittest

from feedback import endpoint_feedback, geometry_key
from joint_adapter import digest


class FeedbackTest(unittest.TestCase):
    def test_losses_do_not_cancel_improvements_or_other_check_types(self):
        original = {"first": (0.1, 0.2), "second": (0.2, 0.3)}
        candidate = {"first": (0.09, 0.4), "second": (0.4, 0.29)}
        result = endpoint_feedback(original, candidate)
        self.assertFalse(result["endpoints_nonworse"])
        for kind in ("setup", "hold"):
            self.assertEqual(result["summary"][kind]["worse"], 1)
            self.assertEqual(result["summary"][kind]["improved"], 1)
            self.assertAlmostEqual(result["summary"][kind]["total_loss_ns"], 0.01)

    def test_zero_tolerance_keeps_tiny_regressions(self):
        original = {"endpoint": (0.5, 1.0)}
        candidate = {"endpoint": (0.5, math.nextafter(1.0, 0.0))}
        result = endpoint_feedback(original, candidate)
        self.assertEqual(result["summary"]["setup"]["worse"], 1)
        self.assertFalse(result["endpoints_nonworse"])

    def test_unconstrained_checks_are_distinct_from_zero_loss(self):
        result = endpoint_feedback(
            {"endpoint": (math.inf, 1.0)}, {"endpoint": (math.inf, 1.0)}
        )
        self.assertEqual(result["summary"]["hold"]["applicable"], 0)
        self.assertEqual(result["summary"]["setup"]["equal"], 1)
        hold = next(row for row in result["rows"] if row["kind"] == "hold")
        self.assertIsNone(hold["delta_ns"])
        self.assertFalse(hold["applicable"])
        for bad in (0.0, -math.inf, math.nan):
            with self.subTest(candidate_hold=bad):
                with self.assertRaisesRegex(ValueError, "nonfinite"):
                    endpoint_feedback(
                        {"endpoint": (math.inf, 1.0)}, {"endpoint": (bad, 1.0)}
                    )

    def test_endpoint_identity_and_columns_are_required(self):
        with self.assertRaisesRegex(ValueError, "identity"):
            endpoint_feedback({"first": (1, 2)}, {"second": (1, 2)})
        with self.assertRaisesRegex(ValueError, "column"):
            endpoint_feedback({"first": (1, 2)}, {"first": (1,)})

    def test_geometry_cache_is_bound_to_original_and_actual_changes(self):
        original = {"first": "old-a", "second": "old-b"}
        transaction = {
            "original_def_sha256": "original",
            "original_signatures": original,
            "expected": {
                "first": {"signature": "new-a"},
                "second": {"signature": "old-b"},
            },
        }
        expected = geometry_key(transaction, digest)
        transaction["expected"] = {"first": {"signature": "new-a"}}
        self.assertEqual(geometry_key(transaction, digest), expected)
        transaction["original_def_sha256"] = "different-original"
        self.assertNotEqual(geometry_key(transaction, digest), expected)


if __name__ == "__main__":
    unittest.main()
