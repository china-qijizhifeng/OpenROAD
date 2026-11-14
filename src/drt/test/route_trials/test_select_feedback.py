# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Read-only selector contract tests; test fixtures stay in diagnosis and use no native results."""

import copy, json, tempfile, unittest
from pathlib import Path
import select_feedback as s


class SelectorTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.pool = [
            {
                "geometry_sha256": f"{i:064x}",
                "features": [i % 5, (i * 7) % 11, 0.0],
                "routes": [[i, i + 1], [100 + i % 7]],
                "metrics": {
                    "vias": 4 + i % 3,
                    "geometry_deviation_dbu": i * 10,
                    "length_dbu": 100 + i,
                },
                "proposal": str(self.root / f"p{i}.json"),
                "proposal_file_sha256": f"{1000+i:064x}",
            }
            for i in range(25)
        ]
        self.path = self.root / "pool.json"
        self.path.write_text(json.dumps({"candidates": self.pool}))
        self.plan = s.make_plan({"regions": [{"region": "R", "pool": str(self.path)}]})
        self.reg = self.plan["regions"][0]
        self.observations = []
        for i, cid in enumerate(self.reg["shared_ids"]):
            self.observations.append(
                self.obs(cid, "shared", 0, setup=(i + 1) * 0.001, hold=(4 - i) * 0.0001)
            )

    def tearDown(self):
        self.tmp.cleanup()

    def obs(self, cid, arm, round, setup=0.001, hold=0.0001):
        return {
            "region": "R",
            "candidate_id": cid,
            "arm": arm,
            "round": round,
            "physical_pass": True,
            "timing_complete": True,
            "setup_max_loss_ns": setup,
            "setup_sum_loss_ns": setup * 3,
            "setup_loss_count": 3,
            "hold_max_loss_ns": hold,
            "hold_sum_loss_ns": hold * 2,
            "hold_loss_count": 2,
        }

    def test_control_and_future_observation_isolation(self):
        expected = s.next_choice(self.plan, self.observations, "R", 1)
        remaining = [
            r["geometry_sha256"]
            for r in self.pool
            if r["geometry_sha256"] not in self.reg["shared_ids"]
        ]
        poison = [
            self.obs(remaining[0], "control", 0, setup=1e6, hold=1e6),
            self.obs(remaining[1], "feedback", 1, setup=0, hold=0),
            self.obs(remaining[2], "feedback", 9, setup=0, hold=0),
        ]
        self.assertEqual(
            expected, s.next_choice(self.plan, self.observations + poison, "R", 1)
        )
        poison[0].update(region="ANOTHER_REGION", arm="shared")
        self.assertEqual(
            expected, s.next_choice(self.plan, self.observations + poison, "R", 1)
        )

    def test_missing_or_negative_labels_never_become_zero(self):
        obs = copy.deepcopy(self.observations)
        for o in obs:
            o["timing_complete"] = False
            o["setup_max_loss_ns"] = None
        r = s.next_choice(self.plan, obs, "R", 1)
        self.assertEqual(r["timing_training_ids"], [])
        self.assertEqual(r["mode"], "no_timing_labels_exploration")
        self.assertTrue(all(p["predicted"] is None for p in r["predictions"]))
        obs = copy.deepcopy(self.observations)
        obs[0]["setup_max_loss_ns"] = -1
        obs[1]["hold_sum_loss_ns"] = float("nan")
        r = s.next_choice(self.plan, obs, "R", 1)
        self.assertNotIn(obs[0]["candidate_id"], r["timing_training_ids"])
        self.assertNotIn(obs[1]["candidate_id"], r["timing_training_ids"])

    def test_route_distance_distinguishes_equal_features(self):
        a = copy.deepcopy(self.pool[0])
        b = copy.deepcopy(a)
        b["routes"] = [[900, 901], [800]]
        self.assertGreater(s.distance(a, b, [1, 1, 1]), 0)
        self.assertEqual(s.distance(a, a, [1, 1, 1]), 0)

    def test_every_third_round_explores_and_no_repeated_feedback(self):
        obs = copy.deepcopy(self.observations)
        chosen = []
        for round in range(1, 11):
            r = s.next_choice(self.plan, obs, "R", round)
            self.assertNotIn(r["candidate_id"], chosen + self.reg["shared_ids"])
            chosen.append(r["candidate_id"])
            self.assertEqual(r["mode"] == "fixed_exploration", round % 3 == 0)
            self.assertEqual(
                set(r["training_candidate_ids"]), {o["candidate_id"] for o in obs}
            )
            obs.append(
                self.obs(
                    r["candidate_id"],
                    "feedback",
                    round,
                    setup=0.0002 * round,
                    hold=0.00002 * (11 - round),
                )
            )
        self.assertEqual(len(set(chosen)), 10)

    def test_setup_hold_predictions_remain_separate_positive_dimensions(self):
        obs = copy.deepcopy(self.observations)
        for i, o in enumerate(obs):
            o.update(
                setup_max_loss_ns=0.001 * (i + 1), hold_max_loss_ns=0.0001 * (4 - i)
            )
        r = s.next_choice(self.plan, obs, "R", 1)
        for p in r["predictions"]:
            self.assertGreaterEqual(p["predicted"]["setup_max_loss_ns"], 0.001)
            self.assertGreaterEqual(p["predicted"]["hold_max_loss_ns"], 0.0001)
            self.assertLessEqual(p["predicted"]["setup_max_loss_ns"], 0.004)
            self.assertLessEqual(p["predicted"]["hold_max_loss_ns"], 0.0004)

    def test_frozen_pool_tamper_and_duplicate_observation_rejected(self):
        with self.assertRaises(AssertionError):
            s.next_choice(self.plan, self.observations + [self.observations[0]], "R", 1)
        self.path.write_text(self.path.read_text() + " ")
        with self.assertRaises(AssertionError):
            s.next_choice(self.plan, self.observations, "R", 1)

    def test_shared_and_arm_budget_single_and_two_region(self):
        self.assertEqual(len(self.reg["shared_ids"]), 4)
        self.assertEqual(len(self.reg["control_ids"]), 10)
        self.assertEqual(self.reg["feedback_rounds"], 10)
        self.assertFalse(set(self.reg["shared_ids"]) & set(self.reg["control_ids"]))
        plan = s.make_plan(
            {
                "regions": [
                    {"region": "R", "pool": str(self.path)},
                    {"region": "S", "pool": str(self.path)},
                ]
            }
        )
        self.assertEqual(sum(len(r["shared_ids"]) for r in plan["regions"]), 4)
        self.assertEqual(sum(len(r["control_ids"]) for r in plan["regions"]), 10)
        self.assertEqual(sum(r["feedback_rounds"] for r in plan["regions"]), 10)
        self.assertEqual(plan["logical_calls_per_arm"], 14)
        self.assertEqual(plan["native_distinct_cap"], 24)

    def test_insufficient_distinct_pool_is_rejected(self):
        self.path.write_text(json.dumps({"candidates": self.pool[:13]}))
        with self.assertRaises(ValueError):
            s.make_plan({"regions": [{"region": "R", "pool": str(self.path)}]})
        duplicated = self.pool[:13] + [copy.deepcopy(self.pool[0])] * 20
        self.path.write_text(json.dumps({"candidates": duplicated}))
        with self.assertRaises(ValueError):
            s.make_plan({"regions": [{"region": "R", "pool": str(self.path)}]})

    def test_missing_past_feedback_round_rejected(self):
        with self.assertRaises(AssertionError):
            s.next_choice(self.plan, self.observations, "R", 2)
        with self.assertRaises(AssertionError):
            s.next_choice(self.plan, self.observations[1:], "R", 1)

    def test_pool_path_relative_to_config_and_proposal_relative_to_pool(self):
        self.path.write_text(
            json.dumps(
                {
                    "candidates": [
                        dict(c, proposal=Path(c["proposal"]).name) for c in self.pool
                    ]
                }
            )
        )
        plan = s.make_plan(
            {"regions": [{"region": "R", "pool": "pool.json"}]}, self.root
        )
        self.assertEqual(plan["regions"][0]["pool"], str(self.path))
        choice = s.next_choice(plan, self.observations, "R", 1)
        self.assertEqual(Path(choice["proposal"]).parent, self.root)


if __name__ == "__main__":
    unittest.main()
