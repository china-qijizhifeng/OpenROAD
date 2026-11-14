# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Frozen sequential route-distance kNN ablation. Never adopts a layout."""

import argparse, hashlib, json, math
from pathlib import Path

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
LABELS = (
    "setup_max_loss_ns",
    "setup_sum_loss_ns",
    "setup_loss_count",
    "hold_max_loss_ns",
    "hold_sum_loss_ns",
    "hold_loss_count",
)


def geometry_rank(c):
    return (
        c["metrics"]["vias"],
        c["metrics"]["geometry_deviation_dbu"],
        c["metrics"]["length_dbu"],
        c["geometry_sha256"],
    )


def stds(pool):
    n = len(pool)
    width = len(pool[0]["features"])
    mean = [sum(c["features"][i] for c in pool) / n for i in range(width)]
    return [
        max(1e-12, math.sqrt(sum((c["features"][i] - mean[i]) ** 2 for c in pool) / n))
        for i in range(width)
    ]


def distance(a, b, scales):
    rms = math.sqrt(
        sum(((x - y) / s) ** 2 for x, y, s in zip(a["features"], b["features"], scales))
        / len(scales)
    )
    bounded = rms / (1 + rms)
    jaccard = sum(
        1 - len(set(x) & set(y)) / max(1, len(set(x) | set(y)))
        for x, y in zip(a["routes"], b["routes"])
    ) / len(a["routes"])
    return 0.5 * bounded + 0.5 * jaccard


def farthest(pool, seen, scales):
    if not seen:
        return min(pool, key=geometry_rank)
    return min(
        pool,
        key=lambda c: (-min(distance(c, x, scales) for x in seen), geometry_rank(c)),
    )


def sequence(pool, scales, count, initial=()):
    selected = list(initial)
    result = []
    for _ in range(count):
        available = [
            c
            for c in pool
            if c["geometry_sha256"] not in {x["geometry_sha256"] for x in selected}
        ]
        if not available:
            raise ValueError("Pool insufficient for fixed shared and per-arm budget")
        chosen = farthest(available, selected, scales)
        result.append(chosen)
        selected.append(chosen)
    return result


def make_plan(config, base_dir=None):
    base_dir = Path.cwd() if base_dir is None else Path(base_dir)
    regions = []
    n = len(config["regions"])
    assert n in (1, 2)
    shared = 4 // n
    follow = 10 // n
    for r in config["regions"]:
        path = (base_dir / r["pool"]).resolve()
        m = json.loads(path.read_text())
        pool = m["candidates"]
        scales = stds(pool)
        initial = sequence(pool, scales, shared)
        control = sequence(pool, scales, follow, initial)
        regions.append(
            {
                "region": r["region"],
                "pool": str(path),
                "pool_sha256": sha(path),
                "candidates": len(pool),
                "feature_scales": scales,
                "shared_ids": [c["geometry_sha256"] for c in initial],
                "control_ids": [c["geometry_sha256"] for c in control],
                "feedback_rounds": follow,
            }
        )
    return {
        "selector_sha256": sha(Path(__file__)),
        "regions": regions,
        "native_distinct_cap": 24,
        "logical_calls_per_arm": 14,
        "old_native_training": False,
        "labels": LABELS,
        "model": "KNN3 inverse(distance+1e-6);distance=.5 mean per-net edge Jaccard+.5 bounded RMS standardized geometry; scale from frozen pool only",
        "selection": "Every third local round exploration=max distance to observed;other rounds Pareto predicted setup/hold maxloss, minimum max of separately scaled positive maxloss;sum/count loss only tie-break. Physical failure modeled separately.",
        "feedback_visibility": "Shared observations plus previous feedback observations only. No control-only or future observations. Missing timing never becomes zero.",
        "final_acceptance": "Independent original strict zero-tolerance all-metric and all-endpoint gates, complete physical and converged power evidence; selection never authorizes acceptance.",
    }


def dominates(a, b):
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def next_choice(plan, observations, region, round_number):
    r = next(r for r in plan["regions"] if r["region"] == region)
    assert 1 <= round_number <= r["feedback_rounds"]
    assert sha(r["pool"]) == r["pool_sha256"]
    m = json.loads(Path(r["pool"]).read_text())
    pool = m["candidates"]
    byid = {c["geometry_sha256"]: c for c in pool}
    scales = r["feature_scales"]
    seen = [
        o
        for o in observations
        if o["region"] == region
        and (
            o["arm"] == "shared" or o["arm"] == "feedback" and o["round"] < round_number
        )
    ]
    assert {o["candidate_id"] for o in seen if o["arm"] == "shared"} == set(
        r["shared_ids"]
    ), "Need exact completed shared seeds"
    assert {o["round"] for o in seen if o["arm"] == "feedback"} == set(
        range(1, round_number)
    ), "Need all previous rounds"
    assert len(seen) == len({o["candidate_id"] for o in seen}), "Duplicate observation"
    for o in seen:
        assert o["candidate_id"] in byid
    available = [
        c for c in pool if c["geometry_sha256"] not in {o["candidate_id"] for o in seen}
    ]
    seen_candidates = [byid[o["candidate_id"]] for o in seen]
    labels = [
        o
        for o in seen
        if o.get("timing_complete") is True
        and all(
            o.get(k) is not None and math.isfinite(o[k]) and o[k] >= 0 for k in LABELS
        )
    ]
    loss_scales = {
        k: max(1e-12, max((o[k] for o in labels), default=0)) for k in LABELS
    }
    predictions = []
    for c in available:
        nearest = sorted(
            labels,
            key=lambda o: (
                distance(c, byid[o["candidate_id"]], scales),
                o["candidate_id"],
            ),
        )[:3]
        weights = [
            1 / (distance(c, byid[o["candidate_id"]], scales) + 1e-6) for o in nearest
        ]
        predicted = (
            {
                k: sum(w * o[k] for o, w in zip(nearest, weights)) / sum(weights)
                for k in LABELS
            }
            if nearest
            else None
        )
        physically_labeled = [
            o for o in seen if o.get("physical_pass") in (True, False)
        ]
        physical_nearest = sorted(
            physically_labeled,
            key=lambda o: (
                distance(c, byid[o["candidate_id"]], scales),
                o["candidate_id"],
            ),
        )[:3]
        pw = [
            1 / (distance(c, byid[o["candidate_id"]], scales) + 1e-6)
            for o in physical_nearest
        ]
        failure = (
            sum(w * (not o["physical_pass"]) for o, w in zip(physical_nearest, pw))
            / sum(pw)
            if pw
            else None
        )
        predictions.append(
            {
                "candidate_id": c["geometry_sha256"],
                "predicted": predicted,
                "physical_failure_probability": failure,
                "nearest_observations": [o["candidate_id"] for o in nearest],
                "distance_to_seen": min(
                    distance(c, x, scales) for x in seen_candidates
                ),
            }
        )
    exploration = round_number % 3 == 0 or not labels
    if exploration:
        chosen = farthest(available, seen_candidates, scales)
        mode = "fixed_exploration" if labels else "no_timing_labels_exploration"
    else:

        def lossvec(p):
            return [
                p["predicted"]["setup_max_loss_ns"],
                p["predicted"]["hold_max_loss_ns"],
            ]

        pareto = [
            p
            for p in predictions
            if not any(dominates(lossvec(q), lossvec(p)) for q in predictions)
        ]

        def rank(p):
            pred = p["predicted"]
            a = pred["setup_max_loss_ns"] / loss_scales["setup_max_loss_ns"]
            b = pred["hold_max_loss_ns"] / loss_scales["hold_max_loss_ns"]
            return (
                p["physical_failure_probability"] or 0,
                max(a, b),
                a + b,
                pred["setup_sum_loss_ns"] / loss_scales["setup_sum_loss_ns"]
                + pred["hold_sum_loss_ns"] / loss_scales["hold_sum_loss_ns"],
                pred["setup_loss_count"] + pred["hold_loss_count"],
                geometry_rank(byid[p["candidate_id"]]),
            )

        selected = min(pareto, key=rank)
        chosen = byid[selected["candidate_id"]]
        mode = "pareto_knn_minimax"
    return {
        "region": region,
        "round": round_number,
        "arm": "feedback",
        "candidate_id": chosen["geometry_sha256"],
        "proposal": str((Path(r["pool"]).parent / chosen["proposal"]).resolve()),
        "proposal_file_sha256": chosen["proposal_file_sha256"],
        "mode": mode,
        "training_candidate_ids": [o["candidate_id"] for o in seen],
        "timing_training_ids": [o["candidate_id"] for o in labels],
        "loss_scales": loss_scales,
        "predictions": predictions,
        "selector_sha256": sha(Path(__file__)),
    }


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("freeze")
    p.add_argument("config", type=Path)
    p.add_argument("output", type=Path)
    p = sub.add_parser("next")
    p.add_argument("plan", type=Path)
    p.add_argument("observations", type=Path)
    p.add_argument("region")
    p.add_argument("round", type=int)
    p.add_argument("output", type=Path)
    a = ap.parse_args()
    assert not a.output.exists()
    if a.command == "freeze":
        result = make_plan(json.loads(a.config.read_text()), a.config.resolve().parent)
        result["config_sha256"] = sha(a.config)
    else:
        plan = json.loads(a.plan.read_text())
        assert plan["selector_sha256"] == sha(Path(__file__)), "Frozen selector changed"
        obs = json.loads(a.observations.read_text())
        result = next_choice(plan, obs["observations"], a.region, a.round)
        result["plan_sha256"] = sha(a.plan)
        result["observations_sha256"] = sha(a.observations)
    a.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in result.items()
                if k
                in ("candidate_id", "region", "round", "mode", "native_distinct_cap")
            }
        )
    )


if __name__ == "__main__":
    main()
