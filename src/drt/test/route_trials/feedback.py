# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Zero-tolerance endpoint feedback; metrics are separate from acceptance."""

import math


def endpoint_feedback(original, candidate):
    if set(original) != set(candidate):
        raise ValueError("endpoint identity mismatch")
    rows = []
    summary = {
        kind: {
            "applicable": 0,
            "worse": 0,
            "improved": 0,
            "equal": 0,
            "worst_delta_ns": 0.0,
            "total_loss_ns": 0.0,
        }
        for kind in ("hold", "setup")
    }
    for name in sorted(original):
        if len(original[name]) != 2 or len(candidate[name]) != 2:
            raise ValueError("endpoint column count")
        for i, kind in enumerate(("hold", "setup")):
            a, b = float(original[name][i]), float(candidate[name][i])
            if a == math.inf and b == math.inf:
                rows.append(
                    {
                        "endpoint": name,
                        "kind": kind,
                        "applicable": False,
                        "original_ns": None,
                        "candidate_ns": None,
                        "delta_ns": None,
                    }
                )
                continue
            if not math.isfinite(a) or not math.isfinite(b):
                raise ValueError(
                    "timing-check applicability/nonfinite changed: " + name
                )
            delta = b - a
            loss = max(0.0, -delta)
            row = {
                "endpoint": name,
                "kind": kind,
                "applicable": True,
                "original_ns": a,
                "candidate_ns": b,
                "delta_ns": delta,
                "loss_ns": loss,
                "worse": delta < 0,
            }
            rows.append(row)
            s = summary[kind]
            s["applicable"] += 1
            s["worse" if delta < 0 else "improved" if delta > 0 else "equal"] += 1
            s["worst_delta_ns"] = min(s["worst_delta_ns"], delta)
            s["total_loss_ns"] += loss
    return {
        "schema": "joint_endpoint_feedback_v1",
        "tolerance": 0,
        "units": "ns",
        "rows": rows,
        "summary": summary,
        "endpoints_nonworse": all(s["worse"] == 0 for s in summary.values()),
    }


def geometry_key(transaction, digest):
    """Identify final changed-net geometry, ignoring proposals' metadata/order."""
    changes = sorted(
        (net, entry["signature"])
        for net, entry in transaction["expected"].items()
        if entry["signature"] != transaction["original_signatures"][net]
    )
    return digest(
        {"original_def_sha256": transaction["original_def_sha256"], "changes": changes}
    )
