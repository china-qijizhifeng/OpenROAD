# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors

"""Generic checkpoint-level multi-net, multi-port route transactions.

Prepares a partial DEF for ONE import into a private OpenDB process. Original
checkpoints are immutable. Publication is a separate fail-closed operation;
DRC/connectivity/RCX/STA evidence is never inferred from topology checks.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import re


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(obj):
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def save(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + "\n")


def fail(condition, message):
    if not condition:
        raise ValueError(message)


def tech_vias(path):
    text = Path(path).read_text()
    result = {}
    for m in re.finditer(r"^VIA (\S+)[^\n]*\n(.*?)^END \1\s*$", text, re.M | re.S):
        metals = re.findall(r"\bLAYER metal(\d+)\s*;", m[2])
        if len(metals) == 2:
            result[m[1]] = tuple(map(int, metals))
    fail(bool(result), "No via definitions")
    return result


def nets(path):
    active = False
    current = []
    name = None
    with Path(path).open() as f:
        for line in f:
            t = line.split()
            if not active:
                if t[:1] == ["NETS"]:
                    active = True
                continue
            if t == ["END", "NETS"]:
                fail(name is None, "Unterminated net")
                return
            if not t:
                continue
            if t[0] == "-" and name is None:
                name = t[1]
                current = t[2:]
            else:
                current.extend(t)
            if current[-1:] == [";"]:
                yield name, current
                name = None
                current = []
    fail(False, "Missing END NETS")


def canonical_segment(layer, a, b, ea=None, eb=None):
    a = tuple(a)
    b = tuple(b)
    fail(
        len(a) == len(b) == 2 and a != b and (a[0] == b[0] or a[1] == b[1]),
        "Non-Manhattan or zero segment",
    )
    if b < a:
        a, b, ea, eb = b, a, eb, ea
    return int(layer), a, b, ea, eb


def parse(tokens, vias, strict=False):
    pos = tokens.index("ROUTED")
    fail(tokens[pos - 1] == "+", "Unsupported route prefix")
    if strict:
        for tag in (
            "CLOCK",
            "NONDEFAULTRULE",
            "TAPER",
            "TAPERRULE",
            "FIXED",
            "COVER",
            "MASK",
            "RECT",
            "STYLE",
            "SHIELD",
            "SUBNET",
        ):
            fail(tag not in tokens, "Unsupported protected route: " + tag)
    segments = []
    vv = []
    layer = None
    point = None
    extension = None
    i = pos
    while i < len(tokens):
        t = tokens[i]
        if t in ("ROUTED", "NEW"):
            fail(re.fullmatch(r"metal\d+", tokens[i + 1]), "Unsupported layer")
            layer = int(tokens[i + 1][5:])
            point = None
            extension = None
            i += 2
        elif t == "(":
            j = tokens.index(")", i)
            coords = tokens[i + 1 : j]
            fail(len(coords) in (2, 3), "Unsupported point")
            q = tuple(
                point[k] if coords[k] == "*" and point is not None else int(coords[k])
                for k in (0, 1)
            )
            ex = int(coords[2]) if len(coords) == 3 else None
            if point is not None and q != point:
                segments.append(canonical_segment(layer, point, q, extension, ex))
            point = q
            extension = ex
            i = j + 1
        elif t in vias:
            fail(point is not None and layer in vias[t], "Bad via connection")
            vv.append((t, point))
            lo, hi = vias[t]
            layer = hi if layer == lo else lo
            i += 1
        elif t == ";":
            i += 1
        elif t == "TAPER" and not strict:
            i += 1
        else:
            raise ValueError("Unsupported route token " + t)
    return tokens[: pos - 1], segments, vv


def contains(a, b, p):
    return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[
        1
    ] <= max(a[1], b[1])


def topology(segments, vv, via_defs, extra_points=()):
    points = defaultdict(set)
    for l, x, y in extra_points:
        points[l].add((x, y))
    for l, a, b, *_ in segments:
        points[l].update((a, b))
    for name, p in vv:
        for l in via_defs[name]:
            points[l].add(p)
    for l, a, b, *_ in segments:
        if a[1] != b[1]:
            continue
        for m, c, d, *_ in segments:
            if (
                l == m
                and c[0] == d[0]
                and contains(a, b, (c[0], a[1]))
                and contains(c, d, (c[0], a[1]))
            ):
                points[l].add((c[0], a[1]))
    edges = []
    adj = defaultdict(list)

    def add(a, b, owner):
        idx = len(edges)
        edges.append((a, b, owner))
        adj[a].append(idx)
        adj[b].append(idx)

    for idx, (l, a, b, *_) in enumerate(segments):
        ps = sorted(p for p in points[l] if contains(a, b, p))
        for p, q in zip(ps, ps[1:]):
            add((l, *p), (l, *q), ("s", idx))
    for idx, (name, p) in enumerate(vv):
        lo, hi = via_defs[name]
        add((lo, *p), (hi, *p), ("v", idx))
    return edges, adj


def connected_tree(segments, vv, ports, via_defs):
    edges, adj = topology(segments, vv, via_defs, ports)
    fail(
        len(ports) >= 2 and len(ports) == len(set(ports)), "Need unique multiple ports"
    )
    fail(set(ports) <= set(adj), "Port not on replacement graph")
    seen = {ports[0]}
    todo = [ports[0]]
    while todo:
        u = todo.pop()
        for ei in adj[u]:
            a, b, _ = edges[ei]
            v = b if u == a else a
            if v not in seen:
                seen.add(v)
                todo.append(v)
    fail(seen == set(adj), "Disconnected replacement graph")
    fail(len(edges) == len(adj) - 1, "Replacement graph contains a cycle or duplicate")


def validate_new_tracks(source, tech, edit, via_defs):
    tracks = defaultdict(list)
    die = None
    with Path(source).open() as f:
        for line in f:
            if line.startswith("COMPONENTS "):
                break
            m = re.fullmatch(
                r"TRACKS ([XY]) (\d+) DO (\d+) STEP (\d+) LAYER metal(\d+) ;\s*", line
            )
            if m:
                axis, start, count, step, layer = m.groups()
                tracks[(int(layer), axis)].append(tuple(map(int, (start, count, step))))
            if line.startswith("DIEAREA "):
                values = list(map(int, re.findall(r"-?\d+", line)))
                fail(len(values) == 4, "Only rectangular die bounds supported")
                die = values
    fail(die is not None, "Missing die boundary")
    text = Path(tech).read_text()
    directions = {}
    for m in re.finditer(
        r"^LAYER metal(\d+)\s*\n(.*?)^END metal\1\s*$", text, re.M | re.S
    ):
        direction = re.search(r"DIRECTION (HORIZONTAL|VERTICAL)", m[2])
        if direction:
            directions[int(m[1])] = direction[1]

    def inside(point):
        return (
            all(type(x) is int for x in point)
            and die[0] <= point[0] <= die[2]
            and die[1] <= point[1] <= die[3]
        )

    def on(layer, axis, value):
        return any(
            step > 0
            and value >= start
            and (value - start) % step == 0
            and (value - start) // step < count
            for start, count, step in tracks[(layer, axis)]
        )

    for layer, x, y in edit["add_atomic_segments"]:
        fail(
            inside(x) and inside(y), "New segment outside die or noninteger coordinate"
        )
        horizontal = x[1] == y[1]
        fail(
            on(layer, "Y" if horizontal else "X", x[1] if horizontal else x[0]),
            "New segment off track",
        )
    for name, point in edit["add_vias"]:
        fail(name in via_defs and inside(point), "Unknown via or outside die")
        for layer in via_defs[name]:
            fail(layer in directions, "Missing layer direction")
            vertical = directions[layer] == "VERTICAL"
            fail(
                on(layer, "X" if vertical else "Y", point[0] if vertical else point[1]),
                "New via off preferred metal track",
            )


def prepare_edit(name, tokens, edit, via_defs):
    prefix, segments, vv = parse(tokens, via_defs, strict=True)
    remove_s = [canonical_segment(*s)[:3] for s in edit["remove_segments"]]
    remove_v = [(v, tuple(p)) for v, p in edit["remove_vias"]]
    fail(
        len(remove_s) == len(set(remove_s)) and len(remove_v) == len(set(remove_v)),
        "Duplicate removal",
    )
    segment_indices = []
    via_indices = []
    for request in remove_s:
        matches = [i for i, s in enumerate(segments) if s[:3] == request]
        fail(len(matches) == 1, "Missing/ambiguous whole segment on " + name)
        segment_indices += matches
    for request in remove_v:
        matches = [i for i, v in enumerate(vv) if v == request]
        fail(len(matches) == 1, "Missing/ambiguous whole via on " + name)
        via_indices += matches
    ports = [tuple(p) for p in edit["ports"]]
    fail(2 <= len(ports) <= 64, "Port budget")
    full_edges, adj = topology(segments, vv, via_defs)
    removed = {("s", i) for i in segment_indices} | {("v", i) for i in via_indices}
    boundary = {
        node
        for node, es in adj.items()
        if any(full_edges[e][2] in removed for e in es)
        and any(full_edges[e][2] not in removed for e in es)
    }
    fail(
        set(ports) == boundary,
        "Declared ports differ from actual retained-route boundary on " + name,
    )
    old_s = [segments[i] for i in segment_indices]
    old_v = [vv[i] for i in via_indices]
    connected_tree(old_s, old_v, ports, via_defs)
    new_s = [canonical_segment(*s) for s in edit["add_atomic_segments"]]
    new_v = [(v, tuple(p)) for v, p in edit["add_vias"]]
    fail(0 < len(new_s) <= 2048 and len(new_v) <= 256, "Addition budget")
    fail(all(v in via_defs for v, p in new_v), "Unknown new via")
    connected_tree(new_s, new_v, ports, via_defs)
    # Every retained centerline attachment must be a declared boundary port.
    kept_s = [s for i, s in enumerate(segments) if i not in segment_indices]
    kept_v = [v for i, v in enumerate(vv) if i not in via_indices]
    joined, joined_adj = topology(kept_s + new_s, kept_v + new_v, via_defs)

    def isnew(owner):
        return owner[1] >= (len(kept_s) if owner[0] == "s" else len(kept_v))

    contacts = {
        node
        for node, es in joined_adj.items()
        if any(isnew(joined[e][2]) for e in es)
        and any(not isnew(joined[e][2]) for e in es)
    }
    fail(contacts == set(ports), "Unexpected retained same-net connection on " + name)
    length = lambda ss: sum(abs(a[0] - b[0]) + abs(a[1] - b[1]) for l, a, b, *_ in ss)
    return (prefix, kept_s + new_s, kept_v + new_v), dict(
        net=name,
        ports=len(ports),
        removed_vias=len(old_v),
        added_vias=len(new_v),
        old_length=length(old_s),
        new_length=length(new_s),
    )


def emit_net(name, route, via_defs):
    prefix, segments, vv = route
    paths = []

    def point(p, ext):
        return (
            "( "
            + " ".join(map(str, p))
            + ((" " + str(ext)) if ext is not None else "")
            + " )"
        )

    for l, a, b, ea, eb in segments:
        paths.append(f"metal{l} {point(a,ea)} {point(b,eb)}")
    for v, p in vv:
        paths.append(f"metal{via_defs[v][0]} {point(p,None)} {v}")
    fail(bool(paths), "Empty routing")
    return (
        "- "
        + name
        + " "
        + " ".join(prefix)
        + "\n  + ROUTED "
        + ("\n    NEW ".join(paths))
        + " ;\n"
    )


def route_signature(segments, vv):
    # Union collinear centerlines permits harmless writer merging/splitting.
    intervals = defaultdict(list)
    for l, a, b, *_ in segments:
        h = a[1] == b[1]
        intervals[(l, "H" if h else "V", a[1] if h else a[0])].append(
            (a[0] if h else a[1], b[0] if h else b[1])
        )
    merged = []
    for key, ranges in sorted(intervals.items()):
        union = []
        for a, b in sorted(ranges):
            if union and a <= union[-1][1]:
                union[-1] = (union[-1][0], max(b, union[-1][1]))
            else:
                union.append((a, b))
        merged.extend([*key, a, b] for a, b in union)
    return digest([merged, sorted(vv)])


def prepare(source, proposal_path, tech, output, roundtrip=()):
    source = Path(source)
    output = Path(output)
    fail(not output.exists(), "Transaction directory already exists")
    via_defs = tech_vias(tech)
    original_sha = sha(source)
    proposal = json.loads(Path(proposal_path).read_text()) if proposal_path else None
    if proposal:
        fail(
            proposal["schema"] == "joint_route_proposal_v1",
            "Unsupported proposal schema",
        )
        bound = {
            k: v
            for k, v in proposal.items()
            if k not in ("proposal_sha256", "acceptance")
        }
        fail(
            digest(bound) == proposal["proposal_sha256"],
            "Proposal content hash mismatch",
        )
        fail(proposal["original_def_sha256"] == original_sha, "Wrong original DEF")
        edits = {e["net"]: e for e in proposal["edits"]}
        fail(
            len(edits) == len(proposal["edits"]) and 1 <= len(edits) <= 16,
            "Repeated net or edit budget",
        )
    else:
        edits = {name: None for name in roundtrip}
    fail(bool(edits), "No nets selected")
    all_tokens = dict(nets(source))
    fail(set(edits) <= set(all_tokens), "Missing net")
    patches = []
    metrics = []
    expected = {}
    original_signatures = {}
    for name, edit in edits.items():
        tokens = all_tokens[name]
        before = parse(tokens, via_defs, strict=True)
        original_signatures[name] = route_signature(before[1], before[2])
        if edit is not None:
            validate_new_tracks(source, tech, edit, via_defs)
            route, row = prepare_edit(name, tokens, edit, via_defs)
            metrics.append(row)
        else:
            route = before
        patches.append(emit_net(name, route, via_defs))
        expected[name] = dict(
            signature=route_signature(route[1], route[2]), prefix=route[0]
        )
    if proposal:
        fail(
            sum(r["added_vias"] - r["removed_vias"] for r in metrics) < 0,
            "Joint via budget is not improved",
        )
        fail(
            sum(r["new_length"] - r["old_length"] for r in metrics) <= 0,
            "Joint wire length increased",
        )
    header = []
    with source.open() as f:
        for line in f:
            if re.match(r"^(VERSION|DIVIDERCHAR|BUSBITCHARS|DESIGN|UNITS)\b", line):
                header.append(line)
            if line.startswith("DIEAREA"):
                break
    fail(any(s.startswith("DESIGN ") for s in header), "Missing design header")
    # Validate all edits first. No partial DEF is emitted after a later-net error.
    output.mkdir(parents=True)
    patch = (
        "".join(header)
        + f"NETS {len(patches)} ;\n"
        + "".join(patches)
        + "END NETS\nEND DESIGN\n"
    )
    temporary = output / "patch.def.tmp"
    temporary.write_text(patch)
    os.replace(temporary, output / "patch.def")
    receipt = dict(
        schema="joint_checkpoint_transaction_v1",
        original_def_sha256=original_sha,
        proposal_sha256=proposal["proposal_sha256"] if proposal else None,
        patch_sha256=sha(output / "patch.def"),
        tech_sha256=sha(tech),
        expected=expected,
        original_signatures=original_signatures,
        metrics=metrics,
        changed_nets=sorted(edits),
        status="PREPARED_NOT_ACCEPTED",
    )
    save(output / "transaction.json", receipt)
    fail(sha(source) == original_sha, "Original checkpoint changed during preparation")
    return receipt


def audit(source, actual, tech, transaction):
    rec = json.loads(Path(transaction).read_text())
    fail(sha(source) == rec["original_def_sha256"], "Original identity changed")
    via_defs = tech_vias(tech)
    a = dict(nets(source))
    b = dict(nets(actual))
    fail(a.keys() == b.keys(), "Net identity changed")
    mismatches = []
    for name, tokens in a.items():
        if name not in rec["expected"]:
            if tokens != b[name]:
                mismatches.append(
                    dict(net=name, reason="Unedited net DEF tokens changed")
                )
            continue
        prefix, ss, vv = parse(b[name], via_defs, strict=True)
        expected = rec["expected"][name]
        if (
            prefix != expected["prefix"]
            or route_signature(ss, vv) != expected["signature"]
        ):
            mismatches.append(
                dict(
                    net=name,
                    reason="Edited net differs from exact proposed centerline/vias or original connections",
                )
            )
    return dict(
        all_equal=not mismatches,
        mismatches=mismatches,
        note="Centerline/via fidelity and untouched-net token identity. Physical metal and RC require native checks.",
        original_unchanged=sha(source) == rec["original_def_sha256"],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("tech", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--proposal", type=Path)
    p.add_argument("--roundtrip", nargs="*", default=[])
    a = p.parse_args()
    print(
        json.dumps(
            prepare(a.source, a.proposal, a.tech, a.output, a.roundtrip), indent=2
        )
    )
