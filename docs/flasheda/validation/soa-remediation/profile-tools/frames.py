"""Shared frame-name simplification for folded stacks."""
import re


def simplify(frame):
    frame = frame.split(' @')[0]
    frame = re.sub(r' \[clone [^\]]*\]', '', frame)
    # Drop parameter lists and template arguments, keep qualified names.
    out, depth = [], 0
    for ch in frame:
        if ch in '(<':
            depth += 1
        elif ch in ')>':
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return ''.join(out).strip()


def load(path):
    stacks = []
    for line in open(path):
        stack, count = line.rstrip('\n').rsplit(' ', 1)
        raw = stack.split(';')
        stacks.append(([simplify(f) for f in raw], raw[-1], int(count)))
    return stacks
