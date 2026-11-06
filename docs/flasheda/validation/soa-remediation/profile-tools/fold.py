"""Symbolize sampler records into folded stacks (outermost first)."""
import argparse
import bisect
import collections
import json
import os
from pathlib import Path
import struct
import subprocess

DEPTH = 64
RECORD = 8 * (DEPTH + 1)


def read_maps(path):
    maps = []
    for line in Path(path).read_text().splitlines():
        parts = line.split(maxsplit=5)
        if len(parts) < 6 or 'x' not in parts[1]:
            continue
        start, end = (int(x, 16) for x in parts[0].split('-'))
        maps.append((start, end, int(parts[2], 16), parts[5]))
    bases = {}
    for line in Path(path).read_text().splitlines():
        parts = line.split(maxsplit=5)
        if len(parts) == 6 and int(parts[2], 16) == 0:
            start = int(parts[0].split('-')[0], 16)
            bases.setdefault(parts[5], start)
    return sorted(maps), bases


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--samples', type=Path, required=True)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    maps, bases = read_maps(str(a.samples) + '.maps')
    starts = [m[0] for m in maps]
    binary = os.path.realpath(a.binary)
    data = a.samples.read_bytes()
    stacks = []
    for offset in range(0, len(data) - RECORD + 1, RECORD):
        words = struct.unpack_from(f'<{DEPTH + 1}Q', data, offset)
        stacks.append(words[1:1 + words[0]])
    # Resolve every pc to (object, address-to-symbolize).
    resolved = {}
    wanted = set()
    for stack in stacks:
        for depth, pc in enumerate(stack):
            key = (pc, depth == 0)
            if key in resolved:
                continue
            i = bisect.bisect_right(starts, pc) - 1
            if i < 0 or pc >= maps[i][1]:
                resolved[key] = ('?', None)
                continue
            path = maps[i][3]
            if os.path.realpath(path) == binary:
                address = pc - bases[path] - (0 if depth == 0 else 1)
                resolved[key] = ('main', address)
                wanted.add(address)
            else:
                resolved[key] = (Path(path).name, None)
    ordered = sorted(wanted)
    names = {}
    if ordered:
        text = subprocess.run(
            ['addr2line', '-a', '-f', '-C', '-i', '-e', str(a.binary)],
            input='\n'.join(hex(x) for x in ordered) + '\n',
            capture_output=True, text=True, check=True).stdout.splitlines()
        current = None
        for line in text:
            if line.startswith('0x') and all(c in '0123456789abcdefx' for c in line):
                current = int(line, 16)
                names[current] = []
            elif current is not None:
                names[current].append(line)
    def frames(key):
        kind, address = resolved[key]
        if kind != 'main':
            return [f'[{kind}]']
        lines = names.get(address, [])
        # addr2line -i: innermost first, pairs of (function, file:line).
        pairs = [(lines[i], lines[i + 1] if i + 1 < len(lines) else '?')
                 for i in range(0, len(lines), 2)]
        out = []
        for i, (function, location) in enumerate(reversed(pairs)):
            label = function.replace(';', ':')
            if key[1] and i == len(pairs) - 1:
                label += ' @' + location.rsplit('/', 1)[-1].split(' ')[0]
            out.append(label)
        return out or ['??']
    folded = collections.Counter()
    for stack in stacks:
        parts = []
        for depth in range(len(stack) - 1, -1, -1):
            parts.extend(frames((stack[depth], depth == 0)))
        folded[';'.join(parts)] += 1
    with a.output.open('w') as f:
        for line, count in folded.most_common():
            f.write(f'{line} {count}\n')
    print(json.dumps({'samples': len(stacks), 'unique_stacks': len(folded),
                      'symbolized': len(ordered)}))


if __name__ == '__main__':
    main()
