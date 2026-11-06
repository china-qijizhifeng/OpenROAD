"""Attribute GPL wall-clock samples to placement phases and compare runs."""
import argparse
import collections
import json

from frames import load

# Ordered: the first phase whose marker appears anywhere in the stack wins.
PHASES = [
    ('timing_feedback', ['gpl::NesterovPlace::runTimingDriven']),
    ('routability_feedback', ['gpl::NesterovPlace::runRoutability']),
    ('initial_placement', ['gpl::Replace::doInitialPlace']),
    ('density_scatter', ['gpl::BinGrid::updateBinsGCellDensityArea']),
    ('fft', ['gpl::FFT::', 'CpuFftBackend::solve']),
    ('density_field_other', ['updateDensityFieldBin']),
    ('density_gather', ['combineGradientsColumns', 'CpuDensityGradientBackend::getCellGradients',
                        'getDensityGradient', 'fillFillerDensityGradients']),
    ('wirelength_force', ['updateWireLengthForceWA']),
    ('wirelength_gather', ['getAllWireLengthGradientsWA',
                           'getWireLengthGradientWA']),
    ('gradient_combine', ['gpl::NesterovBase::updateGradients']),
    ('hpwl', ['getHpwl']),
    ('coordinates', ['nesterovUpdateCoordinates',
                     'updateGCellDensityCenterLocation',
                     'updateDensityCoordiLayoutInside']),
    ('step_length', ['getStepLength', 'nesterovUpdateStepLength']),
    ('update_iteration', ['updateNextIter']),
    ('nesterov_other', ['gpl::NesterovPlace::doNesterovPlace']),
    ('gpl_other', ['gpl::Replace::']),
]


def phase(frames):
    for name, markers in PHASES:
        if any(m in f for f in frames for m in markers):
            return name
    return None


def leaf_name(frames):
    for f in reversed(frames):
        if f and not f.startswith('['):
            return f
    return frames[-1]


def summarize(path, ms_per_sample):
    by_phase = collections.Counter()
    leaves = collections.defaultdict(collections.Counter)
    for frames, raw_leaf, count in load(path):
        p = phase(frames)
        if p is None:
            continue
        by_phase[p] += count
        # Library leaves are labelled with the deepest application frame.
        tail = frames[-1]
        label = leaf_name(frames)
        if tail.startswith('['):
            label = f'{tail} <- {label}'
        leaves[p][label] += count
    return {p: c * ms_per_sample / 1000 for p, c in by_phase.items()}, leaves


def main():
    p = argparse.ArgumentParser()
    p.add_argument('runs', nargs='+', help='label=folded.txt')
    p.add_argument('--leaves', type=int, default=8)
    p.add_argument('--json', default=None)
    a = p.parse_args()
    results = {}
    for item in a.runs:
        label, path = item.split('=', 1)
        results[label] = summarize(path, 1.0)
    labels = list(results)
    names = [n for n, _ in PHASES]
    print('| phase | ' + ' | '.join(labels) + ' |')
    print('|---|' + '---:|' * len(labels))
    for n in names:
        print(f'| {n} | ' + ' | '.join(f'{results[l][0].get(n, 0):.2f}' for l in labels) + ' |')
    print('| **total** | ' + ' | '.join(f'**{sum(results[l][0].values()):.2f}**' for l in labels) + ' |')
    for n in names:
        print(f'\n### {n}')
        for l in labels:
            top = results[l][1][n].most_common(a.leaves)
            if not top:
                continue
            print(f'- {l}: ' + '; '.join(f'{name} {c / 1000:.2f}s' for name, c in top))
    if a.json:
        json.dump({l: {'phases_s': results[l][0],
                       'leaves_s': {n: {k: v / 1000 for k, v in results[l][1][n].most_common(40)}
                                    for n in names}} for l in labels},
                  open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
