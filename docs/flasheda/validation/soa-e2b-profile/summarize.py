#!/usr/bin/env python3
"""Convert the preserved E2B JSON receipts to tables; no benchmark execution."""
import argparse
import csv
import io
import json
from pathlib import Path
import re

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--input', type=Path,
               default=Path(__file__).resolve().parents[4] / 'build/soa/e2b-profile')
p.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
p.add_argument('--require-complete', action='store_true')
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)


def write(name, fields, rows):
    with (a.output / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


phase_runs, phase_detail = [], []
for design in ['aes', 'large01', 'large02']:
    data = json.loads((a.input / f'{design}-summary.json').read_text())
    assert data['complete'] and len(data['runs']) == 4
    for run in data['runs']:
        phase = run['phase']
        totals = phase['totals']
        ompt = run['ompt']
        phase_runs.append({
            'design': run['design'], 'role': run['role'], 'threads': run['threads'],
            'placement_root_s': phase['root_ms'] / 1000,
            'initial_placement_s': totals.get('initial_placement', {}).get('inclusive_ms', 0) / 1000,
            'nesterov_s': totals.get('nesterov', {}).get('inclusive_ms', 0) / 1000,
            'attempts': phase['iteration_attempts'],
            'final_iterations': '/'.join(run['algorithm']['final_iterations']),
            'output_exact': run['profile_output_exact'],
            'algorithm_exact': run['profile_algorithm_trace_exact'],
            'ompt_measured': ompt['measured'],
            'ompt_wait_thread_s': sum(x['thread_wait_ms'] for x in ompt['waits']) / 1000,
            'ompt_initial_thread_wait_s': sum(x['initial_thread_wait_ms'] for x in ompt['waits']) / 1000,
            'ompt_invalid_pairs': sum(x['invalid_pairs'] for x in ompt['finalization']),
            'ompt_parallel_regions': sum(x['parallel_regions'] for x in ompt['finalization']),
        })
        waits = {x['phase']: x for x in ompt['waits']}
        for name, values in totals.items():
            phase_detail.append({
                'design': run['design'], 'role': run['role'], 'threads': run['threads'],
                'phase': name, 'inclusive_s': values['inclusive_ms'] / 1000,
                'exclusive_s': values['exclusive_ms'] / 1000, 'calls': values['calls'],
                'inclusive_us_per_call': values['inclusive_ms'] * 1000 / values['calls'],
                'exclusive_us_per_call': values['exclusive_ms'] * 1000 / values['calls'],
                'ompt_wait_thread_s': waits.get(name, {}).get('thread_wait_ms', 0) / 1000,
                'ompt_initial_thread_wait_s': waits.get(name, {}).get('initial_thread_wait_ms', 0) / 1000,
            })
write('phase-runs.csv', list(phase_runs[0]), phase_runs)
write('phase-detail.csv', list(phase_detail[0]), phase_detail)


def resource(text, name):
    found = re.search(r'^\s*' + re.escape(name) + r':\s*(.+)$', text, re.M)
    return found.group(1) if found else ''


perf_runs, perf_events, perf_top = [], [], []
for design in ['aes', 'large02']:
    path = a.input / f'{design}-perf-summary.json'
    if not path.exists():
        assert not a.require_complete, path
        continue
    data = json.loads(path.read_text())
    if a.require_complete:
        assert data['complete'] and len(data['runs']) == 4
    for run in data['runs']:
        text = run['process_resource']
        user = resource(text, 'User time (seconds)')
        system = resource(text, 'System time (seconds)')
        wall = resource(text, 'Elapsed (wall clock) time (h:mm:ss or m:ss)')
        seconds = sum(float(v) * 60 ** i for i, v in enumerate(reversed(wall.split(':'))))
        rss = resource(text, 'Maximum resident set size (kbytes)')
        perf_runs.append({
            'design': data['design'], 'role': run['role'], 'threads': run['threads'],
            'mode': run['mode'], 'observed_wall_s': seconds,
            'user_s': user, 'sys_s': system,
            'average_occupied_cores': (float(user) + float(system)) / seconds,
            'max_rss_bytes': int(rss) * 1024,
            'stage_placement_ms': '/'.join(run['stage_placement_ms']),
            'output_exact': run['output_exact'], 'algorithm_exact': run['algorithm_exact'],
            'samples': '/'.join(run.get('samples', [])),
            'lost_samples': '/'.join(run.get('lost_samples', [])),
            'gpl_symbols_visible': run.get('gpl_symbols_visible', ''),
            'unknown_symbol_percent_reported_rows': run.get('unknown_symbol_percent_reported_rows', ''),
        })
        for row in csv.reader(io.StringIO(run.get('stat_csv', '')), delimiter=';'):
            if not row or row[0].startswith('#'):
                continue
            assert len(row) >= 5, row
            perf_events.append({
                'design': data['design'], 'role': run['role'], 'threads': run['threads'],
                'event': row[2], 'value': row[0], 'unit': row[1],
                'event_runtime_ns': row[3], 'percent_running': row[4],
                'metric': row[5] if len(row) > 5 else '',
                'metric_unit': row[6] if len(row) > 6 else '',
            })
        flat = run.get('flat_report', '')
        for percent, dso, privilege, symbol in re.findall(
                r'^\s*([\d.]+)%\s+(\S+)\s+\[(.)\]\s+(.*)$', flat, re.M):
            symbol = re.sub(r'\s+-\s+-\s*$', '', symbol).strip()
            perf_top.append({
                'design': data['design'], 'role': run['role'],
                'self_percent': percent, 'dso': dso, 'privilege': privilege, 'symbol': symbol,
            })
if perf_runs:
    write('perf-runs.csv', list(perf_runs[0]), perf_runs)
if perf_events:
    write('perf-events.csv', list(perf_events[0]), perf_events)
if perf_top:
    write('perf-symbols.csv', list(perf_top[0]), perf_top)
print(f'{len(phase_runs)} phase runs; {len(perf_runs)} perf runs; {len(perf_events)} event rows; {len(perf_top)} symbol rows')
