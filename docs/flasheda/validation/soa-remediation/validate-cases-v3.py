import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import shutil
import time

p = argparse.ArgumentParser()
p.add_argument('--baseline', type=Path, required=True)
p.add_argument('--candidate', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--threads',type=int,nargs='+',default=[1,8])
p.add_argument('--baseline-results',type=Path)
p.add_argument('--cases', nargs='+')
p.add_argument('--debug-consistency', action='store_true')
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
cases = ['incremental02', 'incremental03', 'macro01', 'macro02',
         'region01', 'cluster_place01', 'simple01-td-repair',
         'place_ios_td', 'place_ios_rd', 'incremental01']
suffix = '''
set block [ord::get_db_block]
foreach inst [$block getInsts] {
  set box [$inst getBBox]
  puts "SOA_CELL [$inst getName] [$box xMin] [$box yMin] [$inst getOrient] [[$inst getMaster] getName] [$inst getPlacementStatus]"
}
foreach bterm [$block getBTerms] {
  foreach pin [$bterm getBPins] {
    foreach box [$pin getBoxes] {
      puts "SOA_PORT [$bterm getName] [[$box getTechLayer] getName] [$box xMin] [$box yMin] [$box xMax] [$box yMax] [$pin getPlacementStatus]"
    }
  }
}
foreach net [$block getNets] {
  set terms {}
  foreach iterm [$net getITerms] {
    lappend terms "[[$iterm getInst] getName]/[[$iterm getMTerm] getName]"
  }
  set ports {}
  foreach bterm [$net getBTerms] { lappend ports [$bterm getName] }
  puts "SOA_NET [$net getName] [$net getSigType] [lsort $terms] [lsort $ports]"
}
puts SOA_DONE
'''
if a.cases:
    if not set(a.cases) <= set(cases): raise ValueError('Unknown case')
    cases = a.cases
report = {}
for case in cases:
    report[case] = {}
    source = (a.source / 'src/gpl/test' / (case+'.tcl')).read_text()
    # Original Tcl assertions still run. Compare against the current baseline,
    # rather than upstream golden coordinates recorded for a different build.
    source = re.sub(r'^set_thread_count \d+.*$', '', source, flags=re.M)
    source = re.sub(r'^diff_files .*$', '# golden diff replaced by native differential comparison', source, flags=re.M)
    for threads in a.threads:
        outputs = []
        for name, binary in [('baseline', a.baseline), ('soa', a.candidate)]:
            if name == 'baseline' and a.baseline_results and (a.baseline_results / case / str(threads) / name / 'result.json').exists():
                prior = a.baseline_results / case / str(threads) / name
                result = json.loads((prior/'result.json').read_text())
                text = (prior/'run.log').read_text()
                rows = sorted(line for line in text.splitlines() if line.startswith(('SOA_CELL ', 'SOA_PORT ', 'SOA_NET ')))
                digest = hashlib.sha256('\n'.join(rows).encode()).hexdigest()
                if not rows or 'SOA_DONE' not in text or digest != result['hash']:
                    raise RuntimeError(f'Invalid recorded baseline: {prior}')
                cached = a.output / case / str(threads) / name
                cached.mkdir(parents=True)
                for filename in ['run.log', 'result.json']:
                    shutil.copy2(prior/filename, cached/filename)
                outputs.append(dict(result, reused_from=str(prior)))
                continue
            out = a.output / case / str(threads) / name
            out.mkdir(parents=True)
            script = out / 'case.tcl'
            script.write_text(('set_debug_level GPL checkConsistency 1\n' if a.debug_consistency else '') + f'set_thread_count {threads}\nsta::set_thread_count 1\n' + source + suffix)
            env = dict(os.environ, RESULTS_DIR=str(out/'results'), FLASHEDA_GPL_TIMING='0', FLASHEDA_RUDY_MODE='native', FLASHEDA_RUDY_PROFILE='0')
            start = time.monotonic()
            with (out/'run.log').open('w') as log:
                subprocess.run([str(binary), '-no_init', '-exit', str(script)], cwd=a.source/'src/gpl/test', env=env, stdout=log, stderr=subprocess.STDOUT, timeout=1800, check=True)
            text = (out/'run.log').read_text()
            if 'SOA_DONE' not in text or re.search(r'\[ERROR|^Error:', text, re.M):
                raise RuntimeError(f'Case failed: {out}')
            if re.search(r'Using \d+ of \d+ threads for placement', text):
                raise RuntimeError(f'Placement thread cap detected: {out}')
            rows = sorted(line for line in text.splitlines() if line.startswith(('SOA_CELL ', 'SOA_PORT ', 'SOA_NET ')))
            if not rows:
                raise RuntimeError(f'Empty output: {out}')
            digest = hashlib.sha256('\n'.join(rows).encode()).hexdigest()
            result = dict(hash=digest, rows=len(rows), seconds=time.monotonic()-start)
            (out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
            outputs.append(result)
        report[case][str(threads)] = dict(baseline=outputs[0], soa=outputs[1], exact=outputs[0]['hash']==outputs[1]['hash'])
        (a.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
        print(case, threads, report[case][str(threads)], flush=True)
        if not report[case][str(threads)]['exact']:
            raise RuntimeError(f'Placement difference: {case} at {threads} threads; inspect before continuing')
