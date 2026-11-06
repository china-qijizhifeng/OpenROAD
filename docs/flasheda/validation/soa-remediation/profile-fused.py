"""Run sampled/counted AES global placement for each binary configuration."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

CONFIGS = [('soa-fused-openroad', 1, {}), ('soa-fused-openroad', 8, {})]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--old', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--parser', type=Path, required=True)
    p.add_argument('--tcl', type=Path, required=True)
    a = p.parse_args()
    sys.path.insert(0, str(a.parser.parent))
    import validate
    tool = a.task / 'profile-tools'
    output = a.task / 'fork-results' / 'profile-fused'
    output.mkdir(exist_ok=True)
    original = a.tcl.read_text()
    summary = []
    for name, threads, extra in CONFIGS:
        path = output / f'{name}-t{threads}'
        path.mkdir(exist_ok=False)
        script = path / 'run.tcl'
        script.write_text(original.replace('set_thread_count 1\n',
                                           f'set_thread_count {threads}\n', 1))
        env = dict(os.environ, FLASHEDA_OPENROAD_SOURCE='/tmp/flasheda-gpl-timing/OpenROAD',
                   FLASHEDA_DESIGN='aes_cipher_top', FLASHEDA_OUTPUT=str(path / 'placed'),
                   FLASHEDA_RUDY_MODE='native', FLASHEDA_RUDY_PROFILE='0',
                   FLASHEDA_GPL_TIMING='0', **extra)
        # pcount itself must not load the sampler; env applies it to OpenROAD.
        command = [str(tool / 'pcount'), str(path / 'counters.json'), 'env',
                   f"SAMPLER_OUT={path / 'samples.bin'}",
                   f"LD_PRELOAD={tool / 'sampler.so'}",
                   str(a.old / name), '-no_init', '-exit', str(script)]
        with (path / 'run.log').open('w') as f:
            subprocess.run(command, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=1800)
        log = (path / 'run.log').read_text()
        if re.search(r'Using \d+ of \d+ threads for placement', log):
            raise RuntimeError(f'Placement thread cap detected: {path}')
        result = validate.parse(path / 'run.log', False)
        assert result['placement_hash'] == 'bc141af641d8db860dde7e3cee05047e037fb9b43171f0e8ee4b49d2ee903f92'
        result.update(binary=name, threads=threads, env=extra,
                      counters=json.loads((path / 'counters.json').read_text()))
        (path / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
        result['symbolization'] = 'deferred to BP Dev host addr2line'
        summary.append(result)
        (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(name, threads, result['placement_ms'], result['placement_hash'], flush=True)


if __name__ == '__main__':
    main()
