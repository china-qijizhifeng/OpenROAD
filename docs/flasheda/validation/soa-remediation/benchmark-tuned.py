import argparse
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path
import statistics
import subprocess

p=argparse.ArgumentParser()
p.add_argument('--task', type=Path, required=True)
a=p.parse_args(); task=a.task
spec=importlib.util.spec_from_file_location('parser','/gpfs/users/linlizhi/flasheda-gpl-timing-20261105/source/gpl-timing/validate.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
source=Path('/tmp/flasheda-gpl-timing/OpenROAD')
original=Path('/gpfs/users/linlizhi/flasheda-openroad-rudy-20261105/source/rudy/place.tcl').read_text()
output=task/'fork-results/benchmark-tuned';output.mkdir(parents=True,exist_ok=False)
report={'runs':[]}
for repetition in range(3):
    # Alternate execution order to reduce drift bias.
    names=['soa-final-openroad','tuned-openroad']
    names=names[repetition:]+names[:repetition]
    for threads in [1,8]:
        for name in names:
            path=output/f'{repetition}-{threads}-{name}';path.mkdir()
            script=path/'run.tcl';script.write_text(original.replace('set_thread_count 1\n',f'set_thread_count {threads}\n',1))
            env=dict(os.environ, FLASHEDA_OPENROAD_SOURCE=str(source),FLASHEDA_DESIGN='aes_cipher_top',FLASHEDA_OUTPUT=str(path/'placed'),FLASHEDA_RUDY_MODE='native',FLASHEDA_RUDY_PROFILE='0',FLASHEDA_GPL_TIMING='0')
            with (path/'run.log').open('w') as f:
                command=[str(task/name),'-no_init','-exit',str(script)]
                if Path('/usr/bin/time').exists():
                    command=['/usr/bin/time','-v','-o',str(path/'resources.txt'),*command]
                subprocess.run(command,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=1800)
            if re.search(r'Using \d+ of \d+ threads for placement',(path/'run.log').read_text()):
                raise RuntimeError(f'Placement thread cap detected: {path}')
            r=m.parse(path/'run.log',False);r.update(binary=name,threads=threads,repetition=repetition)
            if (path/'resources.txt').exists():
                match=re.search(r'Maximum resident set size \(kbytes\): (\d+)',(path/'resources.txt').read_text())
                if match:r['max_rss_kb']=int(match[1])
            peers=[x for x in report['runs'] if x['threads']==threads]
            if peers:
                for key in ['placement_hash','wns_ns','tns_ns','area_um2','hpwl_um']:
                    if r[key]!=peers[0][key]:
                        (path/'result.json').write_text(json.dumps(r,indent=2)+'\n')
                        raise RuntimeError(f'Nonidentical {key} for {path}; inspect before continuing')
            report['runs'].append(r)
            (path/'result.json').write_text(json.dumps(r,indent=2)+'\n')
            (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
            print(name,threads,repetition,r['placement_ms'],r['placement_hash'],flush=True)
report['summary']={}
for threads in [1,8]:
    for name in ['soa-final-openroad','tuned-openroad']:
        values=[r['placement_ms'] for r in report['runs'] if r['threads']==threads and r['binary']==name]
        report['summary'][f'{name}-{threads}']=dict(median_ms=statistics.median(values),min_ms=min(values),max_ms=max(values),samples=len(values))
report['binary_sha256']={name:hashlib.sha256((task/name).read_bytes()).hexdigest() for name in ['soa-final-openroad','tuned-openroad']}
report['complete']=True
(output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
