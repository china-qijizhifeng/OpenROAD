import argparse, importlib.util, json, os, subprocess
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--binary',required=True);p.add_argument('--label',required=True);a=p.parse_args()
task=Path('/gpfs/users/linlizhi/flasheda-density-cpu-20261106')
spec=importlib.util.spec_from_file_location('parser','/gpfs/users/linlizhi/flasheda-gpl-timing-20261105/source/gpl-timing/validate.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
ref=json.loads((task/'results/placement/summary.json').read_text()); report={}
for design in ['gcd','aes_cipher_top']:
 out=task/'fork-results'/a.label/design;out.mkdir(parents=True,exist_ok=False)
 env=dict(os.environ);env.update(FLASHEDA_OPENROAD_SOURCE='/tmp/flasheda-gpl-timing/OpenROAD',FLASHEDA_DESIGN=design,FLASHEDA_OUTPUT=str(out/'placed'),FLASHEDA_RUDY_MODE='native',FLASHEDA_RUDY_PROFILE='0',FLASHEDA_GPL_TIMING='0')
 with (out/'run.log').open('w') as f:subprocess.run([a.binary,'-no_init','-exit','/gpfs/users/linlizhi/flasheda-openroad-rudy-20261105/source/rudy/place.tcl'],env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=1800)
 r=m.parse(out/'run.log',False);expected=ref['designs'][design]['runs'][0]
 for key in ['placement_hash','wns_ns','tns_ns','area_um2','hpwl_um']:
  if r[key]!=expected[key]:raise RuntimeError(f'{design} {key}: {r[key]} != {expected[key]}')
 report[design]=r;(out/'result.json').write_text(json.dumps(r,indent=2)+'\n'); print(design,json.dumps(r),flush=True)
(task/'fork-results'/a.label/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
