"""Fail-fast native differential acceptance, then three alternating repetitions."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
from statistics import median
import subprocess
import time

CONFIGS = [
 ('large01','gpl','large01','','0'),
 ('large02','gpl','large02','','0'),
 ('large01-rd','gpl','large01','-routability_driven','1'),
]


def parse(path,kind):
 text=path.read_text()
 if re.search(r'^\[ERROR|^Error:',text,re.M) or 'LARGE_DONE 1' not in text:
  raise RuntimeError(f'Error or incomplete state in {path}')
 if re.search(r'Using \d+ of \d+ threads for placement',text):
  raise RuntimeError(f'Placement thread cap in {path}')
 cells=sorted(l for l in text.splitlines() if l.startswith('LARGE_CELL '))
 rows=sorted(l for l in text.splitlines() if l.startswith(('LARGE_CELL ','LARGE_PORT ','LARGE_NET ')))
 if not cells:raise RuntimeError(f'Empty placement {path}')
 r={'placement_ms':float(re.search(r'FLASHEDA_STAGE placement_ms ([0-9.]+)',text)[1]),
    'placement_hash':hashlib.sha256('\n'.join(cells).encode()).hexdigest(),
    'state_hash':hashlib.sha256('\n'.join(rows).encode()).hexdigest(),
    'cells':len(cells),'state_rows':len(rows),
    'area_um2':float(re.search(r'LARGE_AREA (\S+)',text)[1]),
    'hpwl_um':float(re.search(r'LARGE_HPWL (\S+)',text)[1])}
 if kind=='flow':
  q=re.search(r'FLASHEDA_QUALITY (\S+) (\S+)',text)
  if not q or 'FLASHEDA_LEGAL 1' not in text:raise RuntimeError(f'Missing flow quality/legality {path}')
  r.update(wns_ns=float(q[1]),tns_ns=float(q[2]),legal=True)
 for name,pattern in [('instances',r'GPL-0006\] Number of instances:\s+(\d+)'),
                      ('pins',r'GPL-0011\] Number of pins:\s+(\d+)'),
                      ('nets',r'GPL-0010\] Number of nets:\s+(\d+)'),
                      ('utilization_pct',r'GPL-0019\] Utilization:\s+([0-9.]+)')]:
  m=re.search(pattern,text)
  if m:r[name]=float(m[1])
 return r

def check_peers(records,r):
 peers=[x for x in records if x['label']==r['label']]
 for peer in peers:
  for key in ['placement_hash','state_hash','area_um2','hpwl_um','wns_ns','tns_ns','legal']:
   if r.get(key)!=peer.get(key):
    raise RuntimeError(f"Nonidentical {key}: {r['label']} {r['binary']} t{r['threads']} vs {peer['binary']} t{peer['threads']}")

def main():
 p=argparse.ArgumentParser();p.add_argument('--task',type=Path,required=True)
 p.add_argument('--phase',choices=['acceptance','benchmark'],required=True);a=p.parse_args()
 task=a.task;out=task/'results';out.mkdir(exist_ok=True)
 for name,digest in json.loads((task/'input-binary-sha256.json').read_text()).items():
  assert hashlib.file_digest(Path(name).open('rb'),'sha256').hexdigest()==digest,name

 summary=out/'summary.json'
 report=json.loads(summary.read_text()) if summary.exists() else {'runs':[],'complete':False}
 records=report['runs']
 expected=0 if a.phase=='acceptance' else len(CONFIGS)*4
 assert len(records)==expected,(a.phase,len(records))
 if a.phase=='benchmark':
  # Verify the saved acceptance logs rather than trust a partial JSON receipt.
  fixed=[]
  for r in records:
   f=Path(r['log']);assert hashlib.file_digest(f.open('rb'),'sha256').hexdigest()==r['log_sha256']
   current=parse(f,r['kind']);check_peers(fixed,dict(current,label=r['label'],binary=r['binary'],threads=r['threads']))
   for k,v in current.items():assert r[k]==v,(f,k)
   fixed.append(r)
 binary_hash={n:hashlib.file_digest((task/n).open('rb'),'sha256').hexdigest() for n in ['baseline','candidate']}
 report['binary_sha256']=binary_hash
 repetitions=[0] if a.phase=='acceptance' else [1,2]
 for repetition in repetitions:
  for label,kind,design,args,routability in CONFIGS:
   # Acceptance checks native 1/8 first, then candidate; later repetitions
   # alternate complete native/candidate runs for each thread count.
   if repetition==0:sequence=[('baseline',1),('baseline',8),('candidate',1),('candidate',8)]
   else:
    names=['candidate','baseline'] if repetition%2 else ['baseline','candidate']
    sequence=[(n,t) for t in [1,8] for n in names]
   for name,threads in sequence:
    directory=out/f'{repetition}-{label}-{name}-{threads}';directory.mkdir(exist_ok=False)
    script=directory/'run.tcl'
    if kind=='flow':
     script.write_text((task/'flow.tcl').read_text().replace('set_thread_count 1\n',f'set_thread_count {threads}\n',1))
    else:script.write_text((task/'gpl-only.tcl').read_text())
    env=dict(os.environ,FLASHEDA_TASK=str(task),FLASHEDA_OPENROAD_SOURCE='/tmp/flasheda-gpl-timing/OpenROAD',
             FLASHEDA_DESIGN=design,FLASHEDA_OUTPUT=str(directory/'placed'),FLASHEDA_RUDY_MODE='native',
             FLASHEDA_RUDY_PROFILE='0',FLASHEDA_GPL_TIMING='0',FLASHEDA_THREADS=str(threads),
             FLASHEDA_GPL_ARGS=args,FLASHEDA_ROUTABILITY=routability)
    env.pop('LD_PRELOAD',None);env.pop('SAMPLER_OUT',None);env.pop('OMP_PLACES',None)
    command=[str(task/name),'-no_init','-exit',str(script)]
    (directory/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    start=time.monotonic();log=directory/'run.log'
    with log.open('w') as f:subprocess.run(command,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=5400)
    r=parse(log,kind);r.update(label=label,kind=kind,design=design,args=args,binary=name,threads=threads,
                             repetition=repetition,seconds=time.monotonic()-start,log=str(log),
                             log_sha256=hashlib.file_digest(log.open('rb'),'sha256').hexdigest())
    (directory/'result.json').write_text(json.dumps(r,indent=2)+'\n')
    check_peers(records,r)
    records.append(r);summary.write_text(json.dumps(report,indent=2)+'\n')
    print(label,name,threads,repetition,r['placement_ms'],r['placement_hash'],flush=True)
 if a.phase=='acceptance':
  assert len(records)==len(CONFIGS)*4;report['acceptance_complete']=True
 else:
  assert len(records)==len(CONFIGS)*12
  report['summary']={}
  for label,*_ in CONFIGS:
   for threads in [1,8]:
    result={}
    for name in ['baseline','candidate']:
     values=[r['placement_ms'] for r in records if r['label']==label and r['threads']==threads and r['binary']==name]
     assert len(values)==3
     result[name]={'median_ms':median(values),'min_ms':min(values),'max_ms':max(values),'samples':len(values)}
    result['reduction_pct']=100*(1-result['candidate']['median_ms']/result['baseline']['median_ms'])
    report['summary'][f'{label}-{threads}']=result
  report['complete']=True
 summary.write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
