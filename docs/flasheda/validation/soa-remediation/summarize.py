import argparse
import json
from pathlib import Path
from statistics import median
p=argparse.ArgumentParser()
p.add_argument('--task',type=Path,required=True)
a=p.parse_args();task=a.task
result={}
for name in ['benchmark','benchmark-tuned']:
    d=json.loads((task/'fork-results'/name/'summary.json').read_text())
    if not d.get('complete'):raise RuntimeError(f'Incomplete: {name}')
    result[name]={'summary':d['summary'],'binary_sha256':d['binary_sha256'],'resources':{}}
    for binary in sorted({r['binary'] for r in d['runs']}):
        for t in [1,8]:
            values=[r['max_rss_kb'] for r in d['runs'] if r['binary']==binary and r['threads']==t and 'max_rss_kb' in r]
            if values:result[name]['resources'][f'{binary}-{t}']={'median_rss_mib':median(values)/1024,'range_rss_mib':[min(values)/1024,max(values)/1024]}
    result[name]['hashes']=sorted({r['placement_hash'] for r in d['runs']})
result['replay']={}
for call in [2,400,1100]:
    rows=[json.loads(line) for line in (task/f'replay-{call}.jsonl').read_text().splitlines()]
    result['replay'][str(call)]={variant:{stage:median(r[stage] for r in rows if r['variant']==variant) for stage in ['scatter_ms','gather_ms']} for variant in ['previous','separable']}
result['cases']={}
for name in ['cases-final','cases-tuned']:
    d=json.loads((task/'fork-results'/name/'summary.json').read_text())
    assert len(d)==10 and all(len(v)==2 for v in d.values())
    assert all(r['exact'] for v in d.values() for r in v.values())
    result['cases'][name]={'designs':list(d),'comparisons':sum(len(v) for v in d.values()),'exact':True}
(task/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
