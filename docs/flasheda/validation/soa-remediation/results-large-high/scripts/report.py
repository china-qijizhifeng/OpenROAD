import argparse
import json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--task',type=Path,required=True);a=p.parse_args()
task=a.task;d=json.loads((task/'results/summary.json').read_text())
assert d['complete'] and len(d['runs'])==36
out=['# 高利用率配置 v3 最终重复计时','',
     '三个配置的所有 native/v3、1/8 线程、每种三次样本均在本独立 worker 执行。',
     '不复用其他 worker 的计时。正式计时不加载采样器，OpenSTA 固定 1 线程。','',
     '| 配置 | GPL线程 | 原版中位秒（范围） | v3中位秒（范围） | 耗时减少 |',
     '|---|---:|---:|---:|---:|']
for label in ['large01','large02','large01-rd']:
 rows=[r for r in d['runs'] if r['label']==label]
 assert len(rows)==12 and len({r['placement_hash'] for r in rows})==len({r['state_hash'] for r in rows})==1
 for t in [1,8]:
  r=d['summary'][f'{label}-{t}'];n=r['baseline'];c=r['candidate']
  fmt=lambda x:f"{x['median_ms']/1000:.3f}（{x['min_ms']/1000:.3f}–{x['max_ms']/1000:.3f}）"
  out.append(f"| {label} | {t} | {fmt(n)} | {fmt(c)} | {r['reduction_pct']:.1f}% |")
out += ['', '这些配置只有 LEF/DEF：比较完整 cell/port/net 最终状态、面积与 HPWL；不宣称时序或合法化验收。','']
(task/'report.md').write_text('\n'.join(out));print('\n'.join(out))
