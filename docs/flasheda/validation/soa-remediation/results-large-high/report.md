# 高利用率配置 v3 最终重复计时

三个配置的所有 native/v3、1/8 线程、每种三次样本均在本独立 worker 执行。
不复用其他 worker 的计时。正式计时不加载采样器，OpenSTA 固定 1 线程。

| 配置 | GPL线程 | 原版中位秒（范围） | v3中位秒（范围） | 耗时减少 |
|---|---:|---:|---:|---:|
| large01 | 1 | 215.435（212.947–218.313） | 180.010（178.097–180.205） | 16.4% |
| large01 | 8 | 93.309（91.529–97.161） | 51.322（50.725–54.237） | 45.0% |
| large02 | 1 | 243.114（241.486–243.539） | 214.071（210.978–216.140） | 11.9% |
| large02 | 8 | 106.412（100.685–110.652） | 65.464（64.084–68.374） | 38.5% |
| large01-rd | 1 | 468.669（468.176–468.887） | 374.864（373.615–382.126） | 20.0% |
| large01-rd | 8 | 194.173（177.380–202.310） | 105.706（104.200–105.746） | 45.6% |

这些配置只有 LEF/DEF：比较完整 cell/port/net 最终状态、面积与 HPWL；不宣称时序或合法化验收。

作业 `flasheda-soa-large-high-remediation-20261106`（`job-bb6ab222440c-20261106125349`），
CPU Intel Xeon Gold 6330，24 CPU / 64 GiB / 0 GPU，WeaverRSI / proj-32。
36 次累计进程 wall time 7880.863 秒（131.35 分钟；含加载与完整状态导出），平台终态 Completed，
worker 与三个命令均 exit 0；精确 alias 删除成功，回查 `[]`。

每配置 12 次的 cell/port/net 状态、面积、HPWL 完全一致；详见 `exact-comparison.json`。
native 二进制 SHA256：`d3e694dc50bde8d8c848fe153ff7d64e55ce6a96d9d9b4548cc6484ab07f1c0f`。
v3 二进制 SHA256：`4d627682245a2896c22ce73248c5e122ad5b6f01b3f73b387f8bb5f1c07a4252`。

本地归档 [evidence.tar.gz](/Users/shom/Downloads/openroad-gpl-soa/build/soa/remediation-large-high/evidence.tar.gz)：292534034 字节，本地 SHA256 校验通过：
`1afc610b120df0280a6c982b36ea32c0118d39b78d76bf60995c5693d21fec0e`。
归档包含完整状态日志、命令、脚本、源码 manifest 与固定输入字节；二进制和原始大型布局留在 GPFS。

计时原始记录：`summary.json`；环境：`runtime.json` / `lscpu.txt`；
输入/源码：`input-binary-sha256.json` / `source-sha256-v3.json`；
清理回执：`job-completed.json` / `job-delete.json` / `job-after-delete.json`。
