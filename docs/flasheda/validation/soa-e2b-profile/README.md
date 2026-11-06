# HK E2B profile 原始证据

2026-11-06，通过 `hk-linlizhi-jumpserver` 使用 SWEEnvMaker 的 0820 客户端配置。实际发现并使用的模板为 `flasheda-r5-cpu8-t16-20261030` / `vshaqeiaoi96u3cp66t3`，SDK 2.37.1。实例为 Ubuntu 24.04、16 GiB、8 个逻辑 CPU、4 个物理核；GPL 的有效配置为 1/4 线程。Northjob 的 1/8 线程数据不与本实验混合计时。

AES、large01、large02 分到三个独立实例。每个实例内原版与候选串行运行。12 次 phase/OMPT，以及 AES、large02 的 8 次独立 perf stat/record，均核对完整布局、连接、质量和算法轨迹；参考为已验收原版/候选的 1 线程日志。4 线程未触发 placement 线程降档。

`experiment.json` 固定输入、源码、二进制、资源、工具和归档指纹，并记录四个实例的释放结果。原始二进制位于 `build/soa/merge-master-local/e2b-migration.tar.gz`；原始 `perf.data`、日志、Tcl、摘要和调用链报告保存在 `build/soa/e2b-profile/*-evidence.tar.gz`。这些大产物不纳入 Git。HK 的任务目录为 `/root/openroad-soa-e2b-profile-20261106/`，保留第二份归档。

表格通过下面的命令从原始 JSON 重建：

```sh
python3 docs/flasheda/validation/soa-e2b-profile/summarize.py --require-complete
```

| 文件 | 含义 |
| --- | --- |
| `phase-runs.csv` | 12 次归因运行的根阶段、迭代数、精确性、OMPT 状态 |
| `phase-detail.csv` | 阶段 inclusive/exclusive、调用次数、每次耗时、累计等待线程时间 |
| `perf-runs.csv` | 8 次 perf 运行的 wall/user/sys、峰值 RSS 字节、样本与 lost 元数据 |
| `perf-events.csv` | 40 行原始计数、单位、event runtime、percent-running 和派生指标 |
| `perf-symbols.csv` | 四次 record 中 self CPU 百分比 ≥0.1% 的符号；已去除终端列宽填充 |

phase 根耗时为 initial placement 与 Nesterov 之和，不等于整个进程耗时。OMPT wait 为累计线程时间，包含 worker 空闲/停驻；阶段标签在 wait begin 采样，不能从 wall 扣除，不能直接解释为可消除的同步时间。

perf 使用未修改的正式 ELF。独立工具实例安装的 perf 6.8.12 已在实际 kernel 6.1.158 上完成计数和 DWARF 调用链 smoke；工具、动态加载器和依赖均保存 SHA。record 为 99 Hz `cycles:u`、8 KiB DWARF 用户栈，四次 lost 均为 0。self 百分比描述整个进程的 CPU 样本；children 调用链百分比不可相加。部分 libc 地址缺少符号，未冒充源码级定位。

stat 每角色仅运行一次，用于归因和观察；不构成重复性能验收。所有计数的 percent-running 均为 100%。两次候选的 branch-misses 为 0，保留原值但不据此解释收益；guest PMU 也不用于跨宿主机绝对比较。RSS 使用 Linux `/usr/bin/time -v` 的 KiB 转字节。

准备阶段曾将“模板没有 perf”作为组合探测命令的非零退出；上传与完整性检查已先完成。随后集中检查三份环境并改为独立 syscall 能力探测，没有重复业务运行、重传完整输入或重编译 OpenROAD。客户端实际 key 字节检查回执在控制归档中，未向沙箱下发客户端凭证。
