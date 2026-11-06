# GPL SoA 诊断整改（2026-11-06）

依据：[soa-speedup-diagnosis.md](../soa-speedup-diagnosis.md)。
原版基线为 `18216f3506f090aa7e12f83620b47f21ddff4f40`；整改前为 `dd8f2fc`。

## 实现

- **C1**：按 FFT 的 x-major 行号并行计算 density/overflow，存入连续 scratch；
  最终按原 y-major 顺序归约。总量保留原先的 `int64_t += float` 逐次转换，
  使用实际 GPL 线程数，移除无调用方的 BinGrid 线程设置。
- **C2、D1**：中心、pin、坐标候选和梯度预条件直接操作列及行号。
  region 的 pool-row / optimizer-index 映射只在插入、删除、swap-remove、
  handle 重绑定及 filler 恢复时失效；迭代时不打包几何数据。
  pin 同步使用现有 cell→pin CSR，每个 pin 保留整数中心与 offset 的原计算。
- **D2**：单字节 kind/locked 列代替每轮追踪 Instance。
  初始锁从 Instance 读取，GCell::lock 同步两者；增量布局第二阶段解锁后，
  在进入下一轮前刷新 locked 列。filler 快照复制和行复用同时覆盖这些标志。
- **D3**：density view 每个 bulk kernel 构建一次；线长梯度和 overflow 的
  debug 检查在循环外完成。单单元调试/回调入口保留。
- **D4**：rounded bin 边界按轴存为 `nx+1` / `ny+1`；每个 x 重叠宽度只算一次，
  常见小单元的 y 重叠高度用栈数组复用。大宏单元使用无尺寸上限的后备路径。
  scatter 仍逐贡献截断，gather 仍按 x-outer/y-inner 累加。
- **D5（追加）**：CPU-only 构建将密度 gather 与梯度预条件合为一次遍历，
  继续写出密度梯度供回调与原顺序归约使用。GPU 构建保留原 backend 调用链。
  线长 WA 的八个 net 累加值也改为局部标量，按原 pin 顺序累加后一次写回；
  边界和当前 pin 坐标各加载一次。
- **C3**：保持 Ooura 一维变换，按行/列用 OpenMP 批量并行。
  cos/sin 和 bit-reversal 表先串行初始化，列 scratch 每线程独享。
  单线程与小网格沿用原 `ddct2d` / `ddsct2d` / `ddcst2d`。
  不使用固定线程数、每次创建 pthread 的编译开关。
- **C4**：并行计算各坐标差的平方，再按原 x/y 交错顺序串行求和。
- **E**：完成 IPO 与内核 ISA 的定向实测后，因组合稳定退化，最终源码移除
  `GPL_ENABLE_IPO` / `GPL_CPU_ARCH` 两项实验选项；保留 CPU 禁用浮点收缩。
  Bazel 已全局设置 `-ffp-contract=off`。历史实验源码、实际参数、二进制和日志留存。
- **handoff P1**：per-pin WA 中间量按 net→pin CSR 位置存放，几何与 pin 身份行不变；
  `waSlot` 列将 pin 行映射到 WA 行，bulk gather 和单 pin 调试入口使用同一映射。
  每轮 net force 无条件覆盖四个数值和四个 active 标记，清零只在拓扑重建时执行，
  同时覆盖断开 pin 和空洞行。单 net 内累加顺序保持不变。
- **缓存与接口**：optimizer handle 的追加、swap-remove、改绑集中在三个函数置脏；
  `checkConsistency` 调试开关核对派生映射及 kind/locked 与实时 Instance 状态。
  density backend 对子集保留逐单元计算，并在 Release 下检查输出尺寸。

本次不改变 inflation、回退、停止条件、timing repair 次数或 OpenSTA 线程数。
A/B 中会改变 QoR 的算法实验、坐标更新与 scatter 融合、初始布局 P2、PGO 和 GPU 不计入本次整改结果。

## 验收

在 Northjob `flasheda-soa-remediation-20261106`
（identity `job-317d12d0424e-20261106090444`，WeaverRSI / proj-32）运行，
24 CPU、64 GiB、0 GPU，CPU 为 Xeon Gold 6330。
使用固定镜像、同一缓存源码树和 GCC/CMake 构建；本次没有执行 Bazel 构建或 Bazel 测试。
补充尝试确认，本机和现有 BP Dev SSH 主机均无 `bazelisk` 或 `bazel`；
`bazelisk test //src/gpl/test/soa/...` 在启动前即失败，未运行任何 Bazel target。
回执见 [bazel-attempt.json](soa-remediation/results/bazel-attempt.json)。
OpenSTA 固定 1 线程，GPL 分别使用 1/8 线程，`OMP_PROC_BIND=FALSE`、
`OMP_DYNAMIC=FALSE`，不设 `OMP_PLACES`。runner 拒绝静默降低 GPL 线程数。

默认构建已通过 12 个测试案例 × 1/8 线程，共 24 组与原版的差分：

- GCD、AES：完整布局、合法化，以及坐标/朝向/master 哈希、WNS/TNS、面积、HPWL 对照。
- 十个扩展 Tcl 案例：`incremental02`、`incremental03`、`macro01`、`macro02`、
  `region01`、`cluster_place01`、`simple01-td-repair`、`place_ios_td`、
  `place_ios_rd`、`incremental01`。比较全部 cell 坐标、朝向、master、placement status，
  port 几何和 net 连通性；Tcl 中其他断言继续运行。

这里的 12 个案例会复用设计输入，不表示 12 个不同网表。
扩展案例的原版日志可复用，但复用前会重新计算原始输出哈希。
单次功能验收耗时不作为性能结论。

第一阶段默认构建及 IPO/ISA 对照构建已完成上述验收。
追加融合版本已通过同一组 24 个差分对照，三次配对计时也已完成；
最终阶段采样及 handoff 修正版本的小案例验收均已完成，由独立 subagent 执行和监控。
参考 handoff 后另固定 `source-sha256-v3.json`：它包含 CSR WA 布局、缓存变更统一入口、
密度接口修正与实验选项移除。第一阶段、追加融合版本和此最终候选的二进制分别保存，
各自的测量不能互相替代。
最终候选二进制 SHA-256 为
`4d627682245a2896c22ce73248c5e122ad5b6f01b3f73b387f8bb5f1c07a4252`，
已完整构建并复制到独立大设计验收目录；常规 24 组差分已全部通过。
另开启 `checkConsistency` 的十组对照全部通过：GCD/AES × 1/8，
incremental02/03、simple01-td-repair × 1/8。检查覆盖行映射、kind/locked 缓存，
包含 AES routability/filler 回滚和真实时序修复。AES 配对计时和下述大设计差分均已完成。

独立 storage、point、density、metadata 测试已在本机 Clang ASan/UBSan 和
worker GCC 环境执行；FFT 另通过 worker ASan/UBSan，串行原库与 1/2/8 线程
全部 phi/ex/ey 输出逐字节一致，覆盖矩形、2048×2 / 2×2048 窄网格和线程切换。
三个原生 AES 密度快照（调用 2/400/1100）的生产新内核、整改前内核与原始捕获输出一致。

### 第一阶段完整 AES 布局性能

同一 worker、相同输入和运行参数，三种二进制交替顺序运行，各三次。
单位为秒，表中为 `global_placement` 中位数（最小–最大值）。

| GPL 线程 | 原版 | 整改前 SoA | 整改后默认构建 | 比整改前少耗时 | 比原版少耗时 |
|---|---:|---:|---:|---:|---:|
| 1 | 79.870（78.925–80.005） | 75.266（74.729–76.129） | 62.621（62.418–62.915） | 16.8% | 21.6% |
| 8 | 44.033（43.811–49.563） | 41.549（41.527–41.729） | 30.323（30.286–30.337） | 27.0% | 31.1% |

全部 18 次运行的坐标哈希、WNS/TNS、面积和 HPWL 完全一致。
原版 8 线程出现一个偏慢样本，未删样本或据此追加重复，保留完整范围并报告中位数。
本次 worker 无 `/usr/bin/time`，未得到可比较的 peak RSS，不能据此声称内存改善。

### 可选 IPO / ISA 组合

修正架构参数范围后，该构建也通过全部 24 组差分。
另做默认与定向构建各三次交替计时，12 次输出与 QoR 全部一致：

| GPL 线程 | 同批默认构建（秒） | IPO + 内核 icelake-server（秒） | 定向配置多耗时 |
|---|---:|---:|---:|
| 1 | 62.434（62.391–62.442） | 73.683（73.460–73.728） | 18.0% |
| 8 | 30.659（30.284–30.724） | 34.866（34.768–35.134） | 13.7% |

该组合不推荐在此次硬件/编译器配置上开启。
这里测量的是 IPO 与内核 ISA 一起开启的组合，没有分离各自影响，
不能据此断言某一个选项单独导致退化。主结论使用无 IPO、无架构覆盖的默认构建。
参考 [optimization-handoff.md](../optimization-handoff.md) 后，最终源码进一步移除了
这两个仅 CMake 可用、无正向测量收益的实验选项，避免维护没有验收价值的配置分支。

### 大设计验收范围

上面的 21.6% / 31.1% **只适用于此次 AES 配置**，不是全部设计的通用提速结论。
handoff 指出的 large01/large02 首版 SoA 回退已通过下述独立重复计时核实并消除。
另建独立作业 `flasheda-soa-large-remediation-20261106`
（identity `job-388e9a7b6ce7-20261106113510`，WeaverRSI / proj-32，24 CPU / 64 GiB / 0 GPU），
固定追加融合版本（v2）二进制，先对七个配置做原版与候选 1/8 线程单轮差分；
这 28 次运行现已全部精确通过。
配置为 ibex_core-td-rd、medium04、macro03、macro03-rd、large01、large01-rd、large02。
最终 v3 的七配置首轮也已全部精确通过。
ibex_core、medium04、macro03、macro03-rd 的最终重复计时在此 worker 完成，原版首样本来自
同 worker 的 v2 阶段，复用前核对固定输入和原始日志；另两轮按原版/v3 交替运行。
large01、large02、large01-rd 则完整分到独立作业
`flasheda-soa-large-high-remediation-20261106`
（identity `job-bb6ab222440c-20261106125349`），同镜像、CPU 型号和资源，
每个配置的原版/v3 × 1/8 线程 × 各三次全部在新 worker 重新运行，不复用跨 worker 计时。
两作业均固定输入、二进制、线程与运行环境，不加采样器。
首轮完整 cell/port/net 状态、面积和 HPWL 一致；ibex_core 另比较 WNS/TNS 与合法化。
原 worker 的四配置 48 次已完成，脚本生成结果如下。单位为 `global_placement`
秒，均为三次中位数（最小–最大），不包含状态导出、详细布局或时序统计。

| 配置 | GPL 线程 | 原版 | 最终 v3 | 耗时减少 |
|---|---:|---:|---:|---:|
| ibex_core-td-rd | 1 | 175.590（172.270–176.125） | 130.857（130.706–131.060） | 25.5% |
| ibex_core-td-rd | 8 | 92.610（91.218–97.427） | 61.482（61.344–62.515） | 33.6% |
| medium04 | 1 | 221.268（220.708–223.214） | 149.112（149.105–150.757） | 32.6% |
| medium04 | 8 | 87.042（84.673–90.042） | 46.403（45.169–47.421） | 46.7% |
| macro03 | 1 | 189.959（189.556–192.995） | 135.379（135.088–135.940） | 28.7% |
| macro03 | 8 | 72.266（70.912–74.659） | 40.494（39.952–41.651） | 44.0% |
| macro03-rd | 1 | 456.460（454.330–460.040） | 322.789（322.226–322.888） | 29.3% |
| macro03-rd | 8 | 185.937（180.459–192.001） | 91.989（91.353–92.918） | 50.5% |

完整统计和复用原版首样本的来源说明见
[脚本生成表](soa-remediation/results-large/report.md)及
[48 次原始结构化记录](soa-remediation/results-large/v3-final.json)。
高利用率 worker 的 36 次运行也已完成，全部状态一致：

| 配置 | GPL 线程 | 原版 | 最终 v3 | 耗时减少 |
|---|---:|---:|---:|---:|
| large01 | 1 | 215.435（212.947–218.313） | 180.010（178.097–180.205） | 16.4% |
| large01 | 8 | 93.309（91.529–97.161） | 51.322（50.725–54.237） | 45.0% |
| large02 | 1 | 243.114（241.486–243.539） | 214.071（210.978–216.140） | 11.9% |
| large02 | 8 | 106.412（100.685–110.652） | 65.464（64.084–68.374） | 38.5% |
| large01-rd | 1 | 468.669（468.176–468.887） | 374.864（373.615–382.126） | 20.0% |
| large01-rd | 8 | 194.173（177.380–202.310） | 105.706（104.200–105.746） | 45.6% |

完整来源见 [高利用率脚本生成表](soa-remediation/results-large-high/report.md)。
七配置共 84 个最终记录，每配置的两种二进制、1/8 线程、各三次均具有唯一完整
cell/port/net 状态哈希，面积和 HPWL 一致；ibex_core 另验证时序与合法化。
其余六配置只有 LEF/DEF，不宣称完成时序或合法化验收。
本次七配置的中位耗时均下降：单线程 11.9%–32.6%，8 线程 33.6%–50.5%。
该范围描述本次固定设计、参数和硬件，不代表任意流程的提速保证。
WA CSR 相对 v2 的单独贡献没有三次大设计配对实验，不能将整版收益全部归因于 CSR。
跨分片的 [完整计时索引](soa-remediation/results-large/cross-shard-summary.md) 和
[84 条记录与来源](soa-remediation/results-large/cross-shard-summary.json)核对了
配置到 worker 的映射、两个分片源码 manifest 与原版/候选二进制的一致性。

### 追加融合与局部累加

追加版本仅改动 CPU 密度 gather/预条件的遍历组织及线长 net 的累加存放位置。
不改变 `fastExp`、pin 顺序、各轴密度求和顺序、截断规则或梯度总量归约顺序。
回调仍能读取完整 wirelength/density 梯度列。

合成线长对拍覆盖 10,000 个 net、约 124 万 pin，包括空 net、不同 fanout、
零系数和截断条件；全部 net/pin WA 数值列逐位一致，并通过本机 ASan/UBSan。
worker micro 样本混合不同系数与分支条件，耗时范围较宽，因此不将其汇总中位数
当作完整布局提速结论，完整布局采用下述真实 AES 配对计时。

追加融合版本另与第一阶段默认构建交替运行各三次，全部 12 次输出/QoR 一致：

| GPL 线程 | 同批第一阶段（秒） | 追加融合版本（秒） | 少耗时 |
|---|---:|---:|---:|
| 1 | 62.497（62.420–62.857） | 61.886（61.731–62.037） | 0.98% |
| 8 | 30.378（30.084–30.554） | 30.328（30.237–30.417） | 0.16% |

单线程只有约 0.61 秒的额外收益；8 线程差值在样本噪声范围内，不宣称有稳定提升。
该对照同时加入 gather/预条件融合与 net 局部累加，没有分离两者贡献。

最终 CSR WA/缓存/接口候选再与追加融合版各运行三次，全部 12 次输出/QoR 一致：

| GPL 线程 | 同批追加融合版（秒） | 最终候选（秒） | 候选耗时变化 |
|---|---:|---:|---:|
| 1 | 62.080（61.806–62.128） | 62.226（62.006–62.422） | +0.24% |
| 8 | 30.632（30.096–30.722） | 30.501（30.216–30.844） | −0.43% |

两个差值都在样本噪声范围内，不宣称该轮在 AES 上产生额外性能收益。
CSR WA 布局针对大设计的并行写入和清零开销；七配置最终整版测量见上节，
这里的 AES 对照没有显示可归因于 CSR 的额外收益。

### 追加融合版本的阶段采样

沿用诊断中的 1 kHz 主线程墙钟采样方法，1/8 线程各一次；符号化在 BP Dev host
独立目录完成。坐标哈希与原版一致。下表单位秒，包含采样开销，不能代替无采样中位数。
由于 gather 与预条件已融合，两者合并统计，不能再直接比较旧的单独 gather 行。

| 阶段 | 1 线程 | 8 线程 |
|---|---:|---:|
| 时序反馈 | 17.81 | 17.86 |
| 密度 scatter / overflow | 9.78 | 2.49 |
| FFT | 5.02 | 1.34 |
| 密度 gather + 预条件融合 | 10.71 | 1.60 |
| 线长 WA force | 2.49 | 0.51 |
| 线长梯度 gather | 3.38 | 1.96 |
| 坐标 / pin 同步 | 4.76 | 0.73 |
| 全部归类样本合计 | 63.91 | 31.16 |

8 线程的时序反馈占该次归类总时间约 57%，而 FFT 与坐标/pin 的并行化已明显生效。
这是 AES 的剩余热点；高利用率大设计还要以独立采样/计时判断，不能套用同一占比。

### 密度 replay

每个快照、每种内核交替采样 11 次，每次重复 20 遍。下表为每次调用中位数，单位 ms，
单线程；三者均先通过原生捕获输出核对。

| AES 捕获调用 | scatter 整改前 | scatter 整改后 | gather 整改前 | gather 整改后 |
|---|---:|---:|---:|---:|
| 2 | 7.124 | 5.728 | 6.529 | 5.739 |
| 400 | 7.031 | 5.725 | 6.483 | 5.728 |
| 1100 | 6.753 | 5.547 | 6.293 | 5.552 |

scatter 减少 17.9%–19.6%，gather 减少 11.8%–12.1%。这些是同一快照的
scatter/gather 内核时间，不包括 FFT、overflow、时序修复或完整布局。

## 首次验证暴露的构建缓存问题

新增的 `incremental02` 原版成功，候选在解锁后的第二阶段停滞并最终报 GPL-0305。
第一阶段的全部日志数值一致，首次分歧紧随 `Unlocking all instances`。
静态检查发现，第二次源码叠加后 Make 没有重编译任何 GPL object；候选 ELF
也缺少新加入的 `refreshCellLocks` 符号。归档保留的源码 mtime 早于第一轮刚生成
的 object，因此即使磁盘源码哈希正确，二进制仍是缺少解锁同步的旧实现。

发现后停止后续测试和计时。最终构建先核对每个叠加输入的 SHA-256，再刷新
所有输入 mtime，实际重编译 GPL，并核对默认二进制的新增实现符号后才恢复差分。
初次失败的源码、二进制和日志保留；有效验收从命令 08 开始。

## 定向优化暴露的 Eigen 配置问题

首次可选构建对整个 `gpl_lib` 开启 IPO 和 `-march=icelake-server`。
在第一个 GCD 初始布局案例中，`Eigen::internal::set_from_triplets` 释放内存时崩溃。
停止后续差分和计时后，静态核对 worker 实际 Eigen 头文件、预处理宏和反汇编：
默认 ISA 的 `EIGEN_DEFAULT_ALIGN_BYTES=16`、`EIGEN_MALLOC_ALREADY_ALIGNED=1`；
icelake-server 下分别变为 64 和 0。崩溃位置的代码读取 `ptr[-1]` 后调用 `free`，
使用了手工对齐分配的释放契约。将单个库的 Eigen ISA 改掉，会与应用内其他库的
同名内联分配函数形成不兼容定义。

完整排查 GPL 的 Eigen 使用位置后，架构参数改为只作用于 `nesterovBase`、
`fft`、`fftsg`、`fftsg2d`、`hpwl`、`wirelengthGradient`、`wirelengthSoa` 和
`densityGradient`。初始布局、solver、replace 等其他文件使用默认 ISA。
构建检查读取实际 `compile_commands.json` 和八个内核的 `.o.d`，
核对参数范围，并要求这些内核没有 Eigen 传递依赖。
失败二进制和日志另存，修正后的定向构建差分与重复计时均已通过。

命令 15 的构建已成功，随后参数检查误把独立 `fft_test` 的同源编译条目当作
生产库条目，导致检查失败。静态检查确认生产参数正确后，检查器改为只检查
`gpl_lib.dir` 的输出；命令 15b 复用现成编译产物检查并留存，16b/17b 恢复验收和汇总。
保留原失败回执，没有据此重跑失败的动态案例。

第一阶段默认二进制对应 `source-sha256-final.json`；修正架构配置后的第一阶段
定向构建对应 `source-sha256-kernel-isa.json`。两者仅 CMake 可选分支不同。
追加密度遍历融合与线长寄存器累加的阶段工作区对应 `source-sha256-fused.json`，
由命令 22/23 单独构建和差分。命令 24 在启动任何运行前因结果目录误指向
已有 `benchmark-tuned` 目录而退出；修正 runner 路径后使用 24b/25b/26b 继续计时和采样，
保留失败回执，没有重复功能验收。handoff 后的 CSR WA/缓存/接口修正另固定新快照，
由 subagent 编译与验收，小案例与 AES 结果已完成，大设计结果单独生成。

## 复现入口

独立测试：`bazelisk test //src/gpl/test/soa/...`，新增 `fft_test` 比较串行原库和
1/2/8 线程的全部 phi/ex/ey 字节，覆盖矩形、窄网格和动态线程切换。
密度测试使用独立 all-bin oracle，覆盖超过栈缓存容量的高单元、零尺寸、
非方形/舍入网格和宏单元缩放。
新增 `wirelength_test` 用独立标量方程核对 CSR WA 映射及梯度，覆盖几何行置换、
五次拓扑重建、断开 pin/空洞清零及多次 active 状态切换；本机 Release 与
ASan/UBSan 已通过。

密度 backend 接口另复用已有 `TimingRepairIntegrationTest` fixture 做临时 C++ 验收，
实际 OpenROAD 数据结构、Release / `NDEBUG` 构建，明确运行 1 个测试并通过。
完整列表的 bulk 输出与逐单元结果逐位比较；逆序隔项及重复 handle 的子集同样逐位一致；
空列表可用，输出尺寸错误抛出 `invalid_argument`。没有增加业务测试 hook。
该专项通过 CMake worker 已有库手工链接执行，不计为已运行 Bazel target。

专项脚本首次把测试注册对象与已有 `fixture.cpp` 对象输出为同名文件，导致覆盖和
“0 tests”；静态检查 `nm` 后修复独立对象名，并要求输出明确包含运行 1 个测试和 PASS。
另一次复用 runfiles 链接时因路径已存在退出；核对符号和现有链接后使用 32d 运行通过。
这些前置失败均保留，只有 32d 计为该接口的动态验证证据。

锁状态仍在解锁的上层生命周期边界刷新，没有让 `Instance` 反向依赖 `GCell`。
为降低派生状态漏同步风险，optimizer handle 的增删/改绑已集中为三个写入口，
每次使用行映射时可由 `checkConsistency` 开关与实时 Instance 状态核对。
这不是取消缓存，而是明确写入边界并提供真实回调链验证。

初始布局 P2 尚未改写。静态检查确认，`InitialPlace::updatePinInfo`、
`createSparseMatrix` 和 `Instance::setCenterLocation` 的 pin 更新仍经对象访问列。
后续可在一次初始布局求解的拓扑边界建立 net→pin、instance→pin 和 movable-row 映射；
坐标每轮直接读取实时列，不能缓存上一轮坐标。必须保持同坐标 min/max pin 的
首个命中规则、movable IO 排除条件、原始 fanout 的权重分母，以及 triplet 和 RHS
的插入/累加顺序。Eigen 求解器、对齐和 ISA 保持原配置。该方向需要独立源码快照
与矩阵/RHS、最终布局的差分验证，不计入当前候选的性能收益。

[验收脚本](soa-remediation/)保留此次完整布局差分、重复计时、捕获快照 replay 和源码核对入口。
这些 runner 依赖 GPFS 中固定的历史基线、公开设计输入与解析器，脚本明确记录路径。
replay 直接包含本次生产 `densitySoa.h`，与诊断前版本及原生捕获输出同时比较；
它只度量 scatter/gather，不代表完整布局提速。

## 证据与作业状态

完整构建、小案例差分、调试一致性检查、AES 重复计时、采样和专项接口验收的紧凑证据包：
`build/soa/remediation/evidence.tar.gz`（103,319,866 bytes，1001 个归档成员），SHA-256：
`8d93bfef071d72a05aab89869f70486362959477c08acd039d155d0124fa1c9c`。
本地拉取后重新核对一致。原始大布局输出、采样二进制记录和被测 ELF 保留在
GPFS `/gpfs/users/linlizhi/flasheda-soa-remediation-20261106/`。
脚本目录保留源码 manifest、运行入口、结构化结果及完成/清理回执。
清理后追加的 Bazel 工具可用性回执单独保留在脚本目录，不在上述证据包内。

该构建验收作业已写 STOP，`worker.exit=0`，Northjob 终态 Completed；
按精确 alias 和 identity 删除成功，复查返回 `[]`。
大设计原 worker（四配置 48 次最终运行及七配置 28 次 v2 功能运行）也已完成归档：
`build/soa/remediation-large/evidence.tar.gz`，271,107,980 bytes、318 个成员，SHA-256：
`7fb8afe60631c643aaf17754b23b516d0758a7f95406d81ab08f70da90d6c4ef`。
本地核验一致；包含冻结输入、原始日志、结构化结果与命令回执，不包含 ELF。
命令 01–05 和 worker 均退出 0，作业 Completed；精确删除
`job-388e9a7b6ce7-20261106113510` 后，按同 alias 复查 `[]`。
回执位于 `soa-remediation/results-large/`。
独立高利用率作业的证据包为
`build/soa/remediation-large-high/evidence.tar.gz`，
292,534,034 bytes、177 项，SHA-256：
`1afc610b120df0280a6c982b36ea32c0118d39b78d76bf60995c5693d21fec0e`。
本地核验一致；包含固定输入字节、完整状态日志、源码 manifest、脚本和命令回执。
作业 Completed，worker 与三条命令全部退出 0；精确删除
`job-bb6ab222440c-20261106125349` 成功，复查 `[]`，回执位于
`soa-remediation/results-large-high/`。全部本次创建的构建/验收作业均已清理。
