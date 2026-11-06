# GPL SoA 提速有限的原因诊断

2026-11-06。对象：`flasheda/gpl-soa` 分支 `dd8f2fc`（GPL 七类实体全字段 SoA）。
此前测得：AES 单线程 79.3 → 75.0 s（−5.4%），8 线程 44.7 → 41.3 s（−7.5%）。

## 结论

1. **SoA 只覆盖了一半左右的时间。** 时序反馈（17.9 s）、FFT（5.0 s）、初始布局（3.8 s）、
   步长等阶段完全没被改动，合计占最终单线程时间的 **36.9%**。8 线程下比例更高：
   时序反馈一项就占 **42.9%**。
2. **SoA 改了存储，但大部分代码仍按对象访问。** 除密度和线长两个内核外，其余代码
   仍经 `GCellHandle → GCell → SoaRow → shared_ptr<Store> → column[slot]` 逐字段读取。
   每个字段都多出几次相互依赖的加载，同一单元的字段也从一条 cache line 分散到多条。
   只迁字段、不改访问的版本（`soa-fft`）因此变慢到 **108.5 s**。最终版本中，
   坐标更新、线长和梯度合并仍比原版慢 **+4.6 s**，抵消了密度内核约一半的收益。
3. **密度内核的提速来自"去掉对象访问"，与数据是否常驻无关。** 每次调用都打包的
   adapter 和常驻 SoA 的内层循环几乎一样快（scatter 8.76 / 9.31 s，gather 8.51 / 9.58 s）。
   常驻方案省掉了 5.56 s 的打包时间，但上一条的回退把这部分又吃掉了，所以常驻版（77.06 s）
   反而比 adapter（75.03 s）还慢。
4. **新的布局冲突。** bin 改为 FFT 的 x-major 存储，但 overflow 循环为了保留原有归约顺序
   仍按 y-major 遍历，每次访问跨 256 个元素。这个循环从 1.28 s 涨到 2.68 s。
5. **8 线程下原版本来就是并行的。** 原版 scatter 用原子加、gather 和线长本来就是 OpenMP，
   SoA 只减少了已并行部分的耗时。剩下 **75.9%** 的时间串行或基本串行：时序反馈、
   Ooura FFT、overflow 循环、密度中心/pin 同步和步长计算。

总账（单线程，采样秒数）：变快合计 −9.83 s，变慢合计 +5.05 s，净 −4.79 s。

6. **在更真实的高利用率设计上，SoA 是负收益。** large01、large02（27–31 万单元，利用率约 64–66%）
   单线程持平或略慢，8 线程**慢约 12%**。原因是 pin 同步翻倍、线长 force 的并行伪共享和串行清零、
   overflow 循环在 512×512 网格下涨到约 3 倍，以及初始布局 getter 变慢。低利用率设计
   （ibex_core、medium04、macro03）则因密度占比高，单线程有 8–13% 收益。详见
   [更真实用例上的复测](#更真实用例上的复测)。

## 测量方法

- 在 Northjob `flasheda-soa-diagnose-20261106`（identity `job-cd4b4319b377-20261106084135`，
  WeaverRSI / proj-32，CPU 24，Xeon Gold 6330）上运行与验收相同的 AES Tcl
  （seed 42，density 0.5，timing+routability，`sta` 线程 1，`OMP_PROC_BIND=FALSE`）。
- 使用自写的 `LD_PRELOAD` 采样器：对主线程做 1 kHz 墙钟栈采样，用 `addr2line -i` 符号化，
  展开内联函数，再按调用栈归入阶段。多线程运行中，主线程在 OpenMP barrier 里等待的时间
  计入所在阶段，因此结果代表墙钟时间。
- 采样开销约 3%（baseline 82.0 s，而无采样中位数为 79.3 s），各配置相同。
  采样结果与此前插桩计时吻合：scatter 21.3 vs 20.8 s，gather 13.7 vs 13.3 s，
  时序反馈 17.8 vs 17.2 s。
- 硬件计数器 `perf_event_open` 被拒绝（EACCES），因此**没有** IPC 或 cache miss 数据；
  文中对 cache 行为的判断是根据源码和时间分布推断的，未经计数器验证。
- 7 次运行的坐标哈希全部为 `bc141af641d8…`，与验收结果一致。
- 下文所有表格均由 `build/soa/diagnose/tables.py` 从采样证据生成，原样粘贴。

### 采样运行

| 二进制 | 线程 | global_placement ms | 坐标哈希前 12 位 |
|---|---:|---:|---|
| baseline | 1 | 82023 | bc141af641d8 |
| soa-final-openroad | 1 | 77227 | bc141af641d8 |
| soa-fft-openroad | 1 | 108735 | bc141af641d8 |
| optimized-openroad | 1 | 75192 | bc141af641d8 |
| soa-density-openroad | 1 | 78545 | bc141af641d8 |
| baseline | 8 | 45941 | bc141af641d8 |
| soa-final-openroad | 8 | 41715 | bc141af641d8 |

各二进制对应的阶段：`soa-fft` 只把字段迁入列存，访问仍经对象 getter；
`soa-density` 在此基础上把密度 scatter/gather 改为直接读列；
`optimized-openroad` 是 adapter 实验，原版对象不变，每次调用把数据打包成 SoA 后算密度；
`soa-final` 是最终提交 `dd8f2fc`。

## 单线程分阶段

单位为秒，数据来自主线程墙钟采样。

| 阶段 | 原版 | 仅字段 SoA | +直接密度内核 | 密度 adapter | 最终 SoA | 最终−原版 |
|---|---:|---:|---:|---:|---:|---:|
| 时序反馈 | 17.81 | 17.84 | 17.90 | 17.88 | 17.92 | +0.12 |
| 拥塞反馈 | 0.91 | 1.03 | 0.88 | 0.86 | 0.87 | -0.05 |
| 初始布局 | 3.59 | 3.65 | 3.62 | 3.64 | 3.80 | +0.20 |
| 密度 scatter + overflow | 21.30 | 31.42 | 13.41 | 15.30 | 14.32 | -6.98 |
| FFT/DCT | 5.49 | 5.57 | 5.03 | 5.49 | 5.02 | -0.46 |
| 密度场其他 | 0.70 | 1.01 | 0.69 | 0.70 | 0.69 | -0.01 |
| 密度 gather | 13.67 | 23.08 | 11.91 | 12.45 | 11.51 | -2.16 |
| 线长 WA force | 2.85 | 4.60 | 4.55 | 2.84 | 3.17 | +0.32 |
| 线长梯度 gather | 4.33 | 6.25 | 6.18 | 4.40 | 4.94 | +0.61 |
| 梯度合并/预条件 | 2.37 | 2.69 | 2.67 | 2.41 | 3.04 | +0.67 |
| HPWL | 1.07 | 1.43 | 1.40 | 1.15 | 0.89 | -0.17 |
| 坐标更新 + pin 同步 | 5.83 | 8.00 | 7.98 | 5.93 | 8.87 | +3.04 |
| 步长 | 0.81 | 0.83 | 0.82 | 0.80 | 0.81 | +0.00 |
| 迭代状态更新 | 0.80 | 0.75 | 0.72 | 0.87 | 0.84 | +0.04 |
| GPL 其他 | 0.28 | 0.31 | 0.29 | 0.28 | 0.32 | +0.04 |
| **合计** | **81.85** | **108.51** | **78.06** | **75.03** | **77.06** | **-4.79** |

变快合计 -9.83 s；变慢合计 +5.05 s；净 -4.79 s。
最终版中 SoA 未触及的阶段（时序反馈、初始布局、FFT、步长、拥塞反馈）合计 28.42 s，占 36.9%。

## 原因分析

### 1. 存储换成了 SoA，访问方式没变

`SoaRow::get<F>()` 的执行路径是：先读对象里的 `store_`，再取该列 `vector` 的首地址，
最后读 `[slot]`。原版 `GCell` 的字段内联在对象内（约 96 字节，即 1.5 条 cache line），
一次遍历就是一条顺序流。改造后，一次"按对象"访问要先读句柄，再分别访问 k 个列，
即 k 条独立的流，外加 store 头部。列存对向量化有利的前提是循环直接遍历行号并只读少数列；
如果仍按对象逐个调用 getter，就只剩下额外的间接寻址。

"仅字段 SoA"一列直接证明了这一点：算法完全不变，只换存储，GPL 从 81.85 s 涨到 108.51 s。
`getOverlapDensityArea` 每次访问 bin 时，要分别经过 cell 行和 bin 行读 8 个字段，
自身时间从 2.33 s 涨到 13.06 s。

最终版只把密度 scatter/gather 和线长 CSR 改成了直接读列，其余路径仍是对象式访问。
坐标更新的热点（自身时间，单位秒）如下：

| 函数 | 原版 | 仅字段 SoA | 最终 SoA |
|---|---:|---:|---:|
| `gpl::NesterovBase::getDensityCoordiLayoutInsideX const` | 0.71 | 1.68 | 2.57 |
| `gpl::NesterovBase::getDensityCoordiLayoutInsideY const` | 1.13 | 1.93 | 1.53 |
| `gpl::GCell::setDensityCenterLocation` | 0.59 | 1.02 | 0.90 |
| `gpl::GPin::updateDensityLocation` | 0.48 | 1.10 | 0.86 |
| `gpl::NesterovBase::nesterovUpdateCoordinates` | 1.11 | 0.83 | 0.97 |
| `gpl::GCell::isLocked const` | 0.38 | 0.29 | 0.48 |
| `gpl::Instance::isLocked const` | 0.12 | 0.07 | 0.31 |

`nesterovUpdateCoordinates` 对每个单元都要：经句柄检查 `isLocked()`
（`insts` 列 → `Instance*` → Instance 行 → `is_locked` 列），再在 `InsideX/Y` 中经句柄读
`dLx/dUx/dLy/dUy` 四列。`setDensityCenterLocation` 写四列后，对每个 pin 调用
`GPin::updateDensityLocation(this)`，这一步又跨越 pin 行和 cell 行。这些都是"每单元若干
次间接寻址"的路径，原版里它们只访问对象内的字段。

`isLocked` 还多了一层：flag 现在由 `Instance` 行保存，必须先经过 `GCell` 的 `insts` 列
（每个元素是 24 字节的 `std::vector` 头）再跳转。`isInstance()/isMacroInstance()` 也走同一条链；
这条链在 scatter 中每次调用对每个单元都要走一遍（`GCell::isInstance` 0.86 s，`Instance::isMacro` 0.19 s）。

### 2. 密度内核：提速来自去掉对象访问，常驻并非关键

单线程，单位为秒：

| 项目 | adapter | 常驻 SoA |
|---|---:|---:|
| scatter 内层循环 | 8.76 | 9.31 |
| gather 内层循环 | 8.51 | 9.58 |
| 每次调用打包 AoS→SoA | 5.56 | 0.00 |
| overflow/清零循环（updateBinsGCellDensityArea 自身） | 1.36 | 2.68 |

原版 updateBinsGCellDensityArea 自身（清零 + overflow）1.28 s。

- 内层循环速度相同，说明密度上的收益（scatter −6.98 s，gather −2.16 s）来自三方面：
  不再经 `Bin` 对象读写、去掉单线程原子加、去掉逐 visit 的 `getOverlapDensityArea` 调用和句柄解析。
  这些与数据是否常驻无关。
- 常驻方案省下的 5.56 s 打包时间被其余路径的回退抵消：坐标、线长和梯度合并四段，
  最终 SoA 比 adapter（这些路径仍是原版代码）多 4.45 s。
- 两种方案的内层循环仍然偏慢。scatter 每 visit 约 12 ns（9.31 s / 1261 次 / ~58 万 visit），
  而独立微基准（`replay-final.jsonl`，同一快照）的 serial scatter 约 7.3 ns/visit。
  每个 visit 要从 4 个 bin 边界列读数（`bin_lx/ly/ux/uy`），做 min/max 和 float 乘法，
  并按单元可变的 1–4 次循环分支。gather 每处理一个单元，都会在 `getDensityGradient`
  中重新构造一次 `densityView()`（0.61 s）和 `densityRect()`（0.63 s）。

### 3. 新引入的布局冲突：x-major 存储、y-major 遍历

`BinGrid` 按 `x * ny + y` 分配行，以对齐 FFT；但 `updateBinsGCellDensityArea` 的 overflow
循环为保留 float 归约顺序，按 `legacy = y * nx + x` 遍历，再转换为
`b = (legacy % nx) * ny + legacy / nx`。每次相邻迭代在约 11 列上各跨越 `ny` 个元素，
几乎每次访问都落在新的 cache line 上，还要做两次运行期除法。这个循环自身时间
1.28 → 2.68 s。`BinGrid::num_threads_` 恒为 1（源码注释已说明 `setNumThreads` 没有调用方），
因此它在 8 线程下也是串行。

### 4. 8 线程：原版的并行度已经覆盖了 SoA 改动的部分

单位为秒，数据来自主线程墙钟采样：

| 阶段 | 原版 8 线程 | 最终 SoA 8 线程 | 最终−原版 |
|---|---:|---:|---:|
| 时序反馈 | 17.73 | 17.85 | +0.12 |
| 拥塞反馈 | 0.77 | 0.75 | -0.01 |
| 初始布局 | 1.33 | 1.41 | +0.08 |
| 密度 scatter + overflow | 8.31 | 5.36 | -2.94 |
| FFT/DCT | 5.15 | 5.05 | -0.10 |
| 密度场其他 | 0.30 | 0.21 | -0.09 |
| 密度 gather | 2.84 | 1.71 | -1.14 |
| 线长 WA force | 0.54 | 0.79 | +0.25 |
| 线长梯度 gather | 2.28 | 2.19 | -0.10 |
| 梯度合并/预条件 | 0.95 | 0.97 | +0.02 |
| HPWL | 0.21 | 0.11 | -0.09 |
| 坐标更新 + pin 同步 | 4.11 | 3.77 | -0.35 |
| 步长 | 0.83 | 0.86 | +0.03 |
| 迭代状态更新 | 0.18 | 0.21 | +0.03 |
| GPL 其他 | 0.27 | 0.31 | +0.04 |
| **合计** | **45.84** | **41.60** | **-4.24** |

最终 SoA 8 线程中串行或基本串行的部分：

| 部分 | 秒 | 占 8 线程总时间 |
|---|---:|---:|
| 时序反馈（sta 线程 = 1） | 17.85 | 42.9% |
| FFT/DCT（Ooura，未开线程） | 5.05 | 12.1% |
| overflow 循环（BinGrid::num_threads_ 恒为 1） | 3.46 | 8.3% |
| 密度中心/pin 同步循环（无 OpenMP） | 2.95 | 7.1% |
| 步长 getDistance（串行） | 0.86 | 2.1% |
| 初始布局 | 1.41 | 3.4% |
| **合计** | **31.59** | **75.9%** |

8 线程下，已并行的 scatter、gather 和线长合计只剩约 10 s，SoA 能改善的余地很小。
剩余时间几乎全由串行段决定。`fftsg2d.cpp` 没有定义 `USE_FFT2D_PTHREADS`，因此 FFT 是单线程。
`updateGCellDensityCenterLocation` 中设置单元中心和同步 pin 的循环没有 OpenMP。
`getDistance` 是串行归约。

### 5. 不是原因的部分

- **FFT 拷贝。** `solveColumns` 去掉了场的转置拷贝，FFT 只减少 0.46 s；DCT 本身才是主要耗时。
- **线长内核。** CSR 版的 `updateNetForceColumns`（3.04 s）与原版 `updateWireLengthForceWA_native`
  （2.11 s 自身，另加 `clearWaVars`/`updateBox` 等）量级相同。WA 计算以 `fastExp` 和浮点运算为主，
  不受访存限制，数据布局对它影响很小。线长 gather 里，每个单元还要检查一次
  `utl::Logger::debugCheck`，自身耗时 0.37 s。

## 其他加速方向

这些方向按"已有证据支撑的上限"排序。"保持逐位一致"指坐标哈希能否保持不变，
无法保持的都需要另做 QoR 评估。

由上表推出的预算：

- 坐标/线长/梯度合并四段，最终 SoA 比 adapter（这些路径仍是原版代码）多 4.45 s。
- 8 线程下 FFT、overflow、密度中心同步、步长四段串行合计 12.33 s；理想 8 线程并行可省 10.79 s，总时间 41.60 → 30.82 s。
- 时序反馈与初始布局不变、其余 GPL 耗时全部归零的上限：单线程 77.06 / 21.72 = 3.55×；8 线程 41.60 / 19.26 = 2.16×。

迭代构成（来自原版单线程日志）：

- 快照保存于第 327 轮；共 10 次 routability pass，最后一次 pass 前日志最后一行迭代号 910；最终迭代号 1057。
- 快照到最后一次 pass 之间的迭代（含被回退的重跑）约 583 轮，占 1058 轮的 55.1%。
- 最终 SoA 单线程 Nesterov 主循环（扣除时序/初始/拥塞反馈/其他）54.15 s，按轮均摊约 51.2 ms/轮；上述重跑区间按均摊折算约 29.8 s。

时序反馈内部（最终 SoA 单线程，各项 inclusive 秒数，相互有重叠）：

- `rsz::Resizer::findResizeSlacks`：17.61 s
- `rsz::Rebuffer::fullyRebuffer`：7.90 s
- `rsz::RepairDesign::repairDesign`：5.27 s
- `sta::GraphDelayCalc::findDelays`：4.36 s
- `est::EstimateParasitics::estimateParasitics`：2.35 s
- `sta::Search::findRequireds`：1.10 s
- `sta::Search::findAllArrivals`：0.99 s

### A. 算法层：routability 回退重跑（最大单项）

AES 上，10 次 routability pass 每次都从第 327 轮的快照重新收敛约 60 轮。
快照之后到最后一次 pass 的约 583 轮占总轮数的 55.1%，按均摊折算约 29.8 s（单线程）。
它在时间上远大于数据布局的影响，但这段迭代不能全部省掉，因为 inflation 生效后需要重新收敛。
可以尝试：

- 自适应增大每次 inflation 步长，减少 pass 次数（目标是 weighted congestion 1.01；
  从 pass 7 起每次只降 0.01–0.04）；
- 不完全回退，在当前解上直接增加 inflation 并继续迭代（warm start）；
- 改进停止条件：pass 间改善低于阈值时提前结束。

这一项会改变布局结果，必须用 QoR（HPWL、拥塞、后续 GRT/DRT 溢出、WNS/TNS）验收，
不能只看坐标哈希。另有一项开销较小的冗余：回溯使候选计算达 1247 次，比 1058 轮多 18%。

### B. 时序反馈：8 线程下占 42.9%

- 本次测试固定了 `sta::set_thread_count 1`。打开 OpenSTA 的多线程 delay calc 和 search
  后，`findDelays`（4.36 s）、arrival/required 等部分可以并行。真实流程可能本来就不设为 1，
  需要单独测量；OpenSTA 多线程本身应保持确定性，但必须核对哈希。
- `fullyRebuffer`（7.90 s）和 `repairDesign`（5.27 s）是 rsz 的串行算法。可考虑：
  在 GPL 内部只做轻量的 virtual repair、减少 timing-driven 次数（目前为 2 次），
  或复用上一轮的寄生参数和 Steiner 树。这些都会改变 QoR。
- 可以与 FlashEDA 在线 STA 的工作衔接，但收益取决于 repair 本身而不只是 STA；
  此前已注意到 17 s 不能全部算作 STA。

### C. 补齐串行段（8 线程收益最大，且可保持逐位一致）

上限是理想 8 线程并行可省 10.79 s（41.60 → 30.82 s）。

1. **overflow 循环。** 先按 x-major 顺序并行写出每个 bin 的 overflow 到一个连续 float 数组，
   再按原 y-major 顺序串行归约这一个数组。这样既保留归约顺序（逐位一致），又消除了 11 列上的
   跨步访问；同时修复 `BinGrid::num_threads_` 恒为 1 的问题。
2. **密度中心/pin 同步。** 改为按行号的列循环：单元中心计算 `dLx…dUy` 可直接向量化；
   pin 坐标用 DREAMPlace 的 `pin_pos` 方式，即 `pin_x[p] = cell_x[pin2cell[p]] + offset[p]`，
   基于现有的 cell→pin CSR，可以 OpenMP 并行，结果逐位一致。
3. **FFT。** 先尝试启用 Ooura 自带的 `USE_FFT2D_PTHREADS`：它按列/行切分，每列的计算不变，
   有望逐位一致，但需要验证。也可以对 DCT 行/列做批量 OpenMP 并行。换成 FFTW 或 pocketfft
   会改变浮点结果。
4. **`getDistance` 归约。** 可用固定分块的确定性树归约，但会改变 float 求和顺序；
   或者先并行求每元素平方，再串行求和，以保持结果不变。

### D. 把 SoA 做完：让访问方式也变成列式（单线程收益为主）

第一项能找回约 4.45 s 的回退，后几项有机会让密度内核降到微基准水平：

1. **热循环直接遍历行号。** 去掉 `GCellHandle` 的 variant 分发和 `getFillerGCell`，
   改为对 instance 池和 filler 池各循环一次，并直接读写 `PointColumns` 的 x/y 数组。
   坐标更新、梯度合并和步长都按这种方式改写。
2. **单元类型列。** 增加 1 字节的 kind 列（std/macro/filler/skip），在拓扑变化时更新，
   取代每次调用都走一遍 `insts → Instance → is_macro/is_locked` 的链。
3. **把 per-cell 不变量提到循环外。** `densityView()` 每个 kernel 只构造一次；
   `debugCheck` 每个 kernel 只检查一次。
4. **可分离的重叠计算。** bin 边界由 `lround(x * binSize)` 决定，x 方向只依赖 x，
   y 方向只依赖 y。可用 `nx+1` 和 `ny+1` 两个一维边界数组取代 4 个二维边界列。
   每个单元先算出 ≤4 个 x 方向重叠宽度和 ≤4 个 y 方向重叠高度，内层只剩
   `float(w) * float(h) * scale` 和截断累加。原版本来就是 `float(width) * float(height)`，
   所以结果逐位一致。这去掉了每 visit 的 4 次列读取和 4 次 min/max，
   可在现有独立微基准（`integrations/openroad/density-cpu`）上先验证。
5. **算子融合。** 坐标更新后直接 scatter；gather 直接写入梯度合并（包括预条件），
   每轮少遍历 10 万单元 1–2 次。gather 的 x-outer/y-inner 求和顺序不变。

### E. 构建配置（成本最低，建议先测）

测得的二进制使用 `-O3 -DNDEBUG`，未加 `-march`，LTO 关闭。热路径上有不少跨 TU 的
非内联调用，例如 `GCell::isInstance` 定义在 `.cpp` 中，`densityGradient.cpp` 又调用
`NesterovBase::getDensityGradient`。可以对 gpl 单独启用 LTO、`-march=icelake-server`
（AVX-512）和 PGO。注意 `-march` 可能引入 FMA 收缩，从而改变浮点结果：
需要 `-ffp-contract=off` 才能保持逐位一致。

### F. GPU

fork 中已有 GPU 路径，但未验收。按上面的预算，即使时序反馈与初始布局之外的 GPL 耗时全部
归零，上限也只有单线程 3.55×、8 线程 2.16×。只有与 A、B 一起推进，GPU 才值得投入。

## 建议顺序

1. **E**（构建参数）和 **C1/C2**（overflow 循环、中心/pin 同步）：改动小、可保持逐位一致，
   直接针对 8 线程的串行段。
2. **D4** 可分离重叠：先在独立微基准上验证，再接入；同时做 **D1–D3**，消除 4.45 s 的回退。
3. **C3** FFT 线程化，并验证是否逐位一致。
4. **B**：先测 `sta::set_thread_count N` 的实际收益，再评估是否减少 rebuffer/repair。
5. **A** routability 回退策略：收益最大，但会改变 QoR，需要单独做 QoR 实验，
   并覆盖多个设计（不只 AES）。

每一步都用现有的差分验收（GCD/AES 加 8 个扩展案例，1/8 线程）确认；
不能保持逐位一致的方向改用 QoR 对照。

## 更真实用例上的复测

AES 的利用率只有 9.4%，约 8 万个 filler 放大了密度阶段的占比。为此又在 OpenROAD 自带的
更大用例上复测（Northjob `flasheda-soa-designs-20261106`，identity
`job-822bb6b8b2e9-20261106090758`；机器和环境变量同上，采样周期 2 ms）：

- `ibex_core`：沿用 AES 的 timing+routability 脚本；
- `medium04`、`large01`、`large02`：上游测试的默认 `global_placement`；
- `macro03`（swerv_wrapper，含宏单元）：`-density 0.7`；
- `large01`、`macro03` 另跑 routability-driven 变体，routing layer 设置与 AES 相同。

每个配置只跑一次，1/8 线程均为原版与最终 SoA 对照。

| 配置 | 实例 | 可移动 | net | 利用率% | 目标密度 | bin X | 末轮迭代 | RD pass | 线程 | 原版 s | SoA s | 减少 | 哈希一致 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| ibex_core-td-rd | 34184 | 29301 | 33171 | 8.3 | 0.50 | 256 | 1443 | 18 | 1 | 178.8 | 160.8 | 10.1% | True |
| ibex_core-td-rd | 34184 | 29301 | 33171 | 8.3 | 0.50 | 256 | 1443 | 18 | 8 | 95.4 | 83.8 | 12.2% | True |
| medium04 | 135119 | 127555 | 139550 | 10.9 | 0.70 | 512 | 473 | 0 | 1 | 228.5 | 198.3 | 13.2% | True |
| medium04 | 135119 | 127555 | 139550 | 10.9 | 0.70 | 512 | 473 | 0 | 8 | 85.5 | 81.5 | 4.7% | True |
| macro03 | 131271 | 131271 | 143973 | 27.1 | 0.70 | 512 | 466 | 0 | 1 | 196.5 | 175.6 | 10.6% | True |
| macro03 | 131271 | 131271 | 143973 | 27.1 | 0.70 | 512 | 466 | 0 | 8 | 75.0 | 72.1 | 3.9% | True |
| macro03-rd | 131271 | 131271 | 143973 | 27.1 | 0.70 | 512 | 1266 | 7 | 1 | 464.0 | 428.5 | 7.6% | True |
| macro03-rd | 131271 | 131271 | 143973 | 27.1 | 0.70 | 512 | 1266 | 7 | 8 | 195.0 | 176.3 | 9.6% | True |
| large01 | 274700 | 274700 | 290354 | 65.9 | 0.70 | 512 | 572 | 0 | 1 | 220.2 | 219.1 | 0.5% | True |
| large01 | 274700 | 274700 | 290354 | 65.9 | 0.70 | 512 | 572 | 0 | 8 | 92.8 | 104.6 | -12.7% | True |
| large01-rd | 274700 | 274700 | 290354 | 65.9 | 0.70 | 512 | 1279 | 7 | 1 | 471.2 | 467.8 | 0.7% | True |
| large01-rd | 274700 | 274700 | 290354 | 65.9 | 0.70 | 512 | 1279 | 7 | 8 | 202.1 | 227.2 | -12.4% | True |
| large02 | 312529 | 312529 | 406912 | 63.6 | 0.70 | 512 | 530 | 0 | 1 | 251.5 | 257.8 | -2.5% | True |
| large02 | 312529 | 312529 | 406912 | 63.6 | 0.70 | 512 | 530 | 0 | 8 | 113.3 | 126.8 | -11.8% | None |

不带采样器的 large02 8 线程重跑：baseline 114.1 s（哈希 9c1928637a2c）；soa-final-openroad 123.4 s（哈希 9c1928637a2c）
每个配置所有运行（两种二进制、1/8 线程、重跑）的不同坐标哈希数：ibex_core-td-rd 1，medium04 1，macro03 1，macro03-rd 1，large01 1，large01-rd 1，large02 1。

- 坐标哈希：每个配置的所有运行（两种二进制、1/8 线程、重跑）只有一个哈希，SoA 与原版逐位一致。
- `large02` 原版 8 线程首次运行在启动约 2 ms 时以 SIGSEGV 退出，日志为空，采样文件只有 1 条记录。
  同配置带采样重跑一次、不带采样重跑一次都正常，因此判断为采样器在信号处理函数里调用
  `backtrace()`（非异步信号安全）导致的偶发崩溃，不是 OpenROAD 的问题；该次结果已用带采样的重跑替代。
- 不带采样器时，`large02` 8 线程 SoA 也更慢（见上方重跑行），说明变慢不是采样造成的。

### 结论随利用率变化

- **低利用率设计**（ibex_core、medium04、macro03，利用率 8–27%）：filler 多，密度 scatter/gather
  占原版一半左右，SoA 在单线程下减少约 8–13%。这与 AES 的机制相同，只是幅度更大。
- **高利用率设计**（large01、large02，利用率约 64–66%，最接近真实流程）：单线程持平或略慢，
  8 线程反而**变慢约 12%**。密度阶段只占约 20%，主要时间转移到线长、HPWL、pin 同步和初始布局；
  而这些恰好是 SoA 改造中仍走对象式访问、或新引入了并行问题的部分。

### large01，1 线程（秒；括号内为原版占比）

| 阶段 | 原版 | SoA | 差 |
|---|---:|---:|---:|
| 初始布局 | 42.81 (19.5%) | 46.55 | +3.74 |
| 密度 scatter | 29.96 (13.6%) | 28.45 | -1.51 |
| FFT | 12.06 (5.5%) | 10.02 | -2.04 |
| 密度场其他 | 1.97 (0.9%) | 2.20 | +0.23 |
| 密度 gather | 17.98 (8.2%) | 15.25 | -2.73 |
| 线长 force | 24.48 (11.1%) | 24.82 | +0.33 |
| 线长 gather | 54.30 (24.7%) | 47.28 | -7.02 |
| 梯度合并 | 3.38 (1.5%) | 3.92 | +0.54 |
| HPWL | 11.47 (5.2%) | 5.39 | -6.08 |
| 坐标/pin 同步 | 16.67 (7.6%) | 29.85 | +13.18 |
| 步长 | 1.06 (0.5%) | 1.06 | +0.00 |
| 迭代更新 | 2.64 (1.2%) | 2.52 | -0.13 |
| Nesterov 其他 | 0.11 (0.1%) | 0.13 | +0.01 |
| GPL 其他 | 1.07 (0.5%) | 1.51 | +0.44 |
| **合计** | **219.96** | **218.95** | **-1.01** |

### large01，8 线程（秒；括号内为原版占比）

| 阶段 | 原版 | SoA | 差 |
|---|---:|---:|---:|
| 初始布局 | 21.50 (23.2%) | 25.42 | +3.93 |
| 密度 scatter | 18.65 (20.1%) | 17.74 | -0.91 |
| FFT | 10.22 (11.0%) | 9.99 | -0.23 |
| 密度场其他 | 0.52 (0.6%) | 0.72 | +0.20 |
| 密度 gather | 3.15 (3.4%) | 2.67 | -0.48 |
| 线长 force | 5.65 (6.1%) | 10.63 | +4.97 |
| 线长 gather | 9.00 (9.7%) | 9.13 | +0.13 |
| 梯度合并 | 1.22 (1.3%) | 1.29 | +0.07 |
| HPWL | 2.87 (3.1%) | 0.96 | -1.91 |
| 坐标/pin 同步 | 17.46 (18.8%) | 23.08 | +5.62 |
| 步长 | 1.05 (1.1%) | 1.08 | +0.03 |
| 迭代更新 | 0.41 (0.4%) | 0.31 | -0.10 |
| Nesterov 其他 | 0.11 (0.1%) | 0.12 | +0.01 |
| GPL 其他 | 0.89 (1.0%) | 1.33 | +0.44 |
| **合计** | **92.70** | **104.48** | **+11.77** |

### large02，8 线程（秒；括号内为原版占比）

| 阶段 | 原版 | SoA | 差 |
|---|---:|---:|---:|
| 初始布局 | 28.93 (25.6%) | 34.06 | +5.13 |
| 密度 scatter | 20.57 (18.2%) | 18.59 | -1.98 |
| FFT | 9.51 (8.4%) | 9.38 | -0.13 |
| 密度场其他 | 0.54 (0.5%) | 0.68 | +0.14 |
| 密度 gather | 4.63 (4.1%) | 2.97 | -1.65 |
| 线长 force | 7.54 (6.7%) | 14.41 | +6.87 |
| 线长 gather | 12.94 (11.4%) | 11.49 | -1.44 |
| 梯度合并 | 1.46 (1.3%) | 1.53 | +0.07 |
| HPWL | 3.47 (3.1%) | 1.14 | -2.33 |
| 坐标/pin 同步 | 20.39 (18.0%) | 28.84 | +8.45 |
| 步长 | 1.19 (1.1%) | 1.17 | -0.02 |
| 迭代更新 | 0.44 (0.4%) | 0.39 | -0.05 |
| Nesterov 其他 | 0.17 (0.2%) | 0.16 | -0.01 |
| GPL 其他 | 1.29 (1.1%) | 1.67 | +0.37 |
| **合计** | **113.06** | **126.48** | **+13.42** |

### 大设计上的新问题

### large01

| 函数 | large01 1T 原版 | large01 1T SoA | large01 8T 原版 | large01 8T SoA |
|---|---:|---:|---:|---:|
| pin 同步 GPin::updateDensityLocation | 9.43 | 19.31 | 12.19 | 18.81 |
| overflow 循环 updateBinsGCellDensityArea 自身 | 3.19 | 9.64 | 3.66 | 9.94 |
| scatter 内层 scatterDensityAreaInPlace 自身 | 15.41 | 11.17 | 12.43 | 3.37 |
| 线长 force 计算（原版/CSR 内核自身） | 17.39 | 23.56 | 3.87 | 5.34 |
| 线长 force 中 OpenMP 等待（libgomp） | 0.00 | 0.00 | 0.10 | 4.14 |
| 线长 force 中 libc（SoA 为串行清零 8 个 pin 列） | 0.00 | 0.94 | 0.00 | 1.09 |
| 初始布局中 Instance/Pin getter | 0.91 | 1.87 | 0.83 | 1.97 |

### large02

| 函数 | large02 1T 原版 | large02 1T SoA | large02 8T 原版 | large02 8T SoA |
|---|---:|---:|---:|---:|
| pin 同步 GPin::updateDensityLocation | 10.81 | 24.04 | 15.45 | 24.25 |
| overflow 循环 updateBinsGCellDensityArea 自身 | 2.91 | 9.00 | 3.90 | 9.29 |
| scatter 内层 scatterDensityAreaInPlace 自身 | 16.41 | 11.66 | 13.98 | 4.49 |
| 线长 force 计算（原版/CSR 内核自身） | 22.61 | 29.80 | 5.02 | 7.04 |
| 线长 force 中 OpenMP 等待（libgomp） | 0.00 | 0.00 | 0.48 | 5.88 |
| 线长 force 中 libc（SoA 为串行清零 8 个 pin 列） | 0.00 | 1.21 | 0.00 | 1.43 |
| 初始布局中 Instance/Pin getter | 1.02 | 2.50 | 1.05 | 2.49 |

1. **pin 同步翻倍。** `GPin::updateDensityLocation` 在 large01/large02 上单线程从约 10 s 涨到约 20 s。
   这些设计有 29–41 万条 net（large02 共 140 万个 pin），每次坐标更新都要对每个 pin 经过 pin 行、
   cell 行做对象式读写，而且这个循环本来就没有 OpenMP，8 线程下照样是串行。
   这是第 1 条原因（存储换了，访问没换）在 pin 数量很大时的放大版。
2. **线长 force 在 8 线程下变慢。** CSR 内核自身只多了约 1.5–2 s，但 OpenMP 等待从约 0.1–0.5 s
   涨到 4–6 s，另有约 1 s 的串行清零（SoA 版在并行循环前用 `std::fill` 串行清 8 个 pin 列，
   原版是并行清零的）。最可能的原因是**伪共享**：各线程按 net 并行写 pin 列，同一 net 的 pin
   在列中不连续，而 `hasMinExpSumX` 等是 1 字节列（64 个 pin 共享一条 cache line），float 列每条
   cache line 16 个 pin，不同线程写相邻 pin 时互相失效；原版每个 pin 是独立对象，很少共享 cache line。
   这一点没有硬件计数器证实，属于根据源码和等待时间的推断。
3. **overflow 循环涨到约 3 倍。** 512×512 网格下，y-major 遍历 x-major 列的跨步是 512，
   每列 1–2 MiB 放不进 L2；这个循环又是串行的。large01/large02 上从约 3 s 涨到约 9–10 s。
4. **初始布局变慢。** `Instance`/`Pin` 的 getter 在初始布局中从约 1 s 涨到约 2–2.5 s，
   加上 `createSparseMatrix`、OpenDB 访问的零散增加，8 线程下初始布局多出约 4–5 s。
5. **scatter 内层收益仍然成立。** 去掉原子加和私有网格并行后，8 线程 scatter 内层从约 12–14 s
   降到约 3.4–4.5 s；但被同一阶段里变慢的 overflow 循环抵消了大半。

### 时序与 routability 在真实设置下的占比

ibex_core 是唯一带 SDC 的大用例。时序反馈在单线程下约占 19%，8 线程下约占 35%（见下表），
与 AES 的结论一致：8 线程时它是最大的串行段。

### ibex_core-td-rd，8 线程（秒；括号内为原版占比）

| 阶段 | 原版 | SoA | 差 |
|---|---:|---:|---:|
| 时序反馈 | 33.22 (35.0%) | 34.24 | +1.02 |
| 拥塞反馈 | 2.11 (2.2%) | 2.12 | +0.01 |
| 初始布局 | 1.97 (2.1%) | 2.22 | +0.25 |
| 密度 scatter | 20.00 (21.1%) | 9.33 | -10.67 |
| FFT | 7.13 (7.5%) | 6.98 | -0.15 |
| 密度场其他 | 0.37 (0.4%) | 0.29 | -0.08 |
| 密度 gather | 5.29 (5.6%) | 4.03 | -1.26 |
| 线长 force | 1.36 (1.4%) | 1.79 | +0.43 |
| 线长 gather | 5.79 (6.1%) | 6.40 | +0.61 |
| 梯度合并 | 2.05 (2.2%) | 2.36 | +0.31 |
| HPWL | 0.59 (0.6%) | 0.27 | -0.33 |
| 坐标/pin 同步 | 12.12 (12.8%) | 10.20 | -1.92 |
| 步长 | 2.01 (2.1%) | 1.94 | -0.07 |
| 迭代更新 | 0.46 (0.5%) | 0.56 | +0.10 |
| Nesterov 其他 | 0.10 (0.1%) | 0.13 | +0.03 |
| GPL 其他 | 0.44 (0.5%) | 0.50 | +0.06 |
| **合计** | **95.01** | **83.36** | **-11.65** |

开启 routability 后，macro03 和 large01 的末轮迭代号都增加到约 2.2–2.7 倍（D1 表"末轮迭代"列），
总时间也相应翻倍。"routability 回退重跑"这一算法层面的杠杆在大设计上同样存在。

### 修正后的优先级

结合大设计数据，建议把顺序调整为：

1. **pin 同步改为列式并行**（DREAMPlace `pin_pos`：`pin_x[p] = cell_x[pin2cell[p]] + offset[p]`）。
   高利用率设计上这是 SoA 变慢的最大单项，可保持逐位一致。
2. **线长 force 去伪共享**：并行清零；按 net 的 CSR 顺序重排 pin 行，让同一 net 的 pin 连续，
   或在线程内先算到局部缓冲再写回；把 1 字节的 `has*` 标志并入 float 列（例如用 ±inf 表示未激活）。
3. **overflow 循环**：按前文 C1 的方法 x-major 并行加数组归约，同时修复 `num_threads_` 恒为 1。
4. 之后再做前文的 D 类（热循环直接遍历行号、可分离重叠）、FFT 线程化、构建参数，以及时序反馈和
   routability 策略。

只有修完 1–3，SoA 在高利用率设计上才不会是负收益；在此之前不宜把这个分支作为默认实现。

## 局限

- 前半部分只基于 AES；复测补充了 OpenROAD 自带的 6 个更大配置，但它们仍是测试用例。
  large/medium 用例只有 LEF/DEF，没有时序约束；只有 ibex_core 带 SDC。ORFS 从 RTL 综合出的设计
  （如 jpeg、ariane133、ASAP7 设计）尚未测试。
- 每个配置只采样一次；阶段级差值小于约 0.3 s（AES）或约 1 s（大设计）时不应解读。AES 的
  时间差已由此前每种 3 次的计时确认，大设计只有单次运行。
- 没有硬件计数器，cache/分支相关的解释是推断。
- 时序反馈使用 `sta` 1 线程，与真实流程的设置可能不同。

## 证据与复现

- 采样器、符号化、分阶段和表格生成脚本：`build/soa/diagnose/`
  （`sampler.c`、`pcount.c`、`fold.py`、`analyze.py`、`tables.py`、`profile.py`、`commands/`、`job.yaml`）。
- 折叠栈：`build/soa/diagnose/folded/*.txt`（可直接用 flamegraph 渲染）；
  分阶段 JSON：`phases-t1.json`、`phases-t8.json`；表格：`tables.md`。
- 证据包：`build/soa/diagnose/evidence.tar.gz`，SHA-256
  `a8c779d152196195f4d353b50916511c131e17506cbeb038a63ccb15353cdfa7`
  （不含原始 `samples.bin` 和 DEF/Verilog；原始数据保留在 GPFS
  `/gpfs/users/linlizhi/flasheda-soa-diagnose-20261106/`）。
- 被测二进制使用 `/gpfs/users/linlizhi/flasheda-density-cpu-20261106/` 中的原始文件，未重新编译。
- 作业已经 STOP，worker 以 0 退出，状态为 Completed，已按精确 alias 删除，复查返回 `[]`。
- 大设计复测：脚本 `build/soa/designs/`（`run_designs.py`、`rerun.py`、`gpl_only.tcl`、
  `report.py`、`append_doc.py`、`section.md`），表格 `build/soa/designs/tables.md`，折叠栈
  `build/soa/designs/folded/`。证据包 `build/soa/designs/evidence.tar.gz`，SHA-256
  `fc399617ff469cb2b40ce24fdf7c2bea19d2bb5993d1761e8a3ec0969c90d901`；原始采样保留在
  GPFS `/gpfs/users/linlizhi/flasheda-soa-designs-20261106/`。作业同样已 STOP、Completed、
  按精确 alias 删除，复查为 `[]`。
