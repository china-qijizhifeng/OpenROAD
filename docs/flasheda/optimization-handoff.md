# GPL 优化 handoff

> 历史快照（2026-11-06 11:10 CST）。分支合并、构建方式和资源状态以下方时间为准；
> 后续行存储、初始布局与线长 gather 的实现及验收见
> [SOA_ROW_STORAGE.md](validation/SOA_ROW_STORAGE.md)。

更新时间：2026-11-06 11:10 CST。面向下一位接手 OpenROAD GPL 性能工作的 agent 或工程师。
本文只汇总状态和下一步；数据与论证见下列文档，数字以它们为准：

- 诊断：[soa-speedup-diagnosis.md](soa-speedup-diagnosis.md)（AES + 6 个更大配置的采样归因）
- 整改与验收：[validation/SOA_REMEDIATION.md](validation/SOA_REMEDIATION.md)（Codex 会话编写，仍在更新）
- 设计说明：[gpl-soa.md](gpl-soa.md)；首轮验收：[validation/RESULTS.md](validation/RESULTS.md)

## 1. 现状

| 阶段 | 内容 | AES global_placement 中位数（1 / 8 线程） | 状态 |
|---|---|---|---|
| 原版 | OpenROAD `18216f3506` | 约 79.3–79.9 s / 44.0–44.7 s | 基线 |
| 首版 SoA | `dd8f2fc`，七类实体全字段列存 | 75.0 s / 41.3 s | 已提交并推送 |
| 整改版 | 诊断中的 C1–C4、D1–D5、E | 62.6 s / 30.3 s（比原版少 21.6% / 31.1%） | **未提交**，大设计尚未验收 |

- 所有版本在已验收案例上都与原版**坐标逐位一致**。这是硬约束（见第 6 节）。
- 诊断发现，首版 SoA 在高利用率大设计（`large01`/`large02`，利用率约 64–66%）上
  8 线程**慢约 12%**。整改版是否修复了这一点**还没测**，这是当前最重要的未知数。

## 2. 并发注意（先读）

截至本文时间，Codex 会话 `01a1ae5e-ed76-7ee0-bbaf-c27834772ada`（主）和
`01a1aeba-2951-7f43-a5ed-b074f0bb5bde`（subagent，监控构建与验收）**正在**修改
`/Users/shom/Downloads/openroad-gpl-soa` 的工作区，并运行 Northjob 作业
`flasheda-soa-remediation-20261106`（identity `job-317d12d0424e-20261106090444`）。

- 工作区有 19 个已修改文件（`src/gpl/src/*`、`src/gpl/test/soa/*`、BUILD/CMake），以及新增的
  `fft_test.cpp`、`validation/SOA_REMEDIATION.md`、`validation/soa-remediation/`，全部未提交。
- 接手前先确认这两个会话已结束：查看 `~/.codex/sessions/2026/11/06/rollout-*-<id>.jsonl` 的最后一个
  `task_complete`/`turn_aborted` 事件和文件修改时间，并用 `git status` 检查工作区。不要在别人写入时编辑同一文件。
- 不要停止或删除不是自己创建的 Northjob 作业；删除时只按精确 alias 操作。

## 3. 代码与环境

**仓库**

- Fork：<https://github.com/china-qijizhifeng/OpenROAD>，分支 `flasheda/gpl-soa`，最新推送 `dd8f2fc`。
  默认分支仍是上游 `master`，未合并。
- 本地：`/Users/shom/Downloads/openroad-gpl-soa`。这是部分克隆加 sparse checkout，
  只检出 `src/gpl/src`、`src/gpl/include`、`docs/agents`。新增 sparse 范围外的文件需要 `git add --sparse`。
- 提交规则（仓库 AGENTS/CLAUDE.md）：`git commit -s`；C++ 改动先跑 clang-format
  （本机可用 `~/Downloads/flasheda-openroad-rudy-20261105/build/density-cpu/format-tools/clang_format/data/bin/clang-format`，v18.1.8）；
  不要修改 `src/sta`；新测试在 Bazel 注册，新源文件仍需加 CMake 条目。

**远端（BP Dev Northjob）**

- 提交与清理用 `~/.codex/skills/northjob-gpu-test-bp-dev/scripts/northjob.py`（`jobs/submit/delete`）。
  项目 `WeaverRSI`，队列 `proj-32`，namespace `qiji`。
- 作业模板：`build/soa/designs/job.yaml`（24 CPU / 64 GiB / 0 GPU，镜像 digest 已固定，
  `worker.sh` 依次执行 `commands/*.sh`，写入 `STOP` 后退出）。
- 依赖恢复：`/gpfs/users/linlizhi/flasheda-sta-online-20261102/{cmake-deps,sysdeps,flex-header}.tar.gz` 解压到 `/`。
  源码与构建缓存：`/gpfs/users/linlizhi/flasheda-density-cpu-20261106/fork-build-cache.tar.gz`
  （解压到 `/tmp/flasheda-gpl-timing/{OpenROAD,build}`，可增量构建）。
  整改版的构建方式见 `validation/soa-remediation/commands/`。
- 各阶段二进制都在 `/gpfs/users/linlizhi/flasheda-density-cpu-20261106/`：`baseline`（原版）、
  `soa-fft-openroad`、`soa-density-openroad`、`optimized-openroad`（adapter）、`soa-final-openroad`（`dd8f2fc`）等。
- 运行环境：`PATH=/usr/local/bin:$PATH`，
  `LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib`；
  `OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE`，**不要**设置 `OMP_PLACES`；`sta::set_thread_count 1`。

## 4. 待办（按优先级）

### P0：在大设计上验收整改版（只测量，不改代码）

诊断中的负收益全部出现在大设计上，而整改版只在 GCD、AES 和 10 个小案例上验收过。

- 用 `build/soa/designs/run_designs.py` 的 7 个配置（ibex_core td+rd、medium04、macro03、
  macro03-rd、large01、large01-rd、large02）跑原版与整改版，1/8 线程。
- 先做逐位差分：每个配置的所有运行只能有一个坐标哈希。
- 计时不加采样器，每种 3 次、交替顺序；需要归因时另跑一次带采样器的
  （`build/soa/diagnose/sampler.c` + `fold.py` + `analyze.py`）。
- 重点看诊断中大设计的四项回退是否消失：pin 同步（`GPin::updateDensityLocation`）、
  overflow 循环、线长 force 的 OpenMP 等待和串行清零、初始布局 getter。

### P1：线长 WA force 的 8 线程回退（整改版尚未处理）

位置：`src/gpl/src/wirelengthSoa.cpp` 的 `updateWireLengthForceWA_soa`，以及 `wirelengthSoa.h` 的
`updateNetForceColumns` 和 `pinGradientColumns`。

- 现在仍然先串行 `std::fill` 8 个 pin 列，再按 net 并行写 pin 列。large01/large02 8 线程下，
  OpenMP 等待从约 0.1–0.5 s 涨到 4–6 s，另有约 1 s 串行清零。
- 疑似伪共享（未经计数器证实）：同一 net 的 pin 在列中不连续，`has*ExpSum*` 是 1 字节列。
- 推荐做法见第 5.2 节第 2 条：删掉每轮清零（只在拓扑重建时清不属于任何 net 的 pin），
  并按 net→pin CSR 位置存放 per-pin 的 WA 中间量。
- 验收：`validation/soa-remediation/wirelength-replay.cpp` 对拍全部 WA 列，再做 P0 的大设计差分与计时。

### P2：初始布局的对象式访问

`Instance`/`Pin` 的 getter（`placerBase.cpp` 中的 `Instance::cx/cy/isInstance`、`Pin::updateLocation`）
和 `InitialPlace::updatePinInfo`/`createSparseMatrix` 在大设计 8 线程下多花约 4–5 s。
做法与 D1 相同：在拓扑边界建立行号映射，热循环直接读列。Eigen 求解部分不要动 ISA（见第 6 节）。

### P3：收尾整改版（等 Codex 会话结束后接手）

- 可选构建（`GPL_ENABLE_IPO`、`GPL_CPU_ARCH`，只对 8 个不依赖 Eigen 的内核生效）和
  D5 融合版的计时，在 `SOA_REMEDIATION.md` 中标为"生成中"，需要补齐结果。
- 整改通过 P0 后再提交（`git commit -s`）并推送到 `flasheda/gpl-soa`。
  Bazel 构建和测试至今没在 CI 上跑过，提交前要说明这一点或补跑。

### P4：会改变 QoR 的方向（需要另立 QoR 协议）

这些方向不再是逐位一致的优化，要用 HPWL、RUDY 拥塞、GRT/DRT 溢出、WNS/TNS 跨多个设计评估，
不能只看哈希：

- **routability 回退重跑**：AES 上 10 次 pass 的迭代占 55%；macro03/large01 开启 routability 后
  末轮迭代号变为 2.2–2.7 倍。可尝试自适应 inflation 步长、warm start 不回退、提前停止。
- **时序反馈**：AES 8 线程下占 42.9%，ibex_core 8 线程下占 35%。先测 `sta::set_thread_count N`
  （应保持确定性，但需核对）；`fullyRebuffer` 和 `repairDesign` 的削减会改变 QoR。

### P5：更真实的基准（需要用户决策）

ORFS（OpenROAD-flow-scripts）是 OpenROAD 官方的 RTL→GDS 流程脚本集：用 Yosys 综合，再调用 OpenROAD
完成 floorplan、布局、CTS、布线和签核，并自带 Nangate45/ASAP7/sky130 等平台上的一批参考设计和默认参数。
建议选 jpeg、ariane133（含 SRAM 宏单元）和一个 ASAP7 设计，用 ORFS 跑到 GPL 之前，
冻结成 ODB 作为固定输入。需要用户确定 ORFS 版本和镜像来源（现有容器没有 Yosys）。

### P6：GPU

fork 已有 GPU 接口，但没有编译和验收。按诊断中的上限估算，必须先把 P1–P4 推进到位，GPU 才有意义。

## 5. 对整改（未提交工作区）的审阅

审阅对象：2026-11-06 11:10 时工作区相对 `dd8f2fc` 的 diff（19 个文件，+603/−223）。只读审阅，未修改任何文件。

### 5.1 方向正确、实现可信的部分

- 整改范围与诊断一一对应，没有引入改变 QoR 的算法改动。
- 每一项都认真处理了逐位一致：
  - overflow：先按 x-major 写入连续数组，再按原 y-major 顺序做 `int64 += float` 归约；
  - 预条件：仍用物理尺寸 `ux-lx`、`uy-ly` 和 pin 数；
  - 坐标钳制：保留 float 与 int 的比较方式；
  - `getDistance`：先并行求平方，再串行按原 x/y 交错顺序求和；
  - 线长：寄存器累加的顺序与原版相同。
- 可分离重叠（D4）用每轴的 `nx+1`/`ny+1` 边界数组，`float(w) * float(h)` 与原版一致；
  y 方向超过 8 个 bin 的宏单元走无上限的后备路径。
- FFT 线程化在单线程或小网格时退回原库，并有 1/2/8 线程逐字节对拍的 `fft_test`。
- 构建缓存 mtime 和 Eigen ISA 两类问题都是"先停下、静态排查、再恢复"，符合故障修复原则。

### 5.2 需要调整的地方

1. **验收只覆盖了小设计，结论可能过于乐观（最重要）。** 整改只在 GCD、AES 和 10 个小案例上验证；
   AES 利用率 9.4%，正是 SoA 最占优的场景。诊断中的负收益全部出现在 large01/large02 上，
   整改版在那里的表现未知。在 P0 完成前，不要把"比原版少 21.6% / 31.1%"当作一般结论对外报告。

2. **线长 force 的 8 线程回退没有处理。** D5 只把 net 级累加移到寄存器，以下两点保持不变：
   每轮串行 `std::fill` 8 个 pin 列；按 net 并行写散布的 pin 列（疑似伪共享）。具体建议：
   - **删掉每轮清零。** 循环里每个属于某条 net 的 pin，8 个列都会被无条件覆盖
     （`has* = active`，`value = active ? exp : 0`）。需要清零的只有不属于任何 net 的 pin 和空洞行，
     可以在 `rebuildSoaTopology()` 时一次性处理。这一项逐位一致，改动很小。
   - **按 net 顺序存放 per-pin 的 WA 中间量**（DREAMPlace 的做法）：用 net→pin CSR 位置
     `i ∈ [offset_n, offset_n+1)` 作下标，再加一个 `pin → CSR 位置` 映射供 gather 读取。
     这样每个线程写连续区间，读取侧是只读访问，不会有伪共享，也不需要重排 pin 行。

3. **重复的状态源需要收敛。** 新增的 `kind`、`locked` 列分别是 `insts`/`Instance::isMacro()` 和
   `Instance::isLocked()` 的缓存，`kernel_rows_` 是 `nb_gcells_` 的派生索引，靠在约 10 处手工设置
   `kernel_rows_dirty_` 维持一致。`incremental02` 的锁状态 bug 已经证明这类缓存容易漏同步。建议：
   - 把 `nb_gcells_` 的增删和改绑收进少数几个函数，由这些函数统一置脏；
   - 让 `Instance::lock()/unlock()`（或 `PlacerBase::unlockAll`）成为唯一写入口，同步 GCell 列，
     而不是靠 `replace.cpp` 里的 `refreshCellLocks()` 补调用；
   - 加一个调试开关或测试钩子：每轮把 `kernelRows()` 与重新构建的结果、缓存的 `kind/locked`
     与 Instance 的实时状态做比对。覆盖 td-repair、routability 的 filler 删除与恢复、incremental 这几条路径。

4. **接口陷阱。** `CpuDensityGradientBackend::getCellGradients(gCells, out)` 现在忽略 `gCells`，
   直接按 `nb_gcells_` 计算，只有 `assert` 保护，而 Release 构建定义了 `NDEBUG`。目前唯一的调用方确实传入
   `nb_gcells_`，但以后若传入子集，结果会静默出错。建议改签名（不再接收列表），或者改成 Release 下也生效的检查。

5. **E（IPO/ISA）的投入产出要先算清楚。** 这部分已经引出 Eigen 崩溃和检查器误判两轮问题，
   而且只有 CMake 支持。仓库说明里 CMake 已计划在 2026-11-01 移除、改用 Bazel，Bazel 侧还没有对应配置。
   建议等计时结果出来再决定：收益在噪声范围内就从分支里去掉；保留的话，单独成一个提交，并补齐 Bazel 配置。
   `-ffp-contract=off` 对默认构建无影响，可以保留。

6. **初始布局没有处理（P2）。** 大设计 8 线程下多花约 4–5 s，整改范围里没有这一项。

7. **提交粒度。** 当前是一整块未提交的 diff，且有多个 agent 并发。P0 通过后建议拆成几个
   `git commit -s`：缓存列和行索引、坐标/pin/预条件列循环、overflow、D4 重叠、FFT 线程化、
   线长 D5、可选构建。每个提交都能单独做差分，便于定位以后出现的分歧。

### 5.3 已核对、没有发现问题的细节

- 中心更新用 `(dUx-dLx)/2` 的整数除法，pin 更新用 `(dLx+dUx)/2 + offset`，与原版
  `setDensityCenterLocation` 和 `GPin::updateDensityLocation` 一致。
- `locked` 的初值在 GCell 构造时从 Instance 读取；`replace.cpp` 在建立 PlacerBaseCommon 时调用
  `Instance::lock()`，早于 GCell 构造，所以首次读取是正确的。
- scatter 并行版对每个行列表各做一次 `omp for`，用私有整数网格合并，与原版的整数语义一致。

## 6. 验收协议

- **逐位一致**：与同线程数的原版比较坐标、朝向、master 哈希，以及 WNS/TNS、面积、HPWL。
  覆盖 GCD、AES、10 个扩展案例（`validation/soa-remediation/validate-cases.py`），以及 P0 的大设计。
  1 线程和 8 线程都要测，runner 遇到 `Using X of Y threads for placement` 降档必须拒绝。
- **性能**：每种 3 次，原版与候选交替顺序，报告中位数和范围；计时不加采样器；表格由脚本从运行日志生成。
- **按故障修复原则**：动态验证首次暴露一类错误时，先停下做静态排查，批量修复后再恢复。

## 7. 已踩过的坑

- `OMP_PROC_BIND=close` 加 `OMP_PLACES=cores` 会先绑住主线程，导致 OpenROAD 把 GPL 线程数静默降为 1。
- Northjob 容器里 `perf_event_open` 返回 EACCES，也没有 `addr2line`：符号化要在 bp-dev 主机上做。
  作业写出的目录属于 root，主机侧输出要写到新目录。
- 自写采样器在本次约 35 次带采样的运行中出现过一次启动时 SIGSEGV（`backtrace()` 不是异步信号安全的）。
  计时不要加采样器。
- 对整个 `gpl_lib` 加 `-march=icelake-server` 会改变 Eigen 的对齐和分配契约，在 `set_from_triplets` 中崩溃。
  ISA 参数只能作用于不依赖 Eigen 的内核（整改版已这样做，由 `check-kernel-isa.py` 检查）。
- 归档源码的 mtime 早于缓存 object 时，Make 不会重编译。叠加源码后要刷新 mtime（`verify-source.py --touch`），
  并核对二进制里有新增符号。
- `-march` 可能引入 FMA 收缩，导致结果不再逐位一致，必须加 `-ffp-contract=off`。
- macOS 默认忽略规则会把 `build` 匹配到 `BUILD`：新增 Bazel BUILD 文件要 `git add -f --sparse`。
- 浅克隆推送时 pack 很慢：先用 GitHub API 在远端建一个指向基线 SHA 的分支，再正常推送。

## 8. 文件索引

| 用途 | 路径 |
|---|---|
| 诊断报告 | `docs/flasheda/soa-speedup-diagnosis.md` |
| AES 采样工具与证据 | `build/soa/diagnose/`（`sampler.c`、`fold.py`、`analyze.py`、`tables.py`、`evidence.tar.gz`） |
| 大设计复测工具与证据 | `build/soa/designs/`（`run_designs.py`、`rerun.py`、`gpl_only.tcl`、`report.py`、`evidence.tar.gz`） |
| 整改验收脚本 | `docs/flasheda/validation/soa-remediation/` |
| 首轮验收脚本 | `docs/flasheda/validation/` |
| 远端原始数据 | `/gpfs/users/linlizhi/flasheda-soa-diagnose-20261106/`、`/gpfs/users/linlizhi/flasheda-soa-designs-20261106/` |
| 早期密度实验 | `~/Downloads/flasheda-openroad-rudy-20261105/integrations/openroad/density-cpu/` |
