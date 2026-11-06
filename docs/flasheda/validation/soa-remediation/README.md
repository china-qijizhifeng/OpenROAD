# 整改验收入口

结果和适用范围见 [SOA_REMEDIATION.md](../SOA_REMEDIATION.md)。

这些文件记录本次实际执行的 CPU 验收，不是可直接在任意机器运行的通用测试框架。
固定路径中的历史基线、源码缓存、设计输入和解析器均保留在 BP Dev 的 GPFS。
`job.yaml` 固定本次 worker 镜像和资源，项目为 WeaverRSI，队列为 proj-32。

- `validate-resident.py`：GCD/AES 的坐标哈希、WNS/TNS、面积和 HPWL 差分。
- `validate-cases.py`：十个 GPL Tcl 案例，比较全部 cell/port/net 输出；
  复用原版结果时重新核对原始日志哈希。Tcl 中其他断言继续执行。
- `benchmark-resident.py`：原版、整改前 SoA、整改后默认构建各三次，交替顺序。
- `benchmark-tuned.py`：默认与 GPL IPO/架构优化构建各三次，交替顺序。
- `benchmark-fused.py`：第一阶段默认构建与追加融合/寄存器累加的默认构建各三次。
- `density-replay.cpp`：直接包含生产密度头文件，验证三个原生捕获快照，
  与整改前内核比较 scatter/gather；每个变体十一次采样，每次重复二十遍。
- `verify-source.py`：核对 `source-sha256-final.json`，`--touch` 刷新叠加输入
  的 mtime，避免缓存 object 意外复用。远端副本名为 `verify-source-final.py`。
- `summarize.py`：检查计时与扩展差分完整性，汇总中位数、RSS 和 replay。
- `wirelength-replay.cpp`（追加融合阶段）：固定随机种子的 10,000 个 net，覆盖空 net、1–250 fanout、
  力截断和零系数；与保存的上一版生产 `previousWirelength.h` 对拍全部 WA 数值列。
  合成输入的内核耗时只用于筛选优化，不代替真实布局性能。
- `profile-fused.py`、`profile-tools/`：追加融合版（v2）1/8 线程各一次 1 kHz 主线程墙钟采样，
  沿用诊断中的采样器与符号化方式。融合后的密度与预条件应合并统计；
  采样时间不混入三次无采样的性能中位数。

`commands/01` 至 `07` 是首次构建、缓存失效排查和被保护条件阻止执行的历史命令；
**有效默认验收为命令 08–10**。命令 11 的整库 ISA 实验在命令 12 的首个案例中崩溃，
命令 13 被保护条件阻止执行。命令 14 静态审计 Eigen 契约，
命令 15 构建成功后的检查器误把独立测试目标当作生产库；15b 对现成构建复查并留存，
**第一阶段定向构建验收和汇总为命令 15b/16b/17b**。
命令 18 记录静态编译命令。**追加融合版本由命令 22/23 构建、验收，24b/25b 配对计时**。
命令 24 在启动前因输出目录冲突退出，25/26 被前置条件阻止，原回执保留。
26b 在 worker 采样，host 在独立目录符号化，不重跑采样。
**handoff 修正后的 v3 由命令 27–31 构建、正常差分、缓存调试检查、配对计时和汇总**。
`worker.sh` 对每条命令保存 `.log` 和 `.exit`；后续命令检查前置命令成功后才执行。
不要在已有结果目录原地重跑。新一次运行应使用独立任务路径、重新固定输入和二进制哈希。

`source-sha256-final.json` 对应默认二进制，`source-sha256-kernel-isa.json`
对应修正可选架构参数范围后的第一阶段工作区。两者仅 CMake 可选配置分支不同。
`source-sha256-fused.json` 对应追加融合与局部累加阶段。
`source-sha256-v3.json` 对应最终 CSR WA/缓存/接口修正后的工作区。
`check-kernel-isa.py` 检查实际编译命令和传递依赖，防止 ISA 特化内核引入 Eigen。
最终源码移除了无正向收益的可选 IPO/ISA 配置，检查器和旧配置仅作为历史实验保留。
新 `src/gpl/test/soa/wirelength_test.cpp` 以独立标量 oracle 覆盖 CSR WA 映射，
替代几何行与 WA 行一致的旧 replay 进行最终正确性验证。

17b/25b/31 生成阶段完成标志，让 worker 保持待命；最终验收和证据检查完成后才写 `STOP`，
等待 `worker.exit=0`，执行 `archive-evidence.py` 并记录作业终态后清理本次作业。

此处注册的 Bazel 测试可通过 `bazelisk test //src/gpl/test/soa/...` 运行。
本次真实 OpenROAD 验收使用缓存源码树上的 GCC/CMake 构建；不据此声称 Bazel 已执行。
本机和现有 BP Dev 主机均无 Bazel 工具，补充尝试在启动前失败；
`results/bazel-attempt.json` 与 `.log` 单独记录这一限制，未包含在清理前生成的证据包中。
