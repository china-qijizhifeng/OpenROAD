# GPL SoA 合入 fork master

原验收快照为 `fc730887cd`，基于 `18216f3506`；完整原始结果见
[诊断整改报告](SOA_REMEDIATION.md)。该快照已单独提交，历史 manifest 和证据不修改。

整合提交 `8bc1c3e674` 将 fork master `a6e7f507a1` 合入 `flasheda/gpl-soa`。
GitHub 入口为 [PR #2](https://github.com/china-qijizhifeng/OpenROAD/pull/2)，目标是
`china-qijizhifeng/OpenROAD:master`，不是官方上游仓库。
PR 已于 2026-11-06 20:49 CST 合并，合并提交为
`1ba3add7889ce159bdd07a89b767aacde9ba71c7`。合并前完成下述功能验收；
独立 profile 的结果不作为已完成该合并的前提或性能证明。

## 整合处理

- 遵循 master 删除 CMake，使用 Bazel 构建。
- overflow 和 density field 的重叠循环保留 SoA 实现。
- 保留 master 对矩形 FFT 网格的 `bin_size_x / bin_size_y` 修正。
- 保留 master 的时序修复 net 过滤、回调 pin/net 检查和 mirrored IO 约束修复。
- OpenSTA 子模块按 master 的 gitlink 恢复。

master 包含影响计算结果和接口的更新，因此旧版的 84 次运行和性能数字只属于
`fc730887cd` 的验收，不能替代本次整合的验证。

## 整合版本验证

由构建/测评 subagent 在新隔离环境执行；对照为同工具链构建的
`a6e7f507a1`，候选生产源码为 `8bc1c3e674`，测试另包含 `268276e3d1` 的修正。
两侧统一显式使用 `--config=opt`，启用仓库配置的 O3 / ThinLTO。
这与仅指定编译模式 `-c opt` 不同，也与历史 GCC/CMake 工具链不同，
新批耗时不能直接拼接到旧批性能表。
两侧完整 Bazel 构建均成功；六个独立 SoA 测试和三个现有 CPU 数学测试
（FFT、overflow、wirelength gradient）全部通过。

GCD、AES、incremental02/03、place_ios_td/rd、large01、large02，
各原版/候选 × 1/8 线程共 32 次真实运行，完整 cell/port/net 状态精确一致；
GCD/AES 的 WNS/TNS、面积、HPWL 和合法化同时核对。另八次候选缓存一致性检查通过。
large01/02 为 LEF/DEF 布局验证，比较面积与 HPWL，不声称时序或合法化覆盖。
本次仅单遍功能验收，没有补做三次重复性能计时。

未插桩二进制 SHA-256：

- master：`9eedbe284a9330f4269c1c1c557ae869f2d9c8472ae0610b3f7dc0f98b6b3fc4`
- 候选：`bf95389f331b53e91cc3942d3c32a38aa2dcdfcaa57a6d520e43ca7e944a2edd`

紧凑回执与源码/输入/运行来源见 [整合验证目录](soa-master-integration/)。
独立 profile 在通过上述门禁后执行。按用户要求转至 HK SWEEnvMaker 0820 的
三个 E2B 沙箱并行，实际 4 物理核，采用有效的 1/4 线程；12 次 phase/OMPT 和
8 次独立 perf 均精确一致。四份 record 零丢样，事件计数和资源指纹完整保留。
结论见 [剩余瓶颈分析](SOA_NEXT_DIAGNOSIS.md)，
证据见 [E2B profile](soa-e2b-profile/README.md)。
该集成版本未做 GPU 执行验收。

Northjob 的部分 profile 独立归档，没有与 E2B 绝对耗时混合。
本次创建的两个 Northjob 作业（含被替换的准备作业）和四个 E2B 沙箱均已清理；
回执位于本目录的 `soa-master-integration/` 及 E2B 的 `experiment.json`。

## 首次构建暴露的问题

六个独立 SoA target 首次在新 Clang / `-Werror` 下构建时，
`storage_test.cpp` 的直接自赋值表达式触发 `-Wself-assign-overloaded`。
停止后续测试后静态排查同类写法，只发现这一处；提交 `268276e3d1` 改为通过
const 引用别名执行同一个自赋值路径，没有关闭警告或修改业务源码。
恢复后六个 target 全部通过。其生产源码仍对应整合提交 `8bc1c3e674`。

当前 master 原版已通过完整 Bazel 构建。候选第一次构建在公共 CAS 中间件下载上
长时间停止进展；实际 worker 的缓存临时文件大小/时间戳不变，没有编译器子进程，
Bazel 日志确认 `DownloadTimeoutException`。保留日志后仅中断本任务的候选 client，
复用已有 server、分析结果和本地 disk cache，用 `--remote_cache=` 恢复候选构建。
源码、工具链、Werror、O3 和 ThinLTO 均未改变。中断前的 2269 秒是构建流程墙钟，
包含网络等待，不能当作编译 CPU 时间，也不能归因于 GPL 运行时性能。
