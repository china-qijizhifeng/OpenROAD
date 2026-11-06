# GPL SoA 合入 fork master

原验收快照为 `fc730887cd`，基于 `18216f3506`；完整原始结果见
[诊断整改报告](SOA_REMEDIATION.md)。该快照已单独提交，历史 manifest 和证据不修改。

整合提交 `8bc1c3e674` 将 fork master `a6e7f507a1` 合入 `flasheda/gpl-soa`。
GitHub 入口为 [PR #2](https://github.com/china-qijizhifeng/OpenROAD/pull/2)，目标是
`china-qijizhifeng/OpenROAD:master`，不是官方上游仓库。

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
`a6e7f507a1`，候选为 `8bc1c3e674`。当前正在准备 Bazel 构建与真实布局差分。
结果完成后在本节补齐源码、二进制、输入和运行回执。此版本未做 GPU 执行验收。
