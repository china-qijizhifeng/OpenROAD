# GPL 全字段 SoA 验收结果

实现提交：`3bc0628`；Bazel 测试注册：`1910dc2`。基线：
`18216f3506f090aa7e12f83620b47f21ddff4f40` 的既有 native 二进制。

## AES 全局布局计时

Nangate45 `aes_cipher_top`，density 0.5、seed 42、timing-driven +
routability-driven，随后详细布局、合法性检查和 QoR/坐标核对。
每种配置 3 次；原版/候选交错顺序；STA 线程固定为 1。
`OMP_PROC_BIND=FALSE`、未设置 `OMP_PLACES`、`OMP_DYNAMIC=FALSE`。
CPU：Xeon Gold 6330；Northjob 分配 24 CPU / 64 GiB / 0 GPU。

| GPL 线程 | 原版中位数/s | SoA 中位数/s | 时间减少 | 原版范围/s | SoA 范围/s |
|---:|---:|---:|---:|---:|---:|
| 1 | 79.325 | 75.028 | 5.42% | 79.004–80.706 | 74.934–75.434 |
| 8 | 44.702 | 41.336 | 7.53% | 44.491–46.571 | 40.956–42.220 |

所有 12 次运行的坐标、朝向、master 哈希及 WNS/TNS、面积、HPWL 完全一致。
AES 坐标哈希 `bc141af641d8db860dde7e3cee05047e037fb9b43171f0e8ee4b49d2ee903f92`，
WNS 3.881395315 ns、TNS 0、面积 23367.834 µm²、HPWL 440574.9175 µm。
基线和候选都通过详细布局合法性检查。此表只计 global_placement，
不把综合、详细布线或整个 RTL→GDS 时间纳入加速比。

## 峰值内存

单独补测 8 线程配置，Python 子进程使用 `getrusage(RUSAGE_CHILDREN)`
读取 OpenROAD 进程生命周期的最大 RSS；这两次不纳入上述三样本计时。

- 原版：508.2 MiB；用户 CPU 226.8 s，系统 CPU 0.6 s。

- SoA：541.2 MiB；用户 CPU 199.9 s，系统 CPU 0.9 s。

峰值 RSS 增加约 33.0 MiB（6.5%）；不能把本次改造描述为同时省时省内存。

## 差分覆盖

GCD/AES 之外，`macro01`、`macro02`、`region01`、`cluster_place01`、
`simple01-td-repair`、`place_ios_td`、`place_ios_rd`、`incremental01`
均在 1 线程和真正的 8 线程下，与相同线程数的 native 二进制逐条比较
单元位置、朝向、master、placement status、端口形状/层和 net 连接，全部一致。
保留原 Tcl 的功能断言，替换 golden DEF diff 为同环境差分。

首次 requested-8 验收因 `OMP_PROC_BIND=close` 提前绑定主线程，被 OpenROAD
物理核检测降成 1 线程。这批仅作为串行证据；修正配置后完整重跑 8 线程，
新 runner 检查并拒绝任何线程降档。脚本队列 19/20 的依赖顺序失败没有运行
计时/归档；有效计时与归档入口为 21/22。

四个独立 C++ 测试通过本地 Clang ASan/UBSan 和 worker GCC 执行：
row 生命周期、metadata 列复制/复用、PointColumns 回滚/交换、非方形舍入网格
密度 all-bin oracle（含宏单元缩放）。Bazel 目标已注册，当前环境实际使用
固定依赖的 CMake Release 构建与直接 C++ 测试，没有声称跑过 Bazel CI。

## 范围与限制

七类核心实体的全部业务字段（包括标志、关系和关联列表）都驻留列存；
优化器和回滚状态按 x/y 分列。singleton 配置、OpenDB 自身和 IO 约束几何等
辅助结构不属于这次七类实体迁移。详情见 [设计说明](../gpl-soa.md)。
GPU 接口适配了 PointColumns，但未编译/运行 GPU 验收，不报告 GPU 加速。

DREAMPlace 用作连续数组、CSR 拓扑及算子划分的架构参考；保留 OpenROAD 的
逐贡献整数截断、fastExp 和浮点累加顺序，没有移植其不同数值语义。

完整证据、二进制及可恢复构建缓存保留在 GPFS
`/gpfs/users/linlizhi/flasheda-density-cpu-20261106/`。
`source-manifest.json` / `benchmark.json` 记录源码和二进制指纹。
原始失败命令和被纠正的配置结果也保留在归档，避免只展示成功样本。
