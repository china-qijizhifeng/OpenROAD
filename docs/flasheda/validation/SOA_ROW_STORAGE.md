# GPL 行存储、初始布局与 gather 优化

状态：实现与验收完成，最终为 snapshot10b（默认标量、显式 AVX2）。代码已提交为
`eb4c976b51`；44 个源码文件与已验收冻结内容相同。验证构建的输入为 fork master
`2847efec57` 加冻结 overlay，不能把 ELF 内版本标记误写成提交后的版本。
生产 GPL 基线与已验收的 `8bc1c3e674` 相同；旧版内存与 CPU profile 见
[SOA_NEXT_DIAGNOSIS.md](SOA_NEXT_DIAGNOSIS.md)。

本轮证据见 [manifest](soa-row-storage/manifest.json) 和
[E2B 总回执](soa-row-storage/e2b-final-compact.json)。最终 15 个 Linux target、
36 次功能差分、12 项 E2B 双模式门禁和 9 次沙箱选项探针均通过；正式性能、
局部 replay、分配量与未验证范围按版本分别报告。GPU 动态执行仍未验收。

## 本轮存储契约

- Instance/Pin/Net、GCell/GPin/GNet 和 Bin 的常驻对象数组改为 `SoaPool`。
  池仅持有列存储、活跃行号和必要时建立的行号到遍历位置映射。getter 对象是
  按值返回的借用视图，不再为每个实体常驻一个 `shared_ptr` 句柄。
- 跨对象引用是 `SoaRef`：表身份与稳定行号。`ref()` 显式提取身份；不允许把临时
  视图的 C++ 地址作为实体身份。列/对象池扩容和 swap-remove 不会改写稳定行号。
- 同类邻接只存 32 位行号，整张列表共享一个表身份。保留原邻接顺序，移除原先
  另建的 cell/net CSR 和 pin→net 副本；WA 重排在拓扑边界使用临时 permutation。
- pin→Instance/Net/GCell/GNet 关系列共享目标表身份，每行只存 32 位行号。
  并行写入之前必须在串行拓扑边界绑定目标表，避免首次绑定发生并发写。
- physical→GPL 的三张哈希表改为按 physical 行号索引的 32 位映射；OpenDB→physical
  保留以真实 DB 对象为 key 的查找表，但 value 只存行号。
- 已知实例、pin、net、bin、filler 批次先计算规模，统一预留全部字段列和释放表。
  未知增量回调统一扩列，避免各列各自倍增。删除行释放嵌套邻接分配，列容量可复用。
- 释放行表也采用 32 位 ID。filler 快照池提前绑定源列存储，批量预留 clone 行
  与全部坐标/梯度列，覆盖首次快照时尚无实体的情况。
- filler 快照仍是独立值：向快照池插入借用视图会克隆行，不能把同一行交给两个池
  所有。恢复也复制快照值；普通引用不拥有行。

32 位 ID 保留最大值表示空引用；超出支持规模时明确报错，不能截断成有效行。
当前 CPU 和 GPU host 适配接口均随行引用迁移，但 GPU 动态验收需单独标注。

## 初始布局

`updatePinInfo` 直接扫描实时坐标列，严格比较保留同坐标时首个命中规则。
矩阵构建对每个 net 提取一次 pin/instance 状态，在原嵌套顺序中产生 triplet 和
累计 RHS。B2B 只有触及每轴 extrema 的 pin 对能贡献矩阵条目；最多四个 extrema
预先按 pin 顺序列出。非 extrema pin 只访问此前的 extrema，省去零贡献 pin 对，
仍保留原 `(pinIdx1, pinIdx2)` 字典序。保留 movable IO 过滤、原 fanout 权重分母、锁语义、float/double 表达式
和 Eigen 求解顺序。

在 extrema 扫描时按 B2B 边数计算 triplet 容量上界：每轴至多两个 extrema，
每条边最多四个 triplet。这样不增加一次 O(fanout²) 预扫描，也不再固定预留
一百万 triplet 或在大设计里反复扩容。边中存在固定端、相同实例时实际条目较少。

## 线长 gather 候选

借鉴 [issue #3](https://github.com/china-qijizhifeng/OpenROAD/issues/3) 的
inspector–executor 分解：拓扑边界为有效 consumer pin 分配 cell 原 pin 顺序的
输出位置。每个 cell 保存起始位置，每个 pin 保存一个 32 位输出映射；它们是
梯度执行计划，不另复制拓扑邻接。

处理每个 net 的 force 后，在同一并行区复用当前 net/pin 数据，按原公式计算梯度，
逐项乘当前 net 权重后写入唯一输出位置。cell gather 连续读取两列梯度，仍按原
pin 顺序求和；不计算固定/IO 等没有 cell consumer 的 pin 梯度。
net 不变量和 float 分母乘积每个 net 提取一次；保留原公式里的 `1.0` double
提升与除法，不改成倒数乘法或 `1.0f`。

缓存新增有效 consumer pin 约 8 字节，映射每 pin 4 字节、每 cell 4 字节，
需以全流程 RSS 与 force+gather 总耗时验收，不能
单看 gather 阶段减少。调试 force 路径使缓存失效，bulk gather 回退原标量公式；
单 cell getter 继续读取实时列。三种权重 setter 更新表级 revision；bulk 调用发现
与缓存 revision 不同则回退即时求导，不复用旧加权结果。坐标同步也使缓存失效。
拓扑更新时重建 WA 与输出映射并使梯度缓存失效。
`soaMemory` 调试输出分别记录列容量字节、池索引、梯度缓冲和 physical 行映射；
这些是分配量核算，不能直接等同于 RSS，也不包含全部 DB、哈希表与分配器开销。

SELL-C-σ 跨 cell SIMD、pin 工作量均衡、块级缓冲属于独立 replay 对照，不能在未测
情况下全部叠加到默认路径。真实早/中/晚迭代 replay 应逐 pin 加权、逐 cell 比较
bitwise 结果，区分局部求导、写读、累加和映射成本。现有 cycles 采样不能将除法与
访存等待精确分账；不能把 issue 中的候选方案写成已测收益。

## issue #3 的经验与本轮取舍

[issue #3](https://github.com/china-qijizhifeng/OpenROAD/issues/3) 提出的科学计算经验，
本轮按可测的执行方案采用，不能把论文中的结论当成 GPL 已测结果。

| 经验 | 本轮采用方式 | 判据与边界 |
| --- | --- | --- |
| MFEM/libCEED 的算子分解与部分装配 | 分离拓扑映射、net 局部求导和 cell 汇总；默认候选使用完整 consumer 梯度缓冲 | 将 force 增量、scatter、gather、同步与内存合在一起验收；块级缓冲仍待独立比较 |
| inspector–executor | 拓扑变化时建立输出映射，每轮只更新数值 | 拓扑变化重建计划；坐标、系数和权重变化使数值缓存失效，不混淆两类生命周期 |
| SELL-C-σ 的跨 cell 推进 | replay 比较 C=8、σ=128 的可移植布局，一个 lane 对应一个 cell | 单 cell pin 顺序不变；无跨 lane 浮点归约；须检查汇编，不能把交错循环称为已实现 SIMD |
| 粒子模拟的邻接复用与分块 | 保留稳定行号和唯一邻接来源，重排只生成执行映射 | net 连接可能跨越很远的物理距离；没有复用证据前不增加图分块与空间排序 |

本轮优先解决“旧对象/指针结构与新列存并存”的空间浪费，再考虑新增缓冲。
cell 原 pin 顺序是数值契约；允许改变独立 pin 的求导次序，不允许改变 cell 内的
加法次序。AVX gather、parallel scan 或树归约不能直接替代这条约束。

真实布局的早/中/晚期 replay 对照 A（cell 即时求导）、B（net 求导并写缓冲，
cell 连续汇总）和 C（SELL 风格推进），独立 verifier 对照生产导出的逐 pin 加权梯度
及逐 cell 结果。replay 工具使用 Point AoS 和 64 位映射，生产使用 x/y 分列和
32 位输出映射，因此局部耗时/分配量只能比较实验路线，不能冒充生产 E2E 或内存核算。

## AVX2 独立 pin 求导与启用策略

启用方式：在进程启动前设置 `GPL_WIRELENGTH_SIMD=avx2`。snapshot10b 默认保留
标量；只有显式选择并且 CPU/OS 支持 AVX2 时进入四路内核，未知值或其他平台回退
标量。选择在每进程首次使用时缓存，运行中改环境变量不是切换接口。

该策略来自完整计时：snapshot10a 自动启用 AVX2 后，large01 单线程三组配对
均变慢（含反向运行顺序），中位数约增加 1.1%；不能根据局部 replay 更快就
默认上线。snapshot10b 只改变进程选项，intrinsics、公式、net 内阈值与尾处理
完全不变；10a 的正式计时是 AVX2 候选证据，不能冒充 10b ELF 的正式计时。

最终 snapshot10b 的 15 个 Linux target 与 ARM header/storage sanitizer 重验通过。
selector 在三个独立 Linux 进程中验证 unset=false、avx2=true、未知值=false；
ARM 三种均为 false。最终 ELF SHA 为
`2d4ea6220bcad8fc231478f04bff792dadb4657c3b271bce9404036e22e597d3`，
124,458,944 字节。完整源码与 10a 的差异只有选项 helper 和调用选择，见
[策略回执](soa-row-storage/snapshot10b-strategy-delta.json)；两种模式真实门禁均已通过。

用户要求在已有 net 批处理路线上补 SIMD。新增 `wirelengthSimd.h`，以四个 pin
对应四条 lane；先保持原 force 的 net 浮点归约，再向量化独立 pin 求导。
几何行号和最终输出仍可散布，WA 指数及 flag 按 net 原 pin 顺序连续加载。

- `coefficient * coordinate` 和 `coefficient * exp * weighted_sum` 仍用 float，
  原式含 `1.0` 的部分提升到四路 double，以 `divpd` 除法并在原方向梯度位置
  舍入回 float，最后用 float 相减、乘 net 权重。禁止 FMA/倒数近似。
- flag 与有效 consumer 共同决定 mask；inactive lane 用零输入和分母一，避免
  对无消费者或无效方向执行除零。少于四个及尾部 pin 走同一标量公式。
- `target("avx2")` 只约束该函数，运行时检测 CPU/OS 能力；其余平台走标量。
  不对整个 GPL 或 Eigen 添加 ISA 参数，不新增常驻数据列或缓冲。
- cell gather 的原 pin 加法顺序完全不变。随机拓扑、WA 非对齐起点、4pin 边界/
  尾部、无 consumer、inactive flag、独立零系数、正负零权重和大于 2^24 的
  整数坐标通过独立 native oracle 逐 pin 验证，并检查无无效/除零异常。
- 新 replay 的 D 调用真实生产 AVX2，E 调用真实生产 scalar；两者使用相同
  SoA/32 位映射，解码/构造在计时外。只有运行时 AVX2 可用且实际 consumer
  四 pin 批次数大于零，才记录为 SIMD 执行。独立 A/B/C 和捕获原值仍作对照。

首次 Linux 编译发现四路临时数组的 braced 写回可能匹配 `PointColumns::Reference`，
而 ARM 不编译该分支。已静态核查所有同类写回，显式构造 `FloatPoint` 后冻结
snapshot10a；不放宽 Werror，也不把 ARM fallback 通过作为 AVX2 通过的证据。
snapshot10a 的 15 个 Linux target（127 gtest + 6 standalone）通过，实际 AVX2
调用 12,230 次与 native oracle 逐位一致，所有 kernel 检查无 FE_INVALID/DIVBYZERO；
ARM 五项 header 检查与 storage ASan/UBSan 也通过。正式 ELF 为
`89f9d9777b01434c9d7c1a744148a16446bd3ea4490fa61fb1010ca21e6941f1`，
124,458,344 字节。最终 ThinLTO 产物中的 AVX2 helper 包含四条 `vdivpd`，
为四个 pin 的 double 除法；相关编译/汇编证据已归档。
[最终功能差分](soa-row-storage/snapshot10a-small-functional.json) 共 32/32 次通过，
覆盖八个案例的 1/4 线程、非空 filler 删除/恢复与双电压域；44 源码及输入前后不变。
六组真实早/中/晚期 D/E replay 与六组 GPL 门禁也全部精确通过；36 个有效正式计时已完成，见下表。

真实三份 replay 中，D/E 每配置各九次轮换并逐位通过，实际 consumer AVX
批次为每份 145,097 组，即 580,388 个 pin，约占该拓扑 1,400,119 个 pin 的
41.5%；819,731 个短 net/尾部 pin 继续标量。SIMD 覆盖率限制了可获得的增量，
但不是已经证明的性能瓶颈。跨 net 合批需给每 lane 不同 net 的统计量，不能
直接把当前广播公式用于这些尾部。下表是相同生产布局下的局部完整 force+scatter+gather
中位数，单位毫秒，E 为标量，D 为 AVX2；不包含 decode/setup 与其他 GPL 阶段。

| 迭代 | 线程 | E 标量 | D AVX2 |
| --- | ---: | ---: | ---: |
| 10 | 1 | 36.870 | 35.544 |
| 10 | 4 | 17.948 | 13.346 |
| 265 | 1 | 39.342 | 37.668 |
| 265 | 4 | 17.775 | 13.830 |
| 520 | 1 | 49.202 | 46.193 |
| 520 | 4 | 19.371 | 16.292 |

本轮中位数改善，但四线程三份快照的 min/max 全部重叠，不作为稳定 1.2–1.3 倍
收益的证明。独立 pin 求导/scatter 在单线程约快 8–13%，四线程差异更小；
完整 GPL 仍以正式计时为准。完整五路线、独立阶段、范围与指纹见
[SIMD replay 比较](soa-row-storage/snapshot10a-replay-simd-comparison.json)。

snapshot10a 的完整未插桩 GPL 计时，每设计/角色/线程三次交替重复，修正后一共
36 个有效样本。下表为秒；9 是标量候选，10a 自动启用 AVX2，**不是最终 10b
两种模式的正式计时**。两者输出及算法轨迹全部精确一致。

| 设计 | 线程 | 9 标量 | 10a AVX2 | 耗时变化 | 本轮范围重叠 |
| --- | ---: | ---: | ---: | ---: | --- |
| aes | 1 | 35.767 | 35.647 | -0.33% | 是 |
| aes | 4 | 26.461 | 26.242 | -0.83% | 是 |
| large01 | 1 | 89.719 | 90.695 | +1.09% | 否 |
| large01 | 4 | 54.901 | 54.805 | -0.18% | 是 |
| large02 | 1 | 100.299 | 96.229 | -4.06% | 是 |
| large02 | 4 | 57.048 | 56.014 | -1.81% | 是 |

large01 单线程三组配对都变慢，user CPU 时间也增加、sys 基本不变，不能用简单
运行先后解释；未记录周期性频率/steal/cache 数据，不能证明具体根因。新增 target
函数调用、mask 与短 net 尾处理是待证候选。large02 的中位数改善也不构成所有
设计加速的证据，因此最终 10b 选择显式启用。详见 [正式比较](soa-row-storage/snapshot10a-formal-comparison.json)
与 [现有日志归因](soa-row-storage/snapshot10a-order-diagnosis.md)。

一次 large02 正式运行期间，subagent 提前压缩了约 4 MB 证据包。按实际时间
交集唯一定位 `large02/t1-r2-snapshot10a`，保留原结果但排除计时，全部原序列
结束后只补同角色/线程一次；原污染和替补输出均精确。36 有效样本已使用该唯一
替补，没有删除其他慢样本；归档操作新增 idle gate。见 [排除回执](soa-row-storage/snapshot10a-timing-exclusion.json)。

## 另外三处 SIMD 的独立分析

按用户要求，三个 subagent 分别分析坐标更新、`getDistance`、梯度合并与预条件。
以下依据 snapshot9 与 snapshot10a 的实际最终 ELF 反汇编；`nesterovBase.cpp`
冻结版本 SHA 相同，但不声称两份 ELF 字节相同。这三项仅完成分析，尚未新增生产
实现或宣称性能收益。

| 位置 | 当前实际指令与结构 | 最小候选 |
| --- | --- | --- |
| `updateNextCoordinatesColumns` | 每轮一个 row/optimizer index，步长与动量为 `mulss/addss`，四次 `movss` 输出及 locked 分支；现有 `addps` 是同一 cell 的四个坐标裁边界，后接 mask 到分支阶梯 | 按 cell 四/八路推进；先 lower 再 upper 的 ordered compare/blend，最后选择 locked 原值；AVX2 需标量散射，或利用已证明连续的 optimizer-index run |
| `getDistance` 差值/平方 | 已有 SSE 四路 `movups/subps/mulps`，带 alias guard 和 chunk 大小门槛 | 如有收益再做 runtime AVX2 八路独立平方；最终求和严格保持原序 |
| `combineGradientsColumns<true>` | 融合 density gather 后仍是两个 `divss`，没有 packed merge；分子先写回再重读、penalty 重读 | 优先原 cell 内 x/y 两 lane，保留融合，去掉多余中间写读；再比较同循环局部多 cell tile |

坐标边界两次判断都针对**原输入值**；超大 cell 可以同时命中，upper 结果覆盖
lower，不能改为互斥分支或普通 clamp。geometry row 和 optimizer index 都可能
间接/离散；实际内层已无 accessor 调用，不把函数调用误判成当前障碍。

`getDistance` 的求和反汇编仍是一条 `addss` 依赖链，四点展开也依次为
`x0,y0,x1,y1,x2,y2,x3,y3`；无树归约或 FMA。把平方阶段从四路扩成八路，仅对
该阶段有理想 2 倍宽度上限。整个函数还包含约 32 字节/point 的数据流、串行
求和和 OpenMP 开销，须测阶段占比，不能据宽度推断整函数或 GPL 加速。

梯度合并必须保持整数尺寸先减、转 float 后乘面积；`size_t` pin 数到 float 的
舍入；分母 `max(raw,1)` 对 NaN 的原行为；以及 float 乘、加、除的独立舍入。
CPU 原融合路径拆成两遍至少新增每 cell 8 字节 density 读取和 16 字节 row/index
读取及一次调度/同步。只测 merge 算术不足以接受拆循环，应对比完整 fused 阶段。
密度 gather 的单 cell bin 累加序、最终 wire/density/gradSum 的串行序均保持。

原 AES 单线程 cycles 样本中，坐标 outline self 约 3.3%，`getDistance` 平方
outline 约 0.4%；只是已有预算参考，不是新候选的墙钟收益或完整上限。
后续候选应分别做真实状态 replay/逐位 oracle、特殊边界和完整 GPL 配对计时，
不要与线长 AVX2 同时修改，避免无法归因。实际最终产物的
[符号/指令摘要](soa-row-storage/snapshot10a-loop-disassembly.json) 及
[merge 指令摘录](soa-row-storage/snapshot10a-gradient-merge-instructions.json)
保留 ELF、符号与汇编 SHA；源码 action、完整 `.s` 和 vectorizer remarks 通过
[编译回执](soa-row-storage/snapshot10a-assembly-actions.json) 追溯。
三份独立分析分别保留为 [坐标更新](soa-row-storage/nesterov-coordinate-vectorization.md)、
[getDistance](soa-row-storage/getdistance-simd-analysis.json) 和
[梯度合并](soa-row-storage/gradient-merge-simd-analysis.md)。merge 的 pre-link
vectorizer remarks 记录 call/early exit 等阻碍；这是用真实编译 action 额外
`-S -fno-lto` 导出的诊断，最终 ThinLTO 的 scalar 指令结论另由正式 ELF 核实。


## 真实迭代 replay 与分配量

snapshot8a 的临时诊断 ELF 导出 large02 第 10、265、520 次 force/gather 状态。
三份输入各以 1/4 线程验证 A/B/C 路线，每路线九次轮换局部计时；独立 verifier
与实验程序均确认加权逐 pin 梯度、逐 cell 梯度，以及重建 force 的 net sums、
pin exp/flags 逐位一致。下表为局部完整 force+gather 中位数，单位毫秒：

| 迭代 | 线程 | A：即时求导 | B：net 求导＋连续汇总 | C：SELL 风格推进 |
| --- | ---: | ---: | ---: | ---: |
| 10 | 1 | 117.749 | 49.193 | 109.023 |
| 10 | 4 | 41.181 | 23.805 | 36.480 |
| 265 | 1 | 113.844 | 47.866 | 106.134 |
| 265 | 4 | 48.591 | 24.677 | 41.408 |
| 520 | 1 | 120.697 | 56.951 | 112.674 |
| 520 | 4 | 52.152 | 27.029 | 45.919 |

B 在六种配置中均胜过 A/C，支持本轮选用路线。C 比 A 有小幅改善，但编译器报告
含控制流/call 未向量化；pre-link 汇编展开为八次标量 `original()` 调用，
`addps` 仅配对同一 Point 的 x/y，不能宣称跨 cell 的 SIMD。A/C 使用同一 force
函数，四线程该阶段仍有明显散布，因此不从某一个阶段样本推导局部差异的原因。
这组实验的 AoS/64 位映射、完整缓冲和输入复用方式与生产路径存在差别；正式收益
仍由未插桩 GPL E2E 验收。映射构建、padding 和分配量在 replay 回执中单独保留。

同一次诊断运行启用 `set_debug_level GPL soaMemory 1`，已列类别为：
列及其中嵌套邻接/释放表 355,776,651 字节，三个 GPL pool 索引 9,728,356 字节，
consumer 梯度容量 11,198,288 字节，physical→GPL 映射 8,478,240 字节，
合计 385,181,535 字节。该调试口径未报告三个 physical pool 的行号索引，
也未覆盖 Bin、DB/哈希表、region/kernel 索引、optimizer 向量及分配器开销，
不能称为全部存储或 RSS。保留完整缓冲的空间取舍以正式峰值 RSS 为准。
带 snapshot I/O、memory debug 的这次运行不参与任何正式性能表。

replay 的两份源码、Bazel BUILD、实际编译/链接命令、工具链、向量化报告与汇编
已保存于 [复现包回执](soa-row-storage/replay-tools-reproduction.json) 指向的自含归档。

## 验收要求与当前证据

- 编译、测试、性能执行由 subagent 负责，保持源输入/补丁/ELF/工具链指纹。
- snapshot8a 的 15 个 Linux Bazel 测试 target 已全部通过（6 个独立 SoA，
  9 个 GPL 单元测试 target，包含 127 个 gtest case；另 6 个独立 target 为
  shell 包装的 standalone 检查，不能重复计作 gtest）；此前独立 C++20
  storage/metadata/density/point/wirelength 检查也通过。
  storage 行池生命周期和 compact relation copy/clear 另通过 ASan/UBSan。
- 首次 wirelength fixture 失败已定位为删除行主动清理后、fixture 仍在删除前快照
  全列的旧契约；活跃 pin 的 force/gradient/native oracle 均通过。fixture 调整为
  在删除后的真实拓扑边界快照，不能放宽浮点逐位判据。
- 正式 `//:openroad` 已通过 opt 构建，未插桩 ELF SHA-256 为
  `2f2c1364cc30d33cd76a9b7bf8d898d590b16e9bc1d928691fd5311bf7c15c03`，
  大小 124,452,824 字节。
- GCD、AES、incremental02/03、place_ios_td/rd、region01、simple01-td-repair，
  各 native/new × 1/4 线程，共 32 次真实运行通过，无失败或重跑。
  完整 cell/port/net 状态、master/orientation、面积、HPWL、适用案例的
  WNS/TNS/合法化，以及初始布局、逐回溯、gradient 和 routability 轨迹精确一致。
  输入和 42 个冻结生产/测试源码文件在运行前后哈希一致。
  region01 实际进入两个 voltage-domain solver；未观测边界另行标注，
  不用小设计通过代替大设计和性能验收。
- AES 的 1/4 线程分别实际删除 filler 11 次、恢复 1 次，native/new 的次数和
  完整 routability 轨迹一致。place_ios_rd 的相关调用没有实际删除或恢复 filler，
  因而该案例只覆盖 IO/routability，不另计非空快照恢复。
- extrema 配对跳过已由独立枚举对照覆盖 13,444 个输入、1,646,790 个有效配对，
  包括重复 extrema、同实例及 X/Y 贡献顺序，结果一致。
- 行视图迁移后按地址推算遍历索引的调用已静态排查；`initDensity1` 改为显式
  索引循环，不能从按值复制的 `GCellHandle` 地址计算 optimizer 索引。
- 初始布局另做临时插桩的逐迭代矩阵结构/values、RHS、坐标 guess 字节对照；正式
  性能 ELF 不带该插桩。
- 真实矩阵对拍 8/8 对通过：GCD/AES 1/4 线程、incremental02、place_ios_rd、
  region01 和 large02 单线程。每侧共 85 次矩阵构建、1,020 个组件（含 schema/
  endianness、X/Y 稀疏索引和值、RHS、location）逐字节相同，完整状态、QoR 与
  轨迹也一致。诊断后完整源码恢复，42 文件和正式 ELF 哈希复核相同；
  见 [矩阵差分摘要](soa-row-storage/matrix-differential.json)。
- E2B 每设计独立沙箱，native/current/new 三角色在同一沙箱串行，1/4 线程。
  正确性先行；正式计时与 perf 分开。关注峰值 RSS、初始布局、force+gather 总量，
  不跨宿主机汇合绝对计时。沙箱在下载证据后释放。
- AES/large01/large02 的 native/current/new × 1/4 线程，共 18 次 E2B 门禁
  全部通过完整状态、QoR 和选定算法轨迹差分，独立的 54 次轮换正式计时也全部通过。
  每设计每角色每线程三次重复，结果见下节；门禁单轮不混入正式中位数。

## snapshot8a 正式结果与后续修正

下表为同设计沙箱内三次未插桩运行的中位数。current 指已发布 SoA，new 为
snapshot8a；这不是追加句柄修正后的 snapshot9 结果。单位为秒与十进制 GB。

| 设计 | 线程 | current → new GPL | current → new 峰值 RSS |
| --- | ---: | ---: | ---: |
| AES | 1 | 36.030 → 36.774 | 0.577 → 0.553 |
| AES | 4 | 28.005 → 25.604 | 0.580 → 0.562 |
| large01 | 1 | 117.093 → 88.947 | 1.188 → 0.975 |
| large01 | 4 | 64.797 → 55.658 | 1.202 → 0.990 |
| large02 | 1 | 139.580 → 100.269 | 1.481 → 1.175 |
| large02 | 4 | 71.237 → 57.126 | 1.492 → 1.193 |

完整 native/current/new、三次 min/median/max 和 RSS 字节数见
[正式比较](soa-row-storage/snapshot8a-formal-comparison.json)。
AES 单线程耗时增加 2.07%，不能因大设计改善而忽略。

独立 AES phase 单次观察发现，单线程 force+gather 合计 4.484 → 4.318 秒，
四线程 3.722 → 1.890 秒；单线程整体回退出现在其余 Nesterov 工作中。
phase 插桩范围与开销存在差别，不能用单次 phase 推导精确整体百分比。
随后正式 ELF 的两次 99Hz cycles 采样均输出精确、零丢样：current 3,962、
new 4,057 samples。new 的 GCell SoaRow 析构 self 为 2.28%，callgraph 中
1.42% 经 updateNextIter、0.79% 经缓存 gather 的 handle 转换。
cycle period 是近似 CPU 事件权重，不能与 phase 墙钟秒数直接相减。

大设计的独立 phase 观察进一步支持按 net 求导路线。下表只表示 snapshot8a
单次插桩归因，不能代替前表的三次正式计时，也不用于归因 snapshot9：

| 设计 | 线程 | current → new 初始布局（秒） | current → new force+gather（秒） |
| --- | ---: | ---: | ---: |
| large01 | 1 | 29.992 → 27.444 | 49.297 → 23.145 |
| large01 | 4 | 25.889 → 26.048 | 21.183 → 10.585 |
| large02 | 1 | 31.353 → 29.999 | 74.009 → 27.646 |
| large02 | 4 | 31.779 → 29.088 | 25.252 → 10.333 |

large01/02 单线程矩阵构建从 4.822/6.466 秒降到 3.641/4.796 秒；
Eigen solve 仍约 18–19 秒，是初始布局剩余的主要部分。large01 四线程的 solve
在这一次观察中反而增加，抵消了矩阵构建收益，因此不能宣称初始布局在所有配置都加速。
全部 18 次 phase 输出精确一致，OMPT 实际线程数为请求的 1/4，invalid 为零；
OMPT 的累计线程等待时间包含 worker 空闲，不能从墙钟预算中直接相减。
详见 [phase 比较](soa-row-storage/snapshot8a-phase-comparison.json)。

据此追加最小 snapshot9 修正：`SoaPool::ref(index)` 直接返回表身份和稳定行号，
NBC/filler 的按索引 getter 以及 GCellHandle 箭头/转换不再先构造临时 Object
再取 ref；`operator*` 仍返回借用视图。未更改算术、列布局或拓扑计划。
最终速度与功能以 snapshot9 复测为准；不凭调用图直接宣称已消除回退。
snapshot9 的 15 个 Linux target（127 gtest + 6 standalone）以及独立 header/
storage ASan/UBSan 已通过，正式 ELF SHA-256 为
`d7ddd6f2280b9842076ce6f9c38edf40e6351945d1515c8d5621fc8968392e7f`，
大小 124,450,848 字节，源码/实际编译 flags/下载核验回执均已归档。
snapshot9 的 8 案例 × 1/4 线程 × native/new 共 32 次功能运行也全部精确通过，
包含 AES 每线程 11 次非空 filler 删除、1 次恢复，以及 region01 两个电压域。
42 文件及输入验收前后不变；见 [最终功能回执](soa-row-storage/snapshot9-small-functional.json)。
这批运行开启 debug 和状态导出，仅作正确性证据，不参与正式计时。


snapshot9 的独立复测没有使用跨轮基线替换：AES 单线程 current/8a/9 各三次，
中位数 36.579/36.517/36.497 秒，范围重叠，不能归因句柄修正消除了先前 2.07%
回退；四线程为 27.527/25.880/26.150 秒，9 相对 current 范围分离，但 8a/9
之间无可分辨改善。large01 四线程首对曾为 48.453/53.768 秒，补两组后
8a/9 三次中位数 55.784/53.768 秒、范围重叠，没有稳定重现该首对回退。
large01 单线程和 large02 的 1/4 线程只做单对筛查，不能据此宣称稳定加速。
详情见 [snapshot9 复测](soa-row-storage/snapshot9-comparison.json)。


首次完整生产编译以 `keep_going` 收集到 21 个错误，涉及 3 个 TU：
测试仍以 `auto*` 接收行引用，Instance/GCell 仍将 `this` 传给行引用 API，
两个只读 helper 不能绑定临时借用视图，以及对 `GCellHandle` 直接判空。
停止该轮后静态检查全部同类调用，批量改为值绑定、显式 `ref()`、const 引用和
显式行引用，再执行 snapshot8a 的完整 15-target 验证；未关闭诊断或降低 Werror。

旧 native/current ELF 与新构建的 `.bazelrc`、`.bazelversion`、`MODULE.bazel`
和 lock 文件逐字节一致，旧 ELF 的 compiler/linker 标识均为 clang/LLD 23.1.2，
新 hermetic clang 的实际版本为 23.1.2，SHA-256 为
`2667f50297ef526c2a594fae811614d83e55a14ac2ac3ac4414b100b7d0d9642`。
两轮均使用 `--config=opt` 的 O3/ThinLTO 和原浮点契约。

## 最终交付与资源

snapshot10b 的 [36 次功能差分](soa-row-storage/snapshot10b-small-functional.json)
全部精确通过：默认八案例 × native/new × 1/4 线程共 32 次，显式 AVX2 的 AES
另四次。default/explicit 的 native/new AES 在每线程下也跨模式逐项相同，
各覆盖 11 次非空 filler 删除与一次恢复；44 源码及输入前后不变。
[12 项 E2B 门禁](soa-row-storage/snapshot10b-e2b-gates.json) 覆盖 AES、large01、
large02 的 default/explicit × 1/4 线程，均与 snapshot9 的完整状态、质量与
算法轨迹精确相同。每沙箱各三次独立进程探针验证实际 AVX2 能力和选项行为。
这些门禁带状态校验，只证明正确性，不是 snapshot10b 的正式计时。

16 个原始 E2B 证据包在本机与 HK 核验后，三个本任务沙箱已释放，模板保留；
见 [清理回执](soa-row-storage/e2b-cleanup.json)。源码、ELF、测试、真实快照、
编译/链接 action、汇编与原始日志保留在回执指向的本机/HK/GPFS 路径。

本轮 Northjob `flasheda-soa-row-storage-20261106` 编译 worker 也已按精确身份删除，
临时容器及缓存随之释放；删除后回查 absent，并再次核对 GPFS 的正式构建证据与
36 次功能归档 SHA，持久目录仍保留。见 [worker 清理回执](soa-row-storage/northjob-cleanup.json)。
