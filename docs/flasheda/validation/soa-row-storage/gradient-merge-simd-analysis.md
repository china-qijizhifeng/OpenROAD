# 梯度合并与预条件的 SIMD 边界

本次只读生产代码，未编译、未新增 GPL 采样。已读取 watch_build 导出的 snapshot9 正式 ELF 汇编并核验文件 SHA；确认 merge 当前仍为标量。尚无本候选的计时结果，以下没有把“可向量化”写成“已经有收益”。

## 当前路径与必须保持的数值契约

- `src/gpl/src/nesterovBase.cpp:3740` 的 `combineGradientsColumns<GatherDensity>` 在 CPU-only 路径融合 density gather、merge 和 precondition。`3873` 调用 `<true>`；GPU-enabled build 的 CPU fallback 在 `3876` 先走 backend，再调用 `<false>`。优化不能假定所有 CPU 执行都实例化 `<true>`。
- `3767` 遍历存储分组内的 cell。`3768` 的 `r=rows[j]` 是 geometry row，`i=indices[j]` 是 optimizer 输出；两者均不能假定等于 j 或连续。`2924` 的 kernelRows 缓存仅按 Store 分组，`2940` 仍保存任意 row；indices 来自 region 顺序、在每组内递增，却可能有间隔。
- `3775` 面积是 **int 差先完成，再分别转 float，然后 float 乘法**：`float(ux[r]-lx[r]) * float(uy[r]-ly[r])`。不能替换为 float 坐标相减或 int 面积乘积后统一转 float。
- `3777` pin_count 来自 `SoaRows::size()` 的 size_t→float。不能未经范围验证用 signed int32 SIMD conversion 代替；超过 2^24 时也有既定 float 舍入。
- `3778` 分母严格为 `max(pin_count + densityPenalty * area, 1.0f)`，每个操作保持 float 结果；`minPreconditioner` 定义于 `nesterovBase.h:899`。
- `3781` 分子是先 float 乘 `penalty*density`，再 float 加 wire；`3783` 是 float 除法。不能用 FMA 合并、提前倒数后乘、fast reciprocal/Newton 或重排乘加。x/y 可共享分母，但必须分别正确舍入除法。
- `std::max(raw, 1.0f)` 在 raw 为 NaN 时保留 raw；直接 `_mm_max_ps(raw, one)` 会取第二操作数，改变行为。可用 ordered compare+blend，或在已证明有限值的契约下采用正确操作数次序。不能悄然吞 NaN，使后续 invalid-gradient 检测失效。
- `3881` 的全局 reduction 保持串行且不变：wire/density 是 x 后 y 逐项累加；gradSum 是先每 cell 求 `abs(x)+abs(y)`，再加入总和。SIMD merge 不要求改动这些依赖链。
- density gather 的 bin 求和在 `densitySoa.h:65`，当前是 x-outer/y-inner 两条独立的 float 累加链。只把不同 cell 或 x/y 分量并排计算，不改变单 cell 的 bin 求和次序。bin 内横向 reduction、跨 bin 拆分部分和不在本候选范围。

## 最小候选与优先顺序

1. **先做原 cell 内的 x/y SLP，不拆循环。** 正式 snapshot9 汇编确认未做此 SLP：两个 `divss`，零 `divps`；分子计算也全为 `mulss/addss`。可尝试在原 cell gather 后将 `{dx,dy}` 与 `{wx,wy}` 装入两 lane，按原 mul/add/div 序列计算再散射。保持现有 OpenMP、一次 cell 遍历、无新增常驻数组。set/extract 指令成本可能抵消少量算术收益，需实际测量。先把 penalty 显式缓存为只读局部常量、分子留在局部寄存器、最终只存一次结果，也是更小的候选；汇编中现在有中间结果 store/reload 与 penalty 重复读取。若要通过预加载输入或 restrict 解除 alias 约束，必须以实际三个调用点的不同 PointColumns 对象为契约，不能给通用 API 增加未验证的非别名假设。
2. **四或八 cell 的局部 tile，在同一个 fused loop 内计算。** 每个 tile 的 density gather 继续逐 cell 执行原算法，然后将独立 cell 的 area、pin_count、dx/dy、wire 装入 lane；以 divps 计算 merge/precondition，最后按 indices 写回。WA SIMD 中的 runtime dispatch 可复用思路，但不要强制全工程 ISA 或 fast-math。尾部走标量。geometry/indices 的随机 gather 与散射保留，不假设 AVX2 支持 scatter。可检查连续 indices 的 run，仅对真实连续 run 使用连续 load/store。
3. **不建议为 SIMD 直接恢复全局 gather/merge 两遍。** `<false>` 已经具备独立 merge pass，可单独评估；CPU-only `<true>` 当前受益于融合。把 gather 的输出写 density 后再全局读回至少增加每 cell 8 B 的 density 读取、16 B 的 row/index map 再读取，以及一次 OpenMP 调度/同步。缓存可能吸收部分流量，但成本不能算作零。若改成全 NB 顺序 merge，还需要 materialize area/pin_count 映射，会引入新的数组及维护职责。

tile 也并非免费：局部 scratch 若未被编译器保留在寄存器，会增加 density/metadata 的 stack 写读，且跨 cell 的 gather 不规则控制流仍为标量。它的优势仅是保持一次全局 cell 遍历和已有数据的缓存热度，需汇编证明没有严重 spill。

## 针对性验证与净收益判据

- 独立 native oracle，比较 area、raw denominator、clamped denominator、分子、输出 x/y 的逐 bit 结果；从真实 checkpoint replay 捕获 geometry/indices/梯度/penalty，不把生产 helper 再调用一次充当 oracle。
- 覆盖空/1/3/4/5/7/8/9 cell、乱序 geometry、稀疏 indices、多个 Store、tail、filler 的零 pin_count；dimensions 与 pin_count 在 2^24 附近、不同精度量级、阈值 raw 恰好等于/略低/略高 1。
- signed zero 的 wire/density/penalty、正负梯度、有限大/小值与 cancellation；单独校验 NaN 分母行为、无效 lane 的异常标志。非法几何造成原始 int subtraction overflow 不应作为合法 oracle 输入。
- 不改 serial reduction；仍用现有 small32 的完整 state/quality/逐 backtrack trace 差分证明路径一致。
- 性能先测 **完整 fused gather+merge** replay；merge-only kernel 可解释算术收益，但不能单独批准拆循环。使用同一环境串行交替 A/B、1T/4T、小/大与稀疏拓扑；正式计时关闭 debug/export/profile。
- 接受条件是完整 fused phase 的配对耗时收益超过自身噪声，且其他配置无显著回退；总 updateGradients 与完整布局也必须改善或持平。只见 divps 或 merge-only 更快而 fused phase 更慢，应拒绝。

## 实测状态

实际证据：`build/soa/row-storage-20261106/snapshot9-loop-disassembly-v2/receipt.json`，目标为 `combineGradientsColumns<true> [clone .omp_outlined]`，汇编 `0.s` SHA `d7ca028a66ceff289ab49241c33ab6a4b9ceefd4efc3868d2ae4a6e72c9ab059`。对应无 instrumentation 的最终 ThinLTO ELF `new-snapshot9` SHA `d7ddd6f2280b9842076ce6f9c38edf40e6351945d1515c8d5621fc8968392e7f`。

- `0x203b200` 等为 pin_count 的 size_t→float 转换；`0x203b3af` 起含处理 unsigned 高位的分支，不能直接缩成 signed int32 conversion。
- `0x203b233/237` 分别 cvtsi2ss 宽、高，`0x203b23b` mulss 算面积。
- `0x203b243/24b/24f` 读取 penalty、mulss、addss 构成 raw denominator；`0x203b289` 是 `maxss raw,one` 的 AT&T 操作数次序，结果按原 NaN-first 语义保留 raw。
- `0x203b25b/268` 与 `0x203b283/297` 分别是 x/y 标量 mulss/addss。penalty 在 `0x203b243` 和 `0x203b27b` 两次读取。
- `0x203b275` 与 `0x203b2a4` 先把未除的分子写入 sum 数组，随后 `0x203b2aa/2ba` 重新读取，`0x203b2b0/2c0` 分别 divss，再写最终值。不是 packed division，也没有 FMA/reciprocal。
- `0x203b36a` 是每 cell 的 forEachOverlap 调用；density 暂存在栈上，经 `0x203b380/38d` 写出。非内联调用和别名限制是当前编译器难以自动做跨 cell SIMD 的可信障碍，但未做编译器 optimization remarks 验证，不能当作完整归因。

snapshot9 与已冻结 snapshot10a 的 `src/gpl/src/nesterovBase.cpp` SHA 均为 `ba6fe6f877e7013badfbeb1c85a730a1d93b9b5ab9f42ce69dc8c2fa8a4149d5`。这是源码相同的证据；ThinLTO 全程序布局可能变化，不能据此声称两个 ELF 的指令字节相同。尚未编译任何候选、没有候选净收益数据。snapshot10a 的构建与 small32 验收优先，不额外启动 fullGPL profile。

已追加实际 snapshot10a ELF 的独立佐证：`build/soa/row-storage-20261106/snapshot10a-loop-disassembly/receipt.json` 与 `0.s`。ELF SHA `89f9d9777b01434c9d7c1a744148a16446bd3ea4490fa61fb1010ca21e6941f1`；merge assembly SHA `e9ffd1aa95d54642e7d79bda1efb1036614bb528f52f66949a35c739527cb343`。该 clone 仍有 2 divss / 4 mulss / 4 addss / 1 maxss，零 packed 算术：具体除法在 `0x203b720` 和 `0x203b730`。本目录 `gradient-merge-instruction-digest.json` 收录实际 symbol、指纹和精选指令，可直接引用；没有重跑编译或采样。

watch_build 另提供按实际 production CppCompile 参数执行的独立 **pre-link** diagnostic assembly/remarks；诊断额外使用 `-S -fno-lto -gline-tables-only -Rpass*`，不能冒充最终 ThinLTO 的 optimization remarks。`snapshot10a/nesterovBase-vectorization.log` SHA `7b5c537ace049d70591e927541d06f2542e82d19c7e621f5740f51990abae6c2` 中，对 `nesterovBase.cpp:3766` 明确报告 `call instruction cannot be vectorized`、`Cannot vectorize early exit loop`、`instruction cannot be vectorized`；`3748` 的外层存储分组 loop 报 missed。`densitySoa.h:68` 邻近报告 unsupported terminator。这将“调用/控制流阻碍跨 cell 自动向量化”从源码推断提升为实际诊断编译器报告；alias 是否也是独立原因仍未得到该 report 直接确认。最终 ELF 的两条 divss 则由独立 disassembly 确认，二者证据层次保持区分。提取记录在 `gradient-merge-compiler-remarks.json`。
