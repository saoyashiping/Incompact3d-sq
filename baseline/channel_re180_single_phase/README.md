# 单相槽道 Re_tau≈180 基线冻结（2.0）

## A. 目的
本目录用于冻结当前**已验证通过**的单相槽道基线，作为后续刚性/柔性纤维开发前的唯一母体版本。

- 本目录仅包含单相槽道基线归档。
- 本目录**不包含任何纤维功能**、纤维接口或纤维相关参数。
- 本目录的目标是保证：可复现、可回退、可审查、可维护。

## 基线输入来源与选型
冻结输入文件为：`baseline/channel_re180_single_phase/input.i3d`。

其来源是对以下现有文件的**逐字复制**：
- `examples/Channel/input_DNS_Re180_LR_explicittime.i3d`

选用理由（在仓库存在多个 channel 输入文件时）：
1. 该文件文件名明确标注 `DNS_Re180`，与“Re_tau≈180 单相槽道基线”目标直接一致。
2. 该文件是长算设置（`ilast=100000`），更符合“已验证长算基线冻结”场景。
3. `tests/Channel/input_test_*.i3d` 主要用于短步测试（如 `ilast=100`），不适合作为生产级长算冻结基线。
4. `examples/Channel/input.i3d` 与此文件当前内容一致，但 `input_DNS_Re180_LR_explicittime.i3d` 的命名语义更清晰、可审查性更好。

## B. 物理范围
本基线严格限定如下物理范围：

- 单相（single phase）
- 不可压（incompressible）
- 槽道流（channel flow）
- DNS 风格设置（`ilesmod=0`）
- 无 IBM（`iibm=0`）
- 无 scalar（`numscalar=0`）
- 无 particle（本输入未启用粒子模块）
- 无 MHD（本输入未启用 MHD 场景）
- 无纤维（无任何纤维相关配置）

## C. 当前基线关键参数摘要
以下参数来自冻结的 `input.i3d`，用于快速审查与回归对比。

- `itype = 3`：流动类型为槽道流。
- `nx, ny, nz = 128, 65, 64`：计算网格分辨率。
- `xlx, yly, zlz = 8.0, 2.0, 4.0`：计算域尺寸。
- `re = 4200`：当前设置下的雷诺数控制参数。
- `cpg = F`：关闭恒压梯度模式（按该输入解释）。
- `istret = 2`：壁面双侧加密拉伸网格。
- `beta = 0.259065151`：壁面拉伸参数。
- `dt = 0.005`：时间步长。
- `ifirst = 1, ilast = 100000`：时间推进区间（长算）。
- `initstat = 40000`：统计启动步数。
- `ilesmod = 0`：DNS 模式。
- `numscalar = 0`：不求解标量输运。
- `iibm = 0`：禁用浸没边界法。

## D. 编译方式（Ubuntu + CMake + MPI 示例）
以下命令为推荐示例（不要求 root）：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

说明：
- 依赖 MPI 与 2DECOMP&FFT（由工程 CMake 查找）。
- 目标可执行文件为 `build/bin/xcompact3d`（默认输出路径）。

## E. 运行方式（标准示例）
推荐使用独立运行目录，避免污染源码树：

```bash
mkdir run_channel_re180_baseline
cp baseline/channel_re180_single_phase/input.i3d run_channel_re180_baseline/
cd run_channel_re180_baseline
mpirun -np 8 ../build/bin/xcompact3d > run.log 2>&1
```

建议：
- 将 `run.log` 作为主日志归档文件。
- 每次基线复现实验使用独立目录，便于回退与审查。

## F. 输出说明
一次完整长算后，建议重点检查以下输出：

- 日志文件（如 `run.log`）：
  - 用于检查是否正常启动、是否有 NaN、是否发生 MPI abort、是否完整结束。
- `statistics/` 目录：
  - 用于存放时均值与二阶统计量文件，是基线对比核心。
- 主要统计文件（如 `umean.dat*`, `vmean.dat*`, `wmean.dat*`, `uumean.dat*`, `vvmean.dat*`, `wwmean.dat*`, `uvmean.dat*`）：
  - 用于平均剖面与 Reynolds 应力审查。
- restart 文件（如存在）：
  - 用于续算、断点恢复与复现一致性检查。

## G. 基线约束
后续任何新分支（尤其涉及刚性/柔性纤维）在引入新模块前，必须满足：

1. 该单相基线可重复运行；
2. 关键统计输出可与本基线对齐比较；
3. 任何偏离都必须可解释、可审查、可回退；
4. 未通过本基线回归前，不得宣称新物理功能结果有效。

---
本目录仅用于“2.0：冻结单相槽道基线”，不包含 2.1/2.2/2.3 的任何开发内容。
