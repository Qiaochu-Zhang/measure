# SEM CD / LER / LWR — V1_14

基于 V13_modified 改进的独立完整包，支持 PNG 中近竖直的 trench（暗沟槽）和 line（亮线）。不需要旧版目录，也不会导入旧包。

## 主要变化

- 默认修复“阈值交点在梯度峰外侧却没有找到”的漏检，仍保留亮暗拓扑、宽度范围及跟踪跳变约束，不靠强制补点提高检出率。
- 默认识别具有持续成对边界的区域，排除纯色、低对比度平坦背景和不符合线条几何的区域。自动调整未手动指定的定位中心/范围；原图不裁剪、不改写，导出的坐标仍属于原图。
- 总结 Excel 第一张表为 `measurement_summary`：`image`（文件名）、`status`、`method` 在前；每张图依次为 **mixed、V13、V10** 三行。每图一行的完整 24 指标宽表保留在 `image_summary`。
- PSD 分别生成 V10、V13 结果，**不做 group4/其他分组，不做 Welch 平均，不跨结构或跨图片平均频谱**。每条结构的每个足够长的连续有效段独立输出。
- PSD 默认重新以单行、1 px 沿线步长、无横向平滑提取边缘。主 CD/LER/LWR 的 32 行平均等原配方与 PSD 配方独立。
- 单点 LER、单组 LWR 输出空值而不是 0；保留部分引擎结果，附属导出失败不丢主测量结果，默认汇总只纳入 `OK`。
- 合并 V10/V13 重复引擎并移除未调用的历史流水线、旧入口和 via 依赖。运行源码从原包 17 个文件精简为 7 个，另附一个自检脚本。

本版并未实现 SEM 噪声去偏、物理 CD 校准或计量不确定度模型。输出是指定图像、标定和配方下的 measured CD/LER/LWR，不能仅凭软件运行成功称为专业绝对量测。

## 1. 安装和快速运行

建议 Python 3.12；代码使用 Python 3.10+ 语法。本次实际验证环境及版本见 `validation_report.json`，并未逐一验证所有依赖的最低版本或所有操作系统。

在本包目录中运行：

```bash
python -m pip install -r requirements_V1_14.txt
python self_check_V1_14.py --quick
python sem_cd_measure_200k_batch_V1_14.py --root examples/input_trench --pattern trench --pixel-size 1 --trench-reference-nm 60 --max-number 4 --output demo_output --no-auto-machine-comparison
```

`demo_output` 必须不存在或为空。为了避免覆盖、混入旧图/旧表，同版本重跑也需要新的输出目录。

自己的 trench 图：

```bash
python sem_cd_measure_200k_batch_V1_14.py --root /path/to/png_images --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --output /path/to/new_results
```

亮 line：

```bash
python sem_cd_measure_200k_batch_V1_14.py --root /path/to/png_images --pattern line --pixel-size 1.3181 --line-reference-nm 60 --output /path/to/new_line_results
```

这里的像素尺寸和参考 CD 只是命令示例，必须替换为适合你的图片/标定的数据。默认 `--max-number 3` 是每图最多选择的结构数；需要更多可显式增大，默认最少数量为 `min(3,max-number)`。参考宽度用于搜索约束，不是测量真值。

Windows 可以用激活环境后的 Python，也可运行 `./RUN_V1_14.ps1 --root ... --pattern trench ...`；Linux/macOS 可用 `./run_V1_14.sh --root ... --pattern trench ...`。PowerShell 包装脚本未在本次 Linux 环境执行，核心 Python 流程已验证。

输入递归扫描 PNG，支持中文路径、相同文件名位于不同子目录。当前输出目录及带已知版本 settings 的历史输出子目录会被排除。

## 2. 背景排除如何工作

自动选区使用分块代表线扫寻找宽度合理的成对明暗边界，再检查它们在逐行图像中是否持续存在。定位使用的平均/中位数仅用于识别区域，测量读取原始图像；掩码不会把背景像素改造成结构。主测量的整段灰度平均窗口也必须处于有效区域。

输出：

- `ROI/相对路径.png`：背景调暗、选区用绿色框标示。必须检查这些图，尤其是低对比度、短线条、缺陷或特殊图案。
- `roi_summary.csv`：原图坐标下的选区范围、保留比例、候选检查数及阈值。
- 逐点表中的 `background_excluded`：被排除的位置不会进入指标或 PSD；即使 `--edge-continuity 100`，也不会填入背景。
- 全图没有可靠区域时：`SKIPPED_BACKGROUND`，指标为空，不伪造零值；仍保留总结和诊断输出。

默认 `--roi-min-contrast 6`（8bit 灰度尺度）、`--roi-min-length 24`（连续高度 px）。如果真实沟槽较暗弱或较短，可适当降低并人工核查，不应只追求有效率。

关闭自动排除：

```bash
python sem_cd_measure_200k_batch_V1_14.py --root /path/to/images --pattern trench --no-auto-roi --center-x 300 --center-y 250 --meas-area-width 400 --meas-area-height 300 --extend-length 240 --output /path/to/new_manual_results
```

手动中心、ROI 尺寸、测量长度的显式设置优先；背景掩码仍有效，除非指定 `--no-auto-roi`。8bit 输入不再做全动态范围拉伸，以免将接近纯色背景的微小噪声放大为强边缘。

限制：这不是语义分割模型。外观同样是持续平行条纹的非目标区域仍可能通过筛选；严重倾斜、弯曲、极低对比度或小于门槛的真实结构也可能被排除。当前验证包含合成背景/偏置前景，不包含用户真实 SEM 混合背景图片。16bit 输入仍有归一化处理，需要单独验证低动态范围图像。

## 3. CD / LER / LWR 的含义和 Excel 顺序

| 方法 | CD | 左/右 LER | LWR |
|---|---|---|---|
| mixed（第一行） | V13 | V10 | V13 |
| V13 | V13 | V13 | V13 |
| V10 | V10 | V10 | V10 |

每种方法同时输出旋转与未旋转结果。旋转方向来自结构中心线 PCA，不是把整幅原图重新采样。

- CD：有效局部宽度的均值。
- LER：未分组的左右边缘坐标各自 3σ。
- 主表 LWR：分组宽度均值的 3σ，仍由 `--group-size` 控制，默认 4；设为 1 可输出未分组 LWR。
- 原始未分组 LWR、有效采样数、有效组数另外记录在对象/引擎表。
- 单点/单组不足以估计粗糙度时输出空值并触发缺项/REVIEW，不输出假零值。至少两个点/组只是数学门槛，不是可信度保证。
- mixed 是已声明的跨引擎兼容指标，不应解释为同一组边缘的完整自洽统计。

主工作簿：`CD_measurement_200K_V1_14_results.xlsx`。

| 工作表 | 用途 |
|---|---|
| `measurement_summary` | 日常查看，文件名/status/method 在前，每图 mixed、V13、V10 |
| `image_summary` | 每图一行的完整指标、质量、来源和选区信息 |
| `coordinate_results` | 坐标与方法的长表 |
| `condition_summary` | 文件夹汇总，默认仅 OK |
| `trench_objects` / `engine_objects` | 结构级及单引擎结果、真实参与统计标记、支持数 |
| `per_sample_results` | 逐点边缘、有效性、背景及恢复标记 |
| `processing_errors` / `roi_summary` / `settings` | 错误、选区和运行配方 |

原有标注、旋转诊断、Viterbi、ERF、continuity、flyer 操作、任意正整数 group-size、机台对照和统计图功能保留。历史未调用的独立批处理入口不再打包；V13_modified 本身不支持的 via 没有额外引入。

## 4. PSD：两引擎、不平均

默认 PSD 使用独立提取的单行边缘，而不是主指标那份 32 行平均边缘。默认 `--psd-edge-source single-row` 对应：

```text
V10 单行原始边缘 → 每条结构/连续有效段 → periodogram → PSD_V10.xlsx
V13 单行原始边缘 → 每条结构/连续有效段 → periodogram → PSD_V13.xlsx
```

粗定位/选区可以使用代表线扫，但 PSD 最终边缘的 `average_range_px=1`、`smoothing_pixel=1`、`group_size=1`。每个引擎提供左右边缘、width、center 四类信号，旋转/未旋转两坐标。center 是左右边缘的几何中点，不是沿线分组平均。

如只希望取消 PSD 后续分组、仍使用主测量边缘，可显式选择：

```bash
python sem_cd_measure_200k_batch_V1_14.py --root /path/to/images --pattern trench --psd-edge-source measurement --output /path/to/new_results
```

该选项会保留主配方已经存在的 32 行灰度平均/横向平滑，但后续仍不做边缘分组或频谱平均。设置和每条谱均记录实际使用的配方，避免混淆。

PSD 文件位于 `PSD/`：

- `PSD_V10.xlsx`、`PSD_V13.xlsx`：各引擎逐结构、逐连续段的频谱统计。
- `V10_psd_curves.csv`、`V13_psd_curves.csv`：频率—PSD 数值，不是平均曲线。
- `per_structure_psd_summary.csv`、`per_structure_psd_curves.csv`：两引擎合并存放的长表，仍保留 engine/结构/段 ID，不进行合并平均。
- `edge_coordinates.csv`：实际输入 PSD 的边缘轨迹及其原始槽位和有效性，便于核查。
- `measurement_manifest.csv`、`PSD_settings.json`：来源、参数、支持数及非去偏声明。
- 每图每引擎的图：分别叠加各条结构和各连续段，不平均。

默认对每个至少 16 点的连续有效段做 Hann 窗、线性去趋势的 periodogram。去均值/去趋势不是分组平均；Hann 窗也不是噪声去偏。可用 `--psd-window boxcar --psd-detrend constant` 做归一化核查。PSD 使用 SciPy 的 density 定义，频率单位 `1/nm`，密度单位 `nm³`。[SciPy periodogram 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.periodogram.html)

缺口不被压缩，也不跨缺口 FFT；所有足够长的段分别保留。默认排除合成补点。可显式请求短缺口插值或包含合成点，但会记录标记，背景缺口始终禁止插值。旧版 Welch/groupN/跨图平均 PSD 按本次需求移除。

单行边缘对图像噪声更敏感，未平均的谱也会比平均谱更起伏；两者都不代表实现错误。本版没有测量噪声去偏。由于采样、滤波和去趋势不同，PSD 积分所得粗糙度不应直接要求等于主表 group4 LWR。

## 5. 质量与失败策略

- `OK`：当前工程检查通过，不代表经过 SEM 绝对准确度认证。
- `REVIEW`：低支持、部分指标/引擎缺失、人工补点或图像级附属输出异常，需要复核。
- `ERROR`：没有可用主指标或读取/主流程失败，仍保留固定列和错误说明。
- `SKIPPED_BACKGROUND`：没有找到可信的目标区域，保留空指标和选区诊断。

默认条件汇总、机台对照只纳入 OK；需要旧行为时加 `--include-review-in-summary`。逐图及逐点表仍保留 REVIEW 结果。

机台匹配为零会导出未匹配记录；机台文件损坏/缺失、PSD 导出失败等会写 `processing_errors`，不删除已经写出的主 CSV。最终主工作簿因文件占用或磁盘问题写入失败时，程序返回失败，先前 CSV 可用于恢复。`--stop-on-error` 会停止后续图像，但先导出已处理记录。

## 6. 验证和代码结构

```bash
python self_check_V1_14.py --output validation_report.json
```

完整自检涵盖：坐标旋转、共同有效掩码、单点/单组、已知正弦 PSD 积分、交替宽度不被 group4 消除、缺口分段、背景不可插值、阈值窗口约束、Viterbi/ERF、默认 trench/line、旧算法数值复现、偏置前景、内部背景和强制补点、空白图、单/双引擎失败、标注/PSD/机台异常、损坏 PNG、输出目录防覆盖。测试输出使用临时目录。

随包示例的 V13 主测量有效点从旧算法 175/512 提升到 512/512；默认没有开启 continuity 强制补点。默认 V10/V13 PSD 均有 32 条有效谱（4 结构 × 2 坐标 × 4 信号），并非 32 个独立样品。示例 trench/line 是同一几何的反相对照，不是两套独立真实数据。

| 源文件 | 职责 |
|---|---|
| `sem_cd_measure_200k_batch_V1_14.py` | CLI、双引擎编排、坐标统计、Excel/CSV 和机台对照 |
| `cdsem_engine.py` | 共享边缘引擎，参数选择 V10/V13 阈值定义 |
| `cdsem_locator.py` | 合并后的盆地定位和必要亮度工具 |
| `cdsem_refinement.py` | 候选选择、Viterbi、ERF、continuity 和 flyer |
| `cdsem_regions.py` | 图像读取、发现、背景选区和选区诊断 |
| `cdsem_statistics.py` | 统计及对照图所需的公共函数 |
| `cdsem_psd.py` | 原始边缘逐结构/逐段 PSD |

源码可直接阅读修改，没有动态拼接旧模块或运行时读取旧 ZIP。`PACKAGE_SHA256.txt` 列出本包文件摘要（不包含该摘要文件自身）。
