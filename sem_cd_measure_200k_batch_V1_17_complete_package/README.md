# SEM CD / LER / LWR — V1_17

当前版本 `V1_17`（2026-09-21）：line 改用原图中的 trench 搜索方法；相邻两条合格 trench 的内侧边缘就是 line 的两边。保留自动旋转 pitch CD，主 Excel 最后一列为 `旋转_pitch_CD_nm`。详见 [修改记录](CHANGELOG.md) 与 [手册第 12 节](RUN_GUIDE.md#12-line-修复与自动-pitch-cd2026-09-21)。

基于 V1_16（c0ae220）的独立完整包。原 V1_15、V1_16 的代码、文档、报告和 ZIP 保持原样。V1_17 无需旧目录即可运行；背景、端部裁剪、trench 边缘算法和统计/PSD 规则继续沿用。

**完整运行说明及全部 86 个命令行选项见 [RUN_GUIDE.md](RUN_GUIDE.md)。** 其中列出每个参数的含义、默认值、单位、取值范围、冲突规则和示例，也说明所有新算法内部经验常数。

## 运行

建议 Python 3.12；支持 Python 3.10+ 语法。解压后在本目录运行：

```bash
python -m pip install -r requirements_V1_17.txt
python self_check_V1_17.py --quick
python sem_cd_measure_200k_batch_V1_17.py --root examples/input_bright_line --pattern trench --pixel-size 1 --trench-reference-nm 60 --output demo_output --no-auto-machine-comparison
```

自己的图片：

```bash
python sem_cd_measure_200k_batch_V1_17.py --root /path/to/png_images --pattern trench --locator-mode auto --pixel-size 1.3181 --trench-reference-nm 60 --end-trim-fraction 0.05 --output /path/to/new_results
```

像素标定和参考宽度必须按实际图片替换。输出目录须不存在或为空。Linux/macOS 可用 `./run_V1_17.sh`，Windows 可用 `./RUN_V1_17.ps1`，参数原样转发；PowerShell 包装器未在本次 Linux 环境执行。

## 选择定位模式

| 参数值 | 使用场景 |
|---|---|
| `auto`（默认） | 每图观察整个有效区域的 X 灰度周期，判断使用以下哪种模式。 |
| `dark-line` | 保留 V1_14 定位：line 中有假暗条，需要利用灰度及周期关系区分真假 trench。 |
| `bright-line` | line 整体亮，trench 整体暗；两类 K-means 自适应阈值、区域多数像素和周期检查，不强制跳过相邻暗区。 |

这里的 line 是 trench 之间的分隔带。**测暗 trench 时仍用 `--pattern trench`，并可以选择 `--locator-mode bright-line`。** `--pattern line` 先搜索原图中的全部合规 trench，再取相邻合格 trench 之间的 line；不对图像反相。

`--locator-majority 0.70` 控制新模式中暗 trench 核心及两侧亮 line 的最低像素比例。定位阈值自适应计算，不是固定灰度值；`--threshold-left/right` 仍只控制最终边缘交点。

同一个模式供 V10、V13 以及 PSD 使用。后续边缘搜索、跟踪、Viterbi、ERF、continuity、flyer、CD/LER/LWR 和双引擎 PSD 沿用 V1_14 流程。

## 背景与端部

先识别持续成对边界，排除背景。只采用一种长度近似：包含前景的各 X 列首尾 Y 跨度的均值；每列从自身上下端各去掉 `ceil(平均长度 × 0.05)` 像素。平均长 440px 时每端去掉 22px。内部背景缺口继续排除。

主测量的完整灰度平均窗口及 PSD 的单行都必须落在裁后有效区域。补点和 PSD 插值不能跨背景或端部。原图和导出坐标不裁切、不重新编号。默认测量长度仍最多 360px，需要更长观察范围时设置 extend-length。

`--end-trim-fraction 0` 只关闭端部裁剪；`--no-auto-roi` 只关闭自动背景识别。用 `--locator-mode dark-line --no-auto-roi --end-trim-fraction 0`，并保持其他配方一致，可对照 V1_14 无掩码流程。

## 查看结果

- `CD_measurement_200K_V1_17_results.xlsx`：首表 `measurement_summary` 每图 mixed、V13、V10 三行；最后一列 `旋转_pitch_CD_nm`。
- `ROI/`、`roi_summary.csv`：裁后区域、裁前/裁后边界、估计平均长度及每端裁剪像素数。
- `locator_summary.csv`：请求、实际和建议模式，自动判断分数/歧义标记，周期和自适应阈值。
- `locator_candidates.csv`：两个引擎的全部定位候选、原图坐标、选择结果、亮度比例和失败原因。
- `line_trench_pairs.csv` / `line_source_samples.csv`：line 与左右来源 trench 的配对和逐点边缘；与主采样表按来源 ID、sample_index 对照，共享边坐标严格一致。
- `pitch_periods.csv` / `pitch_samples.csv`：相邻周期、逐点边缘和入选记录。按 max-number 个完整周期等权平均；不足时报告实际数和 REVIEW，空值不填 0。mixed pitch 来自 V13。
- `per_sample_results.csv`：主边缘及 background_excluded；非 compact 模式还记录完整平均窗口的 Y 范围。
- `PSD/PSD_V10.xlsx`、`PSD/PSD_V13.xlsx` 及 CSV：保持 V1_14 的单行边缘、逐结构/逐连续段 periodogram，不做分组或频谱平均。
- `settings.json`、`processing_errors.csv`：运行配方和异常。

mixed 的 CD/LWR 来自 V13，LER 来自 V10。主 CD/LER 不分组，LWR 默认 group4；PSD 默认独立提取单行边缘。质量、指标定义及其他输出详见完整手册。

## 测试与限制

```bash
python self_check_V1_17.py --output my_validation_report.json
```

测试涵盖原有数值回归、背景、默认双引擎/双极性、Viterbi/ERF、PSD、故障导出和防覆盖，并新增两种定位、auto/手动一致性、噪声/灰度偏移、窄假暗条、不同长度端部及强制恢复/PSD 插值排除测试。实际结果见 [validation_report.json](validation_report.json)。

`examples/input_bright_line` 和 `examples/input_dark_line` 是已知 60px trench、160px 节距、440px 长度的合成对照，可使用 pattern=trench 测 60px 沟槽，也可使用 pattern=line 测 100px 间隙。input_trench/input_line 是保留的旧版反相对照。这些不是用户真实 SEM 数据。

若真假暗条在宽度、亮度、节距上完全相同，灰度无法唯一判别身份；auto 默认 bright-line，并对有限周期/临界分数给出歧义标记。几何背景筛选也可能误判类似条纹、弱对比度或大倾角结构。请核查 ROI、定位和边缘标注；本次软件回归不代表真实 SEM 精度认证。

## 文件结构

| 文件 | 职责 |
|---|---|
| `sem_cd_measure_200k_batch_V1_17.py` | 参数、区域/模式编排、双引擎、结果导出 |
| `cdsem_lines.py` / `check_shared_lines.py` | 相邻合格 trench 共享边缘构成 line，以及逐点、不跨沟槽、数量和 PSD 回归 |
| `cdsem_pitch.py` / `check_line_pitch.py` | 旋转 pitch 和 line/pitch 回归 |
| `cdsem_localization.py` | 新 bright-line、auto、有效灰度轮廓 |
| `cdsem_locator.py` | 保留的 V1_14 盆地定位及公共亮度工具 |
| `cdsem_regions.py` | 读图、背景排除、端部裁剪、区域图 |
| `cdsem_engine.py` | 两种定位共用的边缘引擎 |
| `cdsem_refinement.py` / `cdsem_statistics.py` / `cdsem_psd.py` | 延用的修正、统计和 PSD |
| `self_check_V1_17.py` / `check_localization.py` | 完整自检及新算法测试 |

内部部分函数名包含 v115 或 V1.15，是继承的历史算法标识；当前发布版本以主程序 SCRIPT_VERSION=V1_17 为准。PACKAGE_SHA256.txt 校验全部交付文件，不包括清单自身。

## 测 line

```bash
python sem_cd_measure_200k_batch_V1_17.py --root examples/input_dark_line --pattern line --pixel-size 1 --line-reference-nm 100 --space-reference-nm 60 --max-number 3 --output results_line_V1_17 --no-auto-machine-comparison
```

line-reference 是 line 参考宽度；space-reference 是来源 trench 的宽度，可不填自动估计。max-number 限制最终 line 条数和 pitch 平均周期数，来源 trench 会先全部搜索。只有两条相邻合格 trench 才能定义一条 line；缺少邻居时保留缺失状态。
