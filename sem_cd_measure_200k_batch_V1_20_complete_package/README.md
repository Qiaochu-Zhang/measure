# SEM CD / LER / LWR — V1_20

V1_20（2026-10-28）基于 V1_19 独立复制，新增面向 wiggling 观察的整图 PSD。旧版本目录、文档、验证结果与 ZIP 保留。

## 运行

建议 Python 3.12；Python 3.10+。在本目录运行：

```bash
python -m pip install -r requirements_V1_20.txt
python self_check_V1_20.py --quick
python sem_cd_measure_200k_batch_V1_20.py --root examples/input_wiggling --pattern trench --pixel-size 1 --trench-reference-nm 60 --output demo_output_V1_20 --no-auto-machine-comparison
```

自己的图片：

```bash
python sem_cd_measure_200k_batch_V1_20.py --root /path/to/images --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --psd-split-wavelength 100 --psd-high-split-wavelength 20 --output /path/to/new_results
```

按图像标定修改 pixel-size 和参考宽度。输出目录必须不存在或为空。Linux/macOS 使用 `./run_V1_20.sh`，Windows 使用 `./RUN_V1_20.ps1`。Windows 包装器未在本 Linux 环境执行。

## 新 PSD 输出

每张样品图片默认生成 **2 张 PSD 图**，不混合 V10 和 V13：

- `PSD/V10/<相对目录>/<原完整文件名>.psd.png`
- `PSD/V13/<相对目录>/<原完整文件名>.psd.png`

每图只有一个坐标面板，红色为 **LER PSD**，蓝色为 **LWR PSD**；只使用中心线旋转修正后的坐标。两条竖直虚线把频率轴分为 LF / MF / HF。默认分界为 0.01 和 0.05 nm⁻¹，对应波长 100 和 20 nm；它们是可调分析配方，不是引用论文的通用阈值。

PSD 纳入 ROI 内全部通过稳定性筛选的结构，不受主表 `max-number` 截断。先在每条结构的连续有效段内计算功率谱，再按实际测量长度汇总到同图同引擎的两条曲线；不拼接 trench 坐标，不跨图片平均。LER 是左右边缘功率谱的均值，LWR 是未分组宽度的功率谱；先平均边缘坐标会抵消不同 trench 的摆动，因此不采用。

默认仍独立提取单行、1px 步长、无 X/Y 平均的边缘；改为旋转后仅去均值（constant），减少额外线性去趋势对低频的抑制。需要旧去趋势策略可显式 `--psd-detrend linear`。原主表 CD/LER/LWR、group4 LWR、旋转与未旋转统计及 pitch 算法保持不变。

没有足够有效连续段时仍生成两张标明缺少支持的占位图，不将缺失数据当作零。`--skip-psd` 关闭全部 PSD；`--skip-statistics-plots` 关闭包括 PSD 在内的统计图，仍输出谱数值和工作簿。

## 查看数值和说明

- [PSD 点数、折线原因与 wiggling 阅读说明](PSD_POINTS_AND_WIGGLING.md)：回应频点数量问题，列明公式、限制与参考来源。
- [完整运行指南](RUN_GUIDE.md)：全部参数；第 4.8 和第 15 节详述新版 PSD。
- `PSD/image_psd_summary.csv`：每图、每引擎、LER/LWR 的结构数、段数、频点数及真实频率分辨率。
- `PSD/image_psd_curves.csv`：最终绘图两条曲线，每频点带贡献结构数、段数和权重；`V10/V13_psd_curves.csv` 分引擎保存同样数据。
- `PSD/image_psd_band_summary.csv`：在原生 FFT 频点积分的 LF/MF/HF 方差、3σ、LER/LWR 功率比和边缘相关性指标；不可测频带留空。
- `PSD/PSD_V10.xlsx`、`PSD/PSD_V13.xlsx`：整图汇总、频带统计、原始分段摘要及来源。
- `PSD/per_structure_psd_summary.csv` / `per_structure_psd_curves.csv`：旋转后的逐结构、逐连续段原生谱，用于审计，不再各画一条曲线。
- `PSD/edge_coordinates.csv` / `measurement_manifest.csv` / `PSD_settings.json`：边缘、实际配方及汇总定义。
- `CD_measurement_200K_V1_20_results.xlsx`：原主测量汇总，mixed 的 CD/LWR 来自 V13，LER 来自 V10。

## 保留功能与验证

支持 TIFF/PNG/JPEG/BMP/WebP 混合输入、中文路径、bright-line/dark-line/auto 定位、背景排除、端部裁剪、trench 全测后筛选、相邻 trench 共享边定义 line、旋转 pitch。测亮间隙用 `--pattern line --line-reference-nm 100 --space-reference-nm 60`。V1_19 既有测量参数说明见本版完整指南。

```bash
python self_check_V1_20.py --output validation_report.json
```

自检涵盖原有测量回归及新增的 PSD 已知低频波形、同相/反相摆动、功率平均、频率支持、旋转间距、引擎隔离、频带标记、空结果、全量结构和同名不同格式导出。发布对照见仓库 `audits/v1_20_release_evidence.json`。`examples/input_wiggling` 为可复现合成图，不是实际 SEM 数据；旧版样例也均为合成图。本版没有噪声去偏或 Palasantzas 拟合，也不自动判定工艺原因。
