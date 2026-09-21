# V1_16 发布说明

V1_16 是独立完整版本，包含 line 测量修复与自动旋转 pitch CD。原 V1_15 的代码、文档、测试报告和 ZIP 均恢复并保留为原始提交 `41e3294` 的内容，逐文件及 ZIP 字节校验一致。

- [V1_16 完整包](sem_cd_measure_200k_batch_V1_16_complete_package.zip)
- [完整运行说明与全部 86 个 CLI 选项](V1_16_RUN_GUIDE.md)
- [快速开始](sem_cd_measure_200k_batch_V1_16_complete_package/README.md)
- [完整自检报告](sem_cd_measure_200k_batch_V1_16_complete_package/validation_report.json)
- [版本保留、算法迁移及独立解压验证](audits/v1_16_release_evidence.json)
- [原 V1_15 发布说明](V1_15_RELEASE_NOTES.md)

## 本版修改

1. 修复候选宽度固定上限约 145nm 导致宽 line 被拒绝的问题，按目标参考宽度扩展范围，line/trench 共用。
2. 修复 line 内部暗条在反相后造成结构分裂、自动估宽错误、ROI 漏选及阈值参考污染的问题。定位使用完整结构，最终边缘仍取自原图。
3. line/trench 均自动测量相邻一个 line + 一个 trench 的完整周期。对相邻同类结构的同侧边缘做 PCA 法向投影，按 `max-number` 个不同周期等权平均，仅输出旋转后的 pitch CD。
4. Excel 首表最后一列为 `旋转_pitch_CD_nm`，mixed/V13 使用 V13，V10 使用 V10；附实际周期数、状态和逐周期/逐点审计。数量不足标记 REVIEW，完全缺失保留空值，不凑数。
5. 为测 N 个完整周期，额外测量一个相邻结构；额外结构不进入主 CD/LER/LWR 或 PSD。背景、裁去的端部及合成补点不进入 pitch。

## 独立版本

- 目录：`sem_cd_measure_200k_batch_V1_16_complete_package/`
- 主程序：`sem_cd_measure_200k_batch_V1_16.py`
- 自检：`self_check_V1_16.py`
- 依赖：`requirements_V1_16.txt`
- 启动器：`run_V1_16.sh`、`RUN_V1_16.ps1`
- 主 Excel：`CD_measurement_200K_V1_16_results.xlsx`
- 默认输出目录：`CD_measure_output_200K_V1_16`
- `settings.json` 的 script_version / patch_version 均为 `V1_16`。

V1_16 不依赖 V1_15 目录；可单独解压运行。递归扫描继续排除 V1_14、V1_15 和 V1_16 的历史结果目录。不同版本的结果建议放在输入目录之外，便于对照。

## 验证与记录

完整自检覆盖 36 个独立端到端/故障场景，其中 13 个 line/pitch 场景，另含解析几何、背景和版本隔离检查。覆盖宽线、内部暗条、同图双模式、反相等价、倾斜、max-number=1/2/5、周期不足、Excel 末列、双引擎 PSD 和 Viterbi/ERF。实际结果见自检及发布验证 JSON。

原始失败复现基于 `41e3294`：60nm trench 能测量，180nm line 被上限拒绝；100nm line 内含 40nm 暗条时定位失败，自动参考宽度还可能误估为约 33nm。记录见 [修复前证据](audits/v1_16_line_pitch_baseline.json)。本版合成样例中 V13 测得宽 line 约 179.953nm、pitch=240.000nm，带暗条 line 约 99.959nm、pitch=160.000nm。

改进最初在提交 `ae11296` 中实现；本次将其独立发布为 V1_16，并恢复 V1_15。保留 Git 提交历史，未改写旧提交。V1_16 发布校验还对迁移前后算法做归一化 AST 比较，确认版本整理未改变这些改进。

测试使用随包和新增合成图，尚未验证用户实际失败的 SEM 图片。灰度完全不可区分的结构仍需要人工核查。

## 运行示例

```bash
python -m pip install -r requirements_V1_16.txt
python self_check_V1_16.py --quick
python sem_cd_measure_200k_batch_V1_16.py --root examples/input_dark_line --pattern line --pixel-size 1 --line-reference-nm 100 --max-number 3 --output results_line_V1_16 --no-auto-machine-comparison
```

将参考宽度和像素尺寸替换为实际标定；同图的 line 与 trench 参考宽度可能不同。
