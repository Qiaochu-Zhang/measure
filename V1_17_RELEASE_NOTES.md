# V1_17 发布说明

V1_17 基于 V1_16，把 line 测量改为：先用原 trench 方法搜索全部合规 trench，相邻两条合格 trench 之间就是 line，左右边缘分别与它们共享。V1_15、V1_16 的文件与 ZIP 保持原样。

- [V1_17 完整包](sem_cd_measure_200k_batch_V1_17_complete_package.zip)
- [完整运行说明](V1_17_RUN_GUIDE.md)
- [修改记录](sem_cd_measure_200k_batch_V1_17_complete_package/CHANGELOG.md)
- [自检报告](sem_cd_measure_200k_batch_V1_17_complete_package/validation_report.json)
- [发布验证](audits/v1_17_release_evidence.json)
- [修复前后对照](audits/v1_17_shared_line_baseline.json)

V1_16 仍然独立检测反相后的 line。内部暗条较宽时，亮线可能被拆成窄片段，导致没有候选；同图 trench 可以测出。V1_17 直接沿用 trench 的实测边界，因此 line 内部的灰度不再触发一遍独立 line 拒绝流程。

搜索全部来源 trench 后，按中心及左右邻居选取最多 max-number 条 line；不会跨过不合格或漏检的 trench 拼接 line。V10/V13、旋转、CD/LER/LWR、分组、PSD、可选边缘修正等规则保留。新增配对和来源边缘表用于逐点核查。

pitch 继续自动测量，每个周期是相邻一个 trench + 一个 line，按 max-number 个不同周期等权平均，只输出旋转结果，放在主 Excel 最后一列。line 流程复用来源 trench，两条 trench 即可测一个 line 和一个 pitch。

```bash
python -m pip install -r requirements_V1_17.txt
python self_check_V1_17.py --quick
python sem_cd_measure_200k_batch_V1_17.py --root examples/input_dark_line --pattern line --pixel-size 1 --line-reference-nm 100 --space-reference-nm 60 --max-number 3 --output results_line_V1_17 --no-auto-machine-comparison
```

`line-reference-nm` 填 line 宽度，`space-reference-nm` 填物理 trench 宽度，未给出时自动估计。左右阈值作用于来源 trench，line 左边继承左 trench 的右阈值，line 右边继承右 trench 的左阈值。请使用实际像素标定与参考宽度。

完整回归覆盖逐点共享边、宽暗条、全部候选搜索、数量边界、缺失/拒绝中间 trench、Viterbi/ERF、非对称阈值和 PSD；发布审计验证原 trench 数值保持一致。具体运行结果见链接报告。测试使用合成图，尚未验证用户实际失败的 SEM 图片。

本次完整自检通过 48 个独立端到端/故障场景，其中 12 个共享边专项场景，另有解析几何、背景和版本隔离检查。宽暗条对照中，V1_16 line 为 ERROR；V1_17 为 OK，line CD=100.1133nm、pitch CD=159.9993nm（合成设定为 100nm、160nm）。
