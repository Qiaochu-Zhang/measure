# V1_18 发布说明

V1_18 基于 V1_17，将 trench 改为先测 ROI 内全部合规候选，再筛选合格 trench，最后按 max-number 输出。line 的共享边算法保留，旧版文件与 ZIP 均保持原样。

- [完整代码 ZIP](sem_cd_measure_200k_batch_V1_18_complete_package.zip)
- [完整参数手册](V1_18_RUN_GUIDE.md)
- [修改记录](sem_cd_measure_200k_batch_V1_18_complete_package/CHANGELOG.md)
- [自检报告](sem_cd_measure_200k_batch_V1_18_complete_package/validation_report.json)
- [发布与旧版对照验证](audits/v1_18_release_evidence.json)

最终选择仍遵守原中心/左右、显式间距先验、max-number 和 min-number 规则；筛选对象必须先通过原稳定性条件。阈值、采样、跟踪、Viterbi/ERF、continuity、飞点、CD/LER/LWR、旋转和 PSD 数学均保持不变。被淘汰或未入选的候选仍保留在新增审计表中。

trench_selection 记录全部已测候选、有效率及入选状态；trench_source_samples 记录所有原始逐点边缘。最终输出通过 trench_source_id 对照来源，筛选不会移动或额外补齐边缘。

pitch 仍自动输出相邻一个 line + 一个 trench 的周期，复用全部合格来源边缘，按 max-number 个不同完整周期等权平均，仅输出旋转结果，位于主 Excel 最后一列。不跨过未合格或漏检的物理 trench，数量不足保留实际值并标记 REVIEW。

```bash
python -m pip install -r requirements_V1_18.txt
python self_check_V1_18.py --quick
python sem_cd_measure_200k_batch_V1_18.py --root examples/input_dark_line --pattern trench --pixel-size 1 --trench-reference-nm 60 --max-number 3 --output results_trench_V1_18 --no-auto-machine-comparison
```

原 V1_15、V1_16、V1_17 的目录、文档、报告及 ZIP 保持不变。V1_18 可独立解压运行，完整自检的 60 个集成场景全部通过，其中 12 个覆盖本次 trench 改动；测试报告随包提供。测试基于合成图，尚未验证用户实际失败的 SEM 图片。

同一张设定 trench=60nm 的合成图，将测量中心放在假暗条处时，V1_17 的 mixed CD 为 53.2607nm；V1_18 排除假中心后为 59.9309nm。对照中的同一真实 trench 原始坐标与有效点判定完全一致，默认及可选拟合下的 line 结果也与 V1_17 完全一致。

参与主汇总或 pitch 平均的对象可能随全量搜索和质量筛选改变，因此最终均值可能变化；边缘坐标算法、统计公式和输出定义保留。
