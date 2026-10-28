# V1_20 发布说明

基于 V1_19（3629670）独立更新，保留所有历史版本。每张样品图片默认生成 V10、V13 两张独立 PSD 图，每张图同时显示红色 LER、蓝色 LWR，只计算旋转修正后的坐标。

PSD 汇总 ROI 内全部合格结构的连续有效边缘功率谱，不再按主表 max-number 截断；先逐段 FFT 再按实际长度加权，不拼接边缘或跨图片平均。原主测量 CD/LER/LWR、group4 和 pitch 算法保持。

- [独立完整包 ZIP](sem_cd_measure_200k_batch_V1_20_complete_package.zip)
- [代码目录与 README](sem_cd_measure_200k_batch_V1_20_complete_package/README.md)
- [完整参数手册](V1_20_RUN_GUIDE.md)
- [折线、频点数量和 wiggling 说明](V1_20_PSD_POINTS_AND_WIGGLING.md)
- [完整自检结果](sem_cd_measure_200k_batch_V1_20_complete_package/validation_report.json)
- [旧版保留、数值对照与独立 ZIP 验证](audits/v1_20_release_evidence.json)

新增 LF/MF/HF 标签和两个竖直虚线边界。默认 100nm / 20nm 波长对应 0.01 / 0.05 nm⁻¹，是可调分析配方，不是论文通用阈值。默认去趋势改为旋转后仅去均值；可用 --psd-detrend linear 保留旧策略。输出整图谱数据、原生分段谱、支持数、频带方差和辅助相关性，缺少支持留空。

```bash
python -m pip install -r requirements_V1_20.txt
python self_check_V1_20.py --quick
python sem_cd_measure_200k_batch_V1_20.py --root examples/input_wiggling --pattern trench --pixel-size 1 --trench-reference-nm 60 --output demo_output_V1_20 --no-auto-machine-comparison
```

请按实际图片修改 pixel-size 和参考宽度。正常启用 PSD 和统计绘图时每图两张；显式 skip-psd / skip-statistics-plots 会关闭相应输出。无有效连续段仍生成说明性占位图。

参考论文 DOI：10.1109/CSTIC.2017.7919761。公开全文暂不可得，本版依据公开摘要及课题组说明借鉴中低频 LER/LWR 对照思想，未声称复现论文全部方法；不含 Palasantzas 拟合或噪声去偏。验证使用合成图，尚未用用户真实 SEM 图片确认测量效果。
