# V1_20 修改记录

## 2026-10-28：整图旋转 PSD 与 wiggling 观察

基于 V1_19（3629670330d358fa028a0dcf384fd286b9055321）创建独立完整包。

- V10、V13 两种边缘引擎分别生成一张图；每图同时显示红色 LER 与蓝色 LWR。
- PSD 仅计算旋转修正坐标，不再计算 unrotated；主测量双坐标统计保留。
- PSD 使用 ROI 内全部合格结构；不再按主表 max-number 截断。保留原稳定性筛选及背景/补点规则。
- 先计算原生连续段 periodogram，再按段的 N×旋转后间距汇总功率。左右边缘功率等权形成 LER，宽度功率形成 LWR；不做坐标拼接、group4 或跨图/引擎平均。
- 不同长度、不同旋转间距的谱在最长物理记录的频率网格内按密度线性对齐，仅在各段真实正频率范围内贡献；保存每频点支持数及原始分段 CSV。
- 默认 psd-detrend 从 linear 改为 constant；保留 linear 可选项。窗函数仍默认 Hann。
- 图上新增 LF/MF/HF、两个竖直虚线边界；沿用 psd-split-wavelength=100nm 作为 LF/MF，新增 psd-high-split-wavelength=20nm 作为 MF/HF，校验 0<高频分界波长<低频分界波长。
- 新增整图汇总、绘图曲线和原生频点频带积分。输出低频 LER/LWR 功率比、中心线方差与归一化边缘相关性，仅供观察，不宣称复现论文拟合/判定公式。
- 无有效段时生成明确占位图，不补零；保留 skip-psd、skip-statistics-plots 的显式关闭行为。
- 新增 PSD_POINTS_AND_WIGGLING.md 回答折线和点数原因，并注明论文公开资料范围、频段阈值来源和频率分辨率限制。
- 新增已知波形及整图 PSD 回归、合成 wiggling 样例；完整报告见 validation_report.json。

原 V1_19 和更早所有交付物逐字节保留；灰度、边缘定位/提取、主 CD/LER/LWR、共享边和 pitch 数学保持。line/trench 模块只为 PSD 额外暴露全部来源；默认结果目录更新为 V1_20，并继续排除 V1_19 历史输出。
