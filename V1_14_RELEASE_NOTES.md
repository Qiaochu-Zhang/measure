# V1_14 发布说明

基于 V13_modified，新增独立 V1_14 完整包；原有两个 ZIP 和原版解压源码没有改动。

- [下载完整代码包](sem_cd_measure_200k_batch_V1_14_complete_package.zip)
- [使用说明与参数](sem_cd_measure_200k_batch_V1_14_complete_package/README.md)
- [完整运行说明与全部命令行参数](V1_14_RUN_GUIDE.md)
- [详细变更记录](sem_cd_measure_200k_batch_V1_14_complete_package/CHANGELOG.md)
- [完整自检结果](sem_cd_measure_200k_batch_V1_14_complete_package/validation_report.json)
- [自检脚本](sem_cd_measure_200k_batch_V1_14_complete_package/self_check_V1_14.py)

## 本次交付

1. 修复阈值交点搜索范围问题，保留原宽度/拓扑/跟踪约束，不用合成补点冒充成功检测。
2. 默认背景排除，自动定位未手动指定的 ROI；保存背景掩码诊断，坐标仍属于原图。无可靠目标的图记 `SKIPPED_BACKGROUND`。
3. 总结 Excel 首表按 `image、status、method` 开头，每图 mixed、V13、V10 三行；另保留每图一行的完整宽表。
4. PSD 分 V10/V13 两套，默认单行、1 px 步长重新提边缘；逐结构、逐连续段 periodogram，不做 group4、Welch 或跨结构/跨图平均。主测量仍保留原 group-size 配方；如需 PSD 复用主测量边缘，可指定 `--psd-edge-source measurement`。
5. 修复单点/单组假零粗糙度、部分引擎结果丢失、附属输出阻断主结果、损坏图片漏出汇总、非空目录覆盖等问题。默认条件汇总只纳入 OK。
6. 运行源码从 17 个文件整合为 7 个；只保留当前 trench/line 流程的实际依赖，继续保留 Viterbi、ERF、continuity、flyer、机台对照和诊断输出。

## 验证要点

- 随包示例 V13 主测量有效边缘：旧搜索 175/512，新版默认 512/512；未开启强制连续性补点。
- 默认 V10、V13 各生成 32 条有效原始谱：4 条结构 × 2 坐标 × 4 信号，每条使用 360 个单行采样位置。
- 关闭自动 ROI、使用 `--threshold-search legacy` 后，旧版示例 mixed 数值可复现到测试容差内。
- 检查了坐标旋转、共同掩码、单点/单组、已知波形 PSD、背景排除、偏置前景、内部断带、各种导出/引擎异常及防覆盖。
- 额外执行了包含机台对照、统计图、双引擎标注/PSD 图的完整示例；没有诊断错误。发布包还经过独立解压检查。

这些是软件回归和随包/合成图验证，**不代表真实 SEM 精度认证**。当前没有用户提供的真实混合背景 SEM 图片；自动背景筛选可能误删低对比度、短线或特殊结构，须检查 `ROI/` 图。单行 PSD 对噪声更敏感，尚未做 SEM 噪声去偏。

运行方式、关闭背景筛选、调整阈值、结果定义和两套 PSD 的读取方法均见包内 README。
