# V1_15 发布说明

基于 V1_14 新增独立完整包，旧版代码及 ZIP 保留。

- [下载完整包](sem_cd_measure_200k_batch_V1_15_complete_package.zip)
- [完整运行说明：全部 86 个 CLI 参数](V1_15_RUN_GUIDE.md)
- [快速开始](sem_cd_measure_200k_batch_V1_15_complete_package/README.md)
- [完整自检报告](sem_cd_measure_200k_batch_V1_15_complete_package/validation_report.json)
- [独立解压和完整导出验证](audits/v1_15_release_evidence.json)

## 交付内容

1. 背景排除后，使用前景列 Y 跨度均值粗估平均 trench 长度，默认上下各排除 1/20；主测量和 PSD 都遵守裁后区域。
2. `--locator-mode dark-line` 保留旧定位；`bright-line` 使用整体亮 line / 整体暗 trench 的两类自适应阈值、多数像素、均匀性及相邻周期规则。
3. 默认 `auto` 分析整个有效观测区域的 X 灰度周期，判断交替深浅/宽窄暗谷特征，每图选择一次并用于两个引擎及 PSD。新增模式、阈值、候选和歧义诊断导出。
4. 两种定位共用后续边缘处理；原定位器、修正、统计和 PSD 模块经 AST 对比确认除版本文字外内容相同。完整参数手册也包含在 ZIP 内。

## 验证结果

- 完整自检通过：22 个端到端/故障场景，加上数学、区域、自动模式和新参数检查。
- 亮 line / 带假暗条 line 两类合成图均选择正确模式，auto 与对应手动模式得到相同边缘；均定位到原图约 X=380、540、700px 的三条真实 trench。
- 平均长度 440px 的合成图每端排除 22px；强制连续性恢复和 PSD 插值仍不能填入端部或内部背景。
- 旧配方数值回归通过；保留示例每引擎 512/512 主采样有效，V10/V13 各生成 32 条独立 PSD 曲线。
- 独立解压 ZIP 后快速自检通过；完整标注、统计图、机台对照、双引擎 PSD 导出通过，0 条处理异常，48 条 PSD 曲线及 36 张导出图。
- 全部 86 个参数覆盖检查、16 条文档示例解析、Python 编译及 Ruff E9/F 检查通过。

所有图片验证均为随包或生成的合成图，尚无用户真实 SEM 混合背景样本。完全相同的真假灰度周期不能唯一判断物理身份，必要时显式指定模式并核查标注。软件回归通过不等于 SEM 计量精度认证。

## 默认运行

```bash
python sem_cd_measure_200k_batch_V1_15.py --root /path/to/images --pattern trench --locator-mode auto --pixel-size 1.3181 --trench-reference-nm 60 --end-trim-fraction 0.05 --output /path/to/new_results
```

请替换实际标定、宽度和路径，输出目录须不存在或为空。
