# V1_19 修改记录

## 2026-09-23：多格式图片输入

基于 V1_18（1fcde1c）创建独立 V1_19 完整包。历史版本文件保持原样。

- 递归扫描 TIF/TIFF、PNG、JPG/JPEG、BMP、WebP，扩展名不区分大小写，可混合输入，支持中文路径。
- 使用现有 OpenCV 解码器；支持 8/16 位灰度和彩色图，多页 TIFF 读取第一页。无需新增依赖。灰度转换、归一化和全部测量数学沿用 V1_18。
- ROI、双引擎标注、旋转诊断和 PSD 图统一为 PNG，输出名保留完整输入扩展名，避免同名不同格式覆盖。
- CLI 提示及图表使用通用图片描述。settings.json 使用 input_image_count、processed_image_count 并记录 supported_image_suffixes。
- 机台未匹配图片导出改为 machine_unmatched_images.csv；保留原有同名歧义排除规则。
- 增加真实格式编解码、16 位/彩色/多页、中文路径、混合批处理、导出防覆盖、非 PNG 机台匹配与损坏文件回归；保留 V1_18 原有自检。

完整验证结果见 validation_report.json；发布与旧版对照见仓库 audits/v1_19_release_evidence.json。JPEG 压缩可能改变测量结果；未使用用户真实 SEM 图验证。Windows PowerShell 包装器未在本次 Linux 环境执行。

验证覆盖 63 个集成场景：首次全量运行已完成全部 60 个原有场景；新增格式检查因测试帮助函数关闭 PSD 绘图而首次断言失败，修正测试开关后，解析检查及 3 个格式专项均通过。随包 JSON 合并两次测试记录，保留执行说明和首轮日志，不将专项重跑表述为第二次全量执行。
