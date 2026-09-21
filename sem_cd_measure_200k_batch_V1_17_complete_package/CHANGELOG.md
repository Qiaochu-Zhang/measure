# V1_17 修改记录

## 2026-09-21：相邻 trench 共享边缘测 line

基于 V1_16（`c0ae220`）创建独立 V1_17 目录、主入口、启动器、依赖清单、文档、自检报告与 ZIP。V1_15、V1_16 保持原样。

- 删除 line 的反相测量路径。ROI、自动参考宽度、定位模式分析均基于原图物理 trench。
- 复用 trench 定位器，搜索 ROI 中全部合规候选并用原 V10/V13 流程测量；来源 trench 不受 max-number 提前截断。
- 仅相邻且均达到原稳定性条件的 trench 可配对。line 左边逐点复制左 trench 右边，右边复制右 trench 左边，采样 Y 必须相同；不再对 line 独立调整边缘。
- 不跨越未合格或漏检的物理 trench；单个 trench 无法构成 line。背景、交叉边和缺失采样无效；合成补点保留标记和原质量限制。
- line 构建完毕后按原中心/左右选择最多 max-number 条。CD/LER/LWR、旋转、分组、mixed、单行 PSD 继续使用原规则。Viterbi/ERF、飞点和 continuity 在来源 trench 上生效。
- line 的 pitch 复用全部来源 trench 同侧边缘。仍只输出旋转 pitch CD，按 max-number 个不同完整周期等权平均，Excel 主表最后一列仍为 `旋转_pitch_CD_nm`。
- 新增 `line_trench_pairs`、`line_source_samples` 表；主采样表记录左右来源 trench ID 和共享坐标。保留全部原 CLI 参数，说明非对称阈值的共享边对应关系。
- 新增共享坐标逐点验证、宽暗条、自动估计、非对称阈值、Viterbi/ERF、max-number=1/2/5、单 trench/双 trench、不合格/缺失中间 trench、单行 PSD 测试；保留既有 trench、pitch、背景、导出和故障回归。

完整自检通过 48 个独立端到端/故障场景，其中 12 个共享边专项场景，另含解析检查。完整结果见 `validation_report.json`；仓库 `audits/v1_17_release_evidence.json` 记录旧版本文件保留、trench 数值对照及独立 ZIP 解压验证。`audits/v1_17_shared_line_baseline.json` 记录同一宽暗条合成图的 V1_16 trench、V1_16 line、V1_17 line 对照。

内部 `v115` / `V1.15` 是继承定位算法的名称；本版输出和版本标识均为 V1_17。验证基于合成图与随包样例，尚无用户实际失败 SEM 图可供对照。
