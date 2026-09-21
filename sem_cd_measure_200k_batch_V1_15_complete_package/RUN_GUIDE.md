# V1_15 完整运行说明与参数手册

适用补丁：`V1_15_line_pitch_1`（2026-09-21）。主程序：`sem_cd_measure_200k_batch_V1_15.py`，基于 V1_14 改进，按本次 V1_15 发布的实际参数解析、校验和执行代码编写。本文覆盖主程序全部 **86 个参数选项（包括帮助选项，`-h` 和 `--help` 为同一个选项）**，另列自检脚本参数。没有在命令行开放的内部常量，不应当作可传入参数。

- [完整代码 ZIP](../sem_cd_measure_200k_batch_V1_15_complete_package.zip)
- [版本概述](CHANGELOG.md)
- [包内 README](README.md)
- [主程序源码](sem_cd_measure_200k_batch_V1_15.py)

## 1. 先了解默认行为

程序用于 PNG 图像中近竖直的暗沟槽 `trench` 或亮线 `line` 的 CD、LER、LWR 测量。一次运行的整个输入目录只能使用一种图案类型；不同类型请分开运行。本版不支持 via，也不自动逐图判断 trench/line。

默认会：

1. 递归读取输入目录及子目录的 PNG，不按文件名或 `200K` 字样筛选；当前输出目录和可识别的历史输出子目录会被排除。
2. 筛除背景，按估计平均 trench 长度从 Y 两端各去掉 5%，自动确定未手动指定的 ROI；每图最多选 3 条结构。根据整个有效区域的 X 灰度周期逐图选择 dark-line 或 bright-line 定位器。
3. 分别运行 V10、V13 边缘算法。总结表中每张图按 `mixed → V13 → V10` 排列。
4. 主测量使用 128 个采样位置、32 行灰度平均、3 像素横向平滑；CD/LER 不分组，LWR 使用 group4。
5. PSD 另外用单行、1 px 沿线步长、无横向平滑重新提取 V10/V13 边缘；逐结构、逐连续段计算，不做 group4、Welch 或跨结构/跨图片频谱平均。
6. 无论 line/trench 都自动测 pitch CD：按 max-number 个相邻完整周期求均值，只输出旋转结果，放在主汇总表最后一列。原 CD/LER/LWR 仍输出旋转与未旋转两套结果、Excel/CSV、标注图和 PSD；若输入根目录存在 `Data.xlsx`，尝试与其机台数据对比。

参考宽度只是检测先验，不是测量真值。像素尺寸必须来自你的图像标定，不能仅凭“200K”倍率推断。

## 2. 安装与第一次运行

### 2.1 安装

解压完整 ZIP，进入包含主程序和 `requirements_V1_15.txt` 的目录。建议 Python 3.12；代码使用 Python 3.10+ 语法，实际测试版本见包内 `validation_report.json`。

```bash
python --version
python -m pip install -r requirements_V1_15.txt
python self_check_V1_15.py --quick
```

如果系统命令叫 `python3`，将命令中的 `python` 替换成 `python3`。建议使用独立虚拟环境安装依赖。`requirements-tested.txt` 记录已测的直接依赖版本，不是完整的跨平台依赖锁定文件。

### 2.2 跑随包示例

```bash
python sem_cd_measure_200k_batch_V1_15.py --root examples/input_trench --pattern trench --pixel-size 1 --trench-reference-nm 60 --max-number 4 --output demo_output --no-auto-machine-comparison
```

`demo_output` 必须不存在或为空。重复运行请换一个输出目录，例如 `demo_output_02`；程序不会覆盖已有结果。

### 2.3 跑自己的图

下面的路径、`1.3181 nm/px` 和 `60 nm` 都只是示例，必须按实际数据替换。命令写成单行，避免 Windows 与 Linux 的续行符差异。

暗沟槽：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root "D:/SEM/input_trench" --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --output "D:/SEM/results_trench_01" --no-auto-machine-comparison
```

亮线：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root "D:/SEM/input_line" --pattern line --pixel-size 1.3181 --line-reference-nm 60 --output "D:/SEM/results_line_01" --no-auto-machine-comparison
```

Linux/macOS 把路径替换为 `/path/to/images` 等实际路径。相对路径基于终端的当前工作目录；含空格的路径加引号。Windows 也可使用 `./RUN_V1_15.ps1`，Linux/macOS 可用 `./run_V1_15.sh` 替代 `python sem_cd_measure_200k_batch_V1_15.py`，参数相同。PowerShell 包装脚本未在本次 Linux 环境执行验证。

输入可以包含背景，程序自动筛选目标区域并排除端部。保留原始像素和正确标定；不要将 trench、line 和历史输出图混放。彩色图会转灰度，8bit 灰度保持原强度；非 8bit 图像会归一化，需要额外检查低动态范围图像。

## 3. 参数阅读规则

- “开关”只写参数名，例如 `--no-auto-roi`；不加 `true`、`false` 或 `1`。不写时保持默认行为。
- `--viterbi`、`--erf-fit` 是需要数值的选项，分别用 `0`/`1`，不是上述开关。
- `px` 表示像素，`nm` 表示纳米；灰度对比度按处理后的 8bit 灰度尺度解释，不是百分比。
- `None/自动` 表示没有固定 CLI 默认值，会按图片或其他参数计算；运行输出的 `settings.json` 记录实际配方。
- 绝大多数边缘参数同时用于 V10 和 V13，不能通过 CLI 分别给两个引擎指定不同值。
- 以下表格记录的是主程序 CLI 默认值，不是内部 `MeasurementParams` 类单独实例化时的默认值。

## 4. 主程序全部参数

### 4.1 输入、输出与运行控制（8 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `-h` / `--help` | 开关 | 显示内置帮助并退出；不需要同时指定 pattern。部分内置帮助较简略，以本手册中的实际执行说明为补充。 |
| `--root` | 路径；`.` | 输入根目录，递归扫描 PNG；必须是存在的文件夹，不是单张图片路径。 |
| `--output` | 路径；自动 | 默认 `root/CD_measure_output_200K_V1_15`。目录必须不存在或为空。不要直接指定到包含输入图的目录。 |
| `--pattern` | **必填**；`trench` / `line` | `trench` 测暗沟槽，`line` 测亮线；对 root 下所有图片统一生效。支持大小写归一化。 |
| `--pixel-size` | 浮点数；`1.3181` nm/px | 像素标定，必须 >0。影响 CD/LER/LWR、宽度先验换算和 PSD 的频率/密度单位。 |
| `--stop-on-error` | 开关；默认不启用 | 逐图主处理发生致命异常或产生 ERROR 时停止后续处理，并先导出已处理结果。**不是对所有错误的全局即时停止**：预扫描读取错误、REVIEW、背景跳过或批次附属导出错误不一定触发停止，详见第 9 节。 |
| `--check-parameters` | 开关；默认不启用 | 解析并校验 CLI 参数，打印 JSON 后退出，不测量图片。仍须指定 pattern；不检查图片实际尺寸、输入目录是否存在、输出目录是否非空或机台表内容。 |
| `--save-debug-masks` | 开关；默认不启用 | 在 `debug/` 额外保存双引擎旋转坐标下的边缘偏差和 CD 诊断曲线。参数沿用历史名称，当前并不是导出二值掩码的开关，也不控制正常 ROI 预览。 |

### 4.2 宽度先验、候选结构及数量（10 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--locator-mode` | `auto` / `dark-line` / `bright-line`；`auto` | trench 粗定位模式。dark-line 保留 V1_14 的盆地定位及跳过 line 内假暗条的周期规则；bright-line 用两个自适应灰度类识别整体亮 line 和整体暗 trench，相邻真实暗区为相邻周期；auto 从整个观测区域的 X 灰度轮廓逐图判断。支持大小写。与测量目标 pattern、边缘引擎 V10/V13 是不同选项，详见第 11 节。 |
| `--locator-majority` | 浮点数；`0.70` | bright-line 候选核心低于自适应阈值、两侧 line 核心高于该阈值的最小像素比例，范围 `(0.5,1]`。同时使用内部均匀性检查。增大更严格；dark-line 原算法不使用此候选门槛。auto 分析也用该值检查有效列支持比例。不是绝对灰度阈值，也不是边缘交点百分比。 |
| `--trench-reference-nm` | 浮点数；自动估计 | 暗沟槽参考宽度，必须 >0。不填时从图像灰度轮廓估计，复杂背景可能使估计不准。line 模式也接受此兼容参数，但数值必须是亮线的参考宽度，不是原 trench 的宽度。 |
| `--line-reference-nm` | 浮点数；自动估计 | 亮线参考宽度，必须 >0，只能用于 line。若与 trench-reference-nm 同时填写，两者必须相同。 |
| `--space-reference-nm` | 浮点数；不显式约束 | 目标之间的间隙宽度，必须 >0，**不是中心距**。显式设置后，候选邻居排序会参考“目标宽度 + 间隙宽度”的节距；不设置时估计值可写入元数据，但不启用此显式排序先验。 |
| `--max-number` | 整数；`3` | 每张图每个引擎最多选取的结构数量，≥1；同时为 pitch 求平均的目标周期数。N 个完整 pitch 需要 N+1 条同类结构的边缘，额外结构仅用于 pitch。不是图片数量，也不是采样点数。 |
| `--min-number` | 整数；`min(3,max-number)` | 质量判定要求的最少稳定结构数，≥1 且不能超过 max-number。不足时仍尽可能保留已有结果，但通常标记 REVIEW。 |
| `--allow-incomplete-triplet` | 开关；默认不启用 | 取消“中心结构及左右邻居齐全”的额外质量要求。不会取消 min-number 要求，也不会凭空增加候选；max-number <3 时原本就不要求完整三条。 |
| `--candidate-dark-tolerance` | 浮点数；`None` | 可选的候选核心灰度容差，≥0：以中心候选核心灰度为参考，过滤比其更亮且超出容差的候选。line 使用反相后的工作灰度，因此不是按原图直接找暗线。不是左右边缘阈值。 |
| `--candidate-min-contrast` | 浮点数；`None` | 可选的候选沟槽/线条最低对比度，≥0；过滤粗定位候选。未指定表示不额外施加此门槛，并不关闭定位器原有检查。 |

只测一条目标可使用 `--max-number 1 --min-number 1`。希望测 4 条但允许只检测到 2 条，可设置 `--max-number 4 --min-number 2 --allow-incomplete-triplet`；这会改变质量准入条件，不能拿来证明图像质量变好了。

### 4.3 自动背景排除与测量区域（9 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--no-auto-roi` | 开关；默认不启用 | 关闭自动背景识别，采用手动矩形或整图作为观测区域；默认仍裁去端部。若要复现完全不使用掩码的旧行为，还需 `--end-trim-fraction 0`。这时手动 ROI / 图像中心默认值沿用 V1_14。 |
| `--end-trim-fraction` | 浮点数；`0.05`（1/20） | 每个 Y 端部各排除“估计平均长度 × 此比例”的像素数，向上取整；总计约排除 10%。范围 `[0,0.5)`；0 只关闭端部裁剪。平均长度只用前景列 Y 跨度均值一种近似方法，见第 11 节。主测量和 PSD 共用裁后掩码，补点不能跨过端部。 |
| `--roi-min-contrast` | 浮点数；`6.0` 灰度级 | 自动背景筛选的最低对比度基准，必须 >0；还会结合局部噪声门槛，并非最终边缘检测的唯一阈值。太大可能漏掉弱结构，太小可能保留背景。关闭自动 ROI 后不参与背景筛选。 |
| `--roi-min-length` | 整数；`24` px | 自动区域保留的最低连通区域高度，至少 4 px；不是最终 PSD 段长。短结构可能因此被排除。关闭自动 ROI 后不参与筛选。 |
| `--center-x` | 浮点数；自动 | 原图坐标中的测量中心 X，向右增大。自动 ROI 开启且未指定时用保留区域包围框中心；无掩码（no-auto-roi 且端部比例为 0）时用 `(图宽−1)/2`。必须在图内。 |
| `--center-y` | 浮点数；自动 | 原图坐标中的测量中心 Y，向下增大。自动 ROI 开启且未指定时用保留区域包围框中心；无掩码时用 `(图高−1)/2`。必须在图内。 |
| `--meas-area-width` | 整数；自动 | 粗定位 ROI 宽度，px，必须 >0 且不超过图宽。自动 ROI 时按保留区域范围及中心到图边的距离调整；无掩码时默认 `min(512,图宽)`。不是最终测出的线宽。 |
| `--meas-area-height` | 整数；自动 | 粗定位 ROI 高度，px，必须 >0 且不超过图高。自动 ROI 时按保留区域范围调整；无掩码时默认 `min(512,图高)`。不是最终沿线采样长度。 |
| `--extend-length` | 整数；自动 | 沿 Y 方向的测量长度，必须在 2 到图像高度之间。自动/手动区域裁剪开启且未显式设置时为 `max(2,min(360,保留区域包围框高度−average-range))`；无掩码时为 `min(360,图高)`。影响主测量和 PSD 测量范围。 |

显式填写的中心、区域尺寸和长度优先于自动值，但**手动 ROI 不会关闭背景掩码**。若想完全关闭自动背景和端部排除，需同时加入 `--no-auto-roi --end-trim-fraction 0`。坐标从左上角 `(0,0)` 开始，始终属于原图；程序不会先裁剪再重新编号。

自动 ROI 是几何和灰度规则，不是能理解图像含义的语义分割。与目标同样呈平行条纹的背景仍可能通过；弱对比度、短线、强弯曲和特殊缺陷可能被排除，必须看 `ROI/` 图。

### 4.4 采样、阈值与边缘搜索（10 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--sample-number` | 整数；`128` | 每条结构主测量的采样位置数，至少 2。不是平均行数。默认 single-row PSD 使用另一套每行采样，不受此数量直接控制。 |
| `--average-range` | 整数；`32` px | 主测量每个采样位置附近沿 Y 平均的图像行数，≥1 且不能超过图高。增大可减小噪声，也可能抑制真实高频粗糙度。默认 single-row PSD 强制为 1。 |
| `--smoothing-pixel` | 整数；`3` px | 主测量横向灰度轮廓的平滑窗口，≥1 且不超过图宽。1 表示无横向窗口平均；通常使用奇数。默认 single-row PSD 强制为 1。 |
| `--search-in` | 整数；`33` px | 非局部跟踪搜索时，从预期边缘向结构内部搜索的范围，≥0。不是从图像中心起算。 |
| `--search-out` | 整数；`30` px | 非局部跟踪搜索时，从预期边缘向结构外部搜索的范围，≥0。局部跟踪/重试主要使用 tracking-radius 等半径，不能只增大此项解决所有漏检。 |
| `--peak-order` | 整数；`1` | 在通过极性、拓扑等检查的梯度峰中，按靠近预期边缘的次序取第几个，≥1；不是按峰高排第几。指定较高次序但候选不足时会失败，不静默退回第一峰。 |
| `--threshold-left` | 浮点数；`50.0` % | 左边缘灰度交点百分比，严格在 0 与 100 之间。阈值为 `暗参考 + 百分比×(亮参考−暗参考)`；不是绝对灰度 50。 |
| `--threshold-right` | 浮点数；`50.0` % | 右边缘灰度交点百分比，严格在 0 与 100 之间。可与左侧不同，但必须明确这会改变 CD 的定义。 |
| `--threshold-search` | `bounded` / `legacy`；`bounded` | bounded 在受限搜索窗口内找正确极性的阈值交点，修复交点在梯度峰外侧导致的漏检；legacy 保留旧搜索逻辑，仅用于复核旧结果。 |
| `--topology-min-contrast` | 浮点数；`2.5` 灰度级 | 边缘内外侧必须满足的最低明暗拓扑对比度，≥0；太高会漏检，太低可能把内部纹理当作边缘。与 ROI 和候选对比度是不同层次的检查。 |

V10 通常使用中心暗参考及峰附近参考，V13 使用侧带均值等参考；当中心受内部条带污染、两侧内部平台一致且更暗时，改用两侧内部平台的中位数均值作为暗参考，具体触发常数见第 12 节。即使同为 50%，交点也不一定相同。line 模式在反相后的工作图上执行这些暗目标规则。

主采样槽位间距约为 `(extend-length−1)/(sample-number−1)` px。增加 sample-number 并不增加图像真实分辨率，32 行平均窗口还可能高度重叠。single-row PSD 改用连续整数行，但自动 extend-length 的推导仍会受到上面的 average-range 设置影响。

### 4.5 跟踪、后过滤与缺点恢复（9 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--tracking-radius` | 浮点数；`13.0` px | 围绕预测/跟踪边缘的初次搜索半径，必须 >0。 |
| `--tracking-retry-radius` | 浮点数；`20.0` px | 初次失败后的重试及恢复搜索半径，必须 >0；放大可能增加找回机会，也可能误锁邻近结构。 |
| `--tracking-max-jump` | 浮点数；`8.0` px | 严格跟踪时单侧边缘相对预测锚点允许的最大跳变，必须 >0；另有内部宽度和中心偏移检查，不是唯一限制。 |
| `--postfilter-window` | 整数；`9` 个采样位置 | 稳健后过滤的局部参考窗口，至少 3；不是 9 nm，也不是灰度平均行数。 |
| `--postfilter-mad-multiplier` | 浮点数；`6.0` | MAD（中位绝对偏差）稳健尺度的容差倍数，≥0。越大通常越宽松；与最低像素容差共同起作用。 |
| `--postfilter-edge-tolerance` | 浮点数；`2.5` px | 后过滤的最低边缘偏差容差，必须 >0；另有内部 CD 容差。即使 flyer-mode 为 reserve，明显几何错误仍可能先被后过滤剔除。 |
| `--minimum-valid-fraction` | 浮点数；`0.70` | 稳定结构要求的真实有效采样比例，在 `(0,1]`。以有资格测量的非背景位置为分母，不把强制合成点当作真实检出点；同时还要求足够的统计点数。 |
| `--no-recovery` | 开关；默认不启用 | 关闭稳健过滤后的常规缺点重测。不会关闭初次跟踪重试，也不会自动关闭另外配置的 edge-continuity 或 Viterbi。 |
| `--edge-continuity` | 整数；`0` | 范围 0–100。0 不额外放宽/强制补点；1–99 逐级放宽并重测缺点；100 还尝试补齐剩余可处理缺点并标记合成来源。背景永不补齐。不能将补点后的完整轨迹当作真实检测率。 |

这些选项在 single-row PSD 的独立重测中也会沿用；“PSD 不平均”不等于关闭所有拓扑约束、离群检查和拟合。默认 `edge-continuity=0`，但普通 recovery 仍然开启。

### 4.6 Flyer、LWR 分组与汇总质量（5 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--flyer-mode` | `reserve` / `remove` / `replace`；`reserve` | reserve 保留通过前置检查的 flyer 并标记；remove 删除被判定为 flyer 的左右成对采样；replace 用局部参考替换异常侧，记录为合成点，并保留几何/路径约束。 |
| `--flyer-left-offset` | 浮点数；`20.0` px | 相对局部参考**向左、负方向**偏移的 flyer 门槛，≥0；同时用于左右两条边缘，不是“左边缘专用参数”。 |
| `--flyer-right-offset` | 浮点数；`20.0` px | 相对局部参考**向右、正方向**偏移的 flyer 门槛，≥0；同样作用于两条边缘。 |
| `--group-size` | 整数；`4` | 主表 LWR 的沿线分组点数，≥1 且不能超过 sample-number。1 表示不分组；CD/LER 不受此分组控制；**两套 PSD 都不使用该分组值**。有效组不足 2 时 LWR 输出空值，不输出假零值。 |
| `--include-review-in-summary` | 开关；默认不启用 | 默认条件汇总和机台对比只使用 OK；加上后也纳入 REVIEW。逐图表原本就保留 REVIEW。此参数不把 REVIEW 改成 OK。 |

Flyer 判断使用局部参考，不是根据“是否偏离某个设计 CD”直接判断。remove/replace 和 group-size 都可能改变粗糙度，比较批次时必须保持配方一致。

### 4.7 Viterbi 全局路径与 ERF 局部拟合（9 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--viterbi` | 整数 `0` / `1`；`0` | 1 启用整条左右边缘对的全局路径选择；用于在多个候选间兼顾沿线连续性。不是保证边缘真实的开关。 |
| `--viterbi-weight` | 浮点数；`5.0` | 路径跳动代价权重，≥0。显式指定必须同时使用 `--viterbi 1`。增大可能更偏好连续路径，也可能抑制真实变化。 |
| `--viterbi-max-jump` | 浮点数；`3.0` px/采样间隔 | 每个预设采样间隔允许的最大横向位移，必须 >0；跨缺点按槽位间距计算。显式指定需 `--viterbi 1`。主测量和单行 PSD 的采样间隔不同，应注意其物理含义。 |
| `--viterbi-gap-cost` | 浮点数；`6.0` | 跳过采样槽位的代价，必须 >0；不是缺口长度上限。显式指定需 `--viterbi 1`。 |
| `--viterbi-candidates` | 整数；`5` | 每个位置最多保留的左右边缘对候选数，≥1；增大会增加搜索开销。显式指定需 `--viterbi 1`。 |
| `--erf-fit` | 整数 `0` / `1`；`0` | 1 启用局部 ERF（误差函数）灰度过渡拟合，细化边缘位置；拟合不合格时保留诊断，不强行接受。 |
| `--erf-window` | 整数；`8` px | ERF 拟合窗口半宽，至少 3 px。显式指定需 `--erf-fit 1`。 |
| `--erf-max-shift` | 浮点数；`2.0` px | 拟合交点相对输入边缘允许的最大位移，必须 >0。显式指定需 `--erf-fit 1`。 |
| `--erf-max-relative-rmse` | 浮点数；`0.2` | 拟合 RMSE / 拟合对比度的上限，必须 >0；0.2 表示比值 0.2，不是绝对灰度误差。显式指定需 `--erf-fit 1`。 |

### 4.8 PSD（13 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--skip-psd` | 开关；默认不启用 | 跳过全部 PSD 计算/输出，主 CD/LER/LWR 照常；不能与任何显式 `--psd-*` 参数一起使用。 |
| `--psd-edge-source` | `single-row` / `measurement`；`single-row` | single-row 分别重提 V10/V13 边缘，沿线每行采样、average-range=1、smoothing-pixel=1；measurement 复用主测量边缘，会继承主配方已经存在的灰度平均/平滑。两者都不做后续 groupN 或频谱平均。 |
| `--psd-axis` | `normal` / `unrotated` / `both`；`both` | normal 为中心线 PCA 旋转后的法向边缘坐标，unrotated 为原图坐标，both 输出两套。只控制 PSD，不改变主表始终输出两套坐标的行为。 |
| `--psd-method` | 仅 `periodogram`；`periodogram` | 每条结构的每个足够长连续段独立计算周期图；本版不能选择 Welch。 |
| `--psd-window` | `hann` / `boxcar` / `blackmanharris`；`hann` | 每段 PSD 的窗函数：Hann；矩形窗；Blackman–Harris 窗。会影响谱泄漏与分辨率，不是分组平均，也不是噪声去偏。 |
| `--psd-detrend` | `linear` / `constant` / `none`；`linear` | 每段计算前去线性趋势、仅去均值或不去趋势。去趋势不是边缘分组平均。none 可能让绝对位置/平均宽度进入直流及低频能量，通常不用于直接解释粗糙度。 |
| `--psd-min-segment` | 整数；`16` 个位置 | 可计算 PSD 的最短连续有效段点数，至少 4。更短的段跳过；不是所有段截为 16 点，也不是 Welch 窗长。 |
| `--psd-gap-mode` | `segments` / `interpolate`；`segments` | segments 在缺点处分段，绝不直接拼接后 FFT；interpolate 允许在两端有效的短内部缺口线性插值，再分段计算。背景缺口禁止插值。 |
| `--psd-max-gap` | 整数；`2` 个缺失位置 | interpolate 模式允许插值的最大连续缺点数，≥0；显式指定必须同时使用 `--psd-gap-mode interpolate`。不是 nm 或频谱平均窗口。 |
| `--psd-include-synthetic` | 开关；默认不启用 | 允许边缘阶段 continuity/flyer 等合成点进入 PSD；默认剔除。不会主动生成合成点，也不会允许背景进入 PSD。PSD 自己的短缺口插值由 gap-mode 另行控制。 |
| `--psd-split-wavelength` | 浮点数；`100.0` nm | 低/高频方差汇总的分界波长，必须 >0；分界频率为 `1/波长`。默认低频包含 `f≤0.01 nm⁻¹`。只划分积分区间，不滤掉任何频率。 |
| `--psd-nperseg` | 遗留解析项；`64` | **本版不要传入**。仅为旧 Welch 接口保留的解析默认值，periodogram 不使用；即使显式填写默认 64，也会因只适用于 Welch 而报错。 |
| `--psd-overlap` | 遗留解析项；`0.5` | **本版不要传入**。旧 Welch 重叠比例，本版不使用；即使显式填写 0.5，也会被校验拒绝。 |

PSD 每个引擎输出 `left`、`right`、`width`、`center` 四类信号；center 是左右边缘的几何中点，不是沿线平均后的边缘。频率单位 `1/nm`，PSD 密度单位 `nm³`。不输出 mixed PSD，因为 mixed 不是一套统一来源的边缘。

默认不做噪声去偏，因此原始 PSD 比平均谱更起伏、单行边缘更容易受噪声影响。改变平均、采样、去趋势或分组定义后，不应要求 PSD 积分粗糙度等于主表 group4 LWR。

### 4.9 机台参考表对比（9 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--machine-excel` | 路径；自动查找 | 指定机台参考 Excel。未指定且允许自动对比时，只检查 `root/Data.xlsx`，不会递归搜索所有 Excel。没有找到自动参考表，不影响普通测量。 |
| `--no-auto-machine-comparison` | 开关；默认不启用 | **只关闭自动查找 root/Data.xlsx**，不关闭测量、PSD 或正常结果 Excel。若同时显式指定 machine-excel，仍会进行对比。 |
| `--machine-sheet` | 字符串；`v2` | 机台表工作表名称。传 `0` 会被当作名称字符串，不是“第一个工作表”索引。 |
| `--machine-start-row` | 整数；`2` | 起始 Excel 行号，按 1 开始计数，包含该行，至少 1。默认跳过第 1 行。 |
| `--machine-end-row` | 整数；`84` | 结束 Excel 行号，包含该行；0 表示到表尾，否则不能小于起始行。超过实际行数时以实际表尾为止。 |
| `--machine-image-column` | Excel 列字母；`B` | 图片名称列，例如 B 或 AA，不填数字列号。建议单元格使用图片文件名。 |
| `--machine-cd-column` | Excel 列字母；`D` | 参考 CD 数值列，单位应为 nm。 |
| `--machine-ler-left-column` | Excel 列字母；`E` | 参考左 LER 数值列，单位应为 nm；应与程序 3σ 定义一致。当前没有右 LER 机台列参数。 |
| `--machine-lwr-column` | Excel 列字母；`F` | 参考 LWR 数值列，单位应为 nm；需核对 3σ、分组、滤波等统计定义是否可比。程序不会自动统一这些定义。 |

对比使用主程序的 mixed 指标（旋转、未旋转分别比较），默认只纳入 OK。匹配键基于文件名去扩展名、忽略大小写，并去掉末尾类似 `(1)` 的数字后缀；不是按完整相对路径匹配。不同子目录存在同名图片可能被判为歧义并排除匹配，应检查 `machine_unmatched_png.csv`。未匹配参考行见 `machine_unmatched_rows.csv`。

### 4.10 图像与表格导出控制（4 项）

| 参数 | 类型 / 默认值 | 含义与注意事项 |
|---|---|---|
| `--no-annotated-images` | 开关；默认不启用 | 不保存双引擎边缘标注图；不关闭 ROI 预览、显式 debug 图或统计/PSD 图。旋转诊断图由 save-debug-masks 独立控制。 |
| `--compact-sample-output` | 开关；默认不启用 | 主逐点结果只保留核心列，省略部分阈值、恢复、Viterbi/ERF 等诊断列；不减少采样数量、不改变算法，也不是关闭独立 PSD 边缘审计表。 |
| `--skip-statistics-plots` | 开关；默认不启用 | 跳过内部趋势图、机台对比统计图及 PSD 图；保留对应数值表。不会关闭边缘标注、旋转诊断、ROI 或显式 debug 图。 |
| `--plot-dpi` | 整数；`180` | 统计和 PSD 图的导出 DPI，必须 >0；不是输入像素标定，不改变测量值。部分标注/旋转诊断使用原图尺寸或自身固定 DPI，不受此项统一控制。 |

## 5. 常用运行配方

以下命令都在代码包目录执行，`input`、`results_*` 是示例路径；请按实际图片填写像素尺寸和参考宽度。每次输出必须是新目录。

### 5.1 只测量图片，不自动对比机台

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --output results_basic --no-auto-machine-comparison
```

仍会输出正常 Excel 和两套 PSD。不要同时填写 machine-excel，否则仍会触发显式对比。

### 5.2 对比自己的机台 Excel

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --output results_machine --machine-excel reference.xlsx --machine-sheet Sheet1 --machine-start-row 2 --machine-end-row 0 --machine-image-column A --machine-cd-column B --machine-ler-left-column C --machine-lwr-column D
```

列位置仅为示例。参考 Excel 不需要放在输入目录里。

### 5.3 手动 ROI（仍裁去端部）

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --no-auto-roi --center-x 300 --center-y 250 --meas-area-width 400 --meas-area-height 300 --extend-length 240 --output results_manual --no-auto-machine-comparison
```

示例中心和尺寸必须适合你的原图。本命令仍从手动区域上下各去掉 5%；若需关闭端部排除，再加 `--end-trim-fraction 0`。只填写 center-x/center-y 不等于关闭背景排除。

### 5.4 只测一条线，并输出未分组主 LWR

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern line --pixel-size 1.3181 --line-reference-nm 60 --max-number 1 --min-number 1 --group-size 1 --output results_single --no-auto-machine-comparison
```

主边缘仍使用默认 32 行平均、3 px 横向平滑；group-size=1 只取消后续 LWR 分组。如需主测量也取消这两种灰度平均，另加 `--average-range 1 --smoothing-pixel 1`，并重新评估噪声影响。

### 5.5 默认单行 PSD，并允许最多 2 个位置的内部缺口插值

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --psd-edge-source single-row --psd-gap-mode interpolate --psd-max-gap 2 --output results_psd_interpolate --no-auto-machine-comparison
```

这是对缺失数据的显式处理，不是实测点；背景缺口不会插值。若不希望任何 PSD 缺口插值，保持默认 segments。

### 5.6 PSD 复用主测量边缘，不进行后续分组平均

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --psd-edge-source measurement --output results_psd_reuse --no-auto-machine-comparison
```

此配方继承主测量的 32 行灰度平均和横向平滑，与默认 single-row 的谱不应直接混作同一配方比较。

### 5.7 启用路径优化与 ERF 拟合

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --viterbi 1 --viterbi-weight 5 --viterbi-max-jump 3 --erf-fit 1 --erf-window 8 --save-debug-masks --output results_path_erf --no-auto-machine-comparison
```

这是可选实验配方，不保证比默认结果更准确。先对照原图、轨迹和已知标准样，再决定是否批量使用。

### 5.8 先只看主指标、减少图表输出

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --pixel-size 1.3181 --trench-reference-nm 60 --skip-psd --skip-statistics-plots --no-annotated-images --compact-sample-output --output results_fast --no-auto-machine-comparison
```

这不是完全不输出图片：自动 ROI 开启时仍保留 ROI 预览。不要在该命令中继续添加 psd-axis 等 PSD 参数。

### 5.9 只检查参数，不跑图片

```bash
python sem_cd_measure_200k_batch_V1_15.py --pattern trench --max-number 4 --group-size 8 --check-parameters
```

用于检查拼写、取值和参数冲突；不能证明图片/参考文件可读或 ROI 尺寸正确。

## 6. 参数之间的限制和优先级

1. 必须指定 pattern；line-reference-nm 只能用于 line，且不能与不同值的 trench-reference-nm 同时使用。
2. max-number ≥ min-number ≥1；sample-number ≥2；1≤group-size≤sample-number。数学上允许的参数，不代表一定有足够有效数据。
3. 阈值百分比严格在 `(0,100)`；minimum-valid-fraction 在 `(0,1]`；edge-continuity 在 `[0,100]`。
4. 显式设置任何 Viterbi 子参数都必须加 `--viterbi 1`；显式设置 ERF 子参数都必须加 `--erf-fit 1`，即使填写的恰好是默认值也一样。
5. 使用 skip-psd 时不能再显式填写任何 psd-* 参数；max-gap 只能显式搭配 interpolate；nperseg、overlap 在当前 periodogram 版本中不能显式传入。
6. 手动中心/尺寸优先，但不会关闭自动背景掩码；no-auto-roi 关闭背景识别，end-trim-fraction=0 关闭端部裁剪，两者独立。auto 在全部裁后观测区域判断模式，不限于默认 360px 测量长度。
7. machine-excel 的显式路径优先于自动查找，也优先于 no-auto-machine-comparison 的“关闭自动”行为。
8. CLI 检查通过后，实际逐图测量仍会检查中心、ROI、平均窗口和长度是否适合该图。不同尺寸图片共用固定手动参数时尤其要注意。

没有 `--auto-roi`、`--auto-machine-comparison`、`--via-reference-nm` 或 `--psd-group-size` 这些可用 CLI 开关；不要把 JSON 中的内部字段名机械转换为命令行选项。

## 7. 输出文件怎么看

先打开 `CD_measurement_200K_V1_15_results.xlsx` 的 `measurement_summary`，再核查 `ROI/` 和双引擎标注图。

| 输出 | 主要用途 |
|---|---|
| `measurement_summary` 工作表 / CSV | 文件名、status、method 在前；每张图 mixed、V13、V10 三行，均含原旋转/未旋转指标；最后一列 `旋转_pitch_CD_nm`，之前为 pitch_source/count/requested_count/status。 |
| `image_summary` | 每图一行的完整宽表、质量状态与诊断字段。 |
| `coordinate_results` | 按坐标/方法展开的长表。 |
| `condition_summary` | 按文件夹汇总，默认仅 OK。 |
| `trench_objects`、`engine_objects` | 每条结构、每个引擎的值、支持数和实际入选信息。 |
| `pitch_periods` / `pitch_samples` | 逐周期和逐采样的 pitch 审计，含入选标记、候选编号、实际支持点、背景/补点排除；最后一列均为 `旋转_pitch_CD_nm`。 |
| `per_sample_results` | 主测量每个采样位置的边缘、有效性、失败原因及补点标记；非 compact 输出还含 average_y0/y1（半开区间）以核查整个平均窗口。 |
| `ROI/`、`roi_summary.csv` | 背景保留区域图和原图坐标范围；绿色框是包围框，暗化部分为排除背景或端部；roi_summary 记录估计长度、每端裁剪像素数及裁前/裁后边界。 |
| `locator_summary.csv` / 同名工作表 | 每图请求/实际/建议定位模式，自动判断分数、歧义标志、周期及阈值。V10/V13_locator_threshold_raw 是定位用阈值，profile_threshold 是模式判断用阈值。 |
| `locator_candidates.csv` / 同名工作表 | 每引擎全部候选的原图 X 坐标、有效性、选择状态、亮度与失败原因。bright-line 额外记录暗核心/亮 line 像素比例及均匀性。 |
| `processing_errors.csv` / 同名工作表 | 读取、引擎及附属输出等异常；即使终端显示完成，也应检查这里。 |
| `settings.json` | 命令行值、逐图实际参数和统计定义，复现实验的重要依据。 |
| `PSD/PSD_V10.xlsx`、`PSD/PSD_V13.xlsx` | 每条结构/连续段的 PSD 汇总与来源清单；不是两个平均谱。 |
| `PSD/V10_psd_curves.csv`、`PSD/V13_psd_curves.csv` | 每个频率点的原始谱数值；Excel 主要是摘要，完整曲线在 CSV。 |
| `PSD/per_structure_psd_summary.csv`、`PSD/per_structure_psd_curves.csv` | 两引擎数据合并存放的长表，未做跨引擎/结构平均。 |
| `PSD/edge_coordinates.csv` | 实际用于 PSD 的单独边缘轨迹；默认不等于主逐点表。 |
| `PSD/measurement_manifest.csv`、`PSD/PSD_settings.json` | PSD 边缘来源、实际平均/平滑值、设置与支持数。 |
| `machine_comparison_detailed.csv`、`V1_15_statistics_and_machine_comparison.xlsx` | 有参考表时的逐图机台对比和统计摘要；另查 unmatched 表定位匹配问题。 |

mixed 的 CD/LWR 来自 V13，左/右 LER 来自 V10；它是兼容性的跨引擎组合，不是一套单一边缘坐标推导的全部统计。旋转是用结构中心线 PCA 做坐标变换，不是重新旋转并插值整幅原图。

CD 为有效宽度均值；LER 为未分组边缘坐标的 3σ；主 LWR 为组均宽度的 3σ。主测量结构内粗糙度使用 `ddof=0`，但不足两个有效点/组时仍输出空值。图像级指标按选中结构分别计算后取算术平均，而不是把所有轨迹拼起来计算一次粗糙度。

## 8. 状态和结果解释

| 状态 | 含义 | 处理建议 |
|---|---|---|
| `OK` | 通过当前工程质量检查 | 仍需核查边缘/标定；不是专业计量精度认证。 |
| `REVIEW` | 部分引擎或指标不足、支持率不足、合成点或部分图像附属输出异常等 | 看 warning、逐点失败原因和图像，不能仅看数值是否存在。默认不参与条件汇总和机台比较。 |
| `ERROR` | 没有可用主指标，或读取/主要处理失败 | 查 processing_errors；保留的空白不是 CD/LER/LWR 等于 0。 |
| `SKIPPED_BACKGROUND` | 没找到可信目标区域 | 核对 pattern、参考宽度、标定和背景筛选门槛；不要直接把背景跳过当成程序崩溃。 |

PSD 表中的 `OK` 只表示该连续段足以计算频谱，**不等同于主测量图像状态 OK**；还应结合 `measurement_stable`、来源和有效点数检查。PSD `SKIPPED` 常表示没有足够长的连续有效段。

本版没有 SEM 噪声去偏、物理 CD 校准或完整不确定度模型。跨图片/跨批次比较时，应固定像素标定、视场/长度、平均与分组、阈值、去趋势等配方。

## 9. 常见问题与排查顺序

### 9.1 检测不到边缘，或整张图被当成背景

依次检查：

1. pattern 是否与实际亮暗目标一致；pixel-size 和参考宽度是否合理。
2. 查看 ROI 图和 roi_summary，目标是否被排除；低对比度或短线可谨慎调整 roi-min-contrast/roi-min-length。
3. 参考宽度自动估计是否受背景影响；已知设计量级时显式设置参考宽度。
4. 检查平均窗口是否大于有效结构长度，手动中心/长度是否超出有效区域。
5. 用 no-auto-roi 加手动 ROI 做对照；保留 threshold-search=bounded，查看逐点 failure_reason，再决定是否调整跟踪/拓扑参数。

不要首先设置 edge-continuity=100。补齐看起来连续的线并不能证明边缘检测正确，也可能改变粗糙度。

### 9.2 明明有数值，为什么是 REVIEW

有可用 CD 不等于所有粗糙度和两引擎都可靠。例如单个有效组无法估计 LWR、稳定结构不足 min-number、左右邻居不齐或某个引擎失败，都可能保留部分数值并标记 REVIEW。

### 9.3 no-auto-machine-comparison 会不会关闭 Excel

不会。它只关闭自动读取 `root/Data.xlsx`；主测量 Excel、CSV、PSD 继续输出。要完全不做机台对比，既加该开关，也不要显式传 machine-excel。

### 9.4 PSD 空、曲线波动大，或与主 LWR 不一致

检查 PSD 表的 reason、edge_coordinates 中的缺口、psd-min-segment 和实际单行检出率。单行提边缘及未平均周期图本来更受噪声影响；分组主 LWR 与单行 PSD 的带宽/滤波/去趋势不同，不应强行要求数值相等。降低段长能输出更多短谱，但频率分辨率也更差。

### 9.5 输出目录非空、Excel 写入失败

换新的输出目录；不要删除旧结果来绕过保护。关闭正被 Excel 占用的目标工作簿，并检查磁盘空间/目录权限。主 Excel 写入失败时，之前已写出的 CSV 可能仍可用于恢复。

### 9.6 使用 stop-on-error 后仍看到后续输出

程序在逐图主处理前有一个发现/读取/尺寸先验扫描阶段，该阶段错误会先记入结果，当前不会因为 stop-on-error 立即结束整个扫描。逐图致命错误触发停止后，还会执行结果导出以保留已完成数据。REVIEW、背景跳过、捕获的附属输出异常也不等同于触发停止的致命主处理错误。

命令正常完成通常返回 0；触发逐图 stop-on-error 或未捕获的批次失败返回 1；参数解析/校验错误返回 2；用户中断返回 130。**返回 0 不保证每张图都是 OK，也不保证 processing_errors 为空**，批处理应同时检查表格状态和诊断。

## 10. 自检脚本与文档核对

自检脚本 `self_check_V1_15.py` 不是主测量入口，它的参数为：

| 参数 | 默认值 | 含义 |
|---|---|---|
| `-h` / `--help` | 不启用 | 查看自检脚本帮助。 |
| `--quick` | 不启用 | 只执行解析/数学/背景等快速检查；不加时还运行端到端与故障场景。 |
| `--output` | 不写报告文件 | 保存本次 JSON 报告的路径；**这个自检选项会覆盖同名报告**，建议使用新文件名，别覆盖随包基准 validation_report.json。 |

```bash
python self_check_V1_15.py --quick
python self_check_V1_15.py --output my_validation_report.json
python sem_cd_measure_200k_batch_V1_15.py --help
```

随包报告包含原有端到端/故障场景，以及新增模式判断、两种定位、端部与 PSD 排除、参数检查；这些是软件回归与随包/合成图验证，不替代真实 SEM 标准样验证。

本手册的覆盖范围以主程序 `build_parser()` 的全部选项为准：编写后逐项检查参数表是否包含每个公开参数，并使用实际解析器检查示例命令及相互冲突的参数。内部 `via_reference_nm` 等历史兼容字段并不是公开接口。本版只在新包内改动；V1_14 包和 ZIP 保留。包内附有同内容的 RUN_GUIDE.md，离线解压即可查阅全部参数。


## 11. V1_15 算法、默认行为和复现

### 11.1 两种定位与 auto

`--pattern trench` 表示测暗沟槽；`--locator-mode bright-line` 表示沟槽之间的 line 整体亮。因此测 trench 时可以且通常应该使用 bright-line，不能把它误写成 `--pattern line`。保留的 `--pattern line` 流程先反相，再应用同一套规则，此时模式名称描述的是反相后工作图的分隔带。

- **dark-line**：调用保留的 V1_14 Basin-First 定位器。代表轮廓采用截尾平均，三类 K-means 分出亮分隔带/盆地，沿 Y 验证盆地，按亮度筛选和周期关系区分真实 trench 与 line 中的假暗条。旧默认周期规则偏好盆地序号间隔 2。新观测掩码会先排除背景和端部。
- **bright-line**：复用代表轮廓和沿 Y 盆地验证；两类 K-means 的暗/亮中心中点给出自适应阈值。暗核心像素低于阈值、两侧 line 核心像素高于阈值的比例均须达到 locator-majority。候选核心的灰度 MAD / 两类灰度差不得超过 0.25，用于排除非均匀区域。所有符合条件的暗区都可成为真实 trench，相邻暗区就是相邻周期，不强制隔一个暗区取一个。检查邻距与估计周期的一致性后，进入原有候选选择和边缘处理。
- **auto（默认）**：在整张图已经识别、裁掉端部的观测区域建立 X 灰度轮廓，再分析是否有周期性重复的深/浅暗谷，或宽/窄暗谷交替。足够证据选 dark-line，否则选 bright-line。判断只做一次，并将同一结果用于 V10、V13 和独立 PSD 提边缘。显式传入模式时始终使用指定模式，同时保留建议值供核查。

模式分析需要观察 line 中心的真实像素：在两条已识别结构之间，仅对长度不超过 2.5 倍参考宽度的内部水平掩码间隔加入 line 灰度上下文。外部背景、整行缺口和裁掉的端部仍排除；该上下文只服务定位，不增加任何可测量采样点。没有有效像素的 X 列只在粗定位轮廓中插值，最终边缘始终取自原图。

auto 至少需 4 个完整暗谷才能形成一组交替证据；使用每连续 4 个暗谷的周期重复性、两种暗谷的灰度差/宽度差，以及同类暗谷重复误差，合成证据的中位数为 `locator_dark_line_score`。分数 ≥0.12 选 dark-line；少于 4 谷或分数在 0.08–0.16 时标记 `locator_auto_ambiguous=True`。这些是可复现的内部经验常数，不是经训练标定的概率，不作为 CLI 参数开放。

若假暗条和真实 trench 的亮度、宽度、节距完全相同，单凭灰度轮廓无法唯一判断物理身份；此时默认 bright-line。严重倾斜、line 太宽、周期不足或两类结构混在一张图内，也可能误判。检查 locator_summary 和候选/标注后用显式模式处理；auto 的歧义标志不自动改变原有 OK/REVIEW 质量策略。

### 11.2 平均长度只用一种近似方法

背景筛选生成前景掩码。对于每一个包含前景的 X 列，记首次有效行 `a_x` 和最后有效行之后的 `b_x`，算 `L_x=b_x-a_x`，再对这些列取算术平均 `L_mean`。这是“近竖直、宽度相近 trench 平均长度”的简单估计，不逐条做精确长度拟合；条宽不等时会按列数加权，同一列有分开的结构时跨度会包含中间缺口。

每端排除 `ceil(L_mean × end-trim-fraction)` 像素。每列只保留 `[a_x+trim, b_x-trim)` 中原本有效的像素；内部背景缺口仍为空，不因裁剪被填上。例如平均长 440px，则上下各排除 22px。若某列太短，裁后该列可全部排除；整图无剩余区域则记 SKIPPED_BACKGROUND。

默认每个主采样的完整 average-range 窗口必须处于裁后有效区域。若手动 extend-length 覆盖端部，被排除的位置标记 background_excluded，不进入 CD/LER/LWR，也不能被 continuity 或 PSD 插值恢复。PSD 单行采样同样受此掩码约束。默认沿线测量长度仍最多 360px，所以长 trench 未必把全部剩余部分都采满；需要覆盖更长区域时按实际尺寸设置 extend-length，整个平均窗口仍受上述限制。

no-auto-roi 时将手动矩形（未给尺寸则整图）视为可观察区域，仍按同一种方法去掉上下端部。这属于手动定义观察对象，不再自动保证矩形内没有背景。

### 11.3 常用模式命令

自动模式、每端默认 5%：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --locator-mode auto --pixel-size 1.3181 --trench-reference-nm 60 --output result_auto --no-auto-machine-comparison
```

已知 line 内有假暗条：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --locator-mode dark-line --pixel-size 1.3181 --trench-reference-nm 60 --output result_dark_line --no-auto-machine-comparison
```

已知 line 整体亮：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --locator-mode bright-line --locator-majority 0.7 --end-trim-fraction 0.05 --pixel-size 1.3181 --trench-reference-nm 60 --output result_bright_line --no-auto-machine-comparison
```

复现 V1_14 关闭自动 ROI 时的行为（仍使用 V1_14 默认 bounded 边缘搜索）：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root input --pattern trench --locator-mode dark-line --no-auto-roi --end-trim-fraction 0 --pixel-size 1.3181 --trench-reference-nm 60 --output result_v114_compatible --no-auto-machine-comparison
```

若复现的是 V1_14 的 threshold-search=legacy 配方，还需同步加上 `--threshold-search legacy`。其余标定、ROI、采样、参考宽度和统计参数也须一致。开启新区域裁剪会改变参与测量的像素/长度，因此不能要求结果与未裁剪版本数值一致。


## 12. line 修复与自动 pitch CD（2026-09-21）

### 12.1 line 失败原因与修复

旧版已在估宽、ROI 和边缘引擎中反相，但定位/灰度模型并非在所有结构上对称：

- 旧盆地宽度限幅实际为 13.5–145 nm，与用户目标宽度无关。180 nm line 被标为 `basin_too_wide`。现使用 `min(30,参考宽度) × 0.45` 至 `max(100,参考宽度) × 1.45` nm，保留既有范围并覆盖目标宽度；双极性采用相同规则。
- 100 nm line 内有 40 nm 假暗条时，反相后成为中间亮条，两侧亮线被拆成约 30 nm 小段。自动估宽可能估成 33 nm；显式给 100 nm 时背景掩码可能只保留外围，最终无候选。
- 新增内部条带识别：对测量极性的灰度轮廓做三类划分，动态范围至少 12，两个类间距均至少占动态范围的 20%；仅在至少两个完整区间中，较高阈值恰好连接两个低灰度平台时启用。每个平台至少 3px、至少占区间 12%，中间间隔至少 3px 和区间 8%；已知参考宽度时完整区间须在其 0.65–1.35 倍内。此识别用于估宽、ROI 和定位，原图像素不会被涂改。
- 边缘阈值的中心参考若比左右 10%–22% 内侧平台均值高出 `max(12,4×两侧差值)`，改用这些平台的中位数均值，保留近边缘真实梯度、跟踪、Viterbi/ERF 和背景约束。逐点 `threshold_reference_mode` 会带 `_edge_interior_reference`。

line 参数填写 line 的参考宽度；同图 trench=60nm、line=100nm 时，不能仅更改 pattern 并仍把参考宽度固定为 60nm。auto/bright-line 的名称按反相后的工作图解释。真实灰度不能区分的结构仍须核查标注。

### 12.2 pitch 定义、旋转和输出

一个 pitch = 一个完整 line + 紧邻的一个完整 trench。使用相邻同类结构的**左边缘到左边缘**定义周期：trench 模式是 trench+line，line 模式是 line+trench。使用实际边缘，不直接相加两个独立平均 CD，也不使用名义参考节距作为测量值。

对周期两边在共同有效 Y 上的边缘中点做 PCA 拟合，沿法向计算 `pitch(y) = (x_next_left(y)-x_left(y)) × |cos(theta)| × pixel_size`。先对一个周期的有效点取均值，再对离测量中心最近的 `max-number` 个不同周期等权平均。仅报告旋转后的 nm 数值。mixed 和 V13 行使用 V13 pitch，V10 行使用 V10 pitch。

主测量选择 N 条结构时通常只有 N−1 个完整周期，所以 pitch 会额外定位并测量一个相邻结构，已有边缘直接复用。此额外结构不计入主 CD/LER/LWR 或 PSD。周期需满足相邻候选序号和观测周期检查：bright-line 的候选序号差为 1；dark-line 允许跨过一条假暗条（序号差最多 2）；距离须为观测周期的 0.60–1.40 倍。此检查使用自动建议的实际形态，避免手动模式把节距翻倍。不跨被漏检结构，不把两倍节距当作一个 pitch。默认要求同一周期至少 2 个共同有效点，且有效比例达到 minimum-valid-fraction。

背景、裁去的端部、无效边缘及 continuity 的合成补点都不进入 pitch。只有 `edge-continuity=100` 补出来的边缘不会被当作 pitch 的实测数据。

`measurement_summary` 的**最后一列**是 `旋转_pitch_CD_nm`；image_summary 的最后一列同样是 mixed/V13 pitch。逐周期和采样审计也提供该末列。`pitch_count` 是实际采用数，`pitch_requested_count` 是 max-number；不足时保留已有周期均值、`pitch_status=PARTIAL`，图片为 REVIEW；完全没有完整周期则为空值、UNAVAILABLE，绝不填 0。图像背景跳过和 ERROR 仍保留列名。主测量原有指标不因 pitch 不足而丢弃。

例：测带内部暗条的亮线，同时自动平均 3 个 pitch：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root examples/input_dark_line --pattern line --pixel-size 1 --line-reference-nm 100 --max-number 3 --output result_line_pitch --no-auto-machine-comparison
```

例：同图测 trench 并自动测相邻周期：

```bash
python sem_cd_measure_200k_batch_V1_15.py --root examples/input_dark_line --pattern trench --pixel-size 1 --trench-reference-nm 60 --max-number 3 --output result_trench_pitch --no-auto-machine-comparison
```

完整自检仍用 `python self_check_V1_15.py --output validation_report.json`；新增测试在 `check_line_pitch.py`，包括同图双模式、宽线、内部暗条、自动估宽、反相对称、倾斜和不同采样支持数、max-number=1/2/5、周期不足、背景/合成点排除和 Excel 末列验证。随包样例及新增图均为合成图，尚无用户实测失败 SEM 图用于验证。
