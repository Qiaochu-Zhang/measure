# V1_19 发布说明

V1_19 基于 V1_18，支持 TIF/TIFF、PNG、JPG/JPEG、BMP、WebP 混合输入，扩展名不区分大小写，支持中文路径。无需增加参数或依赖。

- [完整代码 ZIP](sem_cd_measure_200k_batch_V1_19_complete_package.zip)
- [完整参数手册](V1_19_RUN_GUIDE.md)
- [包内 README](sem_cd_measure_200k_batch_V1_19_complete_package/README.md)
- [修改记录](sem_cd_measure_200k_batch_V1_19_complete_package/CHANGELOG.md)
- [回归验证报告](sem_cd_measure_200k_batch_V1_19_complete_package/validation_report.json)
- [发布与旧版对照验证](audits/v1_19_release_evidence.json)

支持 8/16 位灰度 TIFF/PNG 和彩色图，多页 TIFF 仅读取第一页；像素标定仍通过 pixel-size 指定。沿用 V1_18 的灰度/归一化及 trench、line、pitch、CD/LER/LWR、PSD 算法。

ROI、标注、debug 和 PSD 图统一输出 PNG，并在文件名中保留输入扩展名，避免 sample.tif、sample.png、sample.jpg 的结果互相覆盖。PNG 输入的 ROI 名也改为 sample.png.png。计数字段改为 input_image_count / processed_image_count；机台未匹配图表改为 machine_unmatched_images.csv。机台比较仍按去扩展名的文件名匹配，同名歧义会被排除并记录原因。

```bash
python -m pip install -r requirements_V1_19.txt
python self_check_V1_19.py --quick
python sem_cd_measure_200k_batch_V1_19.py --root /path/to/images --pattern trench --pixel-size 1 --trench-reference-nm 60 --max-number 3 --output results_V1_19 --no-auto-machine-comparison
```

请按实际图像修改 pixel-size 和参考宽度。JPEG 有损压缩可能改变测量结果，优先使用原始无损图。旧版目录、报告及 ZIP 保持原样。测试基于合成图，完整结果见上述报告。
