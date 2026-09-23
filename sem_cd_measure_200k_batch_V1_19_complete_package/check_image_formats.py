"""Real codec, mixed-directory and export checks for V1_19 image inputs."""

import json
from pathlib import Path
import tempfile

import cv2
import numpy as np
import pandas as pd
from openpyxl import Workbook


BASE = Path(__file__).resolve().parent


def analytical_tests(require, app):
    regions = app.regions
    gray = regions.read_gray_image(BASE / "examples/input_dark_line/five/test.png")
    report = {}
    with tempfile.TemporaryDirectory(prefix="v119-formats-") as directory:
        root = Path(directory)
        folder = root / "中文 子目录"
        folder.mkdir()
        for extension in regions.SUPPORTED_IMAGE_SUFFIXES:
            path = folder / ("sample" + extension.upper())
            regions.unicode_imwrite(path, gray)
            decoded = regions.read_gray_image(path)
            require(decoded.shape == gray.shape, f"shape: {extension}")
            require(decoded.dtype == np.uint8, f"dtype: {extension}")
            if extension not in {".jpg", ".jpeg"}:
                require(np.array_equal(decoded, gray), f"lossless pixels: {extension}")
            else:
                require(np.abs(decoded.astype(float) - gray).mean() < 2, "JPEG decode")
            report[extension] = "PASS"

        # Both prescan and edge engine must accept high-bit-depth TIFF and PNG.
        for extension in (".tif", ".png"):
            path = folder / ("gray16" + extension)
            regions.unicode_imwrite(path, gray.astype(np.uint16) * 257)
            require(regions.unicode_imread(path).dtype == np.uint16, "16-bit preserved")
            require(np.array_equal(regions.read_gray_image(path), regions.normalize_gray_uint8(gray.astype(np.uint16) * 257)), "16-bit prescan normalization")
            measured, _ = app.edges.read_image(path)
            expected, _ = app.edges.normalize_to_8bit_float(gray.astype(np.uint16) * 257)
            np.testing.assert_array_equal(measured, expected)
        for extension, channels in ((".tiff", 3), (".png", 4), (".jpeg", 3)):
            color = np.stack([gray] * channels, axis=-1)
            if channels == 4:
                color[..., 3] = 255
            path = folder / ("color" + extension)
            regions.unicode_imwrite(path, color)
            expected = cv2.cvtColor(regions.unicode_imread(path), cv2.COLOR_BGRA2GRAY if channels == 4 else cv2.COLOR_BGR2GRAY)
            np.testing.assert_array_equal(regions.read_gray_image(path), expected)
            np.testing.assert_array_equal(app.edges.read_image(path)[0], expected.astype(float))

        # OpenCV uses the first TIFF page, including for multipage files.
        multipage = root / "multipage.tiff"
        require(cv2.imwritemulti(str(multipage), [gray, 255 - gray]), "write multipage TIFF")
        np.testing.assert_array_equal(regions.read_gray_image(multipage), gray)
        np.testing.assert_array_equal(app.edges.read_image(multipage)[0], gray.astype(float))

        broken = folder / "broken.JpEg"
        broken.write_bytes(b"not an image")
        (folder / "ignored.txt").write_bytes(b"not an image")
        (folder / "directory.tif").mkdir()
        output = root / "current_output"
        output.mkdir()
        (output / "ignored.tif").write_bytes(b"output")
        for version in ("V1_18", "V1_19"):
            previous = root / version
            previous.mkdir()
            (previous / "settings.json").write_text(json.dumps({"script_version": version}))
            (previous / "ignored.jpg").write_bytes(b"output")
        records, errors = regions.discover_images(regions.GeneralConfig(root_dir=root, output_dir=output, force_pattern="trench", pixel_size_nm=1))
        require(len(records) == 13, "all supported files recursively discovered")
        require(len(errors) == 1 and Path(errors[0]["path"]) == broken, "corrupt image reported once")
        require(len({str(regions.output_image_path(str(r.path.relative_to(root)))) for r in records}) == len(records), "unique output paths")
        # Full extension is retained even if a source stem ends in another format.
        keys = ("sample.tif", "sample.png", "sample.tif.png")
        require(len({regions.output_image_path(k) for k in keys}) == len(keys), "compound-name collision")
    report.update(high_bit_depth="PASS", color="PASS", multipage_first_page="PASS", discovery_and_errors="PASS")
    return report


def integration_tests(temp, execute, require, app):
    regions = app.regions
    gray = regions.read_gray_image(BASE / "examples/input_dark_line/five/test.png")
    root = temp / "mixed_format_input"
    folder = root / "中文 子目录"
    names = ["sample.png", "sample.TIF", "sample.tiff", "sample.jpg", "sample.JPEG", "sample.bmp", "sample.webp", "sample.TIF.png", "gray16.tif", "color.tiff"]
    for name in names:
        data = gray
        if name == "gray16.tif":
            data = gray.astype(np.uint16) * 257
        elif name == "color.tiff":
            data = np.stack([gray] * 3, axis=-1)
        regions.unicode_imwrite(folder / name, data)
    book = Workbook()
    book.active.title = "v2"
    book.active.append(["A", "name", "C", "CD", "LER", "LWR"])
    book.active.append([None, "gray16.tif", None, 60, 2, 1])
    book.active.append([None, "sample.tif", None, 60, 2, 1])
    reference = temp / "format_reference.xlsx"
    book.save(reference)
    output, frame = execute("mixed_format_exports", root=root, plots=True, extra=["--max-number", "3", "--save-debug-masks", "--psd-edge-source", "measurement", "--plot-dpi", "60", "--machine-excel", str(reference)])
    require(len(frame) == len(names) and frame.status.eq("OK").all(), "mixed-format measurement statuses")
    require(pd.read_csv(output / "processing_errors.csv").empty, "all format exports succeed")
    settings = json.loads((output / "settings.json").read_text())
    require(settings["input_image_count"] == settings["processed_image_count"] == len(names), "generic image counters")
    require(settings["supported_image_suffixes"] == list(regions.SUPPORTED_IMAGE_SUFFIXES), "record supported formats")
    for name in names:
        require((output / "ROI" / "中文 子目录" / (name + ".png")).is_file(), "ROI export: " + name)
        for engine in ("V10", "V13"):
            for path in (output / f"annotated_未旋转_{engine}" / "中文 子目录" / (name + ".png"), output / "debug" / "中文 子目录" / (name + f"_{engine}_rotated.png"), output / "PSD" / engine / "中文 子目录" / (name + ".psd.png")):
                require(path.is_file(), "independent PNG export: " + str(path))
                require(path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"), "PNG encoding")
    measured = frame.set_index("image")
    columns = [c for c in frame if c.startswith(("旋转_", "未旋转_"))]
    for name in ("sample.TIF", "sample.tiff", "sample.bmp", "sample.webp", "sample.TIF.png", "color.tiff"):
        pd.testing.assert_series_equal(measured.loc["sample.png", columns], measured.loc[name, columns], check_names=False, check_exact=True)
    require((measured["旋转_mixed_CD_nm"] - 60).abs().lt(1).all(), "known trench width across codecs")
    machine = pd.read_csv(output / "machine_comparison_detailed.csv")
    require(list(machine.image) == ["gray16.tif"], "non-PNG machine match and duplicate exclusion")
    unmatched = pd.read_csv(output / "machine_unmatched_images.csv")
    require("duplicate_basename_in_recursive_tree" in set(unmatched.unmatched_reason.dropna()), "ambiguous basenames reported")

    # Exercise independent single-row PSD and the line path on non-PNG input.
    line_root = temp / "tiff_line_input"
    regions.unicode_imwrite(line_root / "line.TIFF", gray.astype(np.uint16) * 257)
    line_output, line = execute("tiff_line_single_row_psd", root=line_root, pattern="line", reference=100, extra=["--no-annotated-images", "--plot-dpi", "60"])
    require(line.status.eq("OK").all(), "16-bit TIFF line status")
    require((line["旋转_mixed_CD_nm"] - 100).abs().lt(1).all(), "known line width")
    psd = pd.read_csv(line_output / "PSD/per_structure_psd_summary.csv")
    require(set(psd.engine) == {"V10", "V13"} and psd.status.eq("OK").all(), "TIFF single-row PSD")

    broken_root = temp / "broken_format_input"
    regions.unicode_imwrite(broken_root / "valid.jpg", gray)
    for extension in ("tif", "png", "jpg", "jpeg"):
        (broken_root / ("broken." + extension)).write_bytes(b"invalid")
    broken_output, broken = execute("mixed_format_errors", root=broken_root, extra=["--skip-psd", "--no-annotated-images", "--max-number", "3"])
    require(sum(broken.status == "ERROR") == 4 and sum(broken.status == "OK") == 1, "bad files recorded; valid input still processed")
    require(len(pd.read_csv(broken_output / "processing_errors.csv")) == 4, "one error per bad image")
    return dict(status="PASS", mixed_images=len(names), exports_per_image=7, lossless_measurements_exact=True, tiff_line_and_psd="PASS", corrupted_images=4)
