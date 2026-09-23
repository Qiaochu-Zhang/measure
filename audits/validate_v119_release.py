#!/usr/bin/env python3
"""Verify V1_19 format support, V1_18 numerics, old artifacts and standalone ZIP."""

import ast
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import zipfile

import cv2
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "sem_cd_measure_200k_batch_V1_19_complete_package"
BASELINE = ROOT / "sem_cd_measure_200k_batch_V1_18_complete_package"
BASELINE_COMMIT = "6c4ea0bf75534835ee187120b68a9eea611ce34f"


def run(args, cwd):
    result = subprocess.run([sys.executable, *map(str, args)], cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return result.stdout


def verify_preserved_files():
    files = subprocess.check_output(
        ["git", "ls-tree", "-rz", "--name-only", BASELINE_COMMIT], cwd=ROOT
    ).decode("utf-8").rstrip("\0").split("\0")
    preserved = [name for name in files if any(version in name.lower() for version in ("v1_14", "v1_15", "v1_16", "v1_17", "v1_18", "v114", "v115", "v116", "v117", "v118"))]
    for name in preserved:
        expected = subprocess.check_output(["git", "show", f"{BASELINE_COMMIT}:{name}"], cwd=ROOT)
        assert (ROOT / name).read_bytes() == expected, name
    return dict(status="PASS", baseline_commit=BASELINE_COMMIT, old_files_verified=len(preserved))


def verify_algorithms():
    unchanged = []
    for name in ("cdsem_engine.py", "cdsem_locator.py", "cdsem_localization.py", "cdsem_refinement.py", "cdsem_lines.py", "cdsem_trenches.py", "cdsem_pitch.py", "cdsem_statistics.py"):
        old = (BASELINE / name).read_text().replace("V1_18", "V1_19")
        assert old == (PACKAGE / name).read_text(), name
        unchanged.append(name)
    regions_before = ast.parse((BASELINE / "cdsem_regions.py").read_text())
    regions_after = ast.parse((PACKAGE / "cdsem_regions.py").read_text())
    for name in ("unicode_imread", "normalize_gray_uint8"):
        before = next(node for node in regions_before.body if isinstance(node, ast.FunctionDef) and node.name == name)
        after = next(node for node in regions_after.body if isinstance(node, ast.FunctionDef) and node.name == name)
        assert ast.dump(before) == ast.dump(after), name
    old_psd = (BASELINE / "cdsem_psd.py").read_text().replace("V1_18", "V1_19")
    expected = old_psd.replace("from pathlib import Path", "from pathlib import Path\nfrom cdsem_regions import output_image_path").replace('Path(key).with_suffix(".psd.png")', 'output_image_path(key, ".psd.png")')
    assert ast.dump(ast.parse(expected)) == ast.dump(ast.parse((PACKAGE / "cdsem_psd.py").read_text()))
    tables = []
    for package, version in ((BASELINE, "V1_18"), (PACKAGE, "V1_19")):
        tree = ast.parse((package / f"sem_cd_measure_200k_batch_{version}.py").read_text())
        tables.append({node.name: ast.dump(node) for node in tree.body if isinstance(node, ast.FunctionDef)})
    preserved = [name for name in tables[0] if name.startswith(("coordinate_metrics", "aggregate", "choose_aggregate", "fit_common"))]
    assert preserved
    for name in preserved:
        assert tables[0][name] == tables[1][name], name
    return dict(unchanged_modules=unchanged, unchanged_main_statistics=preserved, grayscale_and_normalization="unchanged", psd_change="output path only")


def verify_compatibility(temp):
    report = {}
    for pattern, reference in (("trench", 60), ("line", 100)):
        print("PNG compatibility:", pattern, flush=True)
        outputs = []
        for package, version in ((BASELINE, "V1_18"), (PACKAGE, "V1_19")):
            output = temp / (version + pattern)
            run([package / f"sem_cd_measure_200k_batch_{version}.py", "--root", BASELINE / "examples/input_dark_line", "--output", output, "--pattern", pattern, "--pixel-size", "1", f"--{pattern}-reference-nm", str(reference), "--viterbi", "1", "--erf-fit", "1", "--sample-number", "32", "--skip-psd", "--skip-statistics-plots", "--no-auto-machine-comparison", "--no-annotated-images"], ROOT)
            outputs.append(output)
        rows = {}
        for name in ("per_sample_results.csv", "pitch_periods.csv", "pitch_samples.csv", *( ("trench_selection.csv", "trench_source_samples.csv") if pattern == "trench" else ("line_trench_pairs.csv", "line_source_samples.csv") )):
            old, new = [pd.read_csv(out / name) for out in outputs]
            pd.testing.assert_frame_equal(old, new, check_exact=True)
            rows[name] = len(old)
        old, new = [pd.read_csv(out / "measurement_summary.csv") for out in outputs]
        columns = [column for column in old if column.startswith(("旋转", "未旋转", "pitch_")) or column in ("image_key", "status", "method")]
        pd.testing.assert_frame_equal(old[columns], new[columns], check_exact=True)
        assert new.status.eq("OK").all()
        report[pattern] = dict(status="PASS", exact_tables=rows, exact_summary_columns=columns)
    return report


def verify_docs(app):
    guide = (ROOT / "V1_19_RUN_GUIDE.md").read_text()
    actions = app.build_parser()._actions
    for action in actions:
        assert any(f"`{name}`" in guide for name in action.option_strings), action.option_strings
    examples = 0
    for line in guide.splitlines():
        if line.startswith("python sem_cd_measure_200k_batch_V1_19.py ") and "--help" not in line:
            args = shlex.split(line)[2:]
            app.validate_args(app.build_parser().parse_args(args), args)
            examples += 1
    return dict(cli_options=len(actions), validated_examples=examples)


def verify_archive(temp):
    archive = ROOT / (PACKAGE.name + ".zip")
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.testzip() is None
        zipped.extractall(temp)
    unpacked = temp / PACKAGE.name
    manifest = (unpacked / "PACKAGE_SHA256.txt").read_text().splitlines()
    for line in manifest:
        digest, path = line.split("  ", 1)
        assert hashlib.sha256((unpacked / path).read_bytes()).hexdigest() == digest, path
        assert (unpacked / path).read_bytes() == (PACKAGE / path).read_bytes(), path
    quick = run(["self_check_V1_19.py", "--quick"], unpacked).strip()
    fixture = cv2.imread(str(unpacked / "examples/input_dark_line/five/test.png"), cv2.IMREAD_UNCHANGED)
    inputs = temp / "standalone_formats"
    inputs.mkdir()
    for suffix in ("tif", "png", "jpg", "jpeg"):
        assert cv2.imwrite(str(inputs / ("sample." + suffix)), fixture)
    output = temp / "standalone_output"
    run(["sem_cd_measure_200k_batch_V1_19.py", "--root", inputs, "--output", output, "--pattern", "trench", "--pixel-size", "1", "--trench-reference-nm", "60", "--no-auto-machine-comparison", "--psd-edge-source", "measurement", "--plot-dpi", "60"], unpacked)
    images = pd.read_csv(output / "image_summary.csv")
    assert len(images) == 4 and images.status.eq("OK").all()
    assert pd.read_csv(output / "processing_errors.csv").empty
    assert len(list((output / "PSD/V10").glob("*.psd.png"))) == 4
    assert len(list((output / "annotated_未旋转_V13").glob("*.png"))) == 4
    return dict(status="PASS", files_verified=len(manifest), sha256=hashlib.sha256(archive.read_bytes()).hexdigest(), quick_self_check=quick, mixed_images=len(images), mixed_cd_nm=list(images["旋转_mixed_CD_nm"]))


def main():
    sys.path.insert(0, str(PACKAGE))
    import sem_cd_measure_200k_batch_V1_19 as app

    report = dict(version="V1_19", preserved_files=verify_preserved_files(), algorithms=verify_algorithms(), documentation=verify_docs(app))
    with tempfile.TemporaryDirectory(prefix="v119-release-") as directory:
        temp = Path(directory)
        report["png_compatibility"] = verify_compatibility(temp)
        print("Standalone ZIP: mixed formats and exports", flush=True)
        report["archive"] = verify_archive(temp)
    (ROOT / "audits/v1_19_release_evidence.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print("V1_19 RELEASE VALIDATION PASSED", flush=True)


if __name__ == "__main__":
    main()
