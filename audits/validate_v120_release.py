#!/usr/bin/env python3
"""Validate V1_20 isolation, unchanged measurement outputs and standalone ZIP."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "sem_cd_measure_200k_batch_V1_20_complete_package"
BASELINE = ROOT / "sem_cd_measure_200k_batch_V1_19_complete_package"
BASELINE_COMMIT = "3629670330d358fa028a0dcf384fd286b9055321"


def run(args, cwd):
    result = subprocess.run(
        [sys.executable, *map(str, args)], cwd=cwd, capture_output=True, text=True
    )
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return result.stdout


def verify_preserved_repository():
    names = subprocess.check_output(
        ["git", "ls-tree", "-rz", "--name-only", BASELINE_COMMIT], cwd=ROOT
    ).decode("utf-8").rstrip("\0").split("\0")
    for name in names:
        expected = subprocess.check_output(
            ["git", "show", f"{BASELINE_COMMIT}:{name}"], cwd=ROOT
        )
        assert (ROOT / name).read_bytes() == expected, name
    return dict(status="PASS", baseline_commit=BASELINE_COMMIT, files_verified=len(names))


def verify_measurement_compatibility(temp):
    report = {}
    for pattern, reference in (("trench", 60), ("line", 100)):
        outputs = []
        for package, version in ((BASELINE, "V1_19"), (PACKAGE, "V1_20")):
            output = temp / f"{version}_{pattern}"
            args = [
                package / f"sem_cd_measure_200k_batch_{version}.py",
                "--root", BASELINE / "examples/input_dark_line",
                "--output", output,
                "--pattern", pattern,
                "--pixel-size", "1",
                f"--{pattern}-reference-nm", str(reference),
                "--viterbi", "1",
                "--erf-fit", "1",
                "--sample-number", "32",
                "--skip-psd",
                "--skip-statistics-plots",
                "--no-auto-machine-comparison",
                "--no-annotated-images",
            ]
            run(args, ROOT)
            outputs.append(output)
        checked = {}
        names = [
            "measurement_summary.csv", "per_sample_results.csv",
            "pitch_periods.csv", "pitch_samples.csv",
            "trench_selection.csv" if pattern == "trench" else "line_trench_pairs.csv",
            "trench_source_samples.csv" if pattern == "trench" else "line_source_samples.csv",
        ]
        for name in names:
            old, new = [pd.read_csv(output / name) for output in outputs]
            shared = [column for column in old.columns if column in new.columns]
            pd.testing.assert_frame_equal(old[shared], new[shared], check_exact=True)
            checked[name] = len(old)
        report[pattern] = dict(status="PASS", exact_rows=checked)
    return report


def verify_psd_contract(temp):
    output = temp / "psd_contract"
    run(
        [
            PACKAGE / "sem_cd_measure_200k_batch_V1_20.py",
            "--root", PACKAGE / "examples/input_wiggling",
            "--output", output,
            "--pattern", "trench",
            "--pixel-size", "1",
            "--trench-reference-nm", "60",
            "--max-number", "1",
            "--no-auto-machine-comparison",
            "--no-annotated-images",
            "--plot-dpi", "60",
        ],
        ROOT,
    )
    pngs = sorted((output / "PSD").rglob("*.psd.png"))
    assert [path.parent.name for path in pngs] == ["V10", "V13"]
    summary = pd.read_csv(output / "PSD/image_psd_summary.csv")
    assert len(summary) == 4
    assert set(summary.engine) == {"V10", "V13"}
    assert set(summary.signal) == {"LER", "LWR"}
    assert set(summary.axis) == {"normal"}
    assert summary.status.eq("OK").all()
    assert summary.contributing_structures.gt(1).all()
    assert summary.spectrum_averaging.all()
    bands = pd.read_csv(output / "PSD/image_psd_band_summary.csv")
    assert set(bands.band) == {"LF", "MF", "HF"}
    settings = json.loads((output / "PSD/PSD_settings.json").read_text())
    assert settings["axis"] == "normal" and settings["spectrum_averaging"]
    assert settings["low_mid_boundary_per_nm"] == 0.01
    assert settings["mid_high_boundary_per_nm"] == 0.05
    return dict(
        status="PASS", image_count=1, png_count=2,
        engines=["V10", "V13"], signals=["LER", "LWR"],
        all_qualified_structures=int(summary.contributing_structures.min()),
        main_max_number=1,
    )


def verify_archive(temp):
    archive = ROOT / "sem_cd_measure_200k_batch_V1_20_complete_package.zip"
    with zipfile.ZipFile(archive) as source:
        assert source.testzip() is None
        source.extractall(temp / "unpacked")
        members = source.namelist()
    unpacked = temp / "unpacked" / PACKAGE.name
    output = run(["self_check_V1_20.py", "--quick"], unpacked)
    assert "V1_20 SELF-CHECK PASSED" in output
    manifest = {}
    for line in (unpacked / "PACKAGE_SHA256.txt").read_text().splitlines():
        digest, relative = line.split("  ", 1)
        assert hashlib.sha256((unpacked / relative).read_bytes()).hexdigest() == digest
        manifest[relative] = digest
    return dict(
        status="PASS", archive=archive.name, files=len(members),
        sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
        manifest_files=len(manifest), quick_self_check="PASS",
    )


def main():
    report = {"preservation": verify_preserved_repository()}
    with tempfile.TemporaryDirectory(prefix="validate-v120-") as directory:
        temp = Path(directory)
        report["measurement_compatibility"] = verify_measurement_compatibility(temp)
        report["psd_contract"] = verify_psd_contract(temp)
        report["standalone_archive"] = verify_archive(temp)
    target = ROOT / "audits/v1_20_release_evidence.json"
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
