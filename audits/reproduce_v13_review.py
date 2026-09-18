#!/usr/bin/env python3
"""Reproduce review observations without modifying either packaged engine.

Run from repository root after extracting the ZIP packages and installing their
requirements: .venv/bin/python audits/reproduce_v13_review.py --output /tmp/audit.json
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import replace
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy.special import erf

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "sem_cd_measure_200k_batch_V13_modified_complete_package"
sys.path.insert(0, str(PACKAGE))
import sem_cd_measure_200k_batch_V13_modified as app
import cdsem_v13_modified_psd as psd


def serializable(value):
    if isinstance(value, dict):
        return {str(k): serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    if isinstance(value, np.generic):
        return serializable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def make_rows(widths, left=None, right=None):
    w = np.asarray(widths, float)
    left = -w / 2 if left is None else np.asarray(left, float)
    right = w / 2 if right is None else np.asarray(right, float)
    return [dict(space_index=1, space_role="center", sample_index=i,
                 sample_y=float(i), left_edge_x=float(l), right_edge_x=float(r),
                 valid=True) for i, (l, r) in enumerate(zip(left, right))]


def analytic_cases():
    p = app.v5.MeasurementParams(pixel_size_nm=1)
    widths = 60 + (-1.0) ** np.arange(128)
    alternating = {}
    for g in (1, 4, 128):
        alternating[g] = app.coordinate_metrics_by_space(make_rows(widths), p, g)[1]
    single = app.coordinate_metrics_by_space(make_rows([60]), p, 4)[1]
    l = np.array([-30., np.nan, -28., -29., -30., -31., -30., -29.])
    r = np.array([30., 31., np.nan, 32., 30., 31., 32., 33.])
    paired = app.coordinate_metrics_by_space(make_rows(np.zeros(8), l, r), p, 4)[1]
    mask = np.isfinite(l) & np.isfinite(r)
    ref = np.where(mask, r-l, np.nan)
    expected = 3*np.std([np.nanmean(ref[:4]), np.nanmean(ref[4:])])
    assert np.isclose(paired["unrotated_group4_LWR_nm"], expected)

    # All groups have common left/right masks; support differs between groups.
    incomplete = np.array([0., np.nan, np.nan, np.nan, 0., 2., 4., 6.])
    grouped = app.v7.group_mean(incomplete, 4)
    condition_input = []
    for status, value in [("OK", 10.), ("REVIEW", 100.)]:
        condition_input.append(dict(folder_name="same_folder", pattern="trench",
                                    status=status, **{k: value for k in app.RESULT_COLUMNS}))
    condition = app.build_condition_summary(pd.DataFrame(condition_input)).iloc[0]
    a = app.build_parser().parse_args(["--pattern", "trench", "--group-size", "128"])
    app.validate_args(a, ["--pattern", "trench", "--group-size", "128"])
    # A Fourier-bin sinusoid validates PSD scaling independently of the detector.
    a.psd_method = "periodogram"
    a.psd_window = "boxcar"
    a.psd_detrend = "constant"
    signal = 2*np.sin(2*np.pi*np.arange(128)/16)
    f, power, info = psd.spectrum(signal, 1., a)
    spectral = psd.scalar_metrics(f, power, a.psd_split_wavelength)
    assert np.isclose(spectral["variance_psd_nm2"], np.var(signal))
    filled = make_rows([59, 61, 60, 60, 60, 60, 60, 60])
    for i, row in enumerate(filled):
        row["continuity_synthetic"] = i >= 2
    filled_metrics = app.coordinate_metrics_by_space(filled, p, 1)[1]
    return dict(alternating_width_true_3sigma_nm=float(3*np.std(widths)),
                alternating_width_by_group=alternating, single_valid_point=single,
                paired_mask_actual_lwr=paired["unrotated_group4_LWR_nm"],
                paired_mask_expected_lwr=expected,
                incomplete_groups=dict(means=grouped.tolist(), counts=[1, 4]),
                mixed_quality_condition_mean=condition["rotated_V13_CD_nm_mean"],
                group_size_equal_sample_number_accepted=True,
                synthetic_inclusion=dict(detected_only_lwr_nm=3.,
                    including_synthetic_lwr_nm=filled_metrics["rotated_group4_LWR_nm"]),
                psd_parseval=dict(time_variance=float(np.var(signal)), **spectral, **info))


def averaging_cases():
    """Noiseless sinusoidal edge, known analytic edge, common threshold extraction.

    Isolates average_rows_profile, rather than claiming an end-to-end SEM test.
    """
    y = np.arange(640)
    x = np.arange(128)
    outputs = []
    for wavelength in (8, 32, 128):
        true_edge = 64 + .5*np.sin(2*np.pi*y/wavelength)
        img = 20 + 200*.5*(1-erf((x[None, :]-true_edge[:, None])/(np.sqrt(2)*2)))
        for window in (1, 32):
            measured = []
            for center in range(128, 512):
                profile, *_ = app.v13v5.average_rows_profile(img, center, window)
                edge = app.v13v5.interpolate_threshold_crossing(profile, 120., 50, 80., "left")
                measured.append(edge)
            outputs.append(dict(wavelength_px=wavelength, average_rows=window,
                                true_ler_3sigma_px=float(3*np.std(true_edge[128:512])),
                                measured_ler_3sigma_px=float(3*np.std(measured))))
    return outputs


def run_self_checks():
    results = {}
    for name in ("self_check_V13_modified.py", "self_check_partial_results.py"):
        run = subprocess.run([sys.executable, str(PACKAGE/name)], cwd=PACKAGE,
                             text=True, capture_output=True, check=True)
        results[name] = dict(exit_code=run.returncode, stdout=run.stdout.strip())
    return results


def cli_example(kind, output, extra_args=()):
    reference = "--trench-reference-nm" if kind == "trench" else "--line-reference-nm"
    cmd = [sys.executable, str(PACKAGE/"sem_cd_measure_200k_batch_V13_modified.py"),
           "--root", str(PACKAGE/"examples"/f"input_{kind}"), "--output", str(output),
           "--pattern", kind, "--pixel-size", "1", reference, "60", "--max-number", "4",
           "--no-auto-machine-comparison", "--skip-statistics-plots"]
    cmd.extend(extra_args)
    run = subprocess.run(cmd, cwd=PACKAGE, text=True, capture_output=True, check=True)
    primary = pd.read_csv(output/"image_summary.csv").iloc[0]
    samples = pd.read_csv(output/"per_sample_results.csv")
    spectra = pd.read_csv(output/"PSD"/"per_trench_psd_summary.csv")
    wanted = [c for c in primary.index if c.endswith("_nm") and
              c.startswith(("旋转_", "未旋转_"))]
    overview = dict(exit_code=run.returncode, status=primary["status"],
                    warning=primary["warning"], metrics=primary[wanted].to_dict(),
                    engines={}, psd={})
    for engine, group in samples.groupby("engine"):
        valid = group["valid"].astype(bool)
        overview["engines"][engine] = dict(slots=len(group), valid=int(valid.sum()),
            per_space_valid=group.groupby("space_index")["valid"].sum().to_dict(),
            failures=group.loc[~valid, "failure_reason"].value_counts().to_dict())
    for (engine, signal), group in spectra.groupby(["engine", "signal"]):
        overview["psd"][f"{engine}:{signal}"] = group["status"].value_counts().to_dict()
    overview["condition_includes_review"] = int(pd.read_csv(output/"condition_summary.csv")["review_image_count"].iloc[0])
    overview["diagnostic_error_rows"] = len(pd.read_csv(output/"processing_errors.csv"))
    return overview


def summarize_measurement(measurement, params):
    metrics = app.coordinate_metrics_by_space([dict(r) for r in measurement.sample_rows], params, 4)
    selected = app.v7.choose_aggregate_space_indices(measurement.space_rows, params.min_number)
    return dict(valid=measurement.image_row["valid"],
                slots=len(measurement.sample_rows),
                accepted=sum(bool(r["valid"]) for r in measurement.sample_rows),
                stable=measurement.image_row["stable_space_count"],
                failures=dict(Counter(r["failure_reason"] for r in measurement.sample_rows if not r["valid"])),
                synthetic=sum(bool(r["continuity_synthetic"]) for r in measurement.sample_rows),
                metrics=app.aggregate_coordinate_metrics(metrics, selected))


def crossing_experiment():
    args = app.build_parser().parse_args(["--pattern", "trench", "--pixel-size", "1",
        "--trench-reference-nm", "60", "--max-number", "4",
        "--root", str(PACKAGE/"examples/input_trench")])
    app.validate_args(args, [])
    general = app.v17.GeneralConfig(root_dir=args.root, pixel_size_nm=1,
                                    force_pattern="trench", trench_reference_nm=60)
    parsed, errors = app.v17.discover_images(general)
    assert not errors and len(parsed) == 1
    image = app.v17.read_gray_png(parsed[0].path)
    params = app.build_v13_params(app.build_v10_params(parsed[0], image, general, args))
    baseline = app.v13v5.measure_image(parsed[0].path, params)
    original_crossing = app.v13v5.interpolate_threshold_crossing

    def extended(profile, threshold, start_peak_x, center_x, side):
        start = start_peak_x-2 if side == "left" else start_peak_x+2
        return original_crossing(profile, threshold, start, center_x, side)

    # Diagnostic-only in-memory ablation: identical thresholds and quality gates.
    # The package on disk is never changed; original function is always restored.
    try:
        app.v13v5.interpolate_threshold_crossing = extended
        expanded = app.v13v5.measure_image(parsed[0].path, params)
    finally:
        app.v13v5.interpolate_threshold_crossing = original_crossing
    relaxed_params = replace(params, edge_continuity=20)
    relaxed = app.v13v5.measure_image(parsed[0].path, relaxed_params)
    # Baseline exact mask count independently checked against CLI output.
    return dict(baseline=summarize_measurement(baseline, params),
                search_extended_outward_2px_only=summarize_measurement(expanded, params),
                existing_continuity_20=summarize_measurement(relaxed, relaxed_params))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not PACKAGE.is_dir():
        parser.error("Extract sem_cd_measure_200k_batch_V13_modified_complete_package.zip first")
    files = [p for p in ROOT.glob("*.zip")] + sorted(PACKAGE.glob("*.py"))
    original = ROOT/"sem_cd_measure_200k_batch_V1_13_complete_package/cdsem_v13_edge_engine.py"
    if original.exists():
        files.append(original)
    evidence = dict(environment=dict(python=platform.python_version(),
        packages={n: importlib.metadata.version(n) for n in
                  ("numpy", "pandas", "scipy", "matplotlib", "tifffile", "openpyxl", "opencv-python")}),
        sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    print("Running analytic checks and self-checks", flush=True)
    evidence["analytic"] = analytic_cases()
    evidence["averaging"] = averaging_cases()
    evidence["self_checks"] = run_self_checks()
    with tempfile.TemporaryDirectory(prefix="measure-v13-audit-") as temp:
        evidence["examples"] = {}
        for kind in ("trench", "line"):
            print(f"Running full {kind} pipeline, including PSD and annotation", flush=True)
            evidence["examples"][kind] = cli_example(kind, Path(temp)/kind)
        print("Checking single-group LWR with an otherwise valid image", flush=True)
        evidence["examples"]["trench_one_group"] = cli_example("trench", Path(temp)/"one_group",
            ["--group-size", "128", "--edge-continuity", "20"])
        print("Running crossing-search ablation and continuity comparison", flush=True)
        evidence["crossing_experiment"] = crossing_experiment()
    trench = app.v17.unicode_imread(PACKAGE/"examples/input_trench/five/test.png")
    line = app.v17.unicode_imread(PACKAGE/"examples/input_line/five/test.png")
    evidence["polarity_examples"] = dict(shape=list(trench.shape),
        max_abs_line_minus_inverted_trench=float(np.max(np.abs(line.astype(float)-(255-trench.astype(float))))))
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest
               for p, digest in evidence["sha256"].items())
    evidence["package_sources_unchanged_after_audit"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(serializable(evidence), ensure_ascii=False,
                                     indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
