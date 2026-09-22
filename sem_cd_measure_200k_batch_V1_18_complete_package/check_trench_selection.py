"""Regressions for measuring the full trench population before final selection."""

from dataclasses import replace
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import cdsem_regions as regions
from cdsem_refinement import choose_objects
from check_localization import synthetic_pattern
from check_line_pitch import stripe_image


RAW_FIELDS = [
    "sample_index",
    "未旋转_sample_y_px",
    "未旋转_left_edge_x_px",
    "未旋转_right_edge_x_px",
    "valid",
    "background_excluded",
    "continuity_synthetic",
    "left_threshold",
    "right_threshold",
]


def analytical_tests(require, app):
    rows = [
        dict(
            basin_sequence_index=i,
            basin_valid=True,
            brightness_true_trench=real,
            candidate_center_x_px=x,
            basin_center_x_px=x,
            basin_left_x_px=x - 20,
            basin_right_x_px=x + 20,
            estimated_width_px=40,
            basin_quality_score=1,
        )
        for i, (x, real) in enumerate(
            [(100, True), (180, False), (260, True), (340, False), (420, True)]
        )
    ]
    params = app.edges.MeasurementParams(center_x_px=340)
    true = [r for r in rows if r["brightness_true_trench"]]
    legacy = choose_objects(
        true, rows, 340, params, app.edges.choose_center_left_right_v115
    )
    full = choose_objects(
        true,
        rows,
        340,
        replace(params, max_number=1000, min_number=1, require_center_left_right=False),
        app.edges.choose_center_left_right_v115,
    )
    full = [
        (role, r)
        for role, r in full
        if r["basin_valid"] and r["brightness_true_trench"]
    ]
    require(
        any(not r["brightness_true_trench"] for _, r in legacy),
        "reproduce false center in legacy selector",
    )
    require(
        len(full) == 3 and all(r["brightness_true_trench"] for _, r in full),
        "full search rejects false center",
    )
    return dict(
        status="PASS",
        legacy_false_center_reproduced=True,
        full_source_population=[r["candidate_center_x_px"] for _, r in full],
    )


def verify_trench_exports(output, require, maximum=None):
    selection = pd.read_csv(output / "trench_selection.csv")
    sources = pd.read_csv(output / "trench_source_samples.csv")
    chosen = selection[selection.selected]
    require(not chosen.empty, "selected qualified trench outputs")
    require(chosen.qualified.all(), "no unqualified trench in primary output")
    if maximum is not None:
        require(
            chosen.groupby("engine").size().eq(maximum).all(),
            "max-number final outputs",
        )
    samples = pd.read_csv(output / "per_sample_results.csv")
    for (engine, index), group in samples.groupby(["engine", "space_index"]):
        item = chosen[
            (chosen.engine == engine) & (chosen.selected_space_index == index)
        ].iloc[0]
        source = sources[
            (sources.engine == engine)
            & (sources.source_trench_id == item.source_trench_id)
        ]
        require(
            group.trench_source_id.eq(item.source_trench_id).all(),
            "source identity retained",
        )
        pd.testing.assert_frame_equal(
            group[RAW_FIELDS].reset_index(drop=True),
            source[RAW_FIELDS].reset_index(drop=True),
            check_exact=True,
        )
        require(source.selected_for_output.all(), "source selection audit")
    summary = pd.read_csv(output / "measurement_summary.csv")
    require(summary.columns[-1] == "旋转_pitch_CD_nm", "pitch stays last column")
    periods = pd.read_csv(output / "pitch_periods.csv")
    for row in periods[periods.selected].to_dict("records"):
        for side in ("left_basin_index", "right_basin_index"):
            candidate = selection[
                (selection.engine == row["engine"])
                & (selection.source_trench_id == row[side])
            ]
            require(
                len(candidate) == 1 and candidate.qualified.all(),
                "pitch only from qualified sources",
            )
    return dict(
        status="PASS",
        measured_trenches=selection.groupby("engine").size().to_dict(),
        qualified_trenches=selection.groupby("engine").qualified.sum().to_dict(),
        selected_trenches=chosen.groupby("engine").size().to_dict(),
        exact_source_sample_rows=len(samples),
    )


def verify_legacy_quality_gate(output, frame, require):
    selection = pd.read_csv(output / "trench_selection.csv")
    sources = pd.read_csv(output / "trench_source_samples.csv")
    # This legacy/no-ROI recipe formerly included V13 objects with only ~35%
    # valid samples in its fallback mean. V1_18 must exclude them, not preserve
    # a summary which contradicts the requested qualify-before-select policy.
    v13 = selection[selection.engine == "V13"]
    require(
        not v13.empty and not v13.qualified.any() and not v13.selected.any(),
        "legacy low-support V13 objects excluded",
    )
    require(
        frame.status.eq("REVIEW").all()
        and frame["旋转_V13_CD_nm"].isna().all()
        and frame["旋转_V10_CD_nm"].notna().all(),
        "keep qualified engine partial result",
    )
    # Golden per-source coordinates measured independently with V1_17. The new
    # population policy must not move the retained physical trench edges.
    expected = {
        "V10": {
            1: (131.19997034340386, 188.8450712304535, 128),
            2: (231.204325217617, 288.8451092450944, 128),
            3: (331.1997329188556, 388.8465752700666, 128),
        },
        "V13": {
            1: (129.80266205046206, 189.80239956134102, 45),
            2: (229.792714093734, 289.81961938978577, 42),
            3: (329.80480719094055, 389.82005231015626, 44),
        },
    }
    for engine, items in expected.items():
        for tid, (left, right, count) in items.items():
            group = sources[
                (sources.engine == engine) & (sources.source_trench_id == tid)
            ]
            require(
                abs(group["未旋转_left_edge_x_px"].mean() - left) < 1e-9
                and abs(group["未旋转_right_edge_x_px"].mean() - right) < 1e-9
                and group.valid.sum() == count,
                "unchanged legacy per-source edge math",
            )


def integration_tests(temp, execute, require, app):
    report = {}
    quiet = ["--skip-psd", "--no-annotated-images"]
    root = temp / "full_population_input"
    regions.unicode_imwrite(root / "counts.png", stripe_image())
    for maximum in [1, 2, 5]:
        name = f"trench_all_count_{maximum}"
        out, frame = execute(
            name, root=root, extra=[*quiet, "--max-number", str(maximum)]
        )
        require(frame.status.eq("OK").all(), name + " status")
        result = verify_trench_exports(out, require, maximum)
        require(
            min(result["measured_trenches"].values()) >= 7,
            "measure all sources even when max-number=1",
        )
        require(
            pd.read_csv(out / "measurement_summary.csv").pitch_count.eq(maximum).all(),
            "max-number pitch periods",
        )
        report[name] = result

    root = temp / "false_center_input"
    regions.unicode_imwrite(root / "center.png", synthetic_pattern(True))
    out, frame = execute(
        "trench_false_center",
        root=root,
        extra=[
            *quiet,
            "--max-number",
            "3",
            "--center-x",
            "620",
            "--meas-area-width",
            "900",
        ],
    )
    require(
        frame.status.eq("OK").all()
        and abs(frame["旋转_mixed_CD_nm"].iloc[0] - 60) < 1.2,
        "false center does not replace a physical trench",
    )
    selected = pd.read_csv(out / "trench_selection.csv")
    require(
        not selected.source_trench_id.mod(2).any(),
        "interior false dark stripes excluded",
    )
    report["trench_false_center"] = verify_trench_exports(out, require, 3)

    root = temp / "trench_qualification_input"
    regions.unicode_imwrite(root / "quality.png", synthetic_pattern(False))
    out, frame = execute(
        "trench_nearest_bright_neighbors",
        root=root,
        extra=[*quiet, "--max-number", "3"],
    )
    require(frame.status.eq("OK").all(), "bright-line nearest neighbors status")
    selection = pd.read_csv(out / "trench_selection.csv")
    for _, group in selection[selection.selected].groupby("engine"):
        require(
            np.allclose(sorted(group.center_x_px), [380, 540, 700], atol=2),
            "resolved bright-line basin gap is preserved in final selection",
        )
    report["trench_nearest_bright_neighbors"] = verify_trench_exports(out, require, 3)
    original = app.edges.measure_single_space
    for name, allowed, maximum in [
        ("trench_rejected_center", {0, 1, 3, 4}, 3),
        ("trench_insufficient_qualified", {1, 2}, 3),
        ("trench_none_qualified", set(), 3),
    ]:

        def reject(*args, **kwargs):
            result, rows, annotation = original(*args, **kwargs)
            if args[2].v115_basin_sequence_index not in allowed:
                result.stable = False
                result.warning = "injected quality failure to test population selection"
            return result, rows, annotation

        with patch.object(app.edges, "measure_single_space", side_effect=reject):
            out, frame = execute(
                name, root=root, extra=[*quiet, "--max-number", str(maximum)]
            )
        selection = pd.read_csv(out / "trench_selection.csv")
        require(
            selection.groupby("engine").size().eq(5).all(),
            "all five measured before filtering",
        )
        require(
            not selection.loc[~selection.qualified, "selected"].any(),
            "rejected sources never fill quota",
        )
        if allowed:
            result = verify_trench_exports(out, require, min(maximum, len(allowed)))
            require(
                frame.status.eq("REVIEW").all(),
                name + " quality state",
            )
            summary = pd.read_csv(out / "measurement_summary.csv")
            require(
                summary.pitch_count.eq(2 if len(allowed) == 4 else 1).all()
                and summary.pitch_status.eq("PARTIAL").all(),
                "rejected source reduces full pitch count without discarding qualified CD",
            )
            periods = pd.read_csv(out / "pitch_periods.csv")
            for row in periods[periods.selected].to_dict("records"):
                require(
                    row["right_basin_index"] - row["left_basin_index"] == 1,
                    "pitch does not bridge rejected source",
                )
            report[name] = dict(result, injected_rejections=True)
        else:
            require(
                frame.status.eq("ERROR").all()
                and frame["旋转_mixed_CD_nm"].isna().all(),
                "all rejected means no result",
            )
            require(
                not (out / "per_sample_results.csv")
                .read_text(encoding="utf-8-sig")
                .strip(),
                "no unqualified sample output",
            )
            report[name] = dict(
                status="PASS",
                injected_rejections=True,
                measured_per_engine=5,
                selected=0,
            )

    root = temp / "single_physical_trench_input"
    gray = np.full((400, 400), 210.0)
    gray[:, 170:230] = 40
    regions.unicode_imwrite(root / "one.png", gaussian_filter(gray, 1).astype(np.uint8))
    out, frame = execute(
        "trench_single_without_neighbor",
        root=root,
        extra=[
            *quiet,
            "--max-number",
            "1",
            "--no-auto-roi",
            "--end-trim-fraction",
            "0",
        ],
    )
    require(
        frame["旋转_mixed_CD_nm"].notna().all(),
        "a single trench needs no neighboring trench to measure its CD",
    )
    require(frame.status.eq("REVIEW").all(), "one trench has no complete pitch")
    report["trench_single_without_neighbor"] = verify_trench_exports(out, require, 1)

    root = temp / "cross_mode_source_input"
    regions.unicode_imwrite(root / "shared.png", synthetic_pattern(True))
    outputs = {}
    for pattern, reference, space in [("trench", 60, 100), ("line", 100, 60)]:
        out, frame = execute(
            "same_sources_" + pattern,
            root=root,
            pattern=pattern,
            reference=reference,
            extra=[
                *quiet,
                "--max-number",
                "3",
                "--space-reference-nm",
                str(space),
                "--viterbi",
                "1",
                "--erf-fit",
                "1",
                "--threshold-left",
                "40",
                "--threshold-right",
                "60",
                "--sample-number",
                "32",
            ],
        )
        require(frame.status.eq("OK").all(), "cross-mode optional measurement status")
        outputs[pattern] = out
    trench = pd.read_csv(outputs["trench"] / "trench_source_samples.csv")
    line = pd.read_csv(outputs["line"] / "line_source_samples.csv")
    columns = [
        "image_key",
        "engine",
        "source_trench_id",
        "source_trench_qualified",
        *RAW_FIELDS,
    ]
    pd.testing.assert_frame_equal(trench[columns], line[columns], check_exact=True)
    report["same_sources_trench"] = verify_trench_exports(outputs["trench"], require, 3)
    report["same_sources_line"] = dict(
        status="PASS",
        exact_cross_mode_source_rows=len(line),
        viterbi=True,
        erf=True,
        asymmetric_thresholds=True,
    )

    original_image = app.edges.measure_image
    checked = []

    def inspect(*args, **kwargs):
        measurement = original_image(*args, **kwargs)
        p = args[1]
        if (
            p.pattern_kind == "trench"
            and p.average_range_px == 1
            and p.smoothing_pixel == 1
        ):
            payload = measurement.annotation_payload
            selection = payload["trench_selection"]
            require(len(selection) == 5, "PSD also measures full population")
            require(
                all(r["qualified"] for r in selection if r["selected"]),
                "PSD only selected qualified trenches",
            )
            require(len(payload["spaces"]) == 3, "PSD final quota")
            for ann in payload["spaces"]:
                rows = [
                    r
                    for r in measurement.sample_rows
                    if r["space_index"] == ann["space_index"]
                ]
                for side in ["left", "right"]:
                    np.testing.assert_array_equal(
                        ann[side + "_edges"],
                        [r[side + "_edge_x"] if r["valid"] else np.nan for r in rows],
                    )
            checked.append(p.engine_kind)
        return measurement

    with patch.object(app.edges, "measure_image", side_effect=inspect):
        out, frame = execute(
            "trench_all_psd",
            root=root,
            extra=["--max-number", "3", "--no-annotated-images"],
        )
    require(
        frame.status.eq("OK").all() and set(checked) == {"V10", "V13"},
        "both single-row PSD passes use new trench selection",
    )
    psd = pd.read_csv(out / "PSD/per_structure_psd_summary.csv")
    require(
        psd.status.eq("OK").all()
        and psd.average_range_px.eq(1).all()
        and not psd.spectrum_averaging.any(),
        "PSD computation policy retained",
    )
    report["trench_all_psd"] = dict(status="PASS", engines=checked, curves=len(psd))
    return report
