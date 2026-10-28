"""Verify the physical shared-edge contract, including rejected source trenches."""

from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import cdsem_lines as lines
import cdsem_regions as regions
from check_localization import synthetic_pattern
from check_line_pitch import stripe_image


def analytical_tests(require, app):
    params = app.edges.MeasurementParams(pixel_size_nm=2, target_cd_nm=40)
    candidate = SimpleNamespace(
        role="center", center_x_px=30, global_left_edge_px=20, global_right_edge_px=40
    )
    a = dict(
        sample_index=0,
        sample_y=1.0,
        valid=True,
        source_trench_id=0,
        left_edge_x=10.0,
        right_edge_x=20.0,
        continuity_original_valid=True,
        right_threshold=87.0,
        right_is_flyer=False,
        right_peak_fallback=False,
    )
    b = dict(
        a,
        source_trench_id=1,
        left_edge_x=40.0,
        right_edge_x=50.0,
        left_threshold=93.0,
        left_is_flyer=False,
        left_peak_fallback=False,
    )
    row = lines.shared_rows([a], [b], candidate, 1, params)[0]
    require(
        row["valid"]
        and row["left_edge_x"] == 20
        and row["right_edge_x"] == 40
        and row["local_cd_nm"] == 40,
        "line is the gap between trench edges",
    )
    require(
        row["left_threshold"] == 87 and row["right_threshold"] == 93,
        "side-specific thresholds follow physical source edges",
    )
    for change in [
        dict(valid=False),
        dict(background_excluded=True),
        dict(sample_y=2.0),
        dict(left_edge_x=np.nan),
        dict(left_edge_x=19.0),
    ]:
        failed = lines.shared_rows([a], [dict(b, **change)], candidate, 1, params)[0]
        require(
            not failed["valid"] and np.isnan(failed["left_edge_x"]),
            "no line from missing, excluded, mismatched or crossing edges",
        )
    require(not lines.shared_rows([a], [], candidate, 1, params), "missing source slot")
    synthetic = lines.shared_rows(
        [dict(a, continuity_synthetic=True)], [b], candidate, 1, params
    )[0]
    require(
        synthetic["continuity_synthetic"], "synthetic source remains explicitly flagged"
    )
    ca = SimpleNamespace(v115_basin_sequence_index=0, center_x_px=100)
    cb = SimpleNamespace(v115_basin_sequence_index=1, center_x_px=260)
    diag = dict(locator_mode_suggested="bright-line", locator_adjacent_period_px=160.0)
    require(lines.adjacent(ca, cb, diag), "neighbor trenches")
    require(
        not lines.adjacent(
            ca, SimpleNamespace(v115_basin_sequence_index=2, center_x_px=420), diag
        ),
        "missing trench must not create a double-width line",
    )
    dark = dict(diag, locator_mode_suggested="dark-line")
    require(
        lines.adjacent(
            ca, SimpleNamespace(v115_basin_sequence_index=2, center_x_px=260), dark
        ),
        "a rejected interior fake stripe is allowed between physical trenches",
    )
    require(
        not lines.adjacent(
            ca, SimpleNamespace(v115_basin_sequence_index=2, center_x_px=420), dark
        ),
        "period distance check rejects missing physical trench in dark-line mode",
    )
    # The max-number limit applies after constructing line pairs.
    candidates = [
        dict(candidate=SimpleNamespace(center_x_px=x)) for x in [20, 40, 60, 80, 100]
    ]
    selected = lines.select_lines(
        candidates, replace(params, center_x_px=60, max_number=3)
    )
    require(
        [r["candidate"].role for r in selected] == ["center", "left", "right"],
        "original centered triplet selection",
    )
    return dict(
        status="PASS",
        exact_edge_copy=True,
        invalid_or_missing_source_rejected=True,
        synthetic_flag_preserved=True,
        skipped_trench_not_bridged=True,
    )


def verify_shared_exports(output, require):
    samples = pd.read_csv(output / "per_sample_results.csv")
    sources = pd.read_csv(output / "line_source_samples.csv")
    pairs = pd.read_csv(output / "line_trench_pairs.csv")
    selected = pairs[pairs.selected]
    require(
        len(selected) > 0 and selected.status.eq("OK").all(),
        "qualified selected line pairs",
    )
    require(
        selected.left_trench_qualified.all() and selected.right_trench_qualified.all(),
        "both source trenches qualified",
    )
    valid = samples[samples.valid]
    require(set(valid.engine) == {"V10", "V13"}, "both engines have shared edges")
    index = sources.set_index(
        ["image_key", "engine", "source_trench_id", "sample_index"]
    )
    for row in valid.to_dict("records"):
        prefix = (row["image_key"], row["engine"])
        a = index.loc[(*prefix, row["left_trench_id"], row["sample_index"])]
        b = index.loc[(*prefix, row["right_trench_id"], row["sample_index"])]
        require(
            a.valid
            and b.valid
            and a.source_trench_qualified
            and b.source_trench_qualified,
            "valid line comes from qualified measured source trenches",
        )
        require(
            row["未旋转_left_edge_x_px"] == a["未旋转_right_edge_x_px"]
            and row["未旋转_right_edge_x_px"] == b["未旋转_left_edge_x_px"],
            "every exported line coordinate equals its source trench edge exactly",
        )
        require(
            row["未旋转_sample_y_px"]
            == a["未旋转_sample_y_px"]
            == b["未旋转_sample_y_px"],
            "same original Y",
        )
        require(
            row["left_threshold"] == a.right_threshold
            and row["right_threshold"] == b.left_threshold,
            "same source thresholds",
        )
        require(not row["background_excluded"], "background excluded")
    return dict(
        valid_shared_samples=len(valid),
        source_trenches=sources.groupby("engine").source_trench_id.nunique().to_dict(),
        selected_lines=selected.groupby("engine").size().to_dict(),
    )


def integration_tests(temp, execute, require, app):
    report = {}
    quiet = ["--skip-psd", "--no-annotated-images"]

    def verify(name, root, maximum=3, extra=(), reference=100, cd=True):
        out, frame = execute(
            name,
            root=root,
            pattern="line",
            reference=reference,
            extra=[*quiet, "--max-number", str(maximum), *extra],
        )
        require(
            frame.status.eq("OK").all(), name + " status: " + str(frame.warning.iloc[0])
        )
        result = verify_shared_exports(out, require)
        require(
            set(result["selected_lines"].values()) == {maximum}, name + " line count"
        )
        summary = pd.read_csv(out / "measurement_summary.csv")
        require(summary.pitch_count.eq(maximum).all(), name + " pitch count")
        if cd:
            require(abs(frame["旋转_mixed_CD_nm"].iloc[0] - 100) < 1.2, name + " CD")
        result.update(
            status="PASS",
            mixed_cd_nm=float(frame["旋转_mixed_CD_nm"].iloc[0]),
            pitch_nm=summary["旋转_pitch_CD_nm"].tolist(),
        )
        report[name] = result
        return out, frame

    root = temp / "shared_heavy_stain_input"
    regions.unicode_imwrite(root / "heavy.png", synthetic_pattern(True, stain_width=76))
    verify("shared_heavy_stain", root)
    verify("shared_heavy_stain_estimated", root, reference=None)
    verify(
        "shared_heavy_stain_optional",
        root,
        extra=["--viterbi", "1", "--erf-fit", "1", "--group-size", "8"],
    )
    verify(
        "shared_asymmetric_thresholds",
        root,
        extra=["--threshold-left", "40", "--threshold-right", "60"],
        cd=False,
    )
    root = temp / "shared_counts_input"
    regions.unicode_imwrite(root / "counts.png", stripe_image())
    for maximum in [1, 2, 5]:
        out, _ = verify("shared_count_" + str(maximum), root, maximum)
        sources = pd.read_csv(out / "line_source_samples.csv")
        require(
            sources.groupby("engine").source_trench_id.nunique().ge(7).all(),
            "search all sources before max-number cap",
        )
    # Exactly two physical trenches bound one line and one complete pitch.
    root = temp / "shared_two_input"
    gray = np.full((400, 400), 210.0)
    for x in [120, 280]:
        gray[:, x - 30 : x + 30] = 40
    regions.unicode_imwrite(root / "two.png", gaussian_filter(gray, 1).astype(np.uint8))
    verify("shared_two_trenches", root, 1)
    root = temp / "shared_one_input"
    gray[:, 250:310] = 210
    regions.unicode_imwrite(root / "one.png", gaussian_filter(gray, 1).astype(np.uint8))
    out, frame = execute(
        "shared_one_trench",
        root=root,
        pattern="line",
        reference=100,
        extra=[*quiet, "--max-number", "1"],
    )
    require(frame["旋转_mixed_CD_nm"].isna().all(), "one trench cannot bound a line")
    # Preserve the existing empty-table export (BOM/newline without a header).
    require(
        not (out / "per_sample_results.csv").read_text(encoding="utf-8-sig").strip(),
        "no invented exterior line",
    )
    report["shared_one_trench"] = dict(
        status="PASS", image_status=frame.status.tolist()
    )
    # Deliberately reject the middle source after actual edge detection. Neither
    # side may pair across it, even though two farther qualified sources remain.
    root = temp / "shared_rejected_input"
    regions.unicode_imwrite(root / "reject.png", synthetic_pattern(False))
    original = app.edges.measure_single_space

    def reject_middle(*args, **kwargs):
        result, rows, annotation = original(*args, **kwargs)
        if args[2].v115_basin_sequence_index == 2:
            result.stable = False
        return result, rows, annotation

    with patch.object(app.edges, "measure_single_space", side_effect=reject_middle):
        out, _ = execute(
            "shared_rejected_middle",
            root=root,
            pattern="line",
            reference=100,
            extra=[*quiet, "--max-number", "2"],
        )
    pairs = pd.read_csv(out / "line_trench_pairs.csv")
    affected = pairs[(pairs.left_trench_id == 2) | (pairs.right_trench_id == 2)]
    require(
        len(affected) == 4
        and affected.status.eq("UNQUALIFIED_TRENCH").all()
        and not affected.selected.any(),
        "unqualified middle trench disqualifies both adjacent lines",
    )
    require(
        not ((pairs.left_trench_id == 1) & (pairs.right_trench_id == 3)).any(),
        "do not bridge rejected trench",
    )
    report["shared_rejected_middle"] = verify_shared_exports(out, require)
    # Physical middle trench absent: the wide gap must not become a line.
    root = temp / "shared_missing_input"
    gray = synthetic_pattern(False, noise=0)
    gray[80:520, 120 + 380 : 120 + 460] = 210
    regions.unicode_imwrite(root / "missing.png", gray)
    out, _ = execute(
        "shared_missing_middle",
        root=root,
        pattern="line",
        reference=100,
        extra=[*quiet, "--max-number", "2"],
    )
    pairs = pd.read_csv(out / "line_trench_pairs.csv")
    require(
        pairs.status.eq("NONADJACENT").sum() == 2,
        "missing physical trench causes rejected wide gap for each engine",
    )
    report["shared_missing_middle"] = verify_shared_exports(out, require)
    # The PSD pass is separately measured at one-row resolution; its annotation
    # arrays must retain exactly the line sample coordinates as well.
    root = temp / "shared_psd_input"
    regions.unicode_imwrite(root / "psd.png", synthetic_pattern(True))
    real_measure = app.edges.measure_image
    checked = []

    def inspect_measure(*args, **kwargs):
        measured = real_measure(*args, **kwargs)
        p = args[1]
        if (
            p.pattern_kind == "line"
            and p.average_range_px == 1
            and p.smoothing_pixel == 1
        ):
            for ann in measured.annotation_payload["spaces"]:
                rows = [
                    r
                    for r in measured.sample_rows
                    if r["space_index"] == ann["space_index"]
                ]
                np.testing.assert_array_equal(
                    ann["left_edges"],
                    [
                        r["line_left_trench_right_x_px"] if r["valid"] else np.nan
                        for r in rows
                    ],
                )
                np.testing.assert_array_equal(
                    ann["right_edges"],
                    [
                        r["line_right_trench_left_x_px"] if r["valid"] else np.nan
                        for r in rows
                    ],
                )
            checked.append(p.engine_kind)
        return measured

    with patch.object(app.edges, "measure_image", side_effect=inspect_measure):
        out, frame = execute(
            "shared_psd",
            root=root,
            pattern="line",
            reference=100,
            extra=["--max-number", "3", "--no-annotated-images"],
        )
    require(
        frame.status.eq("OK").all() and set(checked) == {"V10", "V13"},
        "both independent PSD passes share source edges",
    )
    psd = pd.read_csv(out / "PSD/per_structure_psd_summary.csv")
    require(
        psd.status.eq("OK").all() and psd.average_range_px.eq(1).all(),
        "PSD output and single-row policy retained",
    )
    report["shared_psd"] = dict(
        status="PASS", shared_psd_engines=checked, curves=len(psd)
    )
    return report
