"""Regressions for asymmetric line failures and measured, rotated full periods."""

from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.ndimage import gaussian_filter

import cdsem_pitch as pitch
import cdsem_regions as regions
from check_localization import synthetic_pattern


PUBLIC_PITCH = "旋转_pitch_CD_nm"


def analytical_tests(require, app):
    params = app.edges.MeasurementParams(
        pixel_size_nm=1.7, max_number=3, center_x_px=340
    )
    candidates, groups = [], {}
    for index, x in enumerate([100, 240, 410, 600]):
        candidates.append(
            SimpleNamespace(center_x_px=x, v115_basin_sequence_index=index)
        )
        groups[index] = [
            dict(
                sample_index=i,
                sample_y=float(i),
                valid=True,
                left_edge_x=x + 0.2 * i,
                right_edge_x=x + 60 + 0.2 * i,
            )
            for i in range(80)
        ]
    summary, rows, samples = pitch.pitch_from_samples(
        groups, candidates, params, "bright-line"
    )
    expected = np.mean([140, 170, 190]) * 1.7 / np.sqrt(1.04)
    require(summary["pitch_count"] == 3, "N+1 boundaries must give N full periods")
    require(
        abs(summary[pitch.PITCH_COLUMN] - expected) < 1e-10,
        "PCA rotated physical pitch",
    )
    require(all(r["selected"] for r in rows), "distinct periods, no double counting")
    one, _, _ = pitch.pitch_from_samples(
        groups, candidates, replace(params, max_number=1), "bright-line"
    )
    require(
        one["pitch_count"] == 1
        and abs(one[pitch.PITCH_COLUMN] - 170 * 1.7 / np.sqrt(1.04)) < 1e-10,
        "max-number=1 must still measure one whole period",
    )
    # Missing Y and synthetic data may not enter pitch, even under forced continuity.
    partial = deepcopy(groups)
    for row in partial[0][:20]:
        row.update(background_excluded=True, left_edge_x=-10000)
    for row in partial[3][-20:]:
        row.update(continuity_synthetic=True, left_edge_x=10000)
    measured, periods, audit = pitch.pitch_from_samples(
        partial, candidates, params, "bright-line"
    )
    require(
        abs(measured[pitch.PITCH_COLUMN] - expected) < 1e-10,
        "equal period means, not sample-weighted; excluded samples cannot bias pitch",
    )
    require(
        not any(
            r["valid"]
            for r in audit
            if r["background_excluded"] or r["synthetic_excluded"]
        ),
        "background and synthetic exclusions",
    )
    skipped, skipped_rows, _ = pitch.pitch_from_samples(
        groups, [candidates[0], candidates[2]], params, "bright-line"
    )
    require(
        skipped["pitch_count"] == 0
        and np.isnan(skipped[pitch.PITCH_COLUMN])
        and skipped_rows[0]["status"] == "NONADJACENT",
        "do not report two periods as one",
    )
    none, _, _ = pitch.pitch_from_samples({}, [], params, "bright-line")
    require(
        none["pitch_status"] == "UNAVAILABLE" and np.isnan(none[pitch.PITCH_COLUMN]),
        "missing is NaN, never zero",
    )
    image = synthetic_pattern(True)
    trench, gap = regions.estimate_trench_references_nm(255 - image, 1)
    require(
        abs(trench - 100) <= 5 and abs(gap - 60) <= 5,
        "automatic line reference must include the interior stain",
    )
    return dict(
        status="PASS",
        rotated_pitch_nm=expected,
        unequal_periods_px=[140, 170, 190],
        rotation_slope=0.2,
        pixel_size_nm=1.7,
        exclusion_and_adjacency="PASS",
    )


def stripe_image(line_width=100, count=7, slope=0, noise=0):
    period = 60 + line_width
    width = period * count + 80
    y, x = np.indices((400, width))
    position = (x - slope * (y - 200) - 40) % period
    image = np.where(position < 60, 40.0, 210.0)
    image = gaussian_filter(image, 1) + np.random.default_rng(92).normal(
        0, noise, image.shape
    )
    padded = np.full((520, width + 120), 155, np.uint8)
    padded[60:460, 60 : width + 60] = np.clip(image, 0, 255).astype(np.uint8)
    return padded


def integration_tests(temp, execute, require):
    report = {}
    quiet = ["--skip-psd", "--no-annotated-images"]

    def verify(
        name,
        root,
        pattern,
        reference,
        expected_pitch,
        expected_cd,
        maximum=3,
        extra=(),
        status="OK",
        count=None,
    ):
        out, image = execute(
            name,
            root=root,
            pattern=pattern,
            reference=reference,
            extra=[*quiet, "--max-number", str(maximum), *extra],
        )
        require(
            image.status.eq(status).all(),
            name + " image status: " + str(image.warning.iloc[0]),
        )
        require(
            abs(image["旋转_mixed_CD_nm"].iloc[0] - expected_cd) < 1.2, name + " CD"
        )
        summary = pd.read_csv(out / "measurement_summary.csv")
        require(summary.columns[-1] == PUBLIC_PITCH, name + " CSV pitch is last column")
        require(
            not any("未旋转" in c and "pitch" in c for c in summary),
            "rotated pitch only",
        )
        require(
            summary.pitch_count.eq(maximum if count is None else count).all(),
            name + " max-number periods",
        )
        require(
            np.allclose(summary[PUBLIC_PITCH], expected_pitch, atol=1.2),
            name + " measured period",
        )
        require(
            summary.loc[summary.method == "mixed", PUBLIC_PITCH].iloc[0]
            == summary.loc[summary.method == "V13", PUBLIC_PITCH].iloc[0],
            "mixed uses V13 pitch",
        )
        book = load_workbook(
            out / "CD_measurement_200K_V1_15_results.xlsx", read_only=True
        )
        for sheet in [
            "measurement_summary",
            "image_summary",
            "pitch_periods",
            "pitch_samples",
        ]:
            header = next(book[sheet].values)
            require(
                header[-1] == PUBLIC_PITCH,
                name + " Excel pitch last column in " + sheet,
            )
        book.close()
        periods = pd.read_csv(out / "pitch_periods.csv")
        for engine in ["V10", "V13"]:
            selected = periods[(periods.engine == engine) & periods.selected]
            value = summary.loc[summary.method == engine, PUBLIC_PITCH].iloc[0]
            require(
                abs(value - selected[PUBLIC_PITCH].mean()) < 1e-10,
                "equal per-period average",
            )
            require(
                not selected.duplicated(
                    ["left_basin_index", "right_basin_index"]
                ).any(),
                "unique periods",
            )
        report[name] = dict(
            status=status,
            pitch_nm=summary[PUBLIC_PITCH].tolist(),
            pitch_count=summary.pitch_count.tolist(),
            mixed_cd_nm=float(image["旋转_mixed_CD_nm"].iloc[0]),
        )
        return out, image

    root = temp / "wide_line_input"
    regions.unicode_imwrite(root / "wide.png", stripe_image(180))
    for pattern, reference in [("trench", 60), ("line", 180)]:
        verify("wide_" + pattern, root, pattern, reference, 240, reference)
    # Same physical image, not an inverted fixture: the original suite missed this case.
    root = temp / "stained_line_input"
    original = synthetic_pattern(True)
    regions.unicode_imwrite(root / "stain.png", original)
    verify("stained_trench_pitch", root, "trench", 60, 160, 60)
    _, line = verify("stained_line_pitch", root, "line", 100, 160, 100)
    verify("stained_line_estimated", root, "line", None, 160, 100)
    verify(
        "stained_line_manual",
        root,
        "line",
        100,
        160,
        100,
        extra=["--locator-mode", "bright-line"],
    )
    inverted = temp / "inverted_stained_input"
    regions.unicode_imwrite(inverted / "stain.png", 255 - original)
    _, trench = verify("inverted_stained_trench", inverted, "trench", 100, 160, 100)
    for metric in app_result_columns():
        require(
            np.allclose(line[metric], trench[metric], atol=1e-9),
            "exact polarity symmetry " + metric,
        )
    root = temp / "tilted_pitch_input"
    slope = 0.08
    regions.unicode_imwrite(root / "tilted.png", stripe_image(slope=slope, noise=2))
    for pattern, reference in [("trench", 60), ("line", 100)]:
        verify(
            "tilted_pitch_" + pattern,
            root,
            pattern,
            reference,
            160 / np.sqrt(1 + slope * slope),
            reference / np.sqrt(1 + slope * slope),
        )
    root = temp / "pitch_count_input"
    regions.unicode_imwrite(root / "counts.png", stripe_image())
    for maximum in [1, 2, 5]:
        verify(f"pitch_count_{maximum}", root, "trench", 60, 160, 60, maximum=maximum)
    # Only three trenches, hence only two *complete* trench+line periods.
    narrow = np.full((400, 560), 210.0, dtype=float)
    for x in [100, 260, 420]:
        narrow[:, x - 30 : x + 30] = 40
    root = temp / "pitch_short_input"
    regions.unicode_imwrite(
        root / "short.png", gaussian_filter(narrow, 1).astype(np.uint8)
    )
    verify("pitch_insufficient", root, "trench", 60, 160, 60, count=2, status="REVIEW")
    return report


def app_result_columns():
    return [
        f"{mode}_{engine}_{metric}_nm"
        for mode in ["旋转", "未旋转"]
        for engine in ["mixed", "V10", "V13"]
        for metric in ["CD", "LER_left", "LER_right", "LWR"]
    ]
