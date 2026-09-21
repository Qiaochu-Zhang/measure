"""Measured full periods: one target band plus its immediately adjacent gap."""

from dataclasses import replace
import math
import numpy as np

import cdsem_engine as edges


PITCH_COLUMN = "rotated_pitch_CD_nm"
PERIOD_COLUMNS = [
    "image_key",
    "engine",
    "pitch_index",
    "left_basin_index",
    "right_basin_index",
    "center_x_px",
    "selected",
    "status",
    "valid_sample_count",
    "eligible_sample_count",
    "rotation_angle_deg",
    PITCH_COLUMN,
]
SAMPLE_COLUMNS = [
    "image_key",
    "engine",
    "pitch_index",
    "sample_index",
    "sample_y",
    "start_edge_x_px",
    "end_edge_x_px",
    "valid",
    "background_excluded",
    "synthetic_excluded",
    PITCH_COLUMN,
]


def summarize_periods(periods, maximum, center_x):
    """Select distinct complete periods nearest the requested center, equally weighted."""
    valid = [r for r in periods if r["status"] == "OK"]
    chosen = sorted(
        valid, key=lambda r: (abs(r["center_x_px"] - center_x), r["center_x_px"])
    )[:maximum]
    for row in periods:
        row["selected"] = any(row is chosen_row for chosen_row in chosen)
    return dict(
        pitch_count=len(chosen),
        pitch_requested_count=maximum,
        pitch_status="OK"
        if len(chosen) == maximum
        else "PARTIAL"
        if chosen
        else "UNAVAILABLE",
        rotated_pitch_CD_nm=float(np.mean([r[PITCH_COLUMN] for r in chosen]))
        if chosen
        else math.nan,
    )


def pitch_from_samples(samples_by_basin, candidates, params, mode, period_px=None):
    """Use same-side edges at common Y, not distances between unlike centers.

    [left edge of band i, left edge of band i+1] contains exactly one band
    and one adjacent gap. PCA of the midpoint fits a shared axis for the
    period; horizontal edge distances are projected onto its normal.
    """
    periods, audit = [], []
    ordered = sorted(candidates, key=lambda c: c.center_x_px)
    for index, (left, right) in enumerate(zip(ordered, ordered[1:]), 1):
        li, ri = left.v115_basin_sequence_index, right.v115_basin_sequence_index
        center = (left.center_x_px + right.center_x_px) / 2
        record = dict(
            pitch_index=index,
            left_basin_index=li,
            right_basin_index=ri,
            center_x_px=center,
            selected=False,
            status="INSUFFICIENT_SUPPORT",
            valid_sample_count=0,
            eligible_sample_count=0,
            rotation_angle_deg=math.nan,
            rotated_pitch_CD_nm=math.nan,
        )
        periods.append(record)
        # Never span a missing structure, even if both distant edges are valid.
        gap = ri - li
        maximum_gap = 2 if mode == "dark-line" else 1
        distance = right.center_x_px - left.center_x_px
        if not 0 < gap <= maximum_gap or (
            period_px is not None
            and np.isfinite(period_px)
            and period_px > 0
            and not 0.60 * period_px <= distance <= 1.40 * period_px
        ):
            record["status"] = "NONADJACENT"
            continue
        a = {r["sample_index"]: r for r in samples_by_basin.get(li, [])}
        b = {r["sample_index"]: r for r in samples_by_basin.get(ri, [])}
        rows = []
        for sample_index in sorted(a.keys() & b.keys()):
            first, second = a[sample_index], b[sample_index]
            y, other_y = float(first["sample_y"]), float(second["sample_y"])
            start = float(first.get("left_edge_x", math.nan))
            end = float(second.get("left_edge_x", math.nan))
            boundary = float(first.get("right_edge_x", math.nan))
            excluded = bool(
                first.get("background_excluded") or second.get("background_excluded")
            )
            synthetic = bool(
                first.get("continuity_synthetic") or second.get("continuity_synthetic")
            )
            valid = (
                bool(first.get("valid"))
                and bool(second.get("valid"))
                and not excluded
                and not synthetic
                and y == other_y
                and np.isfinite([start, end, boundary, y]).all()
                and start < boundary < end
            )
            rows.append(
                dict(
                    pitch_index=index,
                    sample_index=sample_index,
                    sample_y=y,
                    start_edge_x_px=start,
                    end_edge_x_px=end,
                    valid=bool(valid),
                    background_excluded=excluded,
                    synthetic_excluded=synthetic,
                    rotated_pitch_CD_nm=math.nan,
                )
            )
        eligible = sum(not r["background_excluded"] for r in rows)
        valid_rows = [r for r in rows if r["valid"]]
        record.update(
            valid_sample_count=len(valid_rows), eligible_sample_count=eligible
        )
        if (
            len(valid_rows) >= 2
            and len(valid_rows) / max(1, eligible) >= params.minimum_valid_fraction
        ):
            y = np.array([r["sample_y"] for r in valid_rows])
            x0 = np.array([r["start_edge_x_px"] for r in valid_rows])
            x1 = np.array([r["end_edge_x_px"] for r in valid_rows])
            fit = edges.fit_centerline_and_width_factor(y, x0, x1)
            values = (x1 - x0) * fit["normal_width_factor"] * params.pixel_size_nm
            for row, value in zip(valid_rows, values):
                row[PITCH_COLUMN] = float(value)
            record.update(
                status="OK",
                rotation_angle_deg=fit["angle_from_vertical_deg"],
                rotated_pitch_CD_nm=float(np.mean(values)),
            )
        audit.extend(rows)
    return (
        summarize_periods(periods, params.max_number, params.center_x_px),
        periods,
        audit,
    )


def measure_pitch(measurement, params, region_mask=None):
    """Reuse primary edges; measure one extra neighbouring band for N full periods.

    The extra band never enters the primary CD/LER/LWR or PSD statistics.
    All extraction uses the same thresholds, ROI and end/background exclusions.
    """
    if params.pattern_kind == "line":
        data = measurement.annotation_payload["line_pitch_data"]
        return pitch_from_samples(
            data["samples"], data["candidates"], params, data["mode"], data["period_px"]
        )
    extended = replace(
        params,
        max_number=params.max_number + 1,
        min_number=1,
        require_center_left_right=False,
    )
    image = measurement.annotation_payload["image"]
    candidates, diag = edges.detect_space_candidates(image, extended, region_mask)
    existing = {}
    for item in measurement.annotation_payload["spaces"]:
        candidate = item["candidate"]
        existing[candidate.v115_basin_sequence_index] = [
            row
            for row in measurement.sample_rows
            if row["space_index"] == item["space_index"]
        ]
    groups = {}
    for index, candidate in enumerate(candidates, 1):
        key = candidate.v115_basin_sequence_index
        if key in existing:
            groups[key] = existing[key]
        else:
            _, groups[key], _ = edges.measure_single_space(
                measurement.image_row["image_name"],
                image,
                candidate,
                index,
                extended,
                region_mask,
            )
    return pitch_from_samples(
        groups,
        candidates,
        params,
        diag["locator_mode_suggested"],
        diag.get("locator_adjacent_period_px"),
    )
