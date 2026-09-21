#!/usr/bin/env python3


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V1_17: background-aware dual-edge measurement and unaveraged PSD.

Summary: filename/status/method first, mixed before V13/V10. CD/LER are
ungrouped; scalar LWR retains --group-size. PSD never groups edge samples.
"""

from __future__ import annotations
import argparse
import json
import math
import re
import sys
import traceback
from dataclasses import asdict, replace
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Iterable, Optional, Sequence
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
import cdsem_regions as regions
import cdsem_statistics as statistics
import cdsem_engine as edges
import cdsem_localization as localization
from cdsem_psd import PSDBatch
import cdsem_pitch as pitch


SCRIPT_VERSION = "V1_17"
DEFAULT_ROOT = regions.DEFAULT_ROOT
DEFAULT_OUTPUT_NAME = "CD_measure_output_200K_V1_17"
EPS = 1e-12
PUBLIC_GROUP_SIZE = 4
V113_METHOD_CODE = "V1_17_DualCoordinates_V13CD_V10LER_V13LWR"
V113_RESULT_ID = "dual_coordinates__v13_cd__v10_ler__v13_group4_lwr"
V113_STAT_MODE = "rotated_then_unrotated"
V113_AGGREGATION = "mean_per_space"
V10_THRESHOLD_MODE = "v10_center40_dark_median_peak3_bright_median"
V13_THRESHOLD_MODE = "v17_side_band_mean"
V13_ALGORITHM_VERSION = "V13_V17_THRESHOLD_V10_STATISTICS"
METRICS = ("CD", "LER_left", "LER_right", "LWR")
GROUPS = ("mixed", "V13", "V10")
MODES = ("rotated", "unrotated")
OUTPUT_SCHEMA = "metadata_first_methods_raw_psd_v117_shared_lines_pitch"
PATCH_VERSION = "V1_17"
RESULT_COLUMNS = [
    f"{mode}_{group}_{metric}_nm"
    for mode in MODES
    for group in GROUPS
    for metric in METRICS
]
META_COLUMNS = (
    "image_key",
    "image",
    "folder_name",
    "path",
    "pattern",
    "pixel_size_nm",
    "design_trench_nm",
    "design_space_nm",
    "status",
    "warning",
)


def finite_float(value: Any, default: float = math.nan) -> float:
    try:
        result = float(value)
        return result if np.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def safe_mean(values: Iterable[Any]) -> float:
    array = np.asarray([finite_float(value) for value in values], dtype=float)
    array = array[np.isfinite(array)]
    return float(np.mean(array)) if array.size else math.nan


def safe_std(values: Iterable[Any], ddof: int = 1) -> float:
    array = np.asarray([finite_float(value) for value in values], dtype=float)
    array = array[np.isfinite(array)]
    return float(np.std(array, ddof=ddof)) if array.size > ddof else math.nan


def public_text(value):
    if isinstance(value, str):
        return value.replace("group4", f"group{PUBLIC_GROUP_SIZE}").replace(
            "Group4", f"Group{PUBLIC_GROUP_SIZE}"
        )
    if isinstance(value, dict):
        return {public_text(k): public_text(v) for k, v in value.items()}
    if isinstance(value, list):
        return [public_text(v) for v in value]
    return value


def public_frame(frame):
    out = frame.rename(columns={c: public_text(c) for c in frame.columns}).copy()
    for col in out.select_dtypes(include=["object", "string"]).columns:
        out[col] = out[col].map(public_text)
    return label_frame(out)


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame = public_frame(frame)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(
            mode="w",
            encoding="utf-8-sig",
            newline="",
            dir=path.parent,
            suffix=".csv.tmp",
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            frame.to_csv(stream, index=False)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def dataframe_or_empty(
    rows: Sequence[dict[str, Any]], columns: Sequence[str] = ()
) -> pd.DataFrame:
    if rows:
        return pd.DataFrame(rows)
    return pd.DataFrame(columns=list(columns))


def jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def relative_key(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.name


def metadata_for(
    parsed: regions.ParsedImage, general: regions.GeneralConfig
) -> dict[str, Any]:
    return {
        "folder_name": parsed.folder_name,
        "folder_nominal_nm": parsed.folder_nominal_nm,
        "pretreatment_nm": parsed.pretreatment_nm,
        "dose_mj": parsed.dose_mj,
        "sample_id": parsed.sample_id,
        "metadata_source": parsed.metadata_source,
        "image": parsed.path.name,
        "image_key": relative_key(parsed.path, general.root_dir),
        "path": str(parsed.path),
        "pattern": general.force_pattern,
        "pixel_size_nm": general.pixel_size_nm,
        "design_trench_nm": parsed.design_trench_nm,
        "design_space_nm": parsed.design_space_nm,
    }


def build_v10_params(
    parsed: regions.ParsedImage,
    gray: np.ndarray,
    general: regions.GeneralConfig,
    args: argparse.Namespace,
) -> edges.MeasurementParams:
    """Create one V10 parameter set while allowing arbitrary PNG dimensions."""
    height, width = gray.shape
    center_x = (
        float(args.center_x) if args.center_x is not None else (width - 1.0) / 2.0
    )
    center_y = (
        float(args.center_y) if args.center_y is not None else (height - 1.0) / 2.0
    )
    average_range = int(args.average_range)
    smoothing = int(args.smoothing_pixel)
    if not (0 <= center_x < width and 0 <= center_y < height):
        raise ValueError("测量中心必须在图像内")
    if average_range > height or smoothing > width:
        raise ValueError(
            "average-range不能超过图像高度；smoothing-pixel不能超过图像宽度"
        )
    extend_length = (
        int(args.extend_length) if args.extend_length is not None else min(360, height)
    )
    if not 2 <= extend_length <= height:
        raise ValueError("extend-length必须在2至图像高度之间，不再静默限幅")
    measurement_width = (
        args.meas_area_width if args.meas_area_width is not None else min(512, width)
    )
    measurement_height = (
        args.meas_area_height if args.meas_area_height is not None else min(512, height)
    )
    if measurement_width > width or measurement_height > height:
        raise ValueError("测量ROI尺寸不能超过图像；调整meas-area-width/height")

    target_cd = finite_float(parsed.design_trench_nm)
    if not np.isfinite(target_cd) or target_cd <= 0:
        raise ValueError("无法取得有效 Trench 参考宽度")

    return edges.MeasurementParams(
        engine_kind="V10",
        locator_mode=args.locator_mode,
        locator_majority=args.locator_majority,
        threshold_search=args.threshold_search,
        expected_width_px=width,
        expected_height_px=height,
        pixel_size_nm=general.pixel_size_nm,
        target_cd_nm=target_cd,
        center_x_px=center_x,
        center_y_px=center_y,
        search_range_in_px=int(args.search_in),
        search_range_out_px=int(args.search_out),
        extend_length_px=extend_length,
        peak_order=int(args.peak_order),
        edge_threshold_left_pct=float(args.threshold_left),
        edge_threshold_right_pct=float(args.threshold_right),
        smoothing_pixel=smoothing,
        sample_number=int(args.sample_number),
        average_range_px=average_range,
        topology_min_contrast_8bit=float(args.topology_min_contrast),
        tracking_search_radius_px=float(args.tracking_radius),
        tracking_retry_radius_px=float(args.tracking_retry_radius),
        tracking_max_edge_jump_px=float(args.tracking_max_jump),
        postfilter_window=int(args.postfilter_window),
        postfilter_mad_multiplier=float(args.postfilter_mad_multiplier),
        postfilter_min_edge_tolerance_px=float(args.postfilter_edge_tolerance),
        recover_rejected_points=not bool(args.no_recovery),
        edge_continuity=int(args.edge_continuity),
        pattern_kind=args.pattern,
        space_reference_nm=(
            parsed.design_space_nm
            if args.pattern == "line"
            else args.space_reference_nm
        ),
        viterbi_enabled=args.viterbi,
        viterbi_weight=args.viterbi_weight,
        viterbi_max_jump_px=args.viterbi_max_jump,
        viterbi_gap_cost=args.viterbi_gap_cost,
        viterbi_candidates=args.viterbi_candidates,
        erf_enabled=args.erf_fit,
        erf_window_px=args.erf_window,
        erf_max_shift_px=args.erf_max_shift,
        erf_max_relative_rmse=args.erf_max_relative_rmse,
        max_number=int(args.max_number),
        min_number=int(args.min_number),
        candidate_dark_mean_tolerance_8bit=args.candidate_dark_tolerance,
        candidate_min_basin_contrast_8bit=args.candidate_min_contrast,
        require_center_left_right=(
            args.max_number >= 3 and not bool(args.allow_incomplete_triplet)
        ),
        meas_area_width_px=measurement_width,
        meas_area_height_px=measurement_height,
        minimum_valid_fraction=float(args.minimum_valid_fraction),
        flyer_mode=str(args.flyer_mode),
        flyer_left_offset_px=float(args.flyer_left_offset),
        flyer_right_offset_px=float(args.flyer_right_offset),
    )


def build_v13_params(params: edges.MeasurementParams) -> edges.MeasurementParams:
    """Clone every V1.8/V10 public parameter into the V13 edge engine."""
    return replace(params, engine_kind="V13")


def set_measurement_image_key(measurement: Any, image_key: str) -> None:
    measurement.image_row["image_name"] = image_key
    for row in measurement.space_rows:
        row["image_name"] = image_key
    for row in measurement.sample_rows:
        row["image_name"] = image_key
    for row in measurement.annotation_payload.get("candidate_diag", {}).get(
        "candidate_debug_rows", []
    ):
        row["image_name"] = image_key


def run_coordinate_engine(
    engine, params, parsed, key, group_size, meta, errors, region_mask=None
):
    """Keep each engine/structure independent; never require legacy result rows."""
    result = dict(
        measurement=None,
        samples=[],
        metrics={},
        selected=[],
        aggregate={},
        base={},
        issues=[],
        selection_fallback=False,
    )

    def record_failure(stage, exc):
        message = f"{engine} {stage}: {type(exc).__name__}: {exc}"
        result["issues"].append(message)
        errors.append(
            {
                **meta,
                "stage": f"{engine}_{stage}",
                "error_type": type(exc).__name__,
                "error_message": message,
                "traceback": traceback.format_exc(),
            }
        )

    module = edges if engine == "V10" else edges
    try:
        measurement = module.measure_image(parsed.path, params, region_mask=region_mask)
        result["measurement"] = measurement
        set_measurement_image_key(measurement, key)
        result["base"] = dict(measurement.image_row)
        result["samples"] = [dict(row) for row in measurement.sample_rows]
    except Exception as exc:
        record_failure("measurement", exc)
        return result

    groups = {}
    for sample in result["samples"]:
        groups.setdefault(int(sample["space_index"]), []).append(sample)
    for index, samples in groups.items():
        try:
            result["metrics"].update(
                coordinate_metrics_by_space(samples, params, group_size)
            )
        except Exception as exc:
            record_failure(f"structure_{index}_coordinates", exc)
    # Preserve V1.13 selection when it yields usable values. If all chosen
    # structures are unavailable, allow partial measured structures under REVIEW.
    metrics = result["metrics"]
    available = [
        index
        for index, row in metrics.items()
        if any(
            np.isfinite(value)
            for mode in MODES
            for value in engine_values(row, mode).values()
        )
    ]
    selected = statistics.choose_aggregate_space_indices(
        measurement.space_rows, params.min_number
    )
    selected = [index for index in selected if index in available]
    if not selected and available:
        selected = available
        result["selection_fallback"] = True
        result["issues"].append(
            f"{engine}: 原筛选未留下可用结构，使用有有效指标的结构并标为 REVIEW"
        )
    result["selected"] = selected
    result["aggregate"] = aggregate_coordinate_metrics(metrics, selected)
    if not available:
        result["issues"].append(
            f"{engine}: 没有可用 CD/LER/LWR；候选数={len(measurement.space_rows)}，采样数={len(result['samples'])}"
        )
    return result


def primary_from_engines(meta, r10, r13):
    """Output each finite metric; completeness is a quality flag, never a gate."""
    primary = corrected_image_row(meta, r10["aggregate"], r13["aggregate"])
    warnings = []
    quality_ok = True
    for engine, result in (("V10", r10), ("V13", r13)):
        base = result["base"]
        prefix = engine.lower()
        primary[f"{prefix}_measurement_valid"] = bool(base.get("valid", False))
        primary[f"{prefix}_triplet_complete"] = bool(
            base.get("triplet_complete", False)
        )
        primary[f"{prefix}_selected_roles"] = base.get("selected_roles", "")
        primary[f"{prefix}_selected_trench_count"] = base.get("selected_space_count", 0)
        primary[f"{prefix}_stable_trench_count"] = base.get("stable_space_count", 0)
        primary[f"{engine}_selection_fallback"] = result["selection_fallback"]
        missing = [
            column_label(f"{mode}_{engine}_{metric}_nm")
            for mode in MODES
            for metric in METRICS
            if not np.isfinite(primary[f"{mode}_{engine}_{metric}_nm"])
        ]
        primary[f"{engine}_missing_metrics"] = "; ".join(missing)
        primary[f"{engine}_issues"] = " | ".join(result["issues"])
        if base.get("warning"):
            warnings.append(f"{engine}: {base['warning']}")
        warnings.extend(result["issues"])
        if missing:
            warnings.append("缺项: " + "; ".join(missing))
        if not primary[f"{prefix}_measurement_valid"] or result["issues"] or missing:
            quality_ok = False
    primary["rotated_mixed_complete"] = all(
        np.isfinite(primary[f"rotated_mixed_{m}_nm"]) for m in METRICS
    )
    primary["all_result_groups_complete"] = all(
        np.isfinite(primary[c]) for c in RESULT_COLUMNS
    )
    any_result = any(np.isfinite(primary[c]) for c in RESULT_COLUMNS)
    primary["status"] = ("OK" if quality_ok else "REVIEW") if any_result else "ERROR"
    primary["warning"] = " | ".join(warnings)
    primary["patch_version"] = PATCH_VERSION
    return primary


def object_keys_from_metrics(c10, c13, meta):
    """Outer join unique roles; unmatched/ambiguous objects retain their own values."""
    roles10, roles13 = {}, {}
    for metrics, roles in ((c10, roles10), (c13, roles13)):
        for index, item in metrics.items():
            role = str(item.get("space_role", "")).strip()
            roles.setdefault(role, []).append(index)
    rows, used10 = [], set()
    for i13, item in c13.items():
        role = str(item.get("space_role", "")).strip()
        i10 = None
        if role and len(roles10.get(role, [])) == 1 and len(roles13[role]) == 1:
            i10 = roles10[role][0]
        elif (
            not role and i13 in c10 and not str(c10[i13].get("space_role", "")).strip()
        ):
            i10 = i13
        if i10 is not None:
            used10.add(i10)
        rows.append(
            {
                **meta,
                "trench_id": i13,
                "space_role": role,
                "v10_space_index": i10,
                "v13_space_index": i13,
            }
        )
    for i10, item in c10.items():
        if i10 not in used10:
            rows.append(
                {
                    **meta,
                    "trench_id": i10,
                    "space_role": item.get("space_role", ""),
                    "v10_space_index": i10,
                    "v13_space_index": None,
                }
            )
    return rows


def excel_column_index(label: str) -> int:
    text = str(label).strip().upper()
    if not re.fullmatch(r"[A-Z]+", text):
        raise ValueError(f"无效 Excel 列名：{label}")
    value = 0
    for char in text:
        value = value * 26 + ord(char) - ord("A") + 1
    return value - 1


def normalize_match_key(value: Any) -> str:
    name = Path(str(value).strip()).name
    stem = Path(name).stem
    stem = re.sub(r"\(\d+\)$", "", stem).strip()
    return stem.casefold()


def autofit_worksheet(worksheet, max_width: int = 42) -> None:
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True)
    worksheet.freeze_panes = "A2"
    if worksheet.max_row and worksheet.max_column:
        worksheet.auto_filter.ref = worksheet.dimensions
    for cell in worksheet[1]:
        cell.fill = header_fill
        cell.font = header_font
    for column in range(1, worksheet.max_column + 1):
        maximum = 0
        for row in range(1, min(worksheet.max_row, 250) + 1):
            value = worksheet.cell(row=row, column=column).value
            maximum = max(maximum, len(str(value)) if value is not None else 0)
        worksheet.column_dimensions[get_column_letter(column)].width = min(
            max(maximum + 2, 10), max_width
        )
    for row in worksheet.iter_rows():
        for cell in row:
            cell.alignment = Alignment(vertical="top")


def coordinate_metrics_by_space(
    sample_rows: list[dict[str, Any]],
    params: Any,
    group_size: int,
) -> dict[int, dict[str, Any]]:
    """Calculate unrotated and true centerline-rotated metrics per trench.

    The fitted centerline direction is Y'.  Its normal is X'.  Rotating the
    edge coordinates makes LER a 3sigma fluctuation of X', so a perfectly
    straight but tilted edge no longer receives artificial LER from its slope.
    ``sample_rows`` is updated with rotated coordinates for full auditability.
    """
    groups: dict[int, list[dict[str, Any]]] = {}
    for row in sample_rows:
        groups.setdefault(int(row["space_index"]), []).append(row)

    output: dict[int, dict[str, Any]] = {}
    for space_index, rows in groups.items():
        rows.sort(key=lambda item: int(item["sample_index"]))
        y = np.asarray(
            [finite_float(item.get("sample_y")) for item in rows], dtype=float
        )
        left = np.asarray(
            [finite_float(item.get("left_edge_x")) for item in rows], dtype=float
        )
        right = np.asarray(
            [finite_float(item.get("right_edge_x")) for item in rows], dtype=float
        )
        valid_flags = np.asarray([bool(item.get("valid", False)) for item in rows])
        valid = (
            valid_flags
            & np.isfinite(y)
            & np.isfinite(left)
            & np.isfinite(right)
            & (right > left)
        )
        left[~valid] = np.nan
        right[~valid] = np.nan

        fit = edges.fit_centerline_and_width_factor(y, left, right)
        angle_deg = finite_float(fit.get("angle_from_vertical_deg"), 0.0)
        theta = math.radians(angle_deg)
        cos_theta = math.cos(theta)
        sin_theta = math.sin(theta)

        midpoint = 0.5 * (left + right)
        origin_x = safe_mean(midpoint)
        origin_y = safe_mean(y[valid])
        if not np.isfinite(origin_x):
            origin_x = 0.0
        if not np.isfinite(origin_y):
            origin_y = 0.0

        # Rotation about the mean center. Translation does not affect standard
        # deviations or widths, but keeps audit coordinates numerically compact.
        left_prime = cos_theta * (left - origin_x) - sin_theta * (y - origin_y)
        right_prime = cos_theta * (right - origin_x) - sin_theta * (y - origin_y)
        center_y_prime = sin_theta * (midpoint - origin_x) + cos_theta * (y - origin_y)

        for index, row in enumerate(rows):
            row["unrotated_local_cd_nm"] = (
                finite_float((right[index] - left[index]) * params.pixel_size_nm)
                if valid[index]
                else math.nan
            )
            row["rotated_axis_angle_from_vertical_deg"] = angle_deg
            row["rotated_axis_normal_factor"] = abs(cos_theta)
            row["rotated_left_x_prime_px"] = finite_float(left_prime[index])
            row["rotated_right_x_prime_px"] = finite_float(right_prime[index])
            row["rotated_center_y_prime_px"] = finite_float(center_y_prime[index])
            row["rotated_local_cd_nm"] = (
                finite_float(
                    (right_prime[index] - left_prime[index]) * params.pixel_size_nm
                )
                if valid[index]
                else math.nan
            )

        unrotated_raw = statistics.sequence_metrics(
            y,
            left,
            right,
            params.pixel_size_nm,
            1.0,
            params.standard_deviation_ddof,
        )
        rotated_raw = statistics.sequence_metrics(
            center_y_prime,
            left_prime,
            right_prime,
            params.pixel_size_nm,
            1.0,
            params.standard_deviation_ddof,
        )
        unrot_y4, unrot_left4, unrot_right4 = statistics.apply_stat_mode(
            y, left, right, "group4_mean", params, group_size
        )
        rot_y4, rot_left4, rot_right4 = statistics.apply_stat_mode(
            center_y_prime,
            left_prime,
            right_prime,
            "group4_mean",
            params,
            group_size,
        )
        unrotated_group4 = statistics.sequence_metrics(
            unrot_y4,
            unrot_left4,
            unrot_right4,
            params.pixel_size_nm,
            1.0,
            params.standard_deviation_ddof,
        )
        rotated_group4 = statistics.sequence_metrics(
            rot_y4,
            rot_left4,
            rot_right4,
            params.pixel_size_nm,
            1.0,
            params.standard_deviation_ddof,
        )

        output[space_index] = {
            "space_index": space_index,
            "space_role": str(rows[0].get("space_role", "")),
            "rotation_angle_from_vertical_deg": angle_deg,
            "rotation_normal_factor": abs(cos_theta),
            "valid_sample_count": int(valid.sum()),
            "eligible_sample_count": sum(
                not item.get("background_excluded", False) for item in rows
            ),
            "valid_group_count": int(unrotated_group4["valid_count"]),
            "unrotated_raw_LWR_nm": finite_float(unrotated_raw["lwr_x_nm"]),
            "rotated_raw_LWR_nm": finite_float(rotated_raw["lwr_x_nm"]),
            "unrotated_raw_CD_nm": finite_float(unrotated_raw["mean_cd_x_nm"]),
            "unrotated_raw_LER_left_nm": finite_float(unrotated_raw["ler_left_nm"]),
            "unrotated_raw_LER_right_nm": finite_float(unrotated_raw["ler_right_nm"]),
            "unrotated_group4_LWR_nm": finite_float(unrotated_group4["lwr_x_nm"]),
            "rotated_raw_CD_nm": finite_float(rotated_raw["mean_cd_x_nm"]),
            "rotated_raw_LER_left_nm": finite_float(rotated_raw["ler_left_nm"]),
            "rotated_raw_LER_right_nm": finite_float(rotated_raw["ler_right_nm"]),
            "rotated_group4_LWR_nm": finite_float(rotated_group4["lwr_x_nm"]),
        }
    return output


def aggregate_coordinate_metrics(
    metrics_by_space: dict[int, dict[str, Any]],
    selected_indices: Sequence[int],
) -> dict[str, Any]:
    selected = [
        metrics_by_space[index]
        for index in selected_indices
        if index in metrics_by_space
    ]
    if not selected:
        return {}
    metric_names = [
        "unrotated_raw_CD_nm",
        "unrotated_raw_LER_left_nm",
        "unrotated_raw_LER_right_nm",
        "unrotated_group4_LWR_nm",
        "rotated_raw_CD_nm",
        "rotated_raw_LER_left_nm",
        "rotated_raw_LER_right_nm",
        "rotated_group4_LWR_nm",
    ]
    result = {
        name: safe_mean(item.get(name) for item in selected) for name in metric_names
    }
    angles = [
        finite_float(item.get("rotation_angle_from_vertical_deg")) for item in selected
    ]
    result.update(
        {
            "selected_space_count": len(selected),
            "rotation_angle_mean_deg": safe_mean(angles),
            "rotation_angle_abs_mean_deg": safe_mean(abs(value) for value in angles),
            "rotation_angle_min_deg": min(
                (value for value in angles if np.isfinite(value)), default=math.nan
            ),
            "rotation_angle_max_deg": max(
                (value for value in angles if np.isfinite(value)), default=math.nan
            ),
        }
    )
    return result


def configure_statistics_for_v110() -> None:
    """Reuse the proven V9 statistic primitives for both labeled mixed-result coordinate systems."""
    statistics.METHOD_CODE = V113_METHOD_CODE
    statistics.METHOD = {
        "result_id": V113_RESULT_ID,
        "edge_variant": "v13_raw_cd_lwr__v10_raw_ler",
        "stat_mode": V113_STAT_MODE,
        "aggregation_mode": V113_AGGREGATION,
        "cn_name": "V1_17 旋转/未旋转混合结果",
        "physical_meaning": (
            "V13/V1.7侧带50%阈值边缘提供raw CD与group4 LWR；"
            "V10峰值中位数50%阈值边缘提供raw LER；"
            "同时报告图像坐标和以trench中心线为Y'轴的旋转坐标。"
        ),
    }
    statistics.V9_VERSION = SCRIPT_VERSION
    statistics.GEOMETRIES = {
        "unrotated": {
            "folder": "unrotated",
            "title": "Unrotated image X/Y coordinates",
            "metrics": {
                "CD": {
                    "machine_col": "machine_CD_nm",
                    "calc_col": "unrotated_mixed_CD_nm",
                    "signed_pct_col": "unrotated_mixed_CD_signed_diff_pct",
                    "abs_pct_col": "unrotated_mixed_CD_abs_diff_pct",
                    "unit": "nm",
                    "note": "V13 raw CD along the original image X axis.",
                },
                "LER_left": {
                    "machine_col": "machine_LER_left_nm",
                    "calc_col": "unrotated_mixed_LER_left_nm",
                    "signed_pct_col": "unrotated_mixed_LER_left_signed_diff_pct",
                    "abs_pct_col": "unrotated_mixed_LER_left_abs_diff_pct",
                    "unit": "nm",
                    "note": "V10 raw left-edge X-position 3sigma.",
                },
                "LWR": {
                    "machine_col": "machine_LWR_nm",
                    "calc_col": "unrotated_mixed_LWR_nm",
                    "signed_pct_col": "unrotated_mixed_LWR_signed_diff_pct",
                    "abs_pct_col": "unrotated_mixed_LWR_abs_diff_pct",
                    "unit": "nm",
                    "note": "V13 group4 width 3sigma along image X.",
                },
            },
        },
        "rotated": {
            "folder": "rotated",
            "title": "Centerline-rotated X'/Y' coordinates",
            "metrics": {
                "CD": {
                    "machine_col": "machine_CD_nm",
                    "calc_col": "rotated_mixed_CD_nm",
                    "signed_pct_col": "rotated_mixed_CD_signed_diff_pct",
                    "abs_pct_col": "rotated_mixed_CD_abs_diff_pct",
                    "unit": "nm",
                    "note": "V13 raw CD along X', normal to the fitted centerline.",
                },
                "LER_left": {
                    "machine_col": "machine_LER_left_nm",
                    "calc_col": "rotated_mixed_LER_left_nm",
                    "signed_pct_col": "rotated_mixed_LER_left_signed_diff_pct",
                    "abs_pct_col": "rotated_mixed_LER_left_abs_diff_pct",
                    "unit": "nm",
                    "note": "V10 raw left-edge X' 3sigma after centerline rotation.",
                },
                "LWR": {
                    "machine_col": "machine_LWR_nm",
                    "calc_col": "rotated_mixed_LWR_nm",
                    "signed_pct_col": "rotated_mixed_LWR_signed_diff_pct",
                    "abs_pct_col": "rotated_mixed_LWR_abs_diff_pct",
                    "unit": "nm",
                    "note": "V13 group4 width 3sigma along centerline-normal X'.",
                },
            },
        },
    }


def generate_coordinate_statistics(
    detailed: pd.DataFrame,
    output_dir: Path,
    dpi: int,
    make_plots: bool,
) -> pd.DataFrame:
    configure_statistics_for_v110()
    summary_rows: list[dict[str, Any]] = []
    plot_rows: list[dict[str, Any]] = []
    for mode, mode_spec in statistics.GEOMETRIES.items():
        for metric in ("CD", "LER_left", "LWR"):
            spec = mode_spec["metrics"][metric]
            machine = pd.to_numeric(detailed[spec["machine_col"]], errors="coerce")
            calculated = pd.to_numeric(detailed[spec["calc_col"]], errors="coerce")
            usable = (
                np.isfinite(machine) & np.isfinite(calculated) & (machine.abs() > EPS)
            )
            if not usable.any():
                summary_rows.append(
                    {
                        "geometry": mode,
                        "coordinate_mode": mode,
                        "geometry_label": mode_spec["title"],
                        "metric": metric,
                        "n": 0,
                        "status": "SKIPPED",
                        "warning": "没有可比较的有效数值；缺失结果保留为空",
                        "mape_pct": math.nan,
                        "mean_signed_diff_pct": math.nan,
                        "method_code": V113_METHOD_CODE,
                        "method_name": statistics.METHOD["cn_name"],
                    }
                )
                continue
            stats = statistics.metric_statistics(detailed, mode, metric)
            stats["status"] = "OK"
            stats["method_code"] = V113_METHOD_CODE
            stats["method_name"] = statistics.METHOD["cn_name"]
            stats["coordinate_mode"] = mode
            summary_rows.append(stats)

            if not make_plots:
                continue
            metric_dir = output_dir / "plots" / mode_spec["folder"] / metric
            metric_dir.mkdir(parents=True, exist_ok=True)
            plot_specs = [
                (
                    "machine_vs_python",
                    "01_machine_vs_python.png",
                    statistics._save_machine_vs_python,
                ),
                (
                    "per_image_percent_error",
                    "02_percent_error_by_image.png",
                    statistics._save_error_by_image,
                ),
                (
                    "percent_error_histogram",
                    "03_percent_error_histogram.png",
                    statistics._save_error_histogram,
                ),
                ("bland_altman", "04_bland_altman.png", statistics._save_bland_altman),
            ]
            for plot_type, filename, function in plot_specs:
                target = metric_dir / filename
                function(detailed, mode, metric, target, dpi)
                plot_rows.append(
                    {
                        "coordinate_mode": mode,
                        "metric": metric,
                        "plot_type": plot_type,
                        "relative_path": str(target.relative_to(output_dir)),
                    }
                )

        if make_plots:
            mode_summary = pd.DataFrame(
                [
                    row
                    for row in summary_rows
                    if row["coordinate_mode"] == mode and row["n"] > 0
                ]
            )
            if mode_summary.empty:
                continue
            summary_dir = output_dir / "plots" / mode_spec["folder"]
            for record in statistics._save_geometry_summary_bars(
                mode_summary, mode, summary_dir, dpi
            ):
                plot_rows.append(
                    {
                        "coordinate_mode": mode,
                        "metric": record["metric"],
                        "plot_type": record["plot_type"],
                        "relative_path": str(
                            Path(record["relative_path"]).relative_to(output_dir)
                        ),
                    }
                )

    summary = pd.DataFrame(summary_rows)
    write_csv(summary, output_dir / "V1_17_statistics_summary.csv")
    write_csv(pd.DataFrame(plot_rows), output_dir / "V1_17_plot_index.csv")
    write_workbook(
        output_dir / "V1_17_statistics_and_machine_comparison.xlsx",
        {
            "Statistics_Summary": summary,
            "Machine_Per_Image": detailed,
            "Plot_Index": pd.DataFrame(plot_rows),
        },
    )
    return summary


def build_machine_comparison(
    image_df: pd.DataFrame,
    excel_path: Path,
    args: argparse.Namespace,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_excel(excel_path, sheet_name=args.machine_sheet, header=None)
    start = max(0, int(args.machine_start_row) - 1)
    end = min(len(raw), int(args.machine_end_row)) if args.machine_end_row else len(raw)
    image_col = excel_column_index(args.machine_image_column)
    cd_col = excel_column_index(args.machine_cd_column)
    ler_col = excel_column_index(args.machine_ler_left_column)
    lwr_col = excel_column_index(args.machine_lwr_column)

    machine_rows: list[dict[str, Any]] = []
    for index in range(start, end):
        image_value = raw.iat[index, image_col] if image_col < raw.shape[1] else None
        if pd.isna(image_value) or not str(image_value).strip():
            continue
        machine_rows.append(
            {
                "machine_excel_row": index + 1,
                "machine_image_name_raw": str(image_value),
                "match_key": normalize_match_key(image_value),
                "machine_CD_nm": finite_float(raw.iat[index, cd_col])
                if cd_col < raw.shape[1]
                else math.nan,
                "machine_LER_left_nm": finite_float(raw.iat[index, ler_col])
                if ler_col < raw.shape[1]
                else math.nan,
                "machine_LWR_nm": finite_float(raw.iat[index, lwr_col])
                if lwr_col < raw.shape[1]
                else math.nan,
            }
        )
    machine_df = pd.DataFrame(
        machine_rows,
        columns=[
            "machine_excel_row",
            "machine_image_name_raw",
            "match_key",
            "machine_CD_nm",
            "machine_LER_left_nm",
            "machine_LWR_nm",
        ],
    )

    calc = image_df.copy()
    calc["match_key"] = calc["image"].map(normalize_match_key)
    duplicated = calc[calc["match_key"].duplicated(keep=False)].copy()
    unique_calc = calc.drop_duplicates("match_key", keep=False)
    detailed = machine_df.merge(
        unique_calc, on="match_key", how="inner", suffixes=("_machine", "")
    )
    # An unmatched reference is a statistics outcome, not a measurement error.
    detailed["image_name"] = detailed["image_key"]
    detailed["result_id"] = V113_RESULT_ID
    detailed["aggregation_mode"] = V113_AGGREGATION
    detailed["method_code"] = V113_METHOD_CODE

    pairs = [
        (f"{mode}_mixed_{metric}", f"{mode}_mixed_{metric}_nm", f"machine_{metric}_nm")
        for mode in MODES
        for metric in ("CD", "LER_left", "LWR")
    ]
    for prefix, calc_col, machine_col in pairs:
        calc_values = pd.to_numeric(detailed[calc_col], errors="coerce")
        machine_values = pd.to_numeric(detailed[machine_col], errors="coerce")
        signed = (calc_values - machine_values) / machine_values.abs() * 100.0
        signed = signed.where(machine_values.abs() > EPS)
        detailed[f"{prefix}_signed_diff_pct"] = signed
        detailed[f"{prefix}_abs_diff_pct"] = signed.abs()

    matched_machine = set(detailed["machine_excel_row"].astype(int))
    unmatched_machine = machine_df[
        ~machine_df["machine_excel_row"].isin(matched_machine)
    ].copy()
    matched_keys = set(detailed["match_key"].astype(str))
    unmatched_png = calc[~calc["match_key"].isin(matched_keys)].copy()
    if not duplicated.empty:
        duplicated = duplicated.assign(
            unmatched_reason="duplicate_basename_in_recursive_tree"
        )
        unmatched_png = pd.concat(
            [unmatched_png, duplicated], ignore_index=True
        ).drop_duplicates(subset=["image_key"], keep="last")

    detailed_path = output_dir / "machine_comparison_detailed.csv"
    detailed = result_first(detailed)
    write_csv(detailed, detailed_path)
    write_csv(unmatched_machine, output_dir / "machine_unmatched_rows.csv")
    write_csv(unmatched_png, output_dir / "machine_unmatched_png.csv")

    summary = generate_coordinate_statistics(
        detailed,
        output_dir,
        dpi=int(args.plot_dpi),
        make_plots=not bool(args.skip_statistics_plots),
    )
    return (
        detailed,
        summary,
        pd.concat(
            [
                unmatched_machine.assign(unmatched_type="machine_row"),
                unmatched_png.assign(unmatched_type="png"),
            ],
            ignore_index=True,
            sort=False,
        ),
    )


def engine_values(metrics: dict[str, Any], mode: str) -> dict[str, float]:
    """raw denotes ungrouped samples; mode separately identifies the coordinates."""
    return {
        "CD": finite_float(metrics.get(f"{mode}_raw_CD_nm")),
        "LER_left": finite_float(metrics.get(f"{mode}_raw_LER_left_nm")),
        "LER_right": finite_float(metrics.get(f"{mode}_raw_LER_right_nm")),
        "LWR": finite_float(metrics.get(f"{mode}_group4_LWR_nm")),
    }


def metric_groups(
    v10_metrics: dict[str, Any], v13_metrics: dict[str, Any]
) -> dict[str, Any]:
    output = {}
    for mode in MODES:
        engines = {
            "V10": engine_values(v10_metrics, mode),
            "V13": engine_values(v13_metrics, mode),
        }
        mixed = {
            name: engines["V10" if name.startswith("LER") else "V13"][name]
            for name in METRICS
        }
        groups = {"mixed": mixed, **engines}
        output.update(
            {
                f"{mode}_{group}_{name}_nm": groups[group][name]
                for group in GROUPS
                for name in METRICS
            }
        )
    return output


def result_first(frame: pd.DataFrame) -> pd.DataFrame:
    """Guarantee stable columns even if the first image or every image fails."""
    frame = frame.copy()
    for name in RESULT_COLUMNS:
        if name not in frame:
            frame[name] = np.nan
    meta = [name for name in META_COLUMNS if name in frame]
    rest = [name for name in frame if name not in RESULT_COLUMNS and name not in meta]
    leading = [
        name
        for name in ("image", "status", "method", "image_key", "folder_name", "pattern")
        if name in frame
    ]
    meta = [name for name in meta if name not in leading]
    rest = [name for name in rest if name not in leading]
    ordered = frame[leading + meta + RESULT_COLUMNS + rest]
    if pitch.PITCH_COLUMN in ordered:
        ordered = ordered[
            [c for c in ordered if c != pitch.PITCH_COLUMN] + [pitch.PITCH_COLUMN]
        ]
    return ordered


def corrected_image_row(
    raw: dict[str, Any], v10_metrics: dict[str, Any], v13_metrics: dict[str, Any]
) -> dict[str, Any]:
    row = metric_groups(v10_metrics, v13_metrics)
    safe = list(META_COLUMNS) + [
        "v10_measurement_valid",
        "v13_measurement_valid",
        "v10_triplet_complete",
        "v13_triplet_complete",
        "v10_selected_roles",
        "v13_selected_roles",
        "v10_selected_trench_count",
        "v13_selected_trench_count",
        "v10_stable_trench_count",
        "v13_stable_trench_count",
        "selected_trench_count",
        "stable_trench_count",
    ]
    row.update({key: raw[key] for key in safe if key in raw})
    row.update(
        method="mixed / V13 / V10",
        algorithm_version=SCRIPT_VERSION,
        coordinate_mode="both",
        source_CD_engine="V13",
        source_LER_engine="V10",
        source_LWR_engine="V13",
        group_size=PUBLIC_GROUP_SIZE,
        CD_stat_mode="ungrouped_mean",
        LER_stat_mode="ungrouped_3sigma",
        LWR_stat_mode=f"group{PUBLIC_GROUP_SIZE}_3sigma",
    )
    for engine, metrics in (("V13", v13_metrics), ("V10", v10_metrics)):
        row[f"{engine}_rotation_angle_mean_deg"] = finite_float(
            metrics.get("rotation_angle_mean_deg")
        )
        row[f"{engine}_rotation_angle_abs_mean_deg"] = finite_float(
            metrics.get("rotation_angle_abs_mean_deg")
        )
        row[f"{engine}_selected_count"] = metrics.get("selected_space_count", 0)
    return row


def method_summary(image_df):
    columns = [
        "image",
        "status",
        "method",
        "image_key",
        "folder_name",
        "pattern",
        "pixel_size_nm",
        "group_size",
    ]
    columns += [f"{mode}_{metric}_nm" for mode in MODES for metric in METRICS]
    columns += ["method_complete", "CD_source", "LER_source", "LWR_source", "warning"]
    columns += [
        "pitch_source",
        "pitch_count",
        "pitch_requested_count",
        "pitch_status",
        pitch.PITCH_COLUMN,
    ]
    rows = []
    for image in image_df.to_dict("records"):
        for method in GROUPS:
            row = {key: image.get(key) for key in columns[:8]}
            row.update(
                method=method,
                warning=image.get("warning", image.get("error_message", "")),
                CD_source="V13" if method == "mixed" else method,
                LER_source="V10" if method == "mixed" else method,
                LWR_source="V13" if method == "mixed" else method,
            )
            for mode in MODES:
                for metric in METRICS:
                    row[f"{mode}_{metric}_nm"] = finite_float(
                        image.get(f"{mode}_{method}_{metric}_nm")
                    )
            pitch_engine = "V13" if method == "mixed" else method
            row.update(
                pitch_source=pitch_engine,
                pitch_count=image.get(f"{pitch_engine}_pitch_count", 0),
                pitch_requested_count=image.get("requested_max_number"),
                pitch_status=image.get(f"{pitch_engine}_pitch_status", "UNAVAILABLE"),
            )
            row[pitch.PITCH_COLUMN] = finite_float(
                image.get(f"{pitch_engine}_rotated_pitch_CD_nm")
            )
            row["method_complete"] = all(
                np.isfinite(row[f"{mode}_{metric}_nm"])
                for mode in MODES
                for metric in METRICS
            )
            rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def corrected_object_rows(
    raw_rows, v10_metrics, v13_metrics, v10_selected, v13_selected
):
    output = []
    for raw in raw_rows:
        i10, i13 = raw.get("v10_space_index"), raw.get("v13_space_index")
        a, b = v10_metrics.get(i10, {}), v13_metrics.get(i13, {})
        row = metric_groups(a, b)
        keys = list(META_COLUMNS) + [
            "trench_id",
            "space_role",
            "v10_space_index",
            "v13_space_index",
        ]
        row.update({key: raw[key] for key in keys if key in raw})
        row.update(
            coordinate_mode="both",
            algorithm_version=SCRIPT_VERSION,
            group_size=PUBLIC_GROUP_SIZE,
            V10_used_for_image_statistics=i10 in v10_selected,
            V13_used_for_image_statistics=i13 in v13_selected,
            V10_rotation_angle_deg=finite_float(
                a.get("rotation_angle_from_vertical_deg")
            ),
            V13_rotation_angle_deg=finite_float(
                b.get("rotation_angle_from_vertical_deg")
            ),
            V10_valid_sample_count=a.get("valid_sample_count", 0),
            V13_valid_sample_count=b.get("valid_sample_count", 0),
        )
        for engine, metrics in (("V10", a), ("V13", b)):
            for name in (
                "valid_group_count",
                "eligible_sample_count",
                "rotated_raw_LWR_nm",
                "unrotated_raw_LWR_nm",
            ):
                row[f"{engine}_{name}"] = metrics.get(name, math.nan)
        output.append(row)
    return output


def engine_object_rows(metrics_by_space, selected, engine, meta):
    rows = []
    for index, metric in metrics_by_space.items():
        row = {
            f"{mode}_{name}_nm": value
            for mode in MODES
            for name, value in engine_values(metric, mode).items()
        }
        row.update(meta)
        row.update(
            engine=engine,
            coordinate_mode="both",
            space_index=index,
            space_role=metric.get("space_role", ""),
            used_for_image_statistics=index in selected,
            rotation_angle_deg=metric["rotation_angle_from_vertical_deg"],
            valid_sample_count=metric["valid_sample_count"],
            group_size=PUBLIC_GROUP_SIZE,
        )
        for name in (
            "valid_group_count",
            "eligible_sample_count",
            "rotated_raw_LWR_nm",
            "unrotated_raw_LWR_nm",
        ):
            row[name] = metric.get(name, math.nan)
        rows.append(row)
    return rows


def corrected_sample_rows(samples, engine, meta, compact):
    core = [
        "space_index",
        "space_role",
        "left_trench_id",
        "right_trench_id",
        "line_left_trench_right_x_px",
        "line_right_trench_left_x_px",
        "sample_index",
        "valid",
        "failure_reason",
        "background_excluded",
        "rotated_axis_angle_from_vertical_deg",
        "rotated_axis_normal_factor",
        "rotated_left_x_prime_px",
        "rotated_right_x_prime_px",
        "rotated_center_y_prime_px",
        "rotated_local_cd_nm",
        "sample_y",
        "left_edge_x",
        "right_edge_x",
        "unrotated_local_cd_nm",
        "continuity_synthetic",
        "continuity_source",
    ]
    diagnostic = [
        "average_y0",
        "average_y1",
        "actual_average_range",
        "recovered",
        "postfilter_rejected",
        "left_is_flyer",
        "right_is_flyer",
        "left_threshold",
        "right_threshold",
        "dark_level",
        "threshold_reference_mode",
        "continuity_original_valid",
        "continuity_original_failure_reason",
        "continuity_strength_used",
        "continuity_method",
        "pre_path_valid",
        "viterbi_selected",
        "viterbi_reselected",
        "viterbi_skipped",
        "viterbi_candidate_count",
        "viterbi_source",
        "viterbi_path_cost",
        "erf_fitted",
        "erf_failure_reason",
        "erf_left_sigma_px",
        "erf_right_sigma_px",
        "erf_left_rmse",
        "erf_right_rmse",
        "flyer_removed",
        "flyer_replaced",
    ]
    keys = core if compact else core + diagnostic
    return [
        {
            **meta,
            "engine": engine,
            "coordinate_mode": "both",
            **{key: sample[key] for key in keys if key in sample},
        }
        for sample in samples
    ]


def build_condition_summary(image_df):
    rows = []
    for folder, group in image_df.groupby("folder_name", dropna=False, sort=True):
        row = {
            "folder_name": folder,
            "pattern": group["pattern"].iloc[0],
            "image_count": len(group),
            "coordinate_mode": "both",
            "ok_image_count": int((group["status"] == "OK").sum()),
            "review_image_count": int((group["status"] == "REVIEW").sum()),
        }
        for metric in RESULT_COLUMNS:
            values = pd.to_numeric(group[metric], errors="coerce")
            row[f"{metric}_mean"] = safe_mean(values)
            row[f"{metric}_std"] = safe_std(values)
            row[f"{metric}_median"] = (
                float(values.median()) if values.notna().any() else math.nan
            )
            row[f"{metric}_n"] = int(values.notna().sum())
        rows.append(row)
    frame = pd.DataFrame(rows)
    means = [f"{metric}_mean" for metric in RESULT_COLUMNS]
    for column in means:
        if column not in frame:
            frame[column] = pd.Series(dtype=float)
    return frame[means + [c for c in frame if c not in means]]


def build_image_coordinate_results(image_df):
    columns = (
        [
            "coordinate_mode",
            "result_group",
            "CD_nm",
            "LER_left_nm",
            "LER_right_nm",
            "LWR_nm",
        ]
        + list(META_COLUMNS)
        + ["CD_engine", "LER_engine", "LWR_engine", "group_size"]
    )
    rows = []
    for _, item in image_df.iterrows():
        for mode in MODES:
            for group in GROUPS:
                row = {
                    f"{metric}_nm": item[f"{mode}_{group}_{metric}_nm"]
                    for metric in METRICS
                }
                row.update({name: item.get(name) for name in META_COLUMNS})
                row.update(
                    result_group=group,
                    coordinate_mode=mode,
                    group_size=PUBLIC_GROUP_SIZE,
                    CD_engine="V13" if group == "mixed" else group,
                    LER_engine="V10" if group == "mixed" else group,
                    LWR_engine="V13" if group == "mixed" else group,
                )
                rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def save_internal_plots(image_df, output_dir, dpi):
    if image_df.empty:
        return
    target = output_dir / "plots" / "internal"
    target.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 2, figsize=(15, 12), sharex=True)
    x = np.arange(1, len(image_df) + 1)
    for col, mode in enumerate(MODES):
        for row, metric in enumerate(METRICS):
            ax = axes[row, col]
            for group, style in zip(GROUPS, ("-", "--", ":")):
                ax.plot(
                    x,
                    image_df[f"{mode}_{group}_{metric}_nm"],
                    linestyle=style,
                    marker="o",
                    markersize=3,
                    linewidth=1,
                    label=group,
                )
            ax.set_ylabel(metric + " (nm)")
            ax.set_title("Rotated" if mode == "rotated" else "Unrotated")
            ax.grid(True, alpha=0.25)
            ax.legend()
        axes[-1, col].set_xlabel("Processed PNG index")
    fig.suptitle("V1_17: rotated first, unrotated second")
    fig.tight_layout()
    fig.savefig(
        target / "旋转_未旋转_metrics_by_image.png", dpi=dpi, bbox_inches="tight"
    )
    plt.close(fig)


def save_rotated_diagnostic(samples, target, title, pixel_size):
    if not samples or not any("rotated_left_x_prime_px" in row for row in samples):
        return
    fig, axes = plt.subplots(2, 1, figsize=(10, 7))
    frame = pd.DataFrame(samples)
    for index, rows in frame.groupby("space_index", sort=True):
        rows = rows.sort_values("sample_index")
        x = rows["sample_index"]
        for side in ("left", "right"):
            values = rows[f"rotated_{side}_x_prime_px"] * pixel_size
            axes[0].plot(
                x, values - values.mean(), linewidth=0.8, label=f"{index} {side}"
            )
        axes[1].plot(x, rows["rotated_local_cd_nm"], linewidth=0.8, label=str(index))
    axes[0].set_ylabel("Rotated edge residual (nm)")
    axes[1].set_ylabel("Rotated CD (nm)")
    axes[1].set_xlabel("Original sample slot index")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7)
    fig.suptitle(title + " / rotated")
    fig.tight_layout()
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_workbook(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "旋转_mixed_": "176B56",
        "旋转_V13_": "245A81",
        "旋转_V10_": "76539A",
        "未旋转_mixed_": "5C8075",
        "未旋转_V13_": "617E94",
        "未旋转_V10_": "897B96",
    }
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, frame in sheets.items():
            safe = public_frame(frame)
            if safe.empty and len(safe.columns) == 0:
                safe = pd.DataFrame({"note": ["No records"]})
            safe.to_excel(writer, sheet_name=name[:31], index=False)
        for ws in writer.book.worksheets:
            autofit_worksheet(ws)
            ws.freeze_panes = "D2" if ws.title == "measurement_summary" else "A2"
            for cell in ws[1]:
                for prefix, color in colors.items():
                    if str(cell.value).startswith(prefix):
                        cell.fill = PatternFill("solid", fgColor=color)
                        cell.font = Font(color="FFFFFF", bold=True)
                        break


def process_trench(
    general: regions.GeneralConfig, args: argparse.Namespace
) -> dict[str, Any]:
    general.validate()
    output_dir = general.resolved_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)
    # Keep different output schemas in separate folders so obsolete tables cannot be confused.
    settings_path = output_dir / "settings.json"
    if any(output_dir.iterdir()):
        raise ValueError(
            "输出目录非空；请指定新的空目录，避免覆盖已有结果或混入上次运行的文件。"
        )
    parsed_images, discovery_errors = regions.discover_images(general)
    if not parsed_images and not discovery_errors:
        raise RuntimeError("没有找到可处理的 PNG；程序递归扫描 root 并排除输出目录")
    print(
        f"脚本版本：{SCRIPT_VERSION}\n输入目录：{general.root_dir}\n输出目录：{output_dir}"
    )
    print(
        "总结表：文件名、status、method 在前；每图 mixed / V13 / V10。PSD 按引擎分别输出且不平均。"
    )
    image_rows, object_rows, engine_rows, sample_rows, parameter_rows = (
        [],
        [],
        [],
        [],
        [],
    )
    pitch_rows, pitch_samples = [], []
    line_pairs, line_sources = [], []
    roi_rows = []
    locator_rows = []
    locator_candidates = []
    errors = list(discovery_errors)
    for error in discovery_errors:
        key = relative_key(Path(error["path"]), general.root_dir)
        image_rows.append(
            {
                **{c: math.nan for c in RESULT_COLUMNS},
                **error,
                "image_key": key,
                "folder_name": str(Path(key).parent),
                "status": "ERROR",
                "warning": error["error_message"],
                "method": "mixed / V13 / V10",
            }
        )
    stopped_on_error = False
    psd_batch = PSDBatch(args, output_dir) if not args.skip_psd else None
    for number, parsed in enumerate(parsed_images, start=1):
        key = relative_key(parsed.path, general.root_dir)
        print(
            f"[{number}/{len(parsed_images)}] {general.force_pattern.upper()} {key}",
            flush=True,
        )
        meta = metadata_for(parsed, general)
        try:
            gray = regions.read_gray_png(parsed.path)
            search_reference_nm = (
                parsed.design_space_nm
                if args.pattern == "line"
                else parsed.design_trench_nm
            )
            region = None
            image_args = argparse.Namespace(**vars(args))
            if not args.no_auto_roi or args.end_trim_fraction > 0:
                if args.no_auto_roi:
                    region = regions.manual_pattern_region(gray, args)
                else:
                    region = regions.find_pattern_region(
                        gray,
                        search_reference_nm / general.pixel_size_nm,
                        "trench",
                        args.roi_min_contrast,
                        args.roi_min_length,
                    )
                region = regions.trim_region_ends(region, args.end_trim_fraction)
                roi_rows.append({**meta, **region.diagnostics, "status": "SELECTED"})
                x0, x1, y0, y1 = region.bounds
                # Explicit manual ROI settings win; the mask still excludes background.
                if args.center_x is None:
                    image_args.center_x = (x0 + x1) / 2
                if args.center_y is None:
                    image_args.center_y = (y0 + y1) / 2
                if args.meas_area_width is None:
                    image_args.meas_area_width = int(
                        min(
                            x1 - x0,
                            2 * image_args.center_x,
                            2 * (gray.shape[1] - image_args.center_x),
                        )
                    )
                if args.meas_area_height is None:
                    image_args.meas_area_height = int(
                        min(
                            y1 - y0,
                            2 * image_args.center_y,
                            2 * (gray.shape[0] - image_args.center_y),
                        )
                    )
                if args.extend_length is None:
                    image_args.extend_length = max(
                        2, min(360, y1 - y0 - args.average_range)
                    )
                meta.update(region.diagnostics)
            else:
                roi_rows.append(
                    {**meta, "roi_mode": "disabled", "status": "MANUAL_OR_FULL_IMAGE"}
                )
            observation = region or regions.manual_pattern_region(gray, args)
            ox0, ox1, oy0, oy1 = observation.bounds
            observed_image = gray[oy0:oy1, ox0:ox1].astype(float)
            observed_mask = localization.observation_mask(
                observation.mask, search_reference_nm / general.pixel_size_nm
            )
            mode_diag = localization.analyze_mode(
                observed_image,
                observed_mask[oy0:oy1, ox0:ox1],
                args.locator_mode,
                args.locator_majority,
            )
            mode_diag.update(
                locator_observation_x0=ox0,
                locator_observation_x1=ox1,
                locator_observation_y0=oy0,
                locator_observation_y1=oy1,
            )
            mode_diag["locator_line_context_pixels"] = int(
                observed_mask.sum() - observation.mask.sum()
            )
            image_args.locator_mode = mode_diag["locator_mode_selected"]
            print(
                f"  locator: {args.locator_mode} -> {image_args.locator_mode}; "
                f"end trim: {observation.diagnostics.get('roi_end_trim_px', 0)} px/side",
                flush=True,
            )
            meta.update(mode_diag)
            locator_rows.append({**meta, "status": "SELECTED"})
            p10 = build_v10_params(parsed, gray, general, image_args)
            p13 = build_v13_params(p10)
            p10.validate()
            p13.validate()
            parameter_rows.append({"image_key": key, **jsonable(asdict(p10))})
            mask = region.mask if region is not None else None
            r10 = run_coordinate_engine(
                "V10", p10, parsed, key, args.group_size, meta, errors, mask
            )
            r13 = run_coordinate_engine(
                "V13", p13, parsed, key, args.group_size, meta, errors, mask
            )
            m10, m13 = r10["measurement"], r13["measurement"]
            for engine, measurement in (("V10", m10), ("V13", m13)):
                if measurement is not None:
                    payload = measurement.annotation_payload
                    line_pairs.extend(
                        {"image_key": key, "engine": engine, **row}
                        for row in payload.get("line_pairs", [])
                    )
                    line_sources.extend(
                        {"image_key": key, "engine": engine, **row}
                        for row in payload.get("line_source_samples", [])
                    )
                    diag = payload.get("candidate_diag", {})
                    locator_rows[-1][f"{engine}_locator_threshold_raw"] = diag.get(
                        "locator_threshold_raw"
                    )
                    locator_candidates.extend(
                        {
                            "image_key": key,
                            "engine": engine,
                            "locator_mode_selected": image_args.locator_mode,
                            **row,
                        }
                        for row in diag.get("candidate_debug_rows", [])
                    )
            s10, s13 = r10["samples"], r13["samples"]
            c10, c13 = r10["metrics"], r13["metrics"]
            selected10, selected13 = r10["selected"], r13["selected"]
            primary = primary_from_engines(meta, r10, r13)
            primary.update(
                {k: v for k, v in meta.items() if k.startswith(("roi_", "locator_"))}
            )
            raw_objects = object_keys_from_metrics(c10, c13, meta)
            objects = corrected_object_rows(
                raw_objects, c10, c13, selected10, selected13
            )
            primary.update(
                edge_continuity=args.edge_continuity,
                viterbi_enabled=args.viterbi,
                erf_enabled=args.erf_fit,
                requested_max_number=args.max_number,
                required_min_number=args.min_number,
                threshold_left_pct=args.threshold_left,
                threshold_right_pct=args.threshold_right,
            )
            for engine, measurement, params in (("V10", m10, p10), ("V13", m13, p13)):
                summary = {
                    "pitch_count": 0,
                    "pitch_status": "UNAVAILABLE",
                    "rotated_pitch_CD_nm": math.nan,
                }
                if measurement is not None:
                    try:
                        summary, periods, samples = pitch.measure_pitch(
                            measurement, params, mask
                        )
                        pitch_rows.extend(
                            {"image_key": key, "engine": engine, **r} for r in periods
                        )
                        pitch_samples.extend(
                            {"image_key": key, "engine": engine, **r} for r in samples
                        )
                    except Exception as exc:
                        errors.append(
                            {
                                **meta,
                                "stage": f"{engine}_pitch",
                                "error_type": type(exc).__name__,
                                "error_message": str(exc),
                                "traceback": traceback.format_exc(),
                            }
                        )
                primary.update({f"{engine}_{k}": v for k, v in summary.items()})
                if summary["pitch_status"] != "OK" and primary["status"] != "ERROR":
                    primary["status"] = "REVIEW"
                    primary["warning"] += (
                        f" | {engine} pitch: {summary['pitch_count']}/{args.max_number} valid periods"
                    )
            primary[pitch.PITCH_COLUMN] = primary["V13_rotated_pitch_CD_nm"]
            primary["quality_status_before_synthetic_review"] = primary["status"]
            for engine, samples in (("V10", s10), ("V13", s13)):
                synthetic = sum(
                    bool(r.get("continuity_synthetic")) and bool(r.get("valid"))
                    for r in samples
                )
                primary[f"{engine}_synthetic_count"] = synthetic
                primary[f"{engine}_original_valid_count"] = sum(
                    bool(r.get("continuity_original_valid")) for r in samples
                )
                primary[f"{engine}_continuity_redetected_count"] = sum(
                    bool(r.get("valid"))
                    and r.get("continuity_source") == "relaxed_redetect"
                    for r in samples
                )
                eligible_count = sum(
                    not r.get("background_excluded", False) for r in samples
                )
                primary[f"{engine}_eligible_sample_count"] = eligible_count
                primary[f"{engine}_background_excluded_count"] = (
                    len(samples) - eligible_count
                )
                primary[f"{engine}_detected_fraction"] = sum(
                    bool(r.get("valid")) and not bool(r.get("continuity_synthetic"))
                    for r in samples
                ) / max(1, eligible_count)
                primary[f"{engine}_filled_fraction"] = sum(
                    bool(r.get("valid")) for r in samples
                ) / max(1, eligible_count)
                if synthetic and primary["status"] != "ERROR":
                    primary["status"] = "REVIEW"
            for engine, measurement, params, samples, module in (
                ("V10", m10, p10, s10, edges),
                ("V13", m13, p13, s13, edges),
            ):
                if measurement is None:
                    continue
                # Auxiliary exports must not discard successful measurements.
                actions = []
                if not args.no_annotated_images:
                    actions.append(
                        (
                            "annotation",
                            lambda: module.save_annotated_image(
                                parsed.path.name + " / Unrotated image coordinates",
                                measurement.annotation_payload,
                                output_dir / f"annotated_未旋转_{engine}" / Path(key),
                                params,
                            ),
                        )
                    )
                if args.save_debug_masks:
                    target = (
                        output_dir
                        / "debug"
                        / Path(key).parent
                        / f"{parsed.path.stem}_{engine}_rotated.png"
                    )
                    actions.append(
                        (
                            "diagnostic",
                            lambda: save_rotated_diagnostic(
                                samples,
                                target,
                                key + " / " + engine,
                                general.pixel_size_nm,
                            ),
                        )
                    )
                if psd_batch is not None:

                    def add_psd():
                        if args.psd_edge_source == "single-row":
                            # Independent dense integer-row sampling, no X or Y averaging.
                            psd_params = replace(
                                params,
                                average_range_px=1,
                                smoothing_pixel=1,
                                sample_number=params.extend_length_px,
                                integer_sampling=True,
                            )
                            psd_measurement = module.measure_image(
                                parsed.path, psd_params, region_mask=mask
                            )
                        else:
                            psd_params, psd_measurement = params, measurement
                        psd_batch.add_measurement(
                            psd_measurement, key, engine, psd_params
                        )

                    actions.append(("PSD", add_psd))
                for stage, action in actions:
                    try:
                        action()
                    except Exception as exc:
                        message = f"{engine} {stage}: {type(exc).__name__}: {exc}"
                        errors.append(
                            {
                                **meta,
                                "stage": f"{engine}_{stage}",
                                "error_type": type(exc).__name__,
                                "error_message": message,
                                "traceback": traceback.format_exc(),
                            }
                        )
                        primary["warning"] = (
                            primary["warning"] + " | " + message
                        ).strip(" |")
                        if primary["status"] != "ERROR":
                            primary["status"] = "REVIEW"
            if region is not None:
                try:
                    regions.save_region_overlay(
                        gray, region, output_dir / "ROI" / Path(key)
                    )
                except Exception as exc:
                    errors.append(
                        {
                            **meta,
                            "stage": "ROI_annotation",
                            "error_type": type(exc).__name__,
                            "error_message": str(exc),
                        }
                    )
            if primary["status"] == "ERROR":
                errors.append(
                    {
                        **meta,
                        "stage": "metric_availability",
                        "error_type": "NoUsableMetrics",
                        "error_message": primary["warning"],
                        "traceback": "",
                    }
                )
            image_rows.append(primary)
            object_rows.extend(objects)
            engine_rows.extend(engine_object_rows(c13, selected13, "V13", meta))
            engine_rows.extend(engine_object_rows(c10, selected10, "V10", meta))
            sample_rows.extend(
                corrected_sample_rows(s13, "V13", meta, args.compact_sample_output)
            )
            sample_rows.extend(
                corrected_sample_rows(s10, "V10", meta, args.compact_sample_output)
            )
            print(f"  status={primary['status']} objects={len(objects)}")
            for mode in MODES:
                for group in GROUPS:
                    print(
                        "  "
                        + COORDINATE_LABELS[mode]
                        + " "
                        + group
                        + ": "
                        + "  ".join(
                            f"{metric}={primary[f'{mode}_{group}_{metric}_nm']:.4f}"
                            for metric in METRICS
                        )
                        + " nm"
                    )
        except regions.NoPatternRegion as exc:
            roi_rows.append(
                {**meta, "status": "SKIPPED_BACKGROUND", "reason": str(exc)}
            )
            image_rows.append(
                {
                    **{c: math.nan for c in RESULT_COLUMNS},
                    **meta,
                    "status": "SKIPPED_BACKGROUND",
                    "warning": str(exc),
                    "method": "mixed / V13 / V10",
                }
            )
            # Retain an all-excluded overlay to audit a potentially over-strict mask.
            try:
                empty = regions.PatternRegion(
                    np.zeros(gray.shape, bool), (0, gray.shape[1], 0, gray.shape[0]), {}
                )
                regions.save_region_overlay(gray, empty, output_dir / "ROI" / Path(key))
            except Exception as export_exc:
                errors.append(
                    {
                        **meta,
                        "stage": "ROI_annotation",
                        "error_type": type(export_exc).__name__,
                        "error_message": str(export_exc),
                    }
                )
            print(f"  SKIPPED_BACKGROUND: {exc}")
        except Exception as exc:
            errors.append(
                {
                    **meta,
                    "stage": "dual_coordinate_measurement",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "traceback": traceback.format_exc(),
                }
            )
            if image_rows and image_rows[-1].get("image_key") == key:
                image_rows[-1].update(status="REVIEW", error_message=str(exc))
            else:
                image_rows.append(
                    {
                        **{c: math.nan for c in RESULT_COLUMNS},
                        **meta,
                        "status": "ERROR",
                        "coordinate_mode": "both",
                        "algorithm_version": SCRIPT_VERSION,
                        "error_message": str(exc),
                    }
                )
            print(f"  ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
            if not general.continue_on_error:
                stopped_on_error = True
                break
        if image_rows[-1]["status"] == "ERROR" and not general.continue_on_error:
            stopped_on_error = True
            break

    for row in image_rows:
        row.setdefault("requested_max_number", args.max_number)
        row.setdefault(pitch.PITCH_COLUMN, math.nan)
    image_df = result_first(pd.DataFrame(image_rows))
    object_df = result_first(pd.DataFrame(object_rows))
    measured_images = image_df[image_df["status"].isin(["OK", "REVIEW"])].copy()
    valid_images = (
        measured_images
        if args.include_review_in_summary
        else measured_images[measured_images.status == "OK"].copy()
    )
    coordinate_df = build_image_coordinate_results(measured_images)
    condition_df = build_condition_summary(valid_images)
    error_df = pd.DataFrame(
        errors,
        columns=list(META_COLUMNS)
        + ["stage", "error_type", "error_message", "traceback"],
    )
    sample_df = pd.DataFrame(sample_rows)
    engine_df = pd.DataFrame(engine_rows)
    decisions = {row["image_key"]: row for row in locator_rows}
    locator_frame = pd.DataFrame(
        [
            {
                "image_key": row["image_key"],
                "image": row["image"],
                "locator_mode_requested": args.locator_mode,
                "locator_mode_selected": None,
                **decisions.get(row["image_key"], {}),
                "status": row["status"],
                "locator_decision_status": "SELECTED"
                if row["image_key"] in decisions
                else "UNAVAILABLE",
            }
            for row in image_rows
        ]
    )
    frames = {
        "measurement_summary": method_summary(image_df),
        "image_summary": image_df,
        "coordinate_results": coordinate_df,
        "condition_summary": condition_df,
        "trench_objects": object_df,
        "engine_objects": engine_df,
        "line_trench_pairs": dataframe_or_empty(
            line_pairs,
            [
                "image_key",
                "engine",
                "line_index",
                "left_trench_id",
                "right_trench_id",
                "status",
                "selected",
            ],
        ),
        "line_source_samples": dataframe_or_empty(
            line_sources,
            [
                "image_key",
                "engine",
                "source_trench_id",
                "source_trench_qualified",
                "sample_index",
                "sample_y",
                "left_edge_x",
                "right_edge_x",
                "valid",
            ],
        ),
        "pitch_periods": pd.DataFrame(pitch_rows, columns=pitch.PERIOD_COLUMNS),
        "pitch_samples": pd.DataFrame(pitch_samples, columns=pitch.SAMPLE_COLUMNS),
        "per_sample_results": sample_df,
        "processing_errors": error_df,
        "roi_summary": pd.DataFrame(roi_rows),
        "locator_summary": locator_frame,
        "locator_candidates": (
            pd.DataFrame(locator_candidates)
            if locator_candidates
            else pd.DataFrame(
                columns=[
                    "image_key",
                    "engine",
                    "locator_mode_selected",
                    "candidate_center_x_px",
                    "selected",
                ]
            )
        ),
    }
    outputs = {
        "excel": output_dir / "CD_measurement_200K_V1_17_results.xlsx",
        "settings": settings_path,
    }
    for name, frame in frames.items():
        outputs[name] = output_dir / f"{name}.csv"
        write_csv(frame, outputs[name])
    machine_excel = args.machine_excel
    if machine_excel is None and not args.no_auto_machine_comparison:
        candidate = general.root_dir / "Data.xlsx"
        if candidate.exists():
            machine_excel = candidate

    def batch_failure(stage, exc):
        errors.append(
            {
                "stage": stage,
                "error_type": type(exc).__name__,
                "error_message": str(exc),
                "traceback": traceback.format_exc(),
            }
        )
        print(f"  WARNING {stage}: {exc}", file=sys.stderr)

    if machine_excel is not None:
        try:
            if not Path(machine_excel).is_file():
                raise FileNotFoundError(f"Machine 参考表不存在：{machine_excel}")
            if not valid_images.empty:
                detail, summary, unmatched = build_machine_comparison(
                    valid_images, Path(machine_excel), args, output_dir
                )
                frames.update(
                    Machine_Per_Image=detail,
                    Statistics_Summary=summary,
                    Machine_Unmatched=unmatched,
                )
            else:
                print("没有符合汇总质量门槛的结果，跳过机台比较。")
        except Exception as exc:
            batch_failure("machine_comparison", exc)
    if not args.skip_statistics_plots:
        try:
            save_internal_plots(valid_images, output_dir, args.plot_dpi)
        except Exception as exc:
            batch_failure("statistics_plots", exc)
    if psd_batch is not None:
        try:
            outputs["psd_directory"] = psd_batch.finalize()
        except Exception as exc:
            batch_failure("PSD_export", exc)
    error_df = pd.DataFrame(
        errors,
        columns=list(META_COLUMNS)
        + ["stage", "error_type", "error_message", "traceback"],
    )
    frames["processing_errors"] = error_df
    write_csv(error_df, outputs["processing_errors"])
    settings = {
        "script_version": SCRIPT_VERSION,
        "patch_version": PATCH_VERSION,
        "output_schema": OUTPUT_SCHEMA,
        "coordinate_mode": "both",
        "result_column_order": [column_label(c) for c in RESULT_COLUMNS],
        "coordinate_labels": COORDINATE_LABELS,
        "mixed_sources": {
            "CD": "V13",
            "LER_left": "V10",
            "LER_right": "V10",
            "LWR": "V13",
        },
        "statistical_definition": {
            "CD": "mean ungrouped widths, separately for rotated and unrotated coordinates",
            "LER": "3 sigma of ungrouped edge positions, separately for each coordinate system",
            "LWR": f"3 sigma of group{args.group_size} mean widths, separately for each coordinate system",
            "image_aggregation": "arithmetic mean per selected structure, independently for each engine",
            "rotation": "each engine independently fits a PCA centerline per structure",
            "pitch_CD": "rotated only; same-side edges of consecutive bands define one band + one adjacent gap; PCA normal projection; mean per period then equal mean of max-number nearest valid periods; mixed uses V13; partial counts marked REVIEW; synthetic samples excluded",
        },
        "group_size": args.group_size,
        "general": jsonable(asdict(general)),
        "cli": jsonable(vars(args)),
        "per_image_params": parameter_rows,
        "input_png_count": len(parsed_images) + len(discovery_errors),
        "processed_png_count": len(image_rows),
        "stopped_on_error": stopped_on_error,
        "ok_count": int((image_df.status == "OK").sum()),
        "review_count": int((image_df.status == "REVIEW").sum()),
        "error_count": int((image_df.status == "ERROR").sum()),
        "diagnostic_error_count": len(error_df),
        "skipped_background_count": int(
            (image_df.status == "SKIPPED_BACKGROUND").sum()
        ),
        "summary_quality_policy": "OK + REVIEW"
        if args.include_review_in_summary
        else "OK only",
        "psd_policy": "V10/V13 separately; ungrouped periodograms; no spectrum averaging",
        "line_policy": "all qualified physical trenches; line left=left trench right, line right=right trench left; no inversion or independent line edge refit",
        "partial_result_policy": "keep each finite metric; REVIEW for partial/invalid engine; ERROR only if all metrics unavailable",
        "machine_reference": str(machine_excel) if machine_excel is not None else None,
    }
    settings_path.write_text(
        json.dumps(jsonable(settings), ensure_ascii=False, indent=2), encoding="utf-8"
    )
    frames["settings"] = pd.DataFrame(
        [
            ["script_version", SCRIPT_VERSION],
            ["坐标状态", "旋转+未旋转"],
            [
                "summary_layout",
                "文件名、status、method 在前；每图 mixed、V13、V10 三行；完整宽表见 image_summary",
            ],
            ["group_size", args.group_size],
            ["pixel_size_nm", general.pixel_size_nm],
            ["CD_LER_samples", "ungrouped; computed separately for 旋转 / 未旋转"],
            [
                "LWR_samples",
                f"group{args.group_size}; computed separately for 旋转 / 未旋转",
            ],
            [
                "PSDs",
                f"axis={args.psd_axis}; normal=旋转, unrotated=未旋转; independent detrend/window/group policies",
            ],
        ],
        columns=["item", "value"],
    )
    write_workbook(outputs["excel"], frames)
    print(
        f"\n处理完成：{len(image_rows)} 张 PNG，{len(object_df)} 个结构结果，{len(error_df)} 条错误。"
    )
    print(f"Excel：{outputs['excel']}")
    outputs["_stopped_on_error"] = stopped_on_error
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "V1_17：同时输出旋转和未旋转结果；先旋转 mixed/V13/V10，再未旋转 mixed/V13/V10；"
            "V13计算CD/LWR，V10计算LER；递归处理所有PNG，"
            "不按目录名/文件名过滤，由 --pattern 指定整目录图案。"
        )
    )
    # V1.7-compatible public interface and defaults.
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help=f"输入根目录；默认：{DEFAULT_ROOT}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=f"输出目录；默认在 root 下创建 {DEFAULT_OUTPUT_NAME}",
    )
    parser.add_argument(
        "--pixel-size", type=float, default=1.3181, help="nm/pixel；默认 1.3181"
    )
    parser.add_argument(
        "--pattern",
        type=str.casefold,
        choices=("trench", "line"),
        required=True,
        help="必填：整个 root（含全部子文件夹）的统一图案类型",
    )
    parser.set_defaults(via_reference_nm=None)
    parser.add_argument(
        "--trench-reference-nm",
        type=float,
        default=None,
        help="可选：统一 Trench 检测参考宽度",
    )
    parser.add_argument(
        "--space-reference-nm",
        type=float,
        default=None,
        help="可选：目标之间的间隔宽度nm，显式指定时参与候选间距排序",
    )
    parser.add_argument(
        "--save-debug-masks", action="store_true", help="额外保存调试图"
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="任一图出错立即停止；默认记录后继续",
    )

    # V1.8 measurement defaults. Center defaults to the current image center so
    # the V1.7 arbitrary-folder input contract does not depend on one data set.
    parser.add_argument(
        "--center-x", type=float, default=None, help="测量中心 X；默认逐图使用图像中心"
    )
    parser.add_argument(
        "--center-y", type=float, default=None, help="测量中心 Y；默认逐图使用图像中心"
    )
    parser.add_argument("--search-in", type=int, default=33)
    parser.add_argument("--search-out", type=int, default=30)
    parser.add_argument(
        "--extend-length", type=int, default=None, help="默认 min(360, 图像高度)"
    )
    parser.add_argument("--sample-number", type=int, default=128)
    parser.add_argument("--average-range", type=int, default=32)
    parser.add_argument("--smoothing-pixel", type=int, default=3)
    parser.add_argument("--peak-order", type=int, default=1)
    parser.add_argument("--threshold-left", type=float, default=50.0)
    parser.add_argument("--threshold-right", type=float, default=50.0)
    parser.add_argument("--max-number", type=int, default=3)
    parser.add_argument("--min-number", type=int, default=None)
    parser.add_argument("--candidate-dark-tolerance", type=float, default=None)
    parser.add_argument("--candidate-min-contrast", type=float, default=None)
    parser.add_argument("--topology-min-contrast", type=float, default=2.5)
    parser.add_argument("--tracking-radius", type=float, default=13.0)
    parser.add_argument("--tracking-retry-radius", type=float, default=20.0)
    parser.add_argument("--tracking-max-jump", type=float, default=8.0)
    parser.add_argument("--postfilter-window", type=int, default=9)
    parser.add_argument("--postfilter-mad-multiplier", type=float, default=6.0)
    parser.add_argument("--postfilter-edge-tolerance", type=float, default=2.5)
    parser.add_argument("--minimum-valid-fraction", type=float, default=0.70)
    parser.add_argument("--allow-incomplete-triplet", action="store_true")
    parser.add_argument("--no-recovery", action="store_true")
    parser.add_argument(
        "--flyer-mode", choices=("reserve", "remove", "replace"), default="reserve"
    )
    parser.add_argument("--flyer-left-offset", type=float, default=20.0)
    parser.add_argument("--flyer-right-offset", type=float, default=20.0)
    parser.add_argument(
        "--group-size",
        type=int,
        default=4,
        help="LWR分组点数，正整数；1不分组，4/8等均按实际输入计算",
    )
    parser.add_argument("--no-annotated-images", action="store_true")
    parser.add_argument("--compact-sample-output", action="store_true")

    # Optional V1.8-compatible statistics with Data.xlsx.
    parser.add_argument(
        "--machine-excel",
        type=Path,
        default=None,
        help="Machine 参考 Data.xlsx；默认自动检查 root/Data.xlsx",
    )
    parser.add_argument(
        "--no-auto-machine-comparison",
        action="store_true",
        help="即使 root 有 Data.xlsx 也不自动比较",
    )
    parser.add_argument("--machine-sheet", default="v2")
    parser.add_argument("--machine-start-row", type=int, default=2)
    parser.add_argument("--machine-end-row", type=int, default=84)
    parser.add_argument("--machine-image-column", default="B")
    parser.add_argument("--machine-cd-column", default="D")
    parser.add_argument("--machine-ler-left-column", default="E")
    parser.add_argument("--machine-lwr-column", default="F")
    parser.add_argument("--skip-statistics-plots", action="store_true")
    parser.add_argument("--plot-dpi", type=int, default=180)
    parser.add_argument(
        "--edge-continuity",
        type=int,
        default=0,
        help="0不额外放宽/补点；1-99逐级放宽仅重测缺点；100补齐剩余缺点并标记（背景永不补齐）",
    )
    parser.add_argument("--skip-psd", action="store_true", help="不计算PSD")
    parser.add_argument(
        "--psd-axis",
        choices=("normal", "unrotated", "both"),
        default="both",
        help="normal=旋转；unrotated=未旋转；both=两套，默认both",
    )
    parser.add_argument(
        "--psd-method",
        choices=("periodogram",),
        default="periodogram",
        help="V1_17 PSD 不平均：逐结构、逐连续段 periodogram",
    )
    parser.add_argument(
        "--psd-window", choices=("hann", "boxcar", "blackmanharris"), default="hann"
    )
    parser.add_argument(
        "--psd-detrend", choices=("linear", "constant", "none"), default="linear"
    )
    parser.add_argument("--psd-nperseg", type=int, default=64)
    parser.add_argument(
        "--psd-overlap", type=float, default=0.5, help="Welch重叠比例[0,1)"
    )
    parser.add_argument(
        "--psd-min-segment", type=int, default=16, help="PSD最短连续有效段点数"
    )
    parser.add_argument(
        "--psd-gap-mode", choices=("segments", "interpolate"), default="segments"
    )
    parser.add_argument(
        "--psd-max-gap", type=int, default=2, help="PSD插值允许的最大内部缺点数"
    )
    parser.add_argument(
        "--psd-include-synthetic",
        action="store_true",
        help="允许强制补点参与PSD，默认排除",
    )
    parser.add_argument(
        "--psd-split-wavelength",
        type=float,
        default=100.0,
        help="低高频分界波长nm；默认100",
    )
    parser.add_argument(
        "--line-reference-nm", type=float, default=None, help="亮line参考宽度nm"
    )
    parser.add_argument(
        "--meas-area-width",
        type=int,
        default=None,
        help="粗定位ROI宽度px，默认min(512,图宽)",
    )
    parser.add_argument(
        "--meas-area-height",
        type=int,
        default=None,
        help="粗定位ROI高度px，默认min(512,图高)",
    )
    parser.add_argument(
        "--viterbi",
        type=int,
        choices=(0, 1),
        default=0,
        help="1启用整条边缘对的全局路径选择",
    )
    parser.add_argument(
        "--viterbi-weight", type=float, default=5.0, help="路径跳动代价权重"
    )
    parser.add_argument(
        "--viterbi-max-jump",
        type=float,
        default=3.0,
        help="每个预设采样间隔最大横向位移px",
    )
    parser.add_argument(
        "--viterbi-gap-cost", type=float, default=6.0, help="跳过一个采样槽位的代价"
    )
    parser.add_argument(
        "--viterbi-candidates",
        type=int,
        default=5,
        help="每个位置最多保留的左右边缘对候选数",
    )
    parser.add_argument(
        "--erf-fit",
        type=int,
        choices=(0, 1),
        default=0,
        help="1启用局部ERF灰度过渡拟合",
    )
    parser.add_argument("--erf-window", type=int, default=8, help="ERF拟合窗口半宽px")
    parser.add_argument(
        "--erf-max-shift", type=float, default=2.0, help="ERF交点相对输入边缘最大位移px"
    )
    parser.add_argument(
        "--erf-max-relative-rmse",
        type=float,
        default=0.2,
        help="拟合RMSE/拟合对比度上限",
    )
    parser.add_argument(
        "--check-parameters",
        action="store_true",
        help="检查参数并打印解析值，不执行图像测量",
    )
    parser.add_argument(
        "--threshold-search",
        choices=("bounded", "legacy"),
        default="bounded",
        help="bounded 修复交点位于梯度峰外侧的漏检；legacy 仅用于复核旧结果",
    )
    parser.add_argument(
        "--no-auto-roi",
        action="store_true",
        help="关闭自动背景识别，使用手动区域；端部裁剪仍由end-trim-fraction控制",
    )
    parser.add_argument(
        "--locator-mode",
        type=str.casefold,
        choices=("auto", "dark-line", "bright-line"),
        default="auto",
        help="auto逐图判断；dark-line保留旧定位（line内有假暗条）；bright-line假定line整体亮",
    )
    parser.add_argument(
        "--locator-majority",
        type=float,
        default=0.70,
        help="bright-line候选核心中低于自适应阈值的最小像素比例，(0.5,1]",
    )
    parser.add_argument(
        "--end-trim-fraction",
        type=float,
        default=0.05,
        help="上下各裁去估计平均trench长度的比例，默认1/20；[0,0.5)，0关闭",
    )
    parser.add_argument(
        "--roi-min-contrast",
        type=float,
        default=6.0,
        help="背景识别最低灰度对比度（8bit尺度）",
    )
    parser.add_argument(
        "--roi-min-length", type=int, default=24, help="有效结构区域最低连续高度px"
    )
    parser.add_argument(
        "--include-review-in-summary",
        action="store_true",
        help="允许REVIEW进入条件汇总和机台比较；默认只汇总OK",
    )
    parser.add_argument(
        "--psd-edge-source",
        choices=("single-row", "measurement"),
        default="single-row",
        help="默认单行、1px步长、无X/Y平均重提边缘；measurement复用主配方边缘（仍不分组/频谱平均）",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    global PUBLIC_GROUP_SIZE
    parser = build_parser()
    raw_args = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(raw_args)
    try:
        validate_args(args, raw_args)
    except ValueError as exc:
        parser.error(str(exc))
    PUBLIC_GROUP_SIZE = args.group_size
    if args.check_parameters:
        print(json.dumps(jsonable(vars(args)), ensure_ascii=False, indent=2))
        return 0
    output_dir = args.output or (args.root / DEFAULT_OUTPUT_NAME)
    general = regions.GeneralConfig(
        root_dir=args.root,
        output_dir=output_dir,
        pixel_size_nm=args.pixel_size,
        force_pattern=args.pattern,
        via_reference_nm=args.via_reference_nm,
        trench_reference_nm=args.trench_reference_nm,
        space_reference_nm=args.space_reference_nm,
        save_debug_masks=bool(args.save_debug_masks),
        continue_on_error=not bool(args.stop_on_error),
    )
    try:
        outputs = process_trench(general, args)
        return 1 if outputs.get("_stopped_on_error") else 0
    except KeyboardInterrupt:
        print("用户中断。", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"程序失败：{type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 1


"""Explicit coordinate labels for all exported measurement tables."""
COORDINATE_LABELS = {
    "rotated": "旋转",
    "normal": "旋转",
    "unrotated": "未旋转",
    "both": "旋转+未旋转",
}


def column_label(name):
    name = str(name)
    # Match unrotated first; do not perform substring replacement.
    for prefix, label in (("unrotated_", "未旋转_"), ("rotated_", "旋转_")):
        if name.startswith(prefix):
            return label + name[len(prefix) :]
    original_positions = {
        "sample_y": "未旋转_sample_y_px",
        "left_edge_x": "未旋转_left_edge_x_px",
        "right_edge_x": "未旋转_right_edge_x_px",
        "continuity_original_left_x": "未旋转_continuity_original_left_x_px",
        "continuity_original_right_x": "未旋转_continuity_original_right_x_px",
        "pre_path_left_x": "未旋转_pre_path_left_x_px",
        "pre_path_right_x": "未旋转_pre_path_right_x_px",
        "coordinate_mode": "坐标状态",
    }
    return original_positions.get(name, name)


def label_frame(frame):
    """Leave numeric results untouched; label wide columns and long-form rows."""
    out = frame.copy()
    if "coordinate_mode" not in out:
        for axis_field in ("axis", "geometry"):
            if axis_field in out:
                out["coordinate_mode"] = out[axis_field]
                break
    if "coordinate_mode" in out:
        out["coordinate_mode"] = out["coordinate_mode"].map(
            lambda v: COORDINATE_LABELS.get(v, v)
        )
    return out.rename(columns={name: column_label(name) for name in out.columns})


def validate_args(args, explicit=()):
    explicit = {v.split("=")[0] for v in explicit if v.startswith("--")}

    def error(message):
        raise ValueError(message)

    for key, value in vars(args).items():
        if isinstance(value, float) and not math.isfinite(value):
            error(f"--{key.replace('_', '-')} 必须是有限数")
    if args.roi_min_contrast <= 0 or args.roi_min_length < 4:
        error("roi-min-contrast必须大于0；roi-min-length至少4px")
    if not 0 <= args.end_trim_fraction < 0.5:
        error("end-trim-fraction必须在[0,0.5)之间")
    if not 0.5 < args.locator_majority <= 1:
        error("locator-majority必须在(0.5,1]之间")
    positive = [
        "pixel_size",
        "sample_number",
        "average_range",
        "smoothing_pixel",
        "peak_order",
        "max_number",
        "group_size",
        "tracking_radius",
        "tracking_retry_radius",
        "tracking_max_jump",
        "postfilter_edge_tolerance",
        "plot_dpi",
        "psd_nperseg",
        "psd_min_segment",
        "psd_split_wavelength",
        "viterbi_candidates",
        "viterbi_max_jump",
        "viterbi_gap_cost",
        "erf_window",
        "erf_max_shift",
        "erf_max_relative_rmse",
    ]
    for key in positive:
        if getattr(args, key) <= 0:
            error(f"--{key.replace('_', '-')} 必须大于0")
    nonnegative = [
        "search_in",
        "search_out",
        "topology_min_contrast",
        "postfilter_mad_multiplier",
        "flyer_left_offset",
        "flyer_right_offset",
        "psd_max_gap",
        "viterbi_weight",
    ]
    for key in nonnegative:
        if getattr(args, key) < 0:
            error(f"--{key.replace('_', '-')} 不能小于0")
    for key in [
        "via_reference_nm",
        "trench_reference_nm",
        "line_reference_nm",
        "space_reference_nm",
        "extend_length",
        "meas_area_width",
        "meas_area_height",
        "min_number",
    ]:
        value = getattr(args, key)
        if value is not None and value <= 0:
            error(f"--{key.replace('_', '-')} 必须大于0")
    for key in ["candidate_dark_tolerance", "candidate_min_contrast"]:
        value = getattr(args, key)
        if value is not None and value < 0:
            error(f"--{key.replace('_', '-')} 不能小于0")
    if args.min_number is None:
        args.min_number = min(3, args.max_number)
    if args.min_number > args.max_number:
        error("--min-number 不能大于 --max-number")
    if args.sample_number < 2:
        error("--sample-number 至少2")
    if args.group_size > args.sample_number:
        error("--group-size 不能大于 --sample-number")
    if not 0 < args.minimum_valid_fraction <= 1:
        error("--minimum-valid-fraction 必须在(0,1]")
    if not 0 < args.threshold_left < 100 or not 0 < args.threshold_right < 100:
        error("边缘阈值百分比必须在(0,100)")
    if args.postfilter_window < 3:
        error("--postfilter-window 至少3")
    if not 0 <= args.edge_continuity <= 100:
        error("--edge-continuity 必须在0–100")
    if args.psd_nperseg < 4 or args.psd_min_segment < 4:
        error("PSD段长度至少4")
    if (
        not args.skip_psd
        and args.psd_method == "welch"
        and args.psd_nperseg < args.psd_min_segment
    ):
        error("--psd-nperseg 不能小于 --psd-min-segment")
    if not 0 <= args.psd_overlap < 1:
        error("--psd-overlap 必须在[0,1)")
    if args.erf_window < 3:
        error("--erf-window 至少3px")
    if args.machine_start_row < 1 or args.machine_end_row < 0:
        error("Excel起始行至少1，结束行0表示末尾")
    if args.machine_end_row and args.machine_end_row < args.machine_start_row:
        error("Excel结束行不能小于起始行")
    for key in [
        "machine_image_column",
        "machine_cd_column",
        "machine_ler_left_column",
        "machine_lwr_column",
    ]:
        value = getattr(args, key)
        if not isinstance(value, str) or not value.isalpha() or not value.isascii():
            error(f"--{key.replace('_', '-')} 请使用Excel列字母，如B或AA")
    if args.pattern == "line":
        if (
            args.line_reference_nm is not None
            and args.trench_reference_nm is not None
            and args.line_reference_nm != args.trench_reference_nm
        ):
            error("line-reference-nm与兼容别名trench-reference-nm不能指定不同值")
        if args.line_reference_nm is not None:
            args.trench_reference_nm = args.line_reference_nm
    elif args.line_reference_nm is not None:
        error("--line-reference-nm 仅用于line")
    if args.pattern != "via" and args.via_reference_nm is not None:
        error("--via-reference-nm 仅用于via")
    if args.pattern == "via":
        supported = {
            "--root",
            "--output",
            "--pattern",
            "--pixel-size",
            "--via-reference-nm",
            "--save-debug-masks",
            "--stop-on-error",
            "--check-parameters",
            "--help",
            "--skip-psd",
        }
        ignored = explicit - supported
        if ignored:
            error("Via不支持这些line/trench参数：" + ", ".join(sorted(ignored)))

    if not args.viterbi and explicit & {
        "--viterbi-weight",
        "--viterbi-max-jump",
        "--viterbi-gap-cost",
        "--viterbi-candidates",
    }:
        error("设置Viterbi参数时请同时 --viterbi 1")
    if not args.erf_fit and explicit & {
        "--erf-window",
        "--erf-max-shift",
        "--erf-max-relative-rmse",
    }:
        error("设置ERF参数时请同时 --erf-fit 1")
    if args.psd_method == "periodogram" and explicit & {
        "--psd-nperseg",
        "--psd-overlap",
    }:
        error("psd-nperseg和psd-overlap只适用于Welch")
    if args.psd_gap_mode == "segments" and "--psd-max-gap" in explicit:
        error("psd-max-gap只适用于 --psd-gap-mode interpolate")
    if args.skip_psd and any(v.startswith("--psd-") for v in explicit):
        error("skip-psd与显式PSD计算参数不能同时使用")


if __name__ == "__main__":
    raise SystemExit(main())
