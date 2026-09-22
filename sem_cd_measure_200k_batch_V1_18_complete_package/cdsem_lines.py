"""Line measurement from adjacent, independently qualified physical trenches.

No inversion, line-profile classification or independent line edge fitting.
Every line coordinate is copied from the two bounding trench measurements.
"""

from dataclasses import replace
import math
import numpy as np

import cdsem_engine as edges
import cdsem_regions as regions


def adjacent(left, right, diagnostic):
    gap = right.v115_basin_sequence_index - left.v115_basin_sequence_index
    mode = diagnostic.get(
        "locator_mode_suggested", diagnostic.get("locator_mode_selected")
    )
    if not 0 < gap <= (2 if mode == "dark-line" else 1):
        return False
    period = diagnostic.get("locator_adjacent_period_px")
    distance = right.center_x_px - left.center_x_px
    return (
        period is None
        or not np.isfinite(period)
        or period <= 0
        or 0.60 * period <= distance <= 1.40 * period
    )


def shared_rows(left_rows, right_rows, candidate, pair_index, params):
    """Pair only identical sample slots/Y and preserve every shared edge exactly."""
    left_by_slot = {r["sample_index"]: r for r in left_rows}
    right_by_slot = {r["sample_index"]: r for r in right_rows}
    rows = []
    for slot in sorted(left_by_slot.keys() & right_by_slot.keys()):
        a, b = left_by_slot[slot], right_by_slot[slot]
        left = float(a.get("right_edge_x", math.nan))
        right = float(b.get("left_edge_x", math.nan))
        excluded = bool(a.get("background_excluded") or b.get("background_excluded"))
        valid = (
            a.get("valid", False)
            and b.get("valid", False)
            and not excluded
            and a["sample_y"] == b["sample_y"]
            and np.isfinite([left, right]).all()
            and left < right
        )
        row = dict(a)
        # All side-specific diagnostics follow the same physical shared edges.
        for suffix in (
            "peak_x",
            "peak_strength",
            "threshold",
            "peak_level",
            "threshold_bright_level",
            "inside_level",
            "outside_level",
            "is_flyer",
            "peak_fallback",
        ):
            row["left_" + suffix] = a.get("right_" + suffix, math.nan)
            row["right_" + suffix] = b.get("left_" + suffix, math.nan)
        for prefix in ("continuity_original_", "pre_path_"):
            row[prefix + "left_x"] = a.get(prefix + "right_x", math.nan)
            row[prefix + "right_x"] = b.get(prefix + "left_x", math.nan)
        for suffix in ("sigma_px", "rmse"):
            row["erf_left_" + suffix] = a.get("erf_right_" + suffix, math.nan)
            row["erf_right_" + suffix] = b.get("erf_left_" + suffix, math.nan)
        synthetic = bool(a.get("continuity_synthetic") or b.get("continuity_synthetic"))
        cd = (right - left) * params.pixel_size_nm if valid else math.nan
        row.update(
            space_index=pair_index,
            space_role=candidate.role,
            space_center_x=candidate.center_x_px,
            line_pair_index=pair_index,
            left_trench_id=a["source_trench_id"],
            right_trench_id=b["source_trench_id"],
            line_left_trench_right_x_px=left,
            line_right_trench_left_x_px=right,
            left_edge_x=left if valid else math.nan,
            right_edge_x=right if valid else math.nan,
            local_cd_nm=cd,
            local_cd_x_nm=cd,
            local_cd_slant_nm=math.nan,
            valid=bool(valid),
            background_excluded=excluded,
            failure_reason=""
            if valid
            else "background_excluded"
            if excluded
            else "shared_trench_edge_unavailable",
            tracking_anchor_source="shared_trench_edges",
            threshold_reference_mode="shared_trench_edges",
            dark_level=math.nan,
            threshold_dark_level=math.nan,
            anchor_left_x=candidate.global_left_edge_px,
            anchor_right_x=candidate.global_right_edge_px,
            continuity_synthetic=synthetic if valid else False,
            continuity_source="shared_trench_edges_synthetic"
            if synthetic and valid
            else "shared_trench_edges",
            continuity_original_valid=bool(
                a.get("continuity_original_valid")
                and b.get("continuity_original_valid")
            ),
            cd_out_of_recommended_range=bool(
                valid
                and not params.valid_cd_min_factor * params.target_cd_nm
                <= cd
                <= params.valid_cd_max_factor * params.target_cd_nm
            ),
        )
        for flag in (
            "recovered",
            "postfilter_rejected",
            "tracking_retry_used",
            "viterbi_selected",
            "viterbi_reselected",
            "erf_fitted",
            "flyer_removed",
            "flyer_replaced",
        ):
            row[flag] = bool(a.get(flag) or b.get(flag))
        rows.append(row)
    return rows


def summarize_line(image_name, candidate, rows, params):
    """Use the existing metric/rotation functions on immutable shared coordinates."""
    y = np.array([r["sample_y"] for r in rows], float)
    left = np.array([r["left_edge_x"] for r in rows], float)
    right = np.array([r["right_edge_x"] for r in rows], float)
    valid = (
        np.array([r["valid"] for r in rows], bool)
        & np.isfinite(left)
        & np.isfinite(right)
    )
    synthetic = np.array([r["continuity_synthetic"] for r in rows], bool)
    eligible = ~np.array([r["background_excluded"] for r in rows], bool)
    lv, rv = left[valid], right[valid]
    cd = (rv - lv) * params.pixel_size_nm
    lf, lr = edges.mark_flyers(lv, params)
    rf, rr = edges.mark_flyers(rv, params)
    mean, ler_left, ler_right, lwr, count = edges.stats_with_flyer_mode(
        lv, rv, cd, lf, rf, lr, rr, replace(params, flyer_mode="reserve")
    )
    fit = edges.fit_centerline_and_width_factor(y, left, right)
    factor = fit["normal_width_factor"]
    for i in np.flatnonzero(valid):
        rows[i].update(
            local_cd_slant_nm=rows[i]["local_cd_nm"] * factor,
            centerline_angle_from_vertical_deg=fit["angle_from_vertical_deg"],
            normal_width_factor=factor,
        )
    stable = np.sum(valid & ~synthetic) / max(
        1, eligible.sum()
    ) >= params.minimum_valid_fraction and count >= max(
        2, params.standard_deviation_ddof + 1
    )
    warning = "" if stable else "insufficient common measured trench-edge support"
    if (synthetic & valid).any():
        warning = (warning + "; shared synthetic trench edges included; REVIEW").strip(
            "; "
        )
    result = edges.SpaceResult(
        image_name=image_name,
        space_index=rows[0]["space_index"],
        space_role=candidate.role,
        space_center_x=candidate.center_x_px,
        valid_sample_count=int(valid.sum()),
        valid_fraction=float(valid.sum() / max(1, eligible.sum())),
        stable=bool(stable),
        flyer_left_count=sum(bool(r.get("left_is_flyer")) for r in rows if r["valid"]),
        flyer_right_count=sum(
            bool(r.get("right_is_flyer")) for r in rows if r["valid"]
        ),
        cd_out_of_range_count=sum(r["cd_out_of_recommended_range"] for r in rows),
        peak_fallback_count=sum(
            bool(r.get("left_peak_fallback") or r.get("right_peak_fallback"))
            for r in rows
            if r["valid"]
        ),
        hard_rejected_count=sum(r["postfilter_rejected"] for r in rows),
        recovered_count=sum(r["recovered"] for r in rows),
        mean_cd_nm=mean,
        ler_left_nm=ler_left,
        ler_right_nm=ler_right,
        lwr_nm=lwr,
        mean_cd_x_nm=mean,
        lwr_x_nm=lwr,
        mean_cd_slant_nm=mean * factor,
        lwr_slant_nm=lwr * factor,
        centerline_angle_from_vertical_deg=fit["angle_from_vertical_deg"],
        normal_width_factor=factor,
        fit_x_eq_ky_plus_c_k=fit["k"],
        fit_x_eq_ky_plus_c_c_px=fit["c"],
        fit_y_eq_ax_plus_b_a=fit["a"],
        fit_y_eq_ax_plus_b_b_px=fit["b"],
        warning=warning,
    )
    annotation = dict(
        space_index=rows[0]["space_index"],
        candidate=candidate,
        y_values=y,
        left_edges=left,
        right_edges=right,
        initial_left_edges=left.copy(),
        initial_right_edges=right.copy(),
        background_excluded=~eligible,
        stable=bool(stable),
        centerline_fit=fit,
        postfilter_rejected=np.array([r["postfilter_rejected"] for r in rows], bool),
        recovered=np.array([r["recovered"] for r in rows], bool),
        left_flyers=np.array([bool(r.get("left_is_flyer")) for r in rows]),
        right_flyers=np.array([bool(r.get("right_is_flyer")) for r in rows]),
        continuity_sources=[r["continuity_source"] for r in rows],
        continuity_synthetic=synthetic,
        continuity_original_valid=np.array(
            [r["continuity_original_valid"] for r in rows], bool
        ),
        viterbi_selected=np.array([r["viterbi_selected"] for r in rows], bool),
        erf_fitted=np.array([r["erf_fitted"] for r in rows], bool),
    )
    return result, annotation


def select_lines(lines, params):
    if not lines:
        return []
    center = min(
        lines,
        key=lambda item: (
            abs(item["candidate"].center_x_px - params.center_x_px),
            item["candidate"].center_x_px,
        ),
    )
    cx = center["candidate"].center_x_px
    remaining = sorted(
        (line for line in lines if line is not center),
        key=lambda item: abs(item["candidate"].center_x_px - cx),
    )
    chosen = [center]
    if params.require_center_left_right:
        for side in (-1, 1):
            neighbour = next(
                (
                    line
                    for line in remaining
                    if side * (line["candidate"].center_x_px - cx) > 0
                ),
                None,
            )
            if neighbour is not None:
                chosen.append(neighbour)
    for line in remaining:
        if len(chosen) >= params.max_number:
            break
        if not any(line is old for old in chosen):
            chosen.append(line)
    counts = {-1: 0, 1: 0}
    for line in chosen:
        if line is center:
            role = "center"
        else:
            side = -1 if line["candidate"].center_x_px < cx else 1
            counts[side] += 1
            role = ("left" if side == -1 else "right") + (
                str(counts[side]) if counts[side] > 1 else ""
            )
        line["candidate"].role = role
        line["candidate"].is_center_reference = role == "center"
    return chosen


def measure_lines(path, params, image, load_warnings, region_mask=None):
    height, width = image.shape
    if (width, height) != (params.expected_width_px, params.expected_height_px):
        raise ValueError(
            "image dimensions do not match expected measurement dimensions"
        )
    trench_nm = params.space_reference_nm
    if trench_nm is None:
        trench_nm, _ = regions.estimate_trench_references_nm(
            image, params.pixel_size_nm
        )
    source_params = replace(
        params,
        pattern_kind="trench",
        target_cd_nm=trench_nm,
        space_reference_nm=params.target_cd_nm,
        max_number=max(3, params.max_number),
        min_number=1,
        require_center_left_right=False,
    )
    candidates, diag = edges.detect_space_candidates(
        image, source_params, region_mask, all_candidates=True
    )
    ordered = sorted(candidates, key=lambda c: c.center_x_px)
    sources = []
    source_samples = []
    groups = {}
    for index, candidate in enumerate(ordered, 1):
        result, rows, annotation = edges.measure_single_space(
            path.name, image, candidate, index, source_params, region_mask
        )
        tid = candidate.v115_basin_sequence_index
        for row in rows:
            row.update(source_trench_id=tid, source_trench_qualified=result.stable)
        source_samples.extend(dict(row) for row in rows)
        groups[tid] = rows if result.stable else []
        sources.append(
            dict(candidate=candidate, result=result, rows=rows, annotation=annotation)
        )
    pairs = []
    lines = []
    for index, (a, b) in enumerate(zip(sources, sources[1:]), 1):
        ca, cb = a["candidate"], b["candidate"]
        record = dict(
            line_index=index,
            left_trench_id=ca.v115_basin_sequence_index,
            right_trench_id=cb.v115_basin_sequence_index,
            selected=False,
            left_trench_qualified=a["result"].stable,
            right_trench_qualified=b["result"].stable,
        )
        pairs.append(record)
        if not adjacent(ca, cb, diag):
            record["status"] = "NONADJACENT"
            continue
        if not (a["result"].stable and b["result"].stable):
            record["status"] = "UNQUALIFIED_TRENCH"
            continue
        left, right = ca.global_right_edge_px, cb.global_left_edge_px
        if not left < right:
            record["status"] = "INVALID_SHARED_GEOMETRY"
            continue
        candidate = replace(
            ca,
            center_x_px=(left + right) / 2,
            predicted_left_px=left,
            predicted_right_px=right,
            global_left_edge_px=left,
            global_right_edge_px=right,
            global_cd_px=right - left,
            source="adjacent_trench_shared_edges",
            role="candidate",
            is_center_reference=False,
            v115_basin_sequence_index=index - 1,
            v115_basin_left_x_px=left,
            v115_basin_right_x_px=right,
            v115_estimated_width_px=right - left,
        )
        rows = shared_rows(a["rows"], b["rows"], candidate, index, params)
        if not rows:
            record["status"] = "NO_COMMON_SAMPLES"
            continue
        record.update(
            status="OK",
            candidate_center_x_px=candidate.center_x_px,
            shared_left_anchor_px=left,
            shared_right_anchor_px=right,
        )
        lines.append(dict(candidate=candidate, rows=rows, pair=record))
    selected = select_lines(lines, params)
    results = []
    samples = []
    annotations = []
    for index, line in enumerate(selected, 1):
        line["pair"].update(selected=True, line_space_index=index)
        candidate = line["candidate"]
        for row in line["rows"]:
            row.update(space_index=index, space_role=candidate.role)
        result, annotation = summarize_line(path.name, candidate, line["rows"], params)
        results.append(result)
        samples.extend(line["rows"])
        annotations.append(annotation)
    roles = [item["candidate"].role for item in selected]
    diagnostic = dict(
        diag,
        line_method="adjacent_trench_shared_edges",
        physical_trench_reference_nm=trench_nm,
        line_source_trench_count=len(sources),
        line_qualified_trench_count=sum(s["result"].stable for s in sources),
        candidate_seed_count=len(pairs),
        candidate_valid_count=len(lines),
        candidate_selected_count=len(selected),
        triplet_complete={"center", "left", "right"}.issubset(roles),
        selected_roles=",".join(roles),
    )
    diagnostic["candidate_debug_rows"] = [
        dict(row, candidate_kind="source_trench", image_name=path.name)
        for row in diag.get("candidate_debug_rows", [])
    ]
    diagnostic["candidate_debug_rows"].extend(
        dict(row, candidate_kind="line", image_name=path.name) for row in pairs
    )
    x0, x1, y0, y1, _ = edges.measurement_area_bounds(image.shape, params)
    warnings = list(load_warnings)
    if not selected:
        warnings.append("no line bounded by two adjacent qualified trenches")
    result = edges.assemble_measurement(
        path,
        params,
        image,
        image.copy(),
        warnings,
        (x0, x1, y0, y1),
        diagnostic,
        results,
        samples,
        annotations,
    )
    result.annotation_payload.update(
        line_pairs=pairs,
        line_source_samples=source_samples,
        line_pitch_data=dict(
            samples=groups,
            candidates=ordered,
            mode=diag["locator_mode_suggested"],
            period_px=diag.get("locator_adjacent_period_px"),
        ),
    )
    return result
