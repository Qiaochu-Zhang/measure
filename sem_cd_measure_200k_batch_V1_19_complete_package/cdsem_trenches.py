"""Measure all physical trenches, qualify them, then select final trench outputs."""

from dataclasses import replace

import cdsem_engine as edges
from cdsem_refinement import choose_objects
from cdsem_localization import mode_params


def measure_trenches(path, params, image, load_warnings, region_mask=None):
    height, width = image.shape
    if (width, height) != (params.expected_width_px, params.expected_height_px):
        raise ValueError(
            f"image size is {width}x{height}; expected "
            f"{params.expected_width_px}x{params.expected_height_px}"
        )
    # These changes affect only locating the source population. All edge and
    # quality parameters stay unchanged, exactly as in the line source pass.
    source_params = replace(
        params,
        max_number=max(3, params.max_number),
        min_number=1,
        require_center_left_right=False,
    )
    candidates, diagnostic = edges.detect_space_candidates(
        image, source_params, region_mask, all_candidates=True
    )
    ordered = sorted(candidates, key=lambda c: c.center_x_px)
    measured, selection, source_samples, pitch_groups = [], [], [], {}
    for index, candidate in enumerate(ordered, 1):
        result, rows, annotation = edges.measure_single_space(
            path.name, image, candidate, index, source_params, region_mask
        )
        tid = candidate.v115_basin_sequence_index
        for row in rows:
            row.update(
                source_trench_id=tid,
                trench_source_id=tid,
                source_trench_qualified=bool(result.stable),
            )
        source_samples.extend(dict(row) for row in rows)
        record = dict(
            source_trench_id=tid,
            source_space_index=index,
            center_x_px=candidate.center_x_px,
            qualified=bool(result.stable),
            selected=False,
            selected_space_index=None,
            selected_role="",
            valid_sample_count=result.valid_sample_count,
            valid_fraction=result.valid_fraction,
            detected_sample_count=sum(
                bool(row["valid"]) and not row.get("continuity_synthetic", False)
                for row in rows
            ),
            background_excluded_count=sum(
                bool(row["background_excluded"]) for row in rows
            ),
            synthetic_sample_count=sum(
                bool(row.get("continuity_synthetic")) for row in rows
            ),
            status="QUALIFIED_NOT_SELECTED" if result.stable else "UNQUALIFIED",
            warning=result.warning,
        )
        selection.append(record)
        measured.append(
            dict(
                candidate=candidate,
                result=result,
                rows=rows,
                annotation=annotation,
                record=record,
            )
        )
        # Keep the entire candidate order, including rejected trenches, so pitch
        # never bridges a rejected/missing physical structure.
        pitch_groups[tid] = rows if result.stable else []

    qualified = [item for item in measured if item["result"].stable]
    # Retain the existing centered/left/right selection and optional spacing
    # prior, but run it only after actual edge measurement and qualification.
    rows_by_id = {
        int(row.get("basin_sequence_index", row.get("candidate_id", -1))): row
        for row in diagnostic.get("candidate_debug_rows", [])
    }
    items_by_id = {
        item["candidate"].v115_basin_sequence_index: item for item in qualified
    }
    qualified_rows = [rows_by_id[tid] for tid in items_by_id]
    original_selected = [
        rows_by_id[tid]
        for tid, item in items_by_id.items()
        if item["candidate"].v115_original_selected
    ]
    selection_params = mode_params(
        replace(params, locator_mode=diagnostic["locator_mode_selected"])
    )
    chosen = choose_objects(
        original_selected,
        qualified_rows,
        params.center_x_px,
        selection_params,
        edges.choose_center_left_right_v115,
    )
    selected = []
    for role, row in chosen:
        tid = int(row.get("basin_sequence_index", row.get("candidate_id", -1)))
        item = items_by_id[tid]
        item["candidate"].role = role
        item["candidate"].is_center_reference = role == "center"
        selected.append(item)
    results, samples, annotations = [], [], []
    for index, item in enumerate(selected, 1):
        candidate = item["candidate"]
        item["record"].update(
            selected=True,
            selected_space_index=index,
            selected_role=candidate.role,
            status="SELECTED",
        )
        item["result"].space_index = index
        item["result"].space_role = candidate.role
        item["annotation"]["space_index"] = index
        for row in item["rows"]:
            row.update(space_index=index, space_role=candidate.role)
        results.append(item["result"])
        samples.extend(item["rows"])
        annotations.append(item["annotation"])
    selected_by_id = {
        item["candidate"].v115_basin_sequence_index: item["candidate"].role
        for item in selected
    }
    for row in source_samples:
        row["selected_for_output"] = row["source_trench_id"] in selected_by_id
    roles = list(selected_by_id.values())
    diag = dict(
        diagnostic,
        trench_method="measure_all_then_qualify",
        trench_measured_count=len(measured),
        trench_qualified_count=len(qualified),
        trench_unqualified_count=len(measured) - len(qualified),
        candidate_selected_count=len(selected),
        triplet_complete={"center", "left", "right"}.issubset(roles),
        selected_roles=",".join(roles),
    )
    measured_by_id = {r["source_trench_id"]: r for r in selection}
    debug = []
    for row in diagnostic.get("candidate_debug_rows", []):
        tid = int(row.get("basin_sequence_index", row.get("candidate_id", -1)))
        debug.append(
            dict(
                row,
                image_name=path.name,
                source_measured=tid in measured_by_id,
                source_qualified=measured_by_id.get(tid, {}).get("qualified", False),
                selected=tid in selected_by_id,
                selected_role=selected_by_id.get(tid, ""),
            )
        )
    diag["candidate_debug_rows"] = debug
    x0, x1, y0, y1, _ = edges.measurement_area_bounds(image.shape, params)
    warnings = list(load_warnings)
    if not selected:
        warnings.append("no qualified trench after measuring all eligible candidates")
    result = edges.assemble_measurement(
        path,
        params,
        image,
        image.copy(),
        warnings,
        (x0, x1, y0, y1),
        diag,
        results,
        samples,
        annotations,
    )
    # The shared aggregator retains its legacy fallback for existing callers.
    # This path deliberately supplies qualified objects only, even if too few.
    result.image_row["warning"] = result.image_row["warning"].replace(
        "aggregate uses available unstable SPACE results when possible",
        "only qualified trenches are retained; unqualified candidates are not used to fill quota",
    )
    result.annotation_payload.update(
        trench_selection=selection,
        trench_source_samples=source_samples,
        trench_pitch_data=dict(
            samples=pitch_groups,
            candidates=ordered,
            mode=diag["locator_mode_suggested"],
            period_px=diag.get("locator_adjacent_period_px"),
        ),
    )
    return result
