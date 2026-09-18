"""V1_14: consolidated active implementation, derived from V13_modified."""

from __future__ import annotations
import math
import numpy as np

from dataclasses import replace
from scipy.optimize import least_squares
from scipy.special import erf, erfinv


def _anchors(y, left, right, candidate, width):
    good = np.isfinite(left) & np.isfinite(right) & (right > left)
    if good.any():
        return (np.interp(y, y[good], left[good]), np.interp(y, y[good], right[good]))
    l, r = float(candidate.global_left_edge_px), float(candidate.global_right_edge_px)
    if not (np.isfinite(l) and np.isfinite(r) and 0 <= l < r <= width - 1):
        return None
    return np.full(len(y), l), np.full(len(y), r)


def extend_threshold_search(engine, profile, edge, candidate, params, al, ar, t):
    """Recover a same-threshold crossing just OUTSIDE the selected gradient peak.

    V1.10 searches from peak inward only. V13's plateau threshold can lie
    slightly outside that peak. Extend the allowed side progressively without
    changing the requested threshold percentage or its intensity references.
    """
    if getattr(edge, "failure_reason", "") not in (
        "left_threshold_crossing_not_found",
        "right_threshold_crossing_not_found",
    ):
        return None
    peak_l, peak_r = edge.left_peak_x, edge.right_peak_x
    if not np.all(
        np.isfinite([peak_l, peak_r, edge.left_threshold, edge.right_threshold])
    ):
        return None
    extension = min(12.0, 0.25 * (ar - al)) * t

    def crossing(threshold, side, peak, anchor):
        picks = []
        for j in range(len(profile) - 1):
            v0, v1 = profile[j], profile[j + 1]
            correct = (
                (v0 >= threshold >= v1 and v0 > v1)
                if side == "left"
                else (v0 <= threshold <= v1 and v1 > v0)
            )
            if not correct:
                continue
            pos = j + (threshold - v0) / (v1 - v0)
            allowed = (
                (peak - extension <= pos < candidate.center_x_px)
                if side == "left"
                else (candidate.center_x_px < pos <= peak + extension)
            )
            if allowed:
                picks.append(pos)
        return min(picks, key=lambda pos: abs(pos - anchor)) if picks else None

    l = crossing(edge.left_threshold, "left", peak_l, al)
    r = crossing(edge.right_threshold, "right", peak_r, ar)
    if (
        l is None
        or r is None
        or not 0 <= l < candidate.center_x_px < r <= len(profile) - 1
    ):
        return None
    cd = r - l
    target = params.target_cd_nm / params.pixel_size_nm
    if (
        not params.valid_cd_min_factor * target
        <= cd
        <= params.valid_cd_max_factor * target
    ):
        return None
    if (
        abs(l - al) > params.tracking_max_edge_jump_px
        or abs(r - ar) > params.tracking_max_edge_jump_px
        or abs((l + r - al - ar) / 2) > params.tracking_max_center_shift_px
        or abs(cd - (ar - al))
        > params.tracking_max_cd_change_fraction * max(ar - al, target)
    ):
        return None
    lo, li = engine.edge_topology_levels(
        profile, l, "left", params.topology_band_px, params.topology_gap_px
    )
    ro, ri = engine.edge_topology_levels(
        profile, r, "right", params.topology_band_px, params.topology_gap_px
    )
    if (
        not np.all(np.isfinite([lo, li, ro, ri]))
        or min(lo - li, ro - ri) < params.topology_min_contrast_8bit
    ):
        return None
    return replace(
        edge,
        valid=True,
        failure_reason="",
        left_edge_x=l,
        right_edge_x=r,
        local_cd_px=cd,
        left_outside_level=lo,
        left_inside_level=li,
        right_outside_level=ro,
        right_inside_level=ri,
    )


def restore_edges(
    engine,
    image,
    candidate,
    params,
    y,
    left,
    right,
    cd_nm,
    rows,
    rejected,
    left_tol,
    right_tol,
    cd_tol,
):
    """Called after V1.10's FINAL postfilter and before every downstream statistic.

    Integer strengths use the same prefix of 1..100 recovery attempts. Anchors
    are frozen from the baseline so a higher strength cannot discard a recovery
    accepted at a lower strength. Force fill exists ONLY at strength 100.
    """
    strength = int(params.edge_continuity)
    original_valid = np.isfinite(left) & np.isfinite(right) & (right > left)
    for i, row in enumerate(rows):
        row.update(
            continuity_original_valid=bool(original_valid[i]),
            continuity_original_left_x=float(left[i]),
            continuity_original_right_x=float(right[i]),
            continuity_original_failure_reason=row.get("failure_reason", ""),
            continuity_original_postfilter_rejected=bool(rejected[i]),
            continuity_source="baseline" if original_valid[i] else "missing",
            continuity_strength_used=0,
            continuity_synthetic=False,
            continuity_method="original",
        )
    if strength == 0 or original_valid.all():
        return
    anchor = _anchors(y, left, right, candidate, image.shape[1])
    if anchor is None:
        # No valid coarse geometry: never invent a trench outside the image.
        return
    anchor_l, anchor_r = anchor
    profiles = {}
    for i in np.flatnonzero(~original_valid):
        if rows[i].get("background_excluded"):
            continue
        raw, _, _, actual = engine.average_rows_profile(
            image, float(y[i]), params.average_range_px
        )
        if raw is not None:
            profiles[int(i)] = (
                engine.moving_average_reflect(raw, params.smoothing_pixel),
                actual,
            )

    for level in range(1, strength + 1):
        missing = np.flatnonzero(~(np.isfinite(left) & np.isfinite(right)))
        if not len(missing):
            break
        t = level / 100.0
        relaxed = replace(
            params,
            peak_min_contrast_8bit=params.peak_min_contrast_8bit * (1 - 0.9 * t),
            peak_mad_multiplier=params.peak_mad_multiplier * (1 - 0.9 * t),
            topology_min_contrast_8bit=params.topology_min_contrast_8bit
            * (1 - 0.95 * t),
            tracking_max_edge_jump_px=params.tracking_max_edge_jump_px * (1 + 2 * t),
            tracking_max_center_shift_px=params.tracking_max_center_shift_px
            * (1 + 2 * t),
            tracking_max_cd_change_fraction=params.tracking_max_cd_change_fraction
            * (1 + 2 * t),
            valid_cd_min_factor=params.valid_cd_min_factor * (1 - 0.5 * t),
            valid_cd_max_factor=params.valid_cd_max_factor * (1 + 0.5 * t),
        )
        for i in missing:
            item = profiles.get(int(i))
            if item is None:
                continue
            profile, actual = item
            if (
                actual
                < params.minimum_average_range_fraction
                * params.average_range_px
                * (1 - 0.6 * t)
            ):
                continue
            al, ar = float(anchor_l[i]), float(anchor_r[i])
            # Half-CD ceiling prevents the expanded search from freely crossing
            # adjacent trenches; anchor proximity below provides another gate.
            radius = min(
                params.tracking_retry_radius_px * (1 + t), max(3.0, 0.6 * (ar - al))
            )
            edge = engine.detect_edge_pair(
                profile,
                candidate.center_x_px,
                relaxed,
                anchor_left_x=al,
                anchor_right_x=ar,
                search_radius_px=radius,
                strict_tracking=True,
            )
            method = "relaxed_checks"
            if not edge.valid:
                extended = extend_threshold_search(
                    engine, profile, edge, candidate, relaxed, al, ar, t
                )
                if extended is None:
                    continue
                edge = extended
                method = "extended_threshold_search"
            l, r = float(edge.left_edge_x), float(edge.right_edge_x)
            if not (0 <= l < r <= image.shape[1] - 1):
                continue
            # Relax baseline robust tolerances around frozen local anchors.
            if (
                abs(l - al)
                > max(params.postfilter_min_edge_tolerance_px, left_tol) * (1 + 3 * t)
                or abs(r - ar)
                > max(params.postfilter_min_edge_tolerance_px, right_tol) * (1 + 3 * t)
                or abs((r - l) - (ar - al))
                > max(params.postfilter_min_cd_tolerance_px, cd_tol) * (1 + 3 * t)
            ):
                continue
            if getattr(params, "viterbi_enabled", 0):
                good = np.flatnonzero(np.isfinite(left) & np.isfinite(right))
                before, after = good[good < i], good[good > i]
                neighbors = ([before[-1]] if len(before) else []) + (
                    [after[0]] if len(after) else []
                )
                if any(
                    max(abs(l - left[j]), abs(r - right[j]))
                    > params.viterbi_max_jump_px * abs(i - j)
                    for j in neighbors
                ):
                    continue
            left[i], right[i], cd_nm[i] = l, r, (r - l) * params.pixel_size_nm
            engine.update_sample_row_from_edge(rows[i], edge, cd_nm[i])
            rows[i].update(
                valid=True,
                failure_reason="",
                postfilter_rejected=False,
                continuity_source="relaxed_redetect",
                continuity_strength_used=level,
                continuity_method=method,
                tracking_anchor_source="continuity_frozen_baseline",
                cd_out_of_recommended_range=not (
                    params.valid_cd_min_factor * params.target_cd_nm
                    <= cd_nm[i]
                    <= params.valid_cd_max_factor * params.target_cd_nm
                ),
            )
            rejected[i] = False

    if strength == 100:
        remaining = ~(np.isfinite(left) & np.isfinite(right))
        remaining &= ~np.array([row.get("background_excluded", False) for row in rows])
        if remaining.any():
            detected = ~remaining
            filled_l, filled_r = _anchors(y, left, right, candidate, image.shape[1])
            # np.interp: interior linear interpolation, endpoint nearest hold;
            # if nothing detected, use coarse global pair, explicitly labelled.
            for i in np.flatnonzero(remaining):
                left[i], right[i] = float(filled_l[i]), float(filled_r[i])
                cd_nm[i] = (right[i] - left[i]) * params.pixel_size_nm
                source = (
                    "forced_interpolation" if detected.any() else "forced_global_anchor"
                )
                rows[i].update(
                    valid=True,
                    failure_reason="",
                    postfilter_rejected=False,
                    left_edge_x=float(left[i]),
                    right_edge_x=float(right[i]),
                    local_cd_nm=float(cd_nm[i]),
                    continuity_source=source,
                    continuity_strength_used=100,
                    continuity_synthetic=True,
                    continuity_method=source,
                    tracking_anchor_source=source,
                )
                rejected[i] = False


def draw_overlay(ax, annotation, strength):
    """Join consecutive valid samples; NaN gaps remain gaps below strength 100."""
    if not strength:
        return
    y = annotation["y_values"]
    sources = np.asarray(annotation.get("continuity_sources", ["baseline"] * len(y)))
    for values in (annotation["left_edges"], annotation["right_edges"]):
        ax.plot(values, y, "-", color="deepskyblue", alpha=0.45, linewidth=0.65)
        relaxed = (sources == "relaxed_redetect") & np.isfinite(values)
        forced = np.char.startswith(sources.astype(str), "forced_") & np.isfinite(
            values
        )
        if relaxed.any():
            ax.plot(
                values[relaxed],
                y[relaxed],
                "D",
                markersize=3,
                fillstyle="none",
                color="darkorange",
                label="continuity redetected",
            )
        if forced.any():
            ax.plot(
                values[forced],
                y[forced],
                "s",
                markersize=3,
                fillstyle="none",
                color="magenta",
                label="synthetic fill",
            )


def crossings(profile, threshold, side, lo, hi):
    if not np.isfinite(threshold):
        return []
    out = []
    for j in range(max(0, int(np.floor(lo))), min(len(profile) - 1, int(np.ceil(hi)))):
        a, b = profile[j : j + 2]
        ok = (
            (a >= threshold >= b and a > b)
            if side == "left"
            else (a <= threshold <= b and a < b)
        )
        if ok:
            x = j + (threshold - a) / (b - a)
            if lo <= x <= hi:
                out.append(float(x))
    return out


def make_layer(engine, profile, candidate, params, baseline):
    al = float(candidate.global_left_edge_px)
    ar = float(candidate.global_right_edge_px)
    center = candidate.center_x_px
    target = params.target_cd_nm / params.pixel_size_nm
    radius = params.tracking_search_radius_px
    dark = engine.safe_median_segment(
        profile, al + 0.30 * (ar - al), ar - 0.30 * (ar - al)
    )
    gradient = engine.central_difference(profile)
    peaks = {}
    for side, anchor, sign in [("left", al, "negative"), ("right", ar, "positive")]:
        found = []
        start, end = engine.clamp_int_window(
            anchor - radius, anchor + radius, len(profile)
        )
        for order in range(
            params.peak_order, params.peak_order + params.viterbi_candidates
        ):
            pk, _, _ = engine.adaptive_peak_indices(
                profile,
                gradient,
                start,
                end,
                sign,
                order,
                center,
                params.peak_mad_multiplier,
                dark,
                params.peak_min_contrast_8bit,
                anchor,
                radius,
                params.topology_band_px,
                params.topology_gap_px,
                params.topology_min_contrast_8bit,
            )
            if pk is not None and pk not in found:
                found.append(pk)
        peaks[side] = found
    candidates = []
    maxgrad = max(float(np.nanmax(np.abs(gradient))), 1.0)

    def add(l, r, meta):
        if not (0 <= l < center < r < len(profile)):
            return
        if (
            not params.valid_cd_min_factor * target
            <= r - l
            <= params.valid_cd_max_factor * target
        ):
            return
        lo, li = engine.edge_topology_levels(
            profile, l, "left", params.topology_band_px, params.topology_gap_px
        )
        ro, ri = engine.edge_topology_levels(
            profile, r, "right", params.topology_band_px, params.topology_gap_px
        )
        if (
            not np.all(np.isfinite([lo, li, ro, ri]))
            or min(lo - li, ro - ri) < params.topology_min_contrast_8bit
        ):
            return
        evidence = (
            abs(np.interp(l, np.arange(len(profile)), gradient))
            + abs(np.interp(r, np.arange(len(profile)), gradient))
        ) / (2 * maxgrad)
        emission = 2 * (1 - evidence) + 0.1 * (abs(l - al) + abs(r - ar)) / max(
            radius, 1.0
        )
        if any(
            abs(l - c["left"]) < 0.02 and abs(r - c["right"]) < 0.02 for c in candidates
        ):
            return
        candidates.append(
            dict(left=float(l), right=float(r), cost=float(emission), **meta)
        )

    if baseline is not None:
        add(*baseline, dict(source="original"))
    for pl in peaks["left"]:
        for pr in peaks["right"]:
            if params.engine_kind == "V13":
                d, bl, br = engine.v17_threshold_reference_levels(
                    profile, al, ar, pl, pr, target
                )
            else:
                d = dark
                bl = engine.local_peak_level(profile, pl)
                br = engine.local_peak_level(profile, pr)
            if not np.all(np.isfinite([d, bl, br])) or min(bl, br) <= d:
                continue
            tl = d + params.edge_threshold_left_pct / 100 * (bl - d)
            tr = d + params.edge_threshold_right_pct / 100 * (br - d)
            ls = crossings(
                profile,
                tl,
                "left",
                max(al - radius, pl - radius),
                min(center, al + radius),
            )
            rs = crossings(
                profile,
                tr,
                "right",
                max(center, ar - radius),
                min(len(profile) - 1, ar + radius),
            )
            for l in sorted(ls, key=lambda x: abs(x - al))[: params.viterbi_candidates]:
                for r in sorted(rs, key=lambda x: abs(x - ar))[
                    : params.viterbi_candidates
                ]:
                    add(
                        l,
                        r,
                        dict(
                            source="threshold_candidate",
                            left_threshold=tl,
                            right_threshold=tr,
                            dark_level=d,
                            left_bright=bl,
                            right_bright=br,
                        ),
                    )
    return sorted(candidates, key=lambda c: c["cost"])[: params.viterbi_candidates]


def viterbi_path(layers, weight, max_jump, gap_cost):
    """Minimum total emission + transition + missing-slot costs with backtrace.

    Nodes are left/right pairs. Edges may span missing/skipped slots and use
    their true index distance. No data are compressed and no coordinates are
    fabricated by this solver. Leading/trailing skips are charged as well.
    """
    n = len(layers)
    costs = []
    back = []
    for i, layer in enumerate(layers):
        states = np.array([[c["left"], c["right"]] for c in layer], float).reshape(
            -1, 2
        )
        local = np.array([c["cost"] for c in layer])
        best = local + i * gap_cost
        prev = [(-1, -1)] * len(layer)
        for j in range(i):
            if not len(layers[j]) or not len(layer):
                continue
            gap = i - j
            old = np.array([[c["left"], c["right"]] for c in layers[j]])
            delta = states[:, None, :] - old[None, :, :]
            allowed = np.max(np.abs(delta), axis=2) <= max_jump * gap
            # Slower, longer changes cost less than abrupt jumps.
            transition = weight * np.sum((delta / (max_jump * gap)) ** 2, axis=2)
            total = (
                costs[j][None, :] + transition + (gap - 1) * gap_cost + local[:, None]
            )
            total[~allowed] = np.inf
            k = np.argmin(total, axis=1)
            values = total[np.arange(len(layer)), k]
            improve = values < best
            for q in np.flatnonzero(improve):
                prev[q] = (j, int(k[q]))
            best = np.minimum(best, values)
        costs.append(best)
        back.append(prev)
    best_score = n * gap_cost
    end = (-1, -1)
    for i, cost in enumerate(costs):
        if len(cost):
            k = int(np.argmin(cost))
            score = cost[k] + (n - 1 - i) * gap_cost
            if score < best_score:
                best_score = float(score)
                end = (i, k)
    selected = {}
    i, k = end
    while i >= 0:
        selected[i] = layers[i][k]
        i, k = back[i][k]
    return selected, best_score


def fit_erf(profile, seed, side, percent, window, max_shift, max_relative_rmse):
    lo = max(0, int(math.floor(seed - window)))
    hi = min(len(profile), int(math.ceil(seed + window)) + 1)
    if hi - lo < 7:
        return None, "insufficient_window"
    x = np.arange(lo, hi, dtype=float)
    z = np.asarray(profile[lo:hi], float)
    sign = -1.0 if side == "left" else 1.0
    low = float(np.percentile(z, 10))
    high = float(np.percentile(z, 90))
    if high - low < 1e-6:
        return None, "no_contrast"

    def model(p):
        d, b, x0, sigma = p
        return d + (b - d) * 0.5 * (1 + sign * erf((x - x0) / (np.sqrt(2) * sigma)))

    init = np.array(
        [
            np.clip(low, 0.0001, 254.9999),
            np.clip(high, 0.0001, 254.9999),
            seed,
            min(1.5, window / 2),
        ]
    )
    try:
        fit = least_squares(
            lambda p: model(p) - z,
            init,
            bounds=(
                [0, 0, seed - max_shift, 0.15],
                [255, 255, seed + max_shift, float(window)],
            ),
            loss="soft_l1",
            f_scale=max((high - low) * 0.03, 0.5),
            max_nfev=150,
        )
    except (ValueError, FloatingPointError):
        return None, "optimizer_failure"
    d, b, x0, sigma = fit.x
    if not fit.success or b <= d:
        return None, "invalid_fit"
    rmse = float(np.sqrt(np.mean((model(fit.x) - z) ** 2)))
    value = float(x0 + sign * np.sqrt(2) * sigma * erfinv(2 * percent / 100 - 1))
    if not np.isfinite(value) or abs(value - seed) > max_shift:
        return None, "shift_limit"
    if rmse / (b - d) > max_relative_rmse:
        return None, "relative_rmse_limit"
    return dict(
        x=value,
        sigma=float(sigma),
        rmse=rmse,
        contrast=float(b - d),
        dark=float(d),
        bright=float(b),
    ), ""


def refine_path_and_erf(
    engine, image, candidate, params, y, left, right, cd_nm, rows, rejected
):
    valid = np.isfinite(left) & np.isfinite(right)
    for i, row in enumerate(rows):
        row.update(
            pre_path_valid=bool(valid[i]),
            pre_path_left_x=float(left[i]),
            pre_path_right_x=float(right[i]),
            viterbi_selected=False,
            viterbi_reselected=False,
            viterbi_skipped=False,
            viterbi_candidate_count=0,
            erf_fitted=False,
            erf_failure_reason="",
        )
    if not params.viterbi_enabled and not params.erf_enabled:
        return
    profiles = {}
    for i in range(len(y)):
        if rows[i].get("background_excluded"):
            continue
        profile, _, _, actual = engine.average_rows_profile(
            image, float(y[i]), params.average_range_px
        )
        if (
            profile is not None
            and actual
            >= params.minimum_average_range_fraction * params.average_range_px
        ):
            profiles[i] = engine.moving_average_reflect(profile, params.smoothing_pixel)
    if params.viterbi_enabled:
        layers = []
        for i in range(len(y)):
            baseline = (left[i], right[i]) if valid[i] else None
            layer = (
                make_layer(engine, profiles[i], candidate, params, baseline)
                if i in profiles
                else []
            )
            rows[i]["viterbi_candidate_count"] = len(layer)
            layers.append(layer)
        path, score = viterbi_path(
            layers,
            params.viterbi_weight,
            params.viterbi_max_jump_px,
            params.viterbi_gap_cost,
        )
        left[:] = np.nan
        right[:] = np.nan
        cd_nm[:] = np.nan
        for i, row in enumerate(rows):
            row["viterbi_path_cost"] = score
            if i not in path:
                row.update(
                    valid=False,
                    left_edge_x=math.nan,
                    right_edge_x=math.nan,
                    local_cd_nm=math.nan,
                    failure_reason="viterbi_skipped",
                    viterbi_skipped=True,
                )
                continue
            c = path[i]
            left[i] = c["left"]
            right[i] = c["right"]
            cd_nm[i] = (right[i] - left[i]) * params.pixel_size_nm
            row.update(
                valid=True,
                left_edge_x=float(left[i]),
                right_edge_x=float(right[i]),
                local_cd_nm=float(cd_nm[i]),
                failure_reason="",
                postfilter_rejected=False,
                viterbi_selected=True,
                viterbi_reselected=not (
                    valid[i]
                    and abs(left[i] - row["pre_path_left_x"]) < 1e-9
                    and abs(right[i] - row["pre_path_right_x"]) < 1e-9
                ),
                viterbi_source=c["source"],
            )
            for field in ["left_threshold", "right_threshold", "dark_level"]:
                if field in c:
                    row[field] = c[field]
            rejected[i] = False
    if params.erf_enabled:
        for i in np.flatnonzero(np.isfinite(left) & np.isfinite(right)):
            if i not in profiles:
                continue
            lf, le = fit_erf(
                profiles[i],
                left[i],
                "left",
                params.edge_threshold_left_pct,
                params.erf_window_px,
                params.erf_max_shift_px,
                params.erf_max_relative_rmse,
            )
            rf, re = fit_erf(
                profiles[i],
                right[i],
                "right",
                params.edge_threshold_right_pct,
                params.erf_window_px,
                params.erf_max_shift_px,
                params.erf_max_relative_rmse,
            )
            if lf is None or rf is None:
                rows[i]["erf_failure_reason"] = "left:" + le + ";right:" + re
                continue
            l, r = lf["x"], rf["x"]
            if not 0 <= l < candidate.center_x_px < r < image.shape[1]:
                rows[i]["erf_failure_reason"] = "invalid_geometry"
                continue
            target = params.target_cd_nm / params.pixel_size_nm
            if (
                not params.valid_cd_min_factor * target
                <= r - l
                <= params.valid_cd_max_factor * target
            ):
                rows[i]["erf_failure_reason"] = "width_limit"
                continue
            if params.viterbi_enabled:
                good = np.flatnonzero(np.isfinite(left) & np.isfinite(right))
                others = good[good != i]
                neighbors = []
                before = others[others < i]
                after = others[others > i]
                if len(before):
                    neighbors.append(before[-1])
                if len(after):
                    neighbors.append(after[0])
                if any(
                    max(abs(l - left[j]), abs(r - right[j]))
                    > params.viterbi_max_jump_px * abs(i - j)
                    for j in neighbors
                ):
                    rows[i]["erf_failure_reason"] = "viterbi_jump_guard"
                    continue
            left[i], right[i] = l, r
            cd_nm[i] = (r - l) * params.pixel_size_nm
            rows[i].update(
                left_edge_x=l,
                right_edge_x=r,
                local_cd_nm=float(cd_nm[i]),
                erf_fitted=True,
                erf_left_sigma_px=lf["sigma"],
                erf_right_sigma_px=rf["sigma"],
                erf_left_rmse=lf["rmse"],
                erf_right_rmse=rf["rmse"],
            )


def apply_final_flyers(left, right, cd_nm, rows, params, engine):
    valid = np.isfinite(left) & np.isfinite(right)
    indices = np.flatnonzero(valid)
    lf, lref = engine.mark_flyers(left[valid], params)
    rf, rref = engine.mark_flyers(right[valid], params)
    for k, i in enumerate(indices):
        rows[i]["left_is_flyer"] = bool(lf[k])
        rows[i]["right_is_flyer"] = bool(rf[k])
    if params.flyer_mode == "remove":
        for i in indices[lf | rf]:
            left[i] = right[i] = cd_nm[i] = np.nan
            rows[i].update(
                valid=False,
                left_edge_x=math.nan,
                right_edge_x=math.nan,
                local_cd_nm=math.nan,
                failure_reason="flyer_removed",
                flyer_removed=True,
            )
    elif params.flyer_mode == "replace":
        for k, i in enumerate(indices):
            if not (lf[k] or rf[k]):
                continue
            l = lref[k] if lf[k] else left[i]
            r = rref[k] if rf[k] else right[i]
            if r <= l:
                left[i] = right[i] = cd_nm[i] = np.nan
                rows[i].update(
                    valid=False,
                    left_edge_x=math.nan,
                    right_edge_x=math.nan,
                    local_cd_nm=math.nan,
                    failure_reason="flyer_replace_invalid_geometry",
                    flyer_removed=True,
                )
                continue
            if params.viterbi_enabled:
                good = np.flatnonzero(np.isfinite(left) & np.isfinite(right))
                before = good[good < i]
                after = good[good > i]
                neighbors = ([before[-1]] if len(before) else []) + (
                    [after[0]] if len(after) else []
                )
                if any(
                    max(abs(l - left[j]), abs(r - right[j]))
                    > params.viterbi_max_jump_px * abs(i - j)
                    for j in neighbors
                ):
                    left[i] = right[i] = cd_nm[i] = np.nan
                    rows[i].update(
                        valid=False,
                        left_edge_x=math.nan,
                        right_edge_x=math.nan,
                        local_cd_nm=math.nan,
                        failure_reason="flyer_replace_path_guard",
                        flyer_removed=True,
                    )
                    continue
            left[i] = l
            right[i] = r
            cd_nm[i] = (r - l) * params.pixel_size_nm
            rows[i].update(
                left_edge_x=float(l),
                right_edge_x=float(r),
                local_cd_nm=float(cd_nm[i]),
                flyer_replaced=True,
                continuity_synthetic=True,
                continuity_source="flyer_replacement",
            )


def choose_objects(selected_rows, all_rows, center_x, params, legacy_choose):
    # Exact legacy routing for the original default configuration.
    custom = (
        params.max_number != 3
        or params.min_number != 3
        or params.candidate_dark_mean_tolerance_8bit is not None
        or params.candidate_min_basin_contrast_8bit is not None
        or params.space_reference_nm is not None
    )
    if not custom:
        return legacy_choose(selected_rows, all_rows, center_x, params)

    def pos(r):
        return float(
            r.get("candidate_center_x_px", r.get("basin_center_x_px", math.nan))
        )

    pool = [r for r in all_rows if np.isfinite(pos(r)) and r.get("basin_valid", False)]
    if not pool:
        return []
    dark = [r for r in pool if r.get("brightness_true_trench", False)]
    # Brightness classified target set is preferred, never fill quota with a
    # known opposite-polarity region simply to satisfy max-number.
    if dark:
        pool = dark
    center = min(pool, key=lambda r: abs(pos(r) - center_x))
    cp = pos(center)
    if params.candidate_dark_mean_tolerance_8bit is not None:
        ref = float(
            center.get(
                "brightness_core_mean_raw", center.get("basin_core_mean_raw", math.nan)
            )
        )
        if np.isfinite(ref):
            pool = [
                r
                for r in pool
                if float(r.get("brightness_core_mean_raw", ref))
                <= ref + params.candidate_dark_mean_tolerance_8bit
            ]
    chosen = [center]
    remaining = [r for r in pool if r is not center]
    pitch = (
        (params.target_cd_nm + params.space_reference_nm) / params.pixel_size_nm
        if params.space_reference_nm
        else None
    )

    def key(r):
        distance = abs(pos(r) - cp)
        mismatch = abs(distance / pitch - round(distance / pitch)) if pitch else 0.0
        # Prior pitch influences neighbor ranking without moving an edge.
        return (distance + (pitch or 0) * mismatch, distance, pos(r))

    for row in sorted(remaining, key=key):
        if len(chosen) >= params.max_number:
            break
        chosen.append(row)
    left = sorted([r for r in chosen if pos(r) < cp], key=lambda r: cp - pos(r))
    right = sorted([r for r in chosen if pos(r) > cp], key=lambda r: pos(r) - cp)
    named = [("center", center)]
    for side, rows in [("left", left), ("right", right)]:
        for i, row in enumerate(rows, 1):
            named.append((side if i == 1 else f"{side}{i}", row))
    return named
