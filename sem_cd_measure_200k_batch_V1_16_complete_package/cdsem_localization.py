"""V1_16 locator policies; final edge extraction stays in cdsem_engine.

dark-line delegates to the unchanged V1_14 basin selector. bright-line uses
two intensity populations and adjacent periods. auto looks for alternating
deep/shallow or wide/narrow valleys in the complete observation region.
"""

from __future__ import annotations

from dataclasses import replace
import numpy as np
from scipy.ndimage import gaussian_filter1d

import cdsem_locator as legacy


def interior_band_threshold(profile, target_px=None):
    """Locate repeated bands containing an intermediate-intensity interior.

    Accept three classes only if the upper threshold joins two low plateaus
    inside at least two bounded intervals. Pixels used for edges stay intact.
    Operates on measurement polarity, preserving image-inversion symmetry.
    """
    profile = np.asarray(profile, float)
    _, levels, _ = legacy.brightness_kmeans_1d_general(profile, 3, 1, 42)
    dynamic = levels[2] - levels[0]
    if dynamic < 12 or np.min(np.diff(levels)) < 0.20 * dynamic:
        return None
    lower, upper = (levels[:-1] + levels[1:]) / 2
    low_runs = legacy._false_runs(profile >= lower)
    merged = []
    for start, end in legacy._false_runs(profile >= upper):
        if start == 0 or end == profile.size:
            continue
        children = [(a, b) for a, b in low_runs if start <= a < b <= end and b - a >= 3]
        if len(children) != 2:
            continue
        (a, b), (c, d) = children
        width = end - start
        if c - b < max(3, 0.08 * width) or min(b - a, d - c) < 0.12 * width:
            continue
        if target_px is not None and not 0.65 * target_px <= width <= 1.35 * target_px:
            continue
        merged.append(width)
    return float(upper) if len(merged) >= 2 else None


def observation_mask(mask, target_px):
    """Include actual line interiors between nearby supported trenches.

    The measurement mask contains edge context, but can omit a wide line center
    where a fake dark stripe lives. Bridge only bounded horizontal gaps <= 2.5
    target widths for morphology/locating; exterior and whole-row gaps stay out.
    This mask NEVER authorizes extra measurement samples.
    """
    observed = np.asarray(mask, bool).copy()
    maximum_gap = max(1, int(round(2.5 * target_px)))
    for row in observed:
        for start, end in legacy._false_runs(row):
            if start > 0 and end < row.size and end - start <= maximum_gap:
                row[start:end] = True
    return observed


def supported_profile(image, mask=None):
    """Average only observed pixels. Interpolation is for coarse locating only."""
    image = np.asarray(image, dtype=float)
    if mask is None:
        return image.mean(axis=0)
    counts = mask.sum(axis=0)
    valid = counts > 0
    if valid.sum() < 2:
        raise ValueError("有效观测区域不足，无法建立横向灰度轮廓")
    sums = np.where(mask, image, 0.0).sum(axis=0)
    x = np.arange(image.shape[1])
    return np.interp(x, x[valid], sums[valid] / counts[valid])


def masked_block_profile(image, mask, fallback):
    """Keep the original profile function when background selection is disabled."""
    if mask is None:
        return fallback
    reference = supported_profile(image, mask)

    def block(data, y, height):
        start = max(0, int(round(y)) - int(height) // 2)
        end = min(data.shape[0], start + int(height))
        local = mask[start:end]
        counts = local.sum(axis=0)
        sums = np.where(local, data[start:end], 0.0).sum(axis=0)
        return np.divide(sums, counts, out=reference.copy(), where=counts > 0)

    return block


def analyze_mode(image, mask=None, requested="auto", majority=0.70):
    """Deterministic full-region morphology analysis, independent of CD engines.

    Only consecutive observed basins contribute: unknown x gaps cannot become
    evidence for a periodic pattern. At least four valleys support dark-line.
    Equal alternating valleys are intrinsically ambiguous and default bright-line.
    """
    raw = supported_profile(image, mask)
    profile = gaussian_filter1d(raw, 1.5, mode="nearest")
    observed = np.ones(raw.size, bool) if mask is None else mask.any(axis=0)
    values = profile[observed]
    if values.size < 16:
        raise ValueError("横向有效观测区域少于16px")
    _, levels, _ = legacy.brightness_kmeans_1d_general(values, 3, 16, 42)
    threshold = float((levels[1] + levels[2]) / 2)
    dynamic = max(float(levels[2] - levels[0]), 1e-9)
    bright = legacy._remove_short_true_runs(
        legacy._fill_short_false_runs(profile >= threshold, 1), 2
    )
    basins = [
        (a, b)
        for a, b in legacy._false_runs(bright)
        if a > 0 and b < raw.size and b - a >= 3 and observed[a:b].mean() >= majority
    ]
    centers = np.array([(a + b - 1) / 2 for a, b in basins])
    widths = np.array([b - a for a, b in basins], dtype=float)
    cores = np.array(
        [np.median(profile[a + (b - a) // 4 : b - (b - a) // 4]) for a, b in basins]
    )
    scores, periods = [], []
    for i in range(max(0, len(basins) - 3)):
        # A complete pair of periods, with no unobserved background gap.
        a, b = basins[i][0], basins[i + 3][1]
        # Small unsupported line centers are allowed; large gaps are not.
        if any(
            e - s > 0.5 * np.median(np.diff(centers[i : i + 4]))
            for s, e in legacy._false_runs(observed[a:b])
        ):
            continue
        c, w, p = cores[i : i + 4], widths[i : i + 4], centers[i : i + 4]
        pitches = p[2:] - p[:-2]
        repeat_error = abs(pitches[0] - pitches[1]) / max(np.mean(pitches), 1)
        parity_error = (abs(c[0] - c[2]) + abs(c[1] - c[3])) / (2 * dynamic)
        alternating_depth = abs(np.mean(c[::2]) - np.mean(c[1::2])) / dynamic
        # Width evidence helps when a line's stain is as dark as the trench.
        width_gap = abs(np.mean(w[::2]) - np.mean(w[1::2])) / max(w.max(), 1)
        width_error = (abs(w[0] - w[2]) + abs(w[1] - w[3])) / max(2 * w.max(), 1)
        evidence = max(
            alternating_depth - 2 * parity_error,
            0.30 * max(0, width_gap - 2 * width_error),
        )
        scores.append(max(0, evidence) * max(0, 1 - 3 * repeat_error))
        periods.append(float(np.mean(pitches)))
    score = float(np.median(scores)) if scores else 0.0
    suggested = "dark-line" if score >= 0.12 else "bright-line"
    selected = suggested if requested == "auto" else requested
    pitch = (
        float(np.median(periods))
        if selected == "dark-line" and periods
        else float(np.median(np.diff(centers)))
        if len(centers) > 1
        else None
    )
    return dict(
        locator_mode_requested=requested,
        locator_mode_selected=selected,
        locator_mode_suggested=suggested,
        locator_dark_line_score=score,
        locator_auto_decision_threshold=0.12,
        locator_auto_ambiguous=bool(len(basins) < 4 or 0.08 <= score <= 0.16),
        locator_auto_reason=(
            "alternating_valleys"
            if suggested == "dark-line"
            else "uniform_bright_separators_or_insufficient_alternation"
        ),
        locator_profile_valley_count=len(basins),
        locator_period_px=pitch,
        locator_adjacent_period_px=(
            float(np.median(periods))
            if suggested == "dark-line" and periods
            else float(np.median(np.diff(centers)))
            if len(centers) > 1
            else None
        ),
        locator_profile_threshold_raw=threshold,
        locator_profile_low_raw=float(levels[0]),
        locator_profile_high_raw=float(levels[2]),
        locator_profile_observed_columns=int(observed.sum()),
    )


def select_bright_line_trenches(
    image, summary, fixed, block_profile_fn, majority=0.70, mask=None, target_px=None
):
    """Two adaptive intensity classes: every complete dark interval is a trench.

    Reuses V1_14 representative profiles, basin validation, brightness diagnostics
    and candidate geometry. Only intensity classification and period identity differ.
    """
    summary = dict(summary or {})
    raw, smooth, profile_diag = legacy.brightness_build_representative_profile(
        image, fixed, summary, block_profile_fn
    )
    observed = np.ones(raw.size, bool) if mask is None else mask.any(axis=0)
    _, levels, sse = legacy.brightness_kmeans_1d_general(
        smooth[observed], 2, fixed.separator_kmeans_restarts, fixed.random_seed
    )
    threshold = float(np.mean(levels))
    interior_threshold = interior_band_threshold(smooth, target_px)
    if interior_threshold is not None:
        threshold = interior_threshold
    dynamic = max(float(levels[1] - levels[0]), 1e-9)
    bright = legacy._remove_short_true_runs(
        legacy._fill_short_false_runs(smooth >= threshold, 1), 2
    )
    basins = legacy._false_runs(bright)
    bright_runs = legacy._true_runs(bright)
    rows, basin_diag = legacy.build_direct_basin_candidates_v115(
        image, raw, bright, basins, threshold, fixed, block_profile_fn
    )
    for row, (left, right) in zip(rows, basins):
        row.update(
            legacy.brightness_candidate_average_brightness(
                image, row, fixed, summary, block_profile_fn
            )
        )
        inset = max(1, int(0.15 * (right - left)))
        core = image[:, left + inset : right - inset]
        core_mask = (
            np.ones(core.shape, bool)
            if mask is None
            else mask[:, left + inset : right - inset]
        )
        pixels = core[core_mask]
        fraction = float(np.mean(pixels < threshold)) if pixels.size else 0.0
        dispersion = (
            float(np.median(np.abs(pixels - np.median(pixels)))) / dynamic
            if pixels.size
            else 1.0
        )
        # Test real pixels inside the flanking line plateaus as well. These
        # pixels may lie beyond the mask's narrow edge-context margin; use
        # only y rows supported by this trench, never exterior y background.
        valid_y = core_mask.any(axis=1)
        side_fractions = []
        for a, b in bright_runs:
            if b == left or a == right:
                inset_bright = max(1, int(0.15 * (b - a)))
                plateau = image[valid_y, a + inset_bright : b - inset_bright]
                side_fractions.append(
                    float(np.mean(plateau >= threshold)) if plateau.size else 0.0
                )
        line_fraction = min(side_fractions) if len(side_fractions) == 2 else 0.0
        uniform = (
            fraction >= majority and line_fraction >= majority and dispersion <= 0.25
        )
        row.update(
            brightness_threshold_raw=threshold,
            brightness_core_dark_fraction=fraction,
            brightness_line_bright_fraction=line_fraction,
            brightness_core_relative_mad=dispersion,
            brightness_true_trench=bool(row["basin_valid"] and uniform),
            brightness_darkness_rank=1.0,
            locator_mode_selected="bright-line",
            locator_interior_band_merged=interior_threshold is not None,
        )
        if not uniform:
            row["basin_valid"] = False
            row["basin_invalid_reason"] += ";nonuniform_or_not_majority_dark"

    eligible = [r for r in rows if r["basin_valid"] and r["brightness_true_trench"]]
    # Adjacent dark runs define one period. Reject isolated spurious candidates
    # with no neighbour at the dominant pitch (or an integer missing period).
    ordered = sorted(eligible, key=lambda r: r["candidate_center_x_px"])
    positions = np.array([r["candidate_center_x_px"] for r in ordered])
    pitch = float(np.median(np.diff(positions))) if len(positions) > 1 else np.nan
    if len(ordered) >= 3 and pitch > 0:
        for i, row in enumerate(ordered):
            distances = abs(positions - positions[i]) / pitch
            neighbours = (distances >= 0.75) & (
                abs(distances - np.round(distances)) <= 0.25
            )
            if not neighbours.any():
                row["basin_valid"] = row["brightness_true_trench"] = False
                row["basin_invalid_reason"] += ";period_inconsistent"
    eligible = [r for r in ordered if r["basin_valid"]]
    selected = sorted(
        eligible, key=lambda r: abs(r["candidate_center_x_px"] - fixed.center_x_px)
    )[: fixed.preferred_trenches]
    for row in selected:
        row["selected"] = row["brightness_selected"] = True
    return (
        selected,
        rows,
        dict(
            **summary,
            **profile_diag,
            **basin_diag,
            detection_success=len(selected) >= fixed.min_candidate_trenches,
            reason=""
            if len(selected) >= fixed.min_candidate_trenches
            else "bright-line有效周期不足",
            edge_mode="bright_line_adjacent_periods",
            separator_bright_threshold_raw=threshold,
            separator_low_center_raw=float(levels[0]),
            separator_high_center_raw=float(levels[1]),
            separator_kmeans_sse=float(sse),
            pitch_mean_px=pitch,
        ),
    )


def mode_params(params):
    """Use the same downstream object chooser with the appropriate basin gap."""
    return (
        replace(params, locator_pattern_preferred_basin_gap=1)
        if (params.locator_mode == "bright-line")
        else params
    )
