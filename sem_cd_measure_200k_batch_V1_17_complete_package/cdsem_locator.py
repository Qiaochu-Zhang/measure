"""V1_17: consolidated active implementation, derived from V13_modified."""

from __future__ import annotations
import math
import numpy as np

import itertools
from typing import Any, Callable, Sequence
from scipy.ndimage import gaussian_filter1d


brightness_EPS = 1e-12


def brightness__finite(values: Sequence[float] | np.ndarray) -> np.ndarray:
    data = np.asarray(values, dtype=float)
    return data[np.isfinite(data)]


def brightness__mad(
    values: Sequence[float] | np.ndarray,
    floor: float = 1e-9,
) -> float:
    data = brightness__finite(values)
    if data.size == 0:
        return float(floor)

    median = float(np.median(data))
    value = 1.4826 * float(np.median(np.abs(data - median)))

    if not np.isfinite(value) or value < floor:
        q25, q75 = np.percentile(data, [25, 75])
        value = float((q75 - q25) / 1.349)

    if not np.isfinite(value) or value < floor:
        value = float(np.std(data))

    return float(max(value, floor))


def brightness__trimmed_mean(
    values: Sequence[float] | np.ndarray,
    trim_fraction: float,
) -> float:
    data = np.sort(brightness__finite(values))
    if data.size == 0:
        return float("nan")

    trim = int(math.floor(data.size * float(np.clip(trim_fraction, 0.0, 0.40))))

    if 2 * trim >= data.size:
        return float(np.mean(data))

    if trim > 0:
        data = data[trim:-trim]

    return float(np.mean(data))


def brightness__rank01(
    values: Sequence[float],
    higher_is_better: bool,
) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    output = np.zeros(len(x), dtype=float)
    finite_mask = np.isfinite(x)

    if not np.any(finite_mask):
        return output

    finite_indices = np.where(finite_mask)[0]
    finite_values = x[finite_mask]
    order = np.argsort(finite_values)
    sorted_values = finite_values[order]
    ranks = np.empty(len(finite_values), dtype=float)

    start = 0
    while start < len(finite_values):
        end = start + 1
        while (
            end < len(finite_values)
            and abs(sorted_values[end] - sorted_values[start]) <= 1e-12
        ):
            end += 1

        ranks[order[start:end]] = 0.5 * (start + end - 1)
        start = end

    if len(finite_values) == 1:
        scaled = np.ones(1, dtype=float)
    else:
        scaled = ranks / (len(finite_values) - 1)

    if not higher_is_better:
        scaled = 1.0 - scaled

    output[finite_indices] = scaled
    return output


def brightness_build_representative_profile(
    measurement_image: np.ndarray,
    fixed: Any,
    detection_summary: dict[str, Any],
    block_profile_fn: Callable[[np.ndarray, float, int], np.ndarray],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """
    沿y方向取多个block profile，逐x做截尾平均。

    raw_profile:
        原始float灰度代表profile，用于平均亮度诊断。

    smooth_profile:
        仅用于自适应宽度与暗区连通区域，不用于最终CD/LWR。
    """
    image = np.asarray(measurement_image, dtype=np.float64)
    height, _ = image.shape

    center_y = (
        float(fixed.center_y_px)
        if fixed.center_y_px is not None
        else (height - 1) / 2.0
    )

    half_extend = min(
        fixed.extend_length_nm / fixed.pixel_size_nm / 2.0,
        center_y - 3.0,
        height - 4.0 - center_y,
    )

    default_y0 = max(0, int(round(center_y - half_extend)))
    default_y1 = min(height, int(round(center_y + half_extend)))

    y0 = int(detection_summary.get("y_start", default_y0))
    y1 = int(detection_summary.get("y_end", default_y1))
    y0 = int(np.clip(y0, 0, height - 1))
    y1 = int(np.clip(y1, y0 + 1, height))

    sample_count = max(15, int(fixed.width_profile_samples))
    block_height = max(1, int(fixed.width_profile_block_height_px))

    half_block = max(0.5, block_height / 2.0)
    safe_y0 = min(max(y0 + half_block, 0.0), height - 1.0)
    safe_y1 = max(min(y1 - half_block, height - 1.0), safe_y0)

    y_positions = np.linspace(
        safe_y0,
        safe_y1,
        sample_count,
    )

    profiles: list[np.ndarray] = []

    for y in y_positions:
        profile = np.asarray(
            block_profile_fn(
                image,
                float(y),
                block_height,
            ),
            dtype=float,
        )

        if profile.ndim == 1 and np.all(np.isfinite(profile)):
            profiles.append(profile)

    if len(profiles) < max(7, sample_count // 3):
        raise RuntimeError(f"代表profile有效采样不足：{len(profiles)}/{sample_count}")

    stack = np.vstack(profiles)
    sorted_stack = np.sort(stack, axis=0)

    trim_fraction = float(
        np.clip(
            fixed.width_profile_trim_fraction,
            0.0,
            0.30,
        )
    )
    trim = int(math.floor(sorted_stack.shape[0] * trim_fraction))

    if trim > 0 and 2 * trim < sorted_stack.shape[0]:
        trimmed = sorted_stack[trim:-trim]
    else:
        trimmed = sorted_stack

    raw_profile = np.mean(trimmed, axis=0)

    sigma = max(0.0, float(fixed.width_profile_sigma_px))
    smooth_profile = (
        gaussian_filter1d(raw_profile, sigma, mode="nearest")
        if sigma > 0
        else raw_profile.copy()
    )

    diagnostics = {
        "width_profile_requested_samples": sample_count,
        "width_profile_valid_samples": len(profiles),
        "width_profile_y0": y0,
        "width_profile_y1": y1,
        "width_profile_sigma_px": sigma,
    }

    return raw_profile, smooth_profile, diagnostics


def brightness_candidate_average_brightness(
    measurement_image: np.ndarray,
    candidate: dict[str, Any],
    fixed: Any,
    detection_summary: dict[str, Any],
    block_profile_fn: Callable[[np.ndarray, float, int], np.ndarray],
) -> dict[str, Any]:
    image = np.asarray(measurement_image, dtype=np.float64)
    height, width = image.shape

    center_x = float(candidate.get("candidate_center_x_px", np.nan))
    candidate_width_px = float(candidate.get("estimated_width_px", np.nan))

    if (
        not np.isfinite(center_x)
        or not np.isfinite(candidate_width_px)
        or candidate_width_px < 3.0
    ):
        return {
            "brightness_valid": False,
            "brightness_invalid_reason": ("invalid_adaptive_geometry"),
        }

    center_y = (
        float(fixed.center_y_px)
        if fixed.center_y_px is not None
        else (height - 1) / 2.0
    )

    half_extend = min(
        fixed.extend_length_nm / fixed.pixel_size_nm / 2.0,
        center_y - 3.0,
        height - 4.0 - center_y,
    )

    default_y0 = max(0, int(round(center_y - half_extend)))
    default_y1 = min(height, int(round(center_y + half_extend)))

    y0 = int(detection_summary.get("y_start", default_y0))
    y1 = int(detection_summary.get("y_end", default_y1))
    y0 = int(np.clip(y0, 0, height - 1))
    y1 = int(np.clip(y1, y0 + 1, height))

    sample_count = max(9, int(fixed.brightness_samples))
    block_height = max(1, int(fixed.brightness_block_height_px))

    half_block = max(0.5, block_height / 2.0)
    safe_y0 = min(max(y0 + half_block, 0.0), height - 1.0)
    safe_y1 = max(min(y1 - half_block, height - 1.0), safe_y0)

    y_positions = np.linspace(
        safe_y0,
        safe_y1,
        sample_count,
    )

    core_fraction = float(
        np.clip(
            fixed.brightness_core_fraction,
            0.30,
            0.95,
        )
    )
    core_half_width = max(
        1.5,
        0.5 * candidate_width_px * core_fraction,
    )

    core0 = int(
        np.clip(
            math.floor(center_x - core_half_width),
            0,
            width,
        )
    )
    core1 = int(
        np.clip(
            math.ceil(center_x + core_half_width) + 1,
            0,
            width,
        )
    )

    if core1 - core0 < 3:
        return {
            "brightness_valid": False,
            "brightness_invalid_reason": ("brightness_core_too_narrow"),
        }

    per_y_means: list[float] = []
    all_core_pixels: list[np.ndarray] = []

    for y in y_positions:
        profile = np.asarray(
            block_profile_fn(
                image,
                float(y),
                block_height,
            ),
            dtype=float,
        )

        core_values = brightness__finite(profile[core0:core1])

        if core_values.size < 3:
            continue

        per_y_means.append(float(np.mean(core_values)))
        all_core_pixels.append(core_values)

    valid_count = len(per_y_means)

    if valid_count < max(5, sample_count // 3):
        return {
            "brightness_valid": False,
            "brightness_valid_sample_count": valid_count,
            "brightness_sample_count": sample_count,
            "brightness_invalid_reason": ("brightness_insufficient_valid_samples"),
        }

    concatenated = np.concatenate(all_core_pixels)

    trim_fraction = float(
        np.clip(
            fixed.brightness_trim_fraction,
            0.0,
            0.30,
        )
    )

    core_mean_raw = brightness__trimmed_mean(
        concatenated,
        trim_fraction,
    )

    image_values = brightness__finite(image)
    image_p05, image_p95 = np.percentile(
        image_values,
        [5, 95],
    )
    image_dynamic = max(
        float(image_p95 - image_p05),
        1e-9,
    )
    core_mean_normalized = (core_mean_raw - float(image_p05)) / image_dynamic

    per_y_array = np.asarray(per_y_means, dtype=float)
    y_mean_median = float(np.median(per_y_array))
    y_mean_mad = brightness__mad(
        per_y_array,
        floor=max(abs(y_mean_median) * 1e-7, 1e-9),
    )

    return {
        "brightness_valid": True,
        "brightness_valid_sample_count": valid_count,
        "brightness_sample_count": sample_count,
        "brightness_valid_fraction": valid_count / sample_count,
        "brightness_core_x0_px": core0,
        "brightness_core_x1_px": core1,
        "brightness_core_width_px": core1 - core0,
        "brightness_core_mean_raw": float(core_mean_raw),
        "brightness_core_mean_normalized": float(core_mean_normalized),
        "brightness_core_median_raw": float(np.median(concatenated)),
        "brightness_core_p10_raw": float(np.percentile(concatenated, 10)),
        "brightness_core_p90_raw": float(np.percentile(concatenated, 90)),
        "brightness_y_mean_median_raw": y_mean_median,
        "brightness_y_mean_mad_raw": y_mean_mad,
        "brightness_y_relative_variation": float(y_mean_mad / image_dynamic),
        "brightness_invalid_reason": "",
    }


def brightness__dark_cluster_probability(
    brightness: np.ndarray,
    centers: np.ndarray,
    dark_label: int,
) -> np.ndarray:
    values = np.asarray(brightness, dtype=float)
    centers = np.asarray(centers, dtype=float)

    spread = brightness__mad(
        values,
        floor=max(
            abs(float(np.mean(values))) * 1e-7,
            1e-9,
        ),
    )

    distance_squared = (values[:, None] - centers[None, :]) ** 2

    logits = -distance_squared / (2.0 * spread**2)
    logits -= np.max(logits, axis=1, keepdims=True)

    probability = np.exp(logits)
    probability /= np.sum(probability, axis=1, keepdims=True)

    return probability[:, int(dark_label)]


def brightness_kmeans_1d_general(
    values: Sequence[float],
    cluster_count: int,
    restarts: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """通用一维K-means，输出中心按从暗到亮排序。"""
    data = np.asarray(values, dtype=float)
    cluster_count = int(cluster_count)

    if data.ndim != 1:
        raise ValueError("kmeans_1d_general输入必须是一维")
    if len(data) < cluster_count:
        raise ValueError(f"至少需要{cluster_count}个数据，实际{len(data)}")
    if not np.all(np.isfinite(data)):
        raise ValueError("K-means输入包含非有限值")

    if float(np.max(data) - np.min(data)) <= brightness_EPS:
        labels = np.arange(len(data), dtype=int) % cluster_count
        centers = np.full(cluster_count, float(data[0]))
        return labels, centers, 0.0

    rng = np.random.default_rng(seed)
    quantiles = np.linspace(0, 100, cluster_count + 2)[1:-1]

    initial_sets: list[np.ndarray] = [
        np.percentile(data, quantiles).astype(float),
        np.linspace(
            float(np.min(data)),
            float(np.max(data)),
            cluster_count,
        ),
    ]

    for _ in range(max(0, int(restarts) - len(initial_sets))):
        indices = rng.choice(
            len(data),
            size=cluster_count,
            replace=False,
        )
        initial_sets.append(np.sort(data[indices].astype(float)))

    best_labels: np.ndarray | None = None
    best_centers: np.ndarray | None = None
    best_sse = float("inf")

    for initial in initial_sets:
        centers = np.asarray(initial, dtype=float).copy()
        labels = np.zeros(len(data), dtype=int)

        for _ in range(150):
            distances = np.abs(data[:, None] - centers[None, :])
            new_labels = np.argmin(distances, axis=1)

            new_centers = centers.copy()
            valid = True
            for cluster in range(cluster_count):
                members = data[new_labels == cluster]
                if members.size == 0:
                    valid = False
                    break
                new_centers[cluster] = float(np.mean(members))

            if not valid:
                break

            converged = (
                np.array_equal(new_labels, labels)
                and np.max(np.abs(new_centers - centers)) < 1e-9
            )
            labels = new_labels
            centers = new_centers

            if converged:
                break

        counts = np.bincount(labels, minlength=cluster_count)
        if np.any(counts == 0):
            continue

        sse = float(np.sum((data - centers[labels]) ** 2))
        if sse < best_sse:
            best_sse = sse
            best_labels = labels.copy()
            best_centers = centers.copy()

    if best_labels is None or best_centers is None:
        order = np.argsort(data)
        split_groups = np.array_split(order, cluster_count)
        best_labels = np.zeros(len(data), dtype=int)

        for cluster, indices in enumerate(split_groups):
            best_labels[indices] = cluster

        best_centers = np.asarray(
            [
                float(np.mean(data[best_labels == cluster]))
                for cluster in range(cluster_count)
            ]
        )
        best_sse = float(np.sum((data - best_centers[best_labels]) ** 2))

    center_order = np.argsort(best_centers)
    remap = np.empty(cluster_count, dtype=int)
    for new_label, old_label in enumerate(center_order):
        remap[old_label] = new_label

    return (
        remap[best_labels],
        best_centers[center_order],
        best_sse,
    )


EPS = 1e-12


def _finite(values: Sequence[float] | np.ndarray) -> np.ndarray:
    data = np.asarray(values, dtype=float)
    return data[np.isfinite(data)]


def _mad(values: Sequence[float] | np.ndarray, floor: float = 1e-9) -> float:
    data = _finite(values)
    if data.size == 0:
        return float(floor)
    median = float(np.median(data))
    value = 1.4826 * float(np.median(np.abs(data - median)))
    if not np.isfinite(value) or value < floor:
        q25, q75 = np.percentile(data, [25, 75])
        value = float((q75 - q25) / 1.349)
    if not np.isfinite(value) or value < floor:
        value = float(np.std(data))
    return float(max(value, floor))


def _true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    data = np.asarray(mask, dtype=bool)
    runs: list[tuple[int, int]] = []
    i = 0
    while i < len(data):
        if not data[i]:
            i += 1
            continue
        start = i
        while i < len(data) and data[i]:
            i += 1
        runs.append((start, i))
    return runs


def _false_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    data = np.asarray(mask, dtype=bool)
    runs: list[tuple[int, int]] = []
    i = 0
    while i < len(data):
        if data[i]:
            i += 1
            continue
        start = i
        while i < len(data) and not data[i]:
            i += 1
        runs.append((start, i))
    return runs


def _fill_short_false_runs(mask: np.ndarray, maximum_gap: int) -> np.ndarray:
    output = np.asarray(mask, dtype=bool).copy()
    maximum_gap = max(0, int(maximum_gap))
    if maximum_gap <= 0:
        return output
    for start, end in _false_runs(output):
        if (
            end - start <= maximum_gap
            and start > 0
            and end < len(output)
            and output[start - 1]
            and output[end]
        ):
            output[start:end] = True
    return output


def _remove_short_true_runs(mask: np.ndarray, minimum_run: int) -> np.ndarray:
    output = np.asarray(mask, dtype=bool).copy()
    minimum_run = max(1, int(minimum_run))
    for start, end in _true_runs(output):
        if end - start < minimum_run:
            output[start:end] = False
    return output


def _nearest_false_index(mask: np.ndarray, x: int, radius: int) -> int | None:
    """在中心附近寻找最近的非高亮像素。"""
    data = np.asarray(mask, dtype=bool)
    if len(data) == 0:
        return None
    x = int(np.clip(x, 0, len(data) - 1))
    if not data[x]:
        return x
    for distance in range(1, max(1, int(radius)) + 1):
        left = x - distance
        right = x + distance
        if left >= 0 and not data[left]:
            return left
        if right < len(data) and not data[right]:
            return right
    return None


def _false_run_containing(mask: np.ndarray, x: int) -> tuple[int, int] | None:
    data = np.asarray(mask, dtype=bool)
    if x < 0 or x >= len(data) or data[x]:
        return None
    left = x
    while left > 0 and not data[left - 1]:
        left -= 1
    right = x + 1
    while right < len(data) and not data[right]:
        right += 1
    return left, right


def _interval_iou(first: tuple[float, float], second: tuple[float, float]) -> float:
    a0, a1 = first
    b0, b1 = second
    intersection = max(0.0, min(a1, b1) - max(a0, b0))
    union = max(a1, b1) - min(a0, b0)
    return float(intersection / max(union, EPS))


def _measurement_y_range(image: np.ndarray, fixed: Any) -> tuple[int, int]:
    height = int(image.shape[0])
    center_y = (
        float(fixed.center_y_px)
        if fixed.center_y_px is not None
        else (height - 1) / 2.0
    )
    half_extend = min(
        float(fixed.extend_length_nm) / float(fixed.pixel_size_nm) / 2.0,
        center_y - 3.0,
        height - 4.0 - center_y,
    )
    y0 = max(0, int(round(center_y - half_extend)))
    y1 = min(height, int(round(center_y + half_extend)))
    if y1 <= y0:
        y0, y1 = 0, height
    return int(y0), int(y1)


def _default_typical_width_px(fixed: Any) -> float:
    min_px = float(fixed.adaptive_min_width_nm) / float(fixed.pixel_size_nm)
    max_px = float(fixed.adaptive_max_width_nm) / float(fixed.pixel_size_nm)
    return float(math.sqrt(max(min_px, 3.0) * max(max_px, min_px + 1.0)))


def segment_bright_separators_v115(
    smooth_profile: np.ndarray,
    fixed: Any,
) -> tuple[
    np.ndarray,
    list[tuple[int, int]],
    list[tuple[int, int]],
    dict[str, Any],
]:
    """直接由同图代表 profile 得到高亮带和非高亮 Basin。

    与 V1.14 不同：不再需要 Stage-1 候选宽度作为尺度。
    第一次使用物理安全宽度的几何均值做轻量清洗；再从初步 Basin 宽度
    估计本图尺度，执行第二次自适应清洗。
    """
    profile = np.asarray(smooth_profile, dtype=float)
    if profile.ndim != 1 or profile.size < 16:
        raise ValueError("代表profile无效或过短")

    _, centers, sse = brightness_kmeans_1d_general(
        profile,
        cluster_count=3,
        restarts=max(6, int(fixed.separator_kmeans_restarts)),
        seed=int(fixed.random_seed),
    )
    low_center = float(centers[0])
    middle_center = float(centers[1])
    high_center = float(centers[2])
    bright_threshold = 0.5 * (middle_center + high_center)

    raw_bright_mask = profile >= bright_threshold
    default_width = _default_typical_width_px(fixed)

    # 初次轻量清洗，避免用尚未存在的暗候选作为尺度。
    first_min_run = max(
        int(getattr(fixed, "separator_min_run_px", 2)),
        int(round(default_width * float(fixed.separator_min_run_fraction))),
    )
    first_close_gap = max(
        int(getattr(fixed, "separator_close_gap_px", 1)),
        int(round(default_width * float(fixed.separator_close_gap_fraction))),
    )
    first_mask = _fill_short_false_runs(raw_bright_mask, first_close_gap)
    first_mask = _remove_short_true_runs(first_mask, first_min_run)
    first_basins = _false_runs(first_mask)

    min_safe_px = float(fixed.adaptive_min_width_nm) / float(fixed.pixel_size_nm)
    max_safe_px = float(fixed.adaptive_max_width_nm) / float(fixed.pixel_size_nm)
    provisional_widths = np.asarray(
        [
            end - start
            for start, end in first_basins
            if start > 0
            and end < len(profile)
            and 0.45 * min_safe_px <= end - start <= 1.60 * max_safe_px
        ],
        dtype=float,
    )
    typical_width = (
        float(np.median(provisional_widths))
        if provisional_widths.size >= 3
        else default_width
    )

    minimum_run = max(
        int(getattr(fixed, "separator_min_run_px", 2)),
        int(round(typical_width * float(fixed.separator_min_run_fraction))),
    )
    close_gap = max(
        int(getattr(fixed, "separator_close_gap_px", 1)),
        int(round(typical_width * float(fixed.separator_close_gap_fraction))),
    )

    bright_mask = _fill_short_false_runs(raw_bright_mask, close_gap)
    bright_mask = _remove_short_true_runs(bright_mask, minimum_run)
    bright_runs = _true_runs(bright_mask)
    basins = _false_runs(bright_mask)

    diagnostics = {
        "separator_low_center_raw": low_center,
        "separator_middle_center_raw": middle_center,
        "separator_high_center_raw": high_center,
        "separator_bright_threshold_raw": float(bright_threshold),
        "separator_kmeans_sse": float(sse),
        "separator_raw_bright_fraction": float(np.mean(raw_bright_mask)),
        "separator_clean_bright_fraction": float(np.mean(bright_mask)),
        "separator_first_minimum_run_px": int(first_min_run),
        "separator_first_close_gap_px": int(first_close_gap),
        "separator_minimum_run_px": int(minimum_run),
        "separator_close_gap_px": int(close_gap),
        "separator_typical_basin_width_px": float(typical_width),
        "separator_bright_run_count": len(bright_runs),
        "separator_basin_count": len(basins),
    }
    return bright_mask, bright_runs, basins, diagnostics


def validate_basin_across_y_v115(
    measurement_image: np.ndarray,
    basin: tuple[int, int],
    raw_representative_profile: np.ndarray,
    bright_threshold: float,
    fixed: Any,
    block_profile_fn: Callable[[np.ndarray, float, int], np.ndarray],
) -> dict[str, Any]:
    """检查 Basin 是否沿 Y 方向持续存在并且边界稳定。

    代表 profile 只负责生成 Basin；本函数回到多个真实 Y 位置验证：
    - Basin 中心附近是否仍存在非高亮连续区；
    - 左右名义分隔位置附近是否仍有高亮信号；
    - 局部 Basin 与名义 Basin 的重叠；
    - 中心漂移和宽度变化。
    """
    image = np.asarray(measurement_image, dtype=np.float64)
    height, width = image.shape
    left_nominal, right_nominal = map(int, basin)
    nominal_width = float(right_nominal - left_nominal)
    nominal_center = 0.5 * (left_nominal + right_nominal - 1)

    y0, y1 = _measurement_y_range(image, fixed)
    sample_count = max(9, int(getattr(fixed, "basin_validation_samples", 17)))
    block_height = max(1, int(getattr(fixed, "basin_validation_block_height_px", 5)))
    sigma = max(0.0, float(getattr(fixed, "basin_validation_sigma_px", 1.0)))

    half_block = max(0.5, block_height / 2.0)
    safe_y0 = min(max(y0 + half_block, 0.0), height - 1.0)
    safe_y1 = max(min(y1 - half_block, height - 1.0), safe_y0)
    y_positions = np.linspace(safe_y0, safe_y1, sample_count)

    representative_median = float(np.median(raw_representative_profile))
    center_values: list[float] = []
    width_values: list[float] = []
    core_values: list[float] = []
    left_separator_hits = 0
    right_separator_hits = 0
    overlap_hits = 0
    valid_profiles = 0

    separator_window = max(2, int(round(0.12 * nominal_width)))
    center_search_radius = max(2, int(round(0.15 * nominal_width)))
    core_fraction = float(np.clip(fixed.brightness_core_fraction, 0.30, 0.90))

    for y in y_positions:
        profile = np.asarray(
            block_profile_fn(image, float(y), block_height),
            dtype=float,
        )
        if (
            profile.ndim != 1
            or profile.size != width
            or not np.all(np.isfinite(profile))
        ):
            continue
        valid_profiles += 1
        smooth = (
            gaussian_filter1d(profile, sigma, mode="nearest") if sigma > 0 else profile
        )

        # 对每个Y位置只校正整体灰度偏移，不重新自由选择阈值，避免局部K-means抖动。
        local_threshold = float(
            bright_threshold + np.median(smooth) - representative_median
        )
        local_bright = smooth >= local_threshold

        left0 = max(0, left_nominal - separator_window)
        left1 = min(width, left_nominal + separator_window + 1)
        right0 = max(0, right_nominal - separator_window)
        right1 = min(width, right_nominal + separator_window + 1)
        if left1 > left0 and np.mean(local_bright[left0:left1]) >= 0.18:
            left_separator_hits += 1
        if right1 > right0 and np.mean(local_bright[right0:right1]) >= 0.18:
            right_separator_hits += 1

        anchor = _nearest_false_index(
            local_bright,
            int(round(nominal_center)),
            center_search_radius,
        )
        if anchor is None:
            continue
        run = _false_run_containing(local_bright, anchor)
        if run is None:
            continue
        local_left, local_right = run
        local_width = float(local_right - local_left)
        if local_width < 3.0:
            continue

        iou = _interval_iou(
            (float(left_nominal), float(right_nominal)),
            (float(local_left), float(local_right)),
        )
        if iou < float(getattr(fixed, "basin_validation_min_iou", 0.25)):
            continue
        overlap_hits += 1
        center_values.append(0.5 * (local_left + local_right - 1))
        width_values.append(local_width)

        local_center = center_values[-1]
        core_half = max(1.5, 0.5 * local_width * core_fraction)
        core0 = int(np.clip(math.floor(local_center - core_half), 0, width))
        core1 = int(np.clip(math.ceil(local_center + core_half) + 1, 0, width))
        if core1 - core0 >= 3:
            core_values.append(float(np.mean(profile[core0:core1])))

    denominator = max(valid_profiles, 1)
    center_array = np.asarray(center_values, dtype=float)
    width_array = np.asarray(width_values, dtype=float)
    core_array = np.asarray(core_values, dtype=float)

    presence_fraction = len(width_array) / denominator
    left_presence = left_separator_hits / denominator
    right_presence = right_separator_hits / denominator
    overlap_fraction = overlap_hits / denominator
    center_median = (
        float(np.median(center_array)) if center_array.size else nominal_center
    )
    width_median = float(np.median(width_array)) if width_array.size else nominal_width
    center_mad = (
        _mad(center_array, floor=0.0) if center_array.size >= 2 else float("inf")
    )
    width_mad = _mad(width_array, floor=0.0) if width_array.size >= 2 else float("inf")
    width_cv = (
        float(width_mad / max(width_median, EPS))
        if np.isfinite(width_mad)
        else float("inf")
    )

    image_values = _finite(image)
    image_dynamic = max(
        float(np.percentile(image_values, 95) - np.percentile(image_values, 5)), 1e-9
    )
    core_relative_variation = (
        float(_mad(core_array, floor=0.0) / image_dynamic)
        if core_array.size >= 2
        else float("inf")
    )

    return {
        "basin_validation_requested_samples": sample_count,
        "basin_validation_valid_profiles": valid_profiles,
        "basin_vertical_presence_fraction": float(presence_fraction),
        "basin_overlap_presence_fraction": float(overlap_fraction),
        "basin_left_separator_presence_fraction": float(left_presence),
        "basin_right_separator_presence_fraction": float(right_presence),
        "basin_center_median_px": float(center_median),
        "basin_width_median_px": float(width_median),
        "basin_center_mad_px": float(center_mad),
        "basin_width_mad_px": float(width_mad),
        "basin_width_cv_robust": float(width_cv),
        "basin_core_y_relative_variation": float(core_relative_variation),
    }


def build_direct_basin_candidates_v115(
    measurement_image: np.ndarray,
    raw_profile: np.ndarray,
    bright_mask: np.ndarray,
    basins: list[tuple[int, int]],
    bright_threshold: float,
    fixed: Any,
    block_profile_fn: Callable[[np.ndarray, float, int], np.ndarray],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """每个非高亮连续区直接建立一个候选 Basin，不经过暗候选合并。"""
    image_width = int(measurement_image.shape[1])
    min_width_px = (
        float(fixed.adaptive_min_width_nm)
        / float(fixed.pixel_size_nm)
        * float(getattr(fixed, "basin_width_min_factor", 0.45))
    )
    max_width_px = (
        float(fixed.adaptive_max_width_nm)
        / float(fixed.pixel_size_nm)
        * float(getattr(fixed, "basin_width_max_factor", 1.45))
    )
    reject_edge = bool(getattr(fixed, "basin_reject_edge_basins", True))

    rows: list[dict[str, Any]] = []
    valid_count = 0
    invalid_reason_counts: dict[str, int] = {}

    for basin_index, (left, right) in enumerate(basins):
        nominal_width = float(right - left)
        nominal_center = 0.5 * (left + right - 1)
        bounded_left = left > 0 and bool(bright_mask[left - 1])
        bounded_right = right < image_width and bool(bright_mask[right])

        row: dict[str, Any] = {
            "candidate_id": int(basin_index),
            "basin_sequence_index": int(basin_index),
            "basin_left_x_px": float(left),
            "basin_right_x_px": float(right),
            "basin_center_x_px": float(nominal_center),
            "basin_width_px": float(nominal_width),
            "basin_left_bounded_by_separator": bool(bounded_left),
            "basin_right_bounded_by_separator": bool(bounded_right),
            "basin_edge_clipped": bool(not (bounded_left and bounded_right)),
            "basin_direct_candidate": True,
            "basin_group_size": 1,
            "basin_duplicate_merged": False,
            "basin_source_centers_px": f"{nominal_center:.3f}",
            "adaptive_duplicate_merged": False,
            "adaptive_merge_representative": True,
            "adaptive_merge_group_id": int(basin_index),
            "adaptive_merge_group_size": 1,
            "adaptive_merge_source_centers_px": f"{nominal_center:.3f}",
            "adaptive_merge_source_widths_px": f"{nominal_width:.3f}",
            "adaptive_initial_center_x_px": float(nominal_center),
            "adaptive_initial_width_px": float(nominal_width),
            "adaptive_left_x_px": float(left),
            "adaptive_right_x_px": float(right),
            "adaptive_center_x_px": float(nominal_center),
            "adaptive_width_px": float(nominal_width),
            "adaptive_width_ratio": 1.0,
            "adaptive_width_valid": True,
            "candidate_center_x_px": float(nominal_center),
            "estimated_width_px": float(nominal_width),
            "estimated_width_nm": float(nominal_width * fixed.pixel_size_nm),
            "selected": False,
            "brightness_selected": False,
            "brightness_true_trench": False,
        }

        validation = validate_basin_across_y_v115(
            measurement_image,
            (left, right),
            raw_profile,
            bright_threshold,
            fixed,
            block_profile_fn,
        )
        row.update(validation)

        # 使用多Y中位几何作为下游边缘搜索初值，但不越过代表profile盆地边界太远。
        if np.isfinite(float(validation.get("basin_center_median_px", np.nan))):
            robust_center = float(validation["basin_center_median_px"])
            if abs(robust_center - nominal_center) <= max(4.0, 0.25 * nominal_width):
                row["candidate_center_x_px"] = robust_center
                row["adaptive_center_x_px"] = robust_center
        if np.isfinite(float(validation.get("basin_width_median_px", np.nan))):
            robust_width = float(validation["basin_width_median_px"])
            if 0.55 * nominal_width <= robust_width <= 1.55 * nominal_width:
                row["estimated_width_px"] = robust_width
                row["estimated_width_nm"] = robust_width * fixed.pixel_size_nm
                row["adaptive_width_px"] = robust_width

        reasons: list[str] = []
        if reject_edge and row["basin_edge_clipped"]:
            reasons.append("edge_clipped_basin")
        if nominal_width < min_width_px:
            reasons.append("basin_too_narrow")
        if nominal_width > max_width_px:
            reasons.append("basin_too_wide")

        presence = float(validation["basin_vertical_presence_fraction"])
        left_presence = float(validation["basin_left_separator_presence_fraction"])
        right_presence = float(validation["basin_right_separator_presence_fraction"])
        center_mad = float(validation["basin_center_mad_px"])
        width_cv = float(validation["basin_width_cv_robust"])

        min_presence = float(getattr(fixed, "basin_min_vertical_presence", 0.45))
        min_separator_presence = float(
            getattr(fixed, "basin_min_separator_presence", 0.35)
        )
        max_center_jitter = max(
            float(getattr(fixed, "basin_max_center_jitter_px", 8.0)),
            float(getattr(fixed, "basin_max_center_jitter_fraction", 0.22))
            * nominal_width,
        )
        max_width_cv = float(getattr(fixed, "basin_max_width_cv", 0.40))

        if presence < min_presence:
            reasons.append("basin_vertical_presence_low")
        if min(left_presence, right_presence) < min_separator_presence:
            reasons.append("separator_presence_low")
        if not np.isfinite(center_mad) or center_mad > max_center_jitter:
            reasons.append("basin_center_jitter_high")
        if not np.isfinite(width_cv) or width_cv > max_width_cv:
            reasons.append("basin_width_variation_high")

        width_center = math.sqrt(
            max(min_width_px, 1.0) * max(max_width_px, min_width_px + 1.0)
        )
        width_log_distance = abs(math.log(max(nominal_width, 1.0) / width_center))
        width_score = math.exp(-width_log_distance)
        jitter_score = (
            math.exp(-max(center_mad, 0.0) / max(0.18 * nominal_width, 2.0))
            if np.isfinite(center_mad)
            else 0.0
        )
        width_stability_score = (
            math.exp(-max(width_cv, 0.0) / 0.25) if np.isfinite(width_cv) else 0.0
        )
        separator_score = min(left_presence, right_presence)
        quality_score = float(
            np.clip(
                0.30 * presence
                + 0.20 * separator_score
                + 0.18 * width_score
                + 0.16 * jitter_score
                + 0.16 * width_stability_score,
                0.0,
                1.0,
            )
        )

        basin_valid = len(reasons) == 0
        row.update(
            {
                "basin_valid": basin_valid,
                "basin_invalid_reason": ";".join(reasons),
                "basin_quality_score": quality_score,
                "candidate_score": quality_score,
                "precheck_structure_score": quality_score,
                "adaptive_eligible": basin_valid,
                "eligible": basin_valid,
            }
        )
        if basin_valid:
            valid_count += 1
        else:
            for reason in reasons:
                invalid_reason_counts[reason] = invalid_reason_counts.get(reason, 0) + 1
        rows.append(row)

    diagnostics = {
        "basin_direct_candidate_count": len(rows),
        "basin_valid_candidate_count": valid_count,
        "basin_invalid_candidate_count": len(rows) - valid_count,
        "basin_invalid_reason_counts": ";".join(
            f"{key}:{value}" for key, value in sorted(invalid_reason_counts.items())
        ),
        "basin_width_min_allowed_px": float(min_width_px),
        "basin_width_max_allowed_px": float(max_width_px),
        # 兼容V1.14汇总字段；V1.15正常路径不存在合并。
        "basin_unique_candidate_count": len(rows),
        "basin_duplicate_merged_count": 0,
        "eligible_count": valid_count,
    }
    return rows, diagnostics


def _pattern_group_cost_v115(
    combination: Sequence[dict[str, Any]],
    image_center_x: float,
    preferred_basin_gap: int,
) -> float:
    ordered = sorted(combination, key=lambda row: float(row["candidate_center_x_px"]))
    positions = np.asarray([float(row["candidate_center_x_px"]) for row in ordered])
    widths = np.asarray([float(row["estimated_width_px"]) for row in ordered])
    basin_indices = np.asarray([int(row["basin_sequence_index"]) for row in ordered])
    basin_gaps = np.diff(basin_indices)
    if np.any(basin_gaps < 2):
        return float("inf")

    pitches = np.diff(positions)
    if pitches.size == 0 or np.min(pitches) <= 0:
        return float("inf")

    pitch_mean = float(np.mean(pitches))
    pitch_cv = float(np.std(pitches) / max(pitch_mean, EPS))
    pitch_ratio = float(np.max(pitches) / max(np.min(pitches), EPS))
    width_cv = float(np.std(widths) / max(np.mean(widths), EPS))
    center_offset = abs(float(np.median(positions)) - image_center_x)
    pattern_penalty = float(np.mean(np.abs(basin_gaps - int(preferred_basin_gap))))
    darkness_rank = np.asarray(
        [float(row["brightness_darkness_rank"]) for row in ordered]
    )
    quality = np.asarray(
        [float(row.get("basin_quality_score", 0.0)) for row in ordered]
    )
    vertical_presence = np.asarray(
        [float(row.get("basin_vertical_presence_fraction", 0.0)) for row in ordered]
    )

    return float(
        2.25 * pitch_cv
        + 0.75 * max(0.0, pitch_ratio - 1.0)
        + 0.25 * width_cv
        + 0.50 * center_offset / max(pitch_mean, 1.0)
        + 0.42 * pattern_penalty
        - 0.82 * float(np.mean(darkness_rank))
        - 0.42 * float(np.mean(quality))
        - 0.20 * float(np.mean(vertical_presence))
    )


def find_best_pattern_group_v115(
    pool: list[dict[str, Any]],
    preferred_count: int,
    minimum_count: int,
    image_center_x: float,
    preferred_basin_gap: int,
) -> tuple[list[dict[str, Any]] | None, float]:
    best_group: list[dict[str, Any]] | None = None
    best_cost = float("inf")
    maximum_count = min(int(preferred_count), len(pool))
    for count in range(maximum_count, int(minimum_count) - 1, -1):
        for combination in itertools.combinations(pool, count):
            cost = _pattern_group_cost_v115(
                combination,
                image_center_x,
                preferred_basin_gap,
            )
            if cost < best_cost:
                best_cost = cost
                best_group = list(combination)
        if best_group is not None:
            break
    return best_group, best_cost


def select_with_brightness_relaxation_v115(
    candidates: list[dict[str, Any]],
    initial_threshold: float,
    fixed: Any,
    image_center_x: float,
) -> tuple[list[dict[str, Any]] | None, float, int, int, float]:
    """按盆地亮度逐级提高上限，首次形成结构可行组合时停止。"""
    unique_brightness = np.sort(
        np.unique([float(row["brightness_core_mean_raw"]) for row in candidates])
    )
    threshold_sequence = [float(initial_threshold)]
    threshold_sequence.extend(
        [float(value) for value in unique_brightness if value > initial_threshold + EPS]
    )
    threshold_sequence = sorted(set(threshold_sequence))
    initial_pool_count = sum(
        float(row["brightness_core_mean_raw"]) <= initial_threshold + EPS
        for row in candidates
    )

    for step, threshold in enumerate(threshold_sequence):
        pool = [
            row
            for row in candidates
            if float(row["brightness_core_mean_raw"]) <= threshold + EPS
        ]
        if len(pool) < fixed.min_candidate_trenches:
            continue
        pool = sorted(
            pool,
            key=lambda row: (
                float(row["brightness_core_mean_raw"]),
                -float(row.get("basin_quality_score", 0.0)),
            ),
        )[: max(fixed.min_candidate_trenches, int(fixed.brightness_candidate_pool))]
        group, cost = find_best_pattern_group_v115(
            pool,
            preferred_count=fixed.preferred_trenches,
            minimum_count=fixed.min_candidate_trenches,
            image_center_x=image_center_x,
            preferred_basin_gap=int(fixed.pattern_preferred_basin_gap),
        )
        if group is not None:
            return (
                group,
                float(threshold),
                int(step),
                int(initial_pool_count),
                float(cost),
            )

    return (
        None,
        float(threshold_sequence[-1]),
        max(0, len(threshold_sequence) - 1),
        int(initial_pool_count),
        float("inf"),
    )


def select_basin_first_trenches(
    measurement_image: np.ndarray,
    detection_summary: dict[str, Any] | None,
    fixed: Any,
    block_profile_fn: Callable[[np.ndarray, float, int], np.ndarray],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """V1.15 Basin-First 主入口。"""
    image = np.asarray(measurement_image, dtype=np.float64)
    if image.ndim != 2:
        raise ValueError("V1.15只支持二维灰度SEM图像")

    summary = dict(detection_summary or {})
    y0, y1 = _measurement_y_range(image, fixed)
    summary.setdefault("y_start", y0)
    summary.setdefault("y_end", y1)
    summary["basin_first_enabled"] = True
    summary["basin_first_fallback_used"] = False
    summary["edge_mode"] = "v115_basin_first"

    if not bool(getattr(fixed, "basin_first_enabled", True)):
        return (
            [],
            [],
            {
                **summary,
                "detection_success": False,
                "reason": "basin_first_disabled",
            },
        )

    raw_profile, smooth_profile, profile_diagnostics = (
        brightness_build_representative_profile(
            image,
            fixed,
            summary,
            block_profile_fn,
        )
    )
    bright_mask, bright_runs, basins, separator_diagnostics = (
        segment_bright_separators_v115(
            smooth_profile,
            fixed,
        )
    )
    bright_threshold = float(separator_diagnostics["separator_bright_threshold_raw"])

    rows, basin_diagnostics = build_direct_basin_candidates_v115(
        image,
        raw_profile,
        bright_mask,
        basins,
        bright_threshold,
        fixed,
        block_profile_fn,
    )

    eligible = [row for row in rows if bool(row.get("basin_valid", False))]
    if len(eligible) < fixed.min_candidate_trenches:
        summary.update(
            {
                "detection_success": False,
                "reason": (
                    f"Basin-First有效盆地不足{fixed.min_candidate_trenches}条"
                    f"（实际{len(eligible)}条）"
                ),
                "brightness_selected_count": 0,
                **profile_diagnostics,
                **separator_diagnostics,
                **basin_diagnostics,
            }
        )
        return [], rows, summary

    brightness_eligible: list[dict[str, Any]] = []
    for row in eligible:
        metrics = brightness_candidate_average_brightness(
            image,
            row,
            fixed,
            summary,
            block_profile_fn,
        )
        row.update(metrics)
        if bool(row.get("brightness_valid", False)):
            brightness_eligible.append(row)

    if len(brightness_eligible) < fixed.min_candidate_trenches:
        summary.update(
            {
                "detection_success": False,
                "reason": (
                    f"Basin-First亮度有效盆地不足{fixed.min_candidate_trenches}条"
                    f"（实际{len(brightness_eligible)}条）"
                ),
                "brightness_selected_count": 0,
                **profile_diagnostics,
                **separator_diagnostics,
                **basin_diagnostics,
            }
        )
        return [], rows, summary

    brightness_values = np.asarray(
        [float(row["brightness_core_mean_raw"]) for row in brightness_eligible]
    )
    labels, centers, brightness_sse = brightness_kmeans_1d_general(
        brightness_values,
        cluster_count=2,
        restarts=max(6, int(fixed.brightness_kmeans_restarts)),
        seed=int(fixed.random_seed),
    )
    dark_center = float(centers[0])
    bright_center = float(centers[1])
    initial_threshold = 0.5 * (dark_center + bright_center)
    probabilities = brightness__dark_cluster_probability(
        brightness_values,
        centers,
        dark_label=0,
    )
    darkness_ranks = brightness__rank01(
        brightness_values,
        higher_is_better=False,
    )

    for index, row in enumerate(brightness_eligible):
        row.update(
            {
                "brightness_cluster_label": int(labels[index]),
                "brightness_dark_cluster_label": 0,
                "brightness_bright_cluster_label": 1,
                "brightness_initial_true_trench": bool(labels[index] == 0),
                "brightness_dark_probability": float(probabilities[index]),
                "brightness_darkness_rank": float(darkness_ranks[index]),
                "brightness_initial_threshold_raw": float(initial_threshold),
                "brightness_dark_cluster_center_raw": dark_center,
                "brightness_bright_cluster_center_raw": bright_center,
                "brightness_cluster_gap_raw": bright_center - dark_center,
                "brightness_kmeans_sse": float(brightness_sse),
            }
        )

    image_center_x = (
        float(fixed.center_x_px)
        if fixed.center_x_px is not None
        else (image.shape[1] - 1) / 2.0
    )
    best_group, final_threshold, relaxation_steps, initial_pool_count, group_cost = (
        select_with_brightness_relaxation_v115(
            brightness_eligible,
            initial_threshold,
            fixed,
            image_center_x,
        )
    )

    if best_group is None:
        summary.update(
            {
                "detection_success": False,
                "reason": "Basin-First逐级放宽后仍没有满足拓扑/周期约束的组合",
                "brightness_initial_threshold_raw": initial_threshold,
                "brightness_final_threshold_raw": final_threshold,
                "brightness_relaxation_steps": relaxation_steps,
                "brightness_initial_pool_count": initial_pool_count,
                "brightness_selected_count": 0,
                **profile_diagnostics,
                **separator_diagnostics,
                **basin_diagnostics,
            }
        )
        return [], rows, summary

    selected_ids = {id(row) for row in best_group}
    relaxed = bool(final_threshold > initial_threshold + EPS)
    for row in brightness_eligible:
        row["brightness_final_threshold_raw"] = float(final_threshold)
        row["brightness_threshold_raw"] = float(final_threshold)
        row["brightness_threshold_relaxed"] = relaxed
        row["brightness_relaxation_steps"] = int(relaxation_steps)
        row["brightness_true_trench"] = bool(
            float(row["brightness_core_mean_raw"]) <= final_threshold + EPS
        )
        row["brightness_selected"] = id(row) in selected_ids
        row["selected"] = row["brightness_selected"]

    selected_rows = sorted(
        best_group, key=lambda row: float(row["candidate_center_x_px"])
    )
    positions = np.asarray(
        [float(row["candidate_center_x_px"]) for row in selected_rows]
    )
    selected_widths = np.asarray(
        [float(row["estimated_width_px"]) for row in selected_rows]
    )
    selected_brightness = np.asarray(
        [float(row["brightness_core_mean_raw"]) for row in selected_rows]
    )
    selected_basin_indices = np.asarray(
        [int(row["basin_sequence_index"]) for row in selected_rows]
    )
    pitches = np.diff(positions)

    image_values = _finite(image)
    image_p05, image_p95 = np.percentile(image_values, [5, 95])
    image_dynamic = max(float(image_p95 - image_p05), 1e-9)

    summary.update(
        {
            "detection_success": True,
            "reason": "",
            "basin_first_enabled": True,
            "basin_first_fallback_used": False,
            "brightness_filter_enabled": True,
            "brightness_valid_candidate_count": len(brightness_eligible),
            "brightness_dark_cluster_center_raw": dark_center,
            "brightness_bright_cluster_center_raw": bright_center,
            "brightness_cluster_gap_raw": bright_center - dark_center,
            "brightness_initial_threshold_raw": float(initial_threshold),
            "brightness_final_threshold_raw": float(final_threshold),
            "brightness_initial_threshold_normalized": (
                initial_threshold - float(image_p05)
            )
            / image_dynamic,
            "brightness_final_threshold_normalized": (
                final_threshold - float(image_p05)
            )
            / image_dynamic,
            "brightness_threshold_relaxed": relaxed,
            "brightness_relaxation_steps": int(relaxation_steps),
            "brightness_initial_pool_count": int(initial_pool_count),
            "brightness_final_pool_count": sum(
                float(row["brightness_core_mean_raw"]) <= final_threshold + EPS
                for row in brightness_eligible
            ),
            "brightness_selected_count": len(selected_rows),
            "brightness_selected_mean_raw": float(np.mean(selected_brightness)),
            "brightness_selected_median_raw": float(np.median(selected_brightness)),
            "pattern_selected_basin_indices": ",".join(
                str(v) for v in selected_basin_indices
            ),
            "pattern_selected_basin_gaps": ",".join(
                str(v) for v in np.diff(selected_basin_indices)
            ),
            "pattern_group_cost": float(group_cost),
            "selected_count": len(selected_rows),
            "selected_target_count": len(selected_rows),
            "selected_centers_px": positions.tolist(),
            "selected_estimated_widths_nm": (
                selected_widths * float(fixed.pixel_size_nm)
            ).tolist(),
            "estimated_width_median_nm": float(
                np.median(selected_widths) * fixed.pixel_size_nm
            ),
            "pitch_mean_px": float(np.mean(pitches)) if pitches.size else np.nan,
            "pitch_cv": (
                float(np.std(pitches) / max(np.mean(pitches), EPS))
                if pitches.size
                else np.nan
            ),
            "pitch_ratio": (
                float(np.max(pitches) / max(np.min(pitches), EPS))
                if pitches.size
                else np.nan
            ),
            "edge_mode": "v115_basin_first_direct_regions",
            **profile_diagnostics,
            **separator_diagnostics,
            **basin_diagnostics,
        }
    )
    return selected_rows, rows, summary
