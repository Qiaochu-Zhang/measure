"""V1_15: consolidated active implementation, derived from V13_modified."""

from __future__ import annotations
import math
import numpy as np

from pathlib import Path
from typing import Any, Iterable, Optional
from dataclasses import asdict, dataclass
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks
import cdsem_locator as basin_v115


EPS = 1e-9


@dataclass
class MeasurementParams:
    engine_kind: str = "V13"
    threshold_search: str = "bounded"
    locator_mode: str = "dark-line"
    locator_majority: float = 0.70
    integer_sampling: bool = False
    # Image
    expected_width_px: int = 1024
    expected_height_px: int = 1024
    pixel_size_nm: float = 0.8789

    # Measurement
    target_cd_nm: float = 60.0
    center_x_px: float = 586.0
    center_y_px: float = 512.0
    search_range_in_px: int = 33
    search_range_out_px: int = 30
    angle_deg: float = 0.0
    extend_length_px: int = 360

    # Edge algorithm
    peak_order: int = 1
    edge_threshold_left_pct: float = 50.0
    edge_threshold_right_pct: float = 50.0
    smoothing_pixel: int = 3
    sample_number: int = 128
    average_range_px: int = 32
    peak_mad_multiplier: float = 1.0
    peak_min_contrast_8bit: float = 2.0

    # V2: edge topology and tracking.  These limits are deliberately much
    # tighter than the original Search In/Out window, so an edge point cannot
    # silently jump into the middle of a trench or onto a dark mark inside a line.
    topology_band_px: int = 7
    topology_gap_px: int = 1
    topology_min_contrast_8bit: float = 2.5
    tracking_search_radius_px: float = 13.0
    tracking_retry_radius_px: float = 20.0
    tracking_max_edge_jump_px: float = 8.0
    tracking_max_cd_change_fraction: float = 0.28
    tracking_max_center_shift_px: float = 7.0

    # Robust post-filter. Gross geometric errors are removed before the
    # documented Flyer=Reserve logic; ordinary roughness points are still kept.
    postfilter_window: int = 9
    postfilter_mad_multiplier: float = 6.0
    postfilter_min_edge_tolerance_px: float = 2.5
    postfilter_min_cd_tolerance_px: float = 4.0
    recover_rejected_points: bool = True

    pattern_kind: str = "trench"
    space_reference_nm: Optional[float] = None
    viterbi_enabled: int = 0
    viterbi_weight: float = 5.0
    viterbi_max_jump_px: float = 3.0
    viterbi_gap_cost: float = 6.0
    viterbi_candidates: int = 5
    erf_enabled: int = 0
    erf_window_px: int = 8
    erf_max_shift_px: float = 2.0
    erf_max_relative_rmse: float = 0.2
    edge_continuity: int = 0  # No extra relaxation or synthetic recovery by default.

    # Flyer
    flyer_mode: str = "reserve"  # reserve/remove/replace
    flyer_left_offset_px: float = 20.0
    flyer_right_offset_px: float = 20.0
    flyer_local_reference_window: int = 5

    # MultiPoint
    max_number: int = 3
    min_number: int = 3
    forbidden_zone: int = 1
    gray_tolerance: float = 10.0
    grad_tolerance: float = 1.0
    meas_area_width_px: int = 512
    meas_area_height_px: int = 512

    # V2: candidate generation uses an averaged basin profile. A narrow dark
    # stain inside a bright line is diluted by this averaging, whereas a true
    # trench remains dark across most of its CD.
    candidate_dark_window_factor: float = 0.55
    candidate_min_distance_factor: float = 0.55
    candidate_min_prominence_8bit: float = 0.35
    candidate_center_dedup_factor: float = 0.45
    candidate_min_basin_contrast_8bit: Optional[float] = None
    candidate_dark_mean_tolerance_8bit: Optional[float] = None
    candidate_min_width_factor: float = 0.55
    candidate_max_width_factor: float = 1.55
    require_center_left_right: bool = True

    # V1.15 Basin-First：只用于粗定位，不参与128点最终边缘算法
    locator_adaptive_min_width_nm: float = 30.0
    locator_adaptive_max_width_nm: float = 100.0
    locator_width_profile_samples: int = 41
    locator_width_profile_block_height_px: int = 5
    locator_width_profile_trim_fraction: float = 0.10
    locator_width_profile_sigma_px: float = 1.50
    locator_brightness_samples: int = 19
    locator_brightness_block_height_px: int = 5
    locator_brightness_core_fraction: float = 0.76
    locator_brightness_trim_fraction: float = 0.10
    locator_brightness_candidate_pool: int = 24
    locator_brightness_kmeans_restarts: int = 16
    locator_separator_kmeans_restarts: int = 16
    locator_separator_min_run_fraction: float = 0.05
    locator_separator_close_gap_fraction: float = 0.03
    locator_separator_min_run_px: int = 2
    locator_separator_close_gap_px: int = 1
    locator_pattern_preferred_basin_gap: int = 2
    locator_basin_width_min_factor: float = 0.45
    locator_basin_width_max_factor: float = 1.45
    locator_basin_reject_edge_basins: bool = True
    locator_basin_validation_samples: int = 17
    locator_basin_validation_block_height_px: int = 5
    locator_basin_validation_sigma_px: float = 1.0
    locator_basin_validation_min_iou: float = 0.25
    locator_basin_min_vertical_presence: float = 0.45
    locator_basin_min_separator_presence: float = 0.35
    locator_basin_max_center_jitter_px: float = 8.0
    locator_basin_max_center_jitter_fraction: float = 0.22
    locator_basin_max_width_cv: float = 0.40
    random_seed: int = 42

    # Statistics and validity
    roughness_multiplier: float = 3.0
    standard_deviation_ddof: int = 0
    minimum_valid_fraction: float = 0.70
    minimum_average_range_fraction: float = 0.75
    valid_cd_min_factor: float = 0.55
    valid_cd_max_factor: float = 1.55

    def validate(self) -> None:
        if self.locator_mode not in {"auto", "dark-line", "bright-line"}:
            raise ValueError("unsupported locator mode")
        if not 0.5 < self.locator_majority <= 1:
            raise ValueError("locator_majority must be in (0.5,1]")
        if self.engine_kind not in {"V10", "V13"} or self.threshold_search not in {
            "bounded",
            "legacy",
        }:
            raise ValueError("unsupported edge engine or threshold search")
        if not 0 <= self.edge_continuity <= 100:
            raise ValueError("edge_continuity must be in 0..100")
        if self.expected_width_px <= 0 or self.expected_height_px <= 0:
            raise ValueError("expected image size must be positive")
        if self.pixel_size_nm <= 0:
            raise ValueError("pixel_size_nm must be positive")
        if self.target_cd_nm <= 0:
            raise ValueError("target_cd_nm must be positive")
        if self.sample_number < 2:
            raise ValueError("sample_number must be >= 2")
        if self.average_range_px < 1:
            raise ValueError("average_range_px must be >= 1")
        if self.smoothing_pixel < 1:
            raise ValueError("smoothing_pixel must be a positive integer")
        if self.peak_order < 1:
            raise ValueError("peak_order must be >= 1")
        if self.max_number < 1 or self.min_number < 1:
            raise ValueError("max_number and min_number must be >= 1")
        if self.min_number > self.max_number:
            raise ValueError("min_number cannot exceed max_number")
        if self.require_center_left_right and self.max_number < 3:
            raise ValueError("require_center_left_right needs max_number >= 3")
        if not (0.0 < self.minimum_valid_fraction <= 1.0):
            raise ValueError("minimum_valid_fraction must be in (0, 1]")
        if self.standard_deviation_ddof < 0:
            raise ValueError("standard_deviation_ddof must be >= 0")
        if self.flyer_mode.lower() not in {"reserve", "remove", "replace"}:
            raise ValueError("flyer_mode must be reserve, remove, or replace")
        if self.topology_band_px < 2:
            raise ValueError("topology_band_px must be >= 2")
        if self.tracking_search_radius_px <= 0 or self.tracking_retry_radius_px <= 0:
            raise ValueError("tracking search radii must be positive")
        if self.postfilter_window < 3:
            raise ValueError("postfilter_window must be >= 3")
        if abs(self.angle_deg) > EPS:
            raise ValueError(
                "当前 V2 脚本只实现 Angle=0°（竖直 SPACE）。"
                "非零角度需要先围绕量测中心旋转图像。"
            )


@dataclass
class Candidate:
    center_x_px: float
    predicted_left_px: float
    predicted_right_px: float
    global_left_edge_px: float
    global_right_edge_px: float
    global_cd_px: float
    dark_level: float
    left_peak_strength: float
    right_peak_strength: float
    mean_grad_strength: float
    contrast: float
    score: float
    is_center_reference: bool
    source: str
    role: str = "candidate"  # center / left / right
    basin_mean: float = math.nan
    left_shoulder_mean: float = math.nan
    right_shoulder_mean: float = math.nan
    v115_basin_sequence_index: int = -1
    v115_basin_left_x_px: float = math.nan
    v115_basin_right_x_px: float = math.nan
    v115_estimated_width_px: float = math.nan
    v115_quality_score: float = math.nan
    v115_brightness_core_mean: float = math.nan
    v115_original_selected: bool = False


@dataclass
class EdgeDetection:
    valid: bool
    failure_reason: str = ""
    left_peak_x: float = math.nan
    right_peak_x: float = math.nan
    left_peak_strength: float = math.nan
    right_peak_strength: float = math.nan
    left_threshold: float = math.nan
    right_threshold: float = math.nan
    dark_level: float = math.nan
    left_peak_level: float = math.nan
    right_peak_level: float = math.nan
    left_edge_x: float = math.nan
    right_edge_x: float = math.nan
    local_cd_px: float = math.nan
    left_peak_fallback: bool = False
    right_peak_fallback: bool = False
    left_outside_level: float = math.nan
    left_inside_level: float = math.nan
    right_inside_level: float = math.nan
    right_outside_level: float = math.nan
    anchor_left_x: float = math.nan
    anchor_right_x: float = math.nan
    # V13 audit fields. They do not participate in candidate selection or
    # tracking; they record the intensity references used only for the final
    # 50% threshold crossing.
    threshold_reference_mode: str = ""
    threshold_dark_level: float = math.nan
    left_threshold_bright_level: float = math.nan
    right_threshold_bright_level: float = math.nan


@dataclass
class SpaceResult:
    image_name: str
    space_index: int
    space_role: str
    space_center_x: float
    valid_sample_count: int
    valid_fraction: float
    stable: bool
    flyer_left_count: int
    flyer_right_count: int
    cd_out_of_range_count: int
    peak_fallback_count: int
    hard_rejected_count: int
    recovered_count: int
    mean_cd_nm: float
    ler_left_nm: float
    ler_right_nm: float
    lwr_nm: float
    # Width direction diagnostics. Primary V2 result remains x-axis Δx.
    mean_cd_x_nm: float
    lwr_x_nm: float
    mean_cd_slant_nm: float
    lwr_slant_nm: float
    centerline_angle_from_vertical_deg: float
    normal_width_factor: float
    fit_x_eq_ky_plus_c_k: float
    fit_x_eq_ky_plus_c_c_px: float
    fit_y_eq_ax_plus_b_a: float
    fit_y_eq_ax_plus_b_b_px: float
    warning: str


@dataclass
class ImageMeasurement:
    image_row: dict[str, Any]
    space_rows: list[dict[str, Any]]
    sample_rows: list[dict[str, Any]]
    annotation_payload: dict[str, Any]


def robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0
    med = float(np.median(values))
    return float(np.median(np.abs(values - med)))


def normalize_to_8bit_float(image: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """将8/16-bit等输入线性归一化到0~255 float64。"""
    warnings: list[str] = []
    arr = np.asarray(image)
    if arr.ndim == 2 and arr.dtype == np.uint8:
        # Do not amplify a nearly uniform 8-bit background to full contrast.
        return arr.astype(np.float64), warnings

    # Multi-page TIFF: use first page. RGB/RGBA: convert to luminance.
    if arr.ndim == 3:
        if arr.shape[-1] in (3, 4):
            rgb = arr[..., :3].astype(np.float64)
            arr = 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]
            warnings.append("RGB/RGBA TIFF converted to grayscale luminance")
        else:
            arr = arr[0]
            warnings.append("multi-page TIFF detected; only the first page was used")
    elif arr.ndim > 3:
        arr = np.squeeze(arr)
        if arr.ndim > 2:
            arr = arr.reshape((-1,) + arr.shape[-2:])[0]
            warnings.append(
                "high-dimensional TIFF detected; only the first 2D plane was used"
            )

    if arr.ndim != 2:
        raise ValueError(f"expected a 2D grayscale image, got shape={arr.shape}")

    arr = arr.astype(np.float64, copy=False)
    finite = np.isfinite(arr)
    if not np.all(finite):
        if not np.any(finite):
            raise ValueError("image contains no finite pixels")
        fill = float(np.median(arr[finite]))
        arr = np.where(finite, arr, fill)
        warnings.append("non-finite pixels replaced with the finite-pixel median")

    lo = float(np.min(arr))
    hi = float(np.max(arr))
    if hi - lo < EPS:
        raise ValueError("image is constant; edge measurement is impossible")

    normalized = (arr - lo) * (255.0 / (hi - lo))
    return normalized, warnings


def read_image(path: Path) -> tuple[np.ndarray, list[str]]:
    from cdsem_regions import unicode_imread, cv2

    raw = unicode_imread(path)
    if raw.ndim == 3:
        raw = cv2.cvtColor(
            raw, cv2.COLOR_BGRA2GRAY if raw.shape[2] == 4 else cv2.COLOR_BGR2GRAY
        )
    return normalize_to_8bit_float(raw)


def moving_average_reflect(profile: np.ndarray, window: int) -> np.ndarray:
    profile = np.asarray(profile, dtype=np.float64)
    if window <= 1:
        return profile.copy()
    before = window // 2
    after = window - 1 - before
    padded = np.pad(profile, (before, after), mode="reflect")
    kernel = np.ones(window, dtype=np.float64) / float(window)
    return np.convolve(padded, kernel, mode="valid")


def central_difference(profile: np.ndarray) -> np.ndarray:
    profile = np.asarray(profile, dtype=np.float64)
    grad = np.empty_like(profile)
    if profile.size < 3:
        grad[:] = np.nan
        return grad
    grad[1:-1] = (profile[2:] - profile[:-2]) / 2.0
    grad[0] = profile[1] - profile[0]
    grad[-1] = profile[-1] - profile[-2]
    return grad


def clamp_int_window(start: float, end: float, size: int) -> tuple[int, int]:
    s = max(0, min(size, int(math.floor(start))))
    e = max(0, min(size, int(math.ceil(end))))
    if e < s:
        s, e = e, s
    return s, e


def safe_median_segment(profile: np.ndarray, start: float, end: float) -> float:
    s, e = clamp_int_window(start, end, profile.size)
    if e <= s:
        return math.nan
    return float(np.median(profile[s:e]))


def safe_mean_segment(profile: np.ndarray, start: float, end: float) -> float:
    s, e = clamp_int_window(start, end, profile.size)
    if e <= s:
        return math.nan
    return float(np.mean(profile[s:e]))


def local_peak_level(profile: np.ndarray, peak_x: int) -> float:
    s = max(0, peak_x - 1)
    e = min(profile.size, peak_x + 2)
    if e <= s:
        return math.nan
    return float(np.median(profile[s:e]))


def v17_threshold_reference_levels(
    profile: np.ndarray,
    predicted_left: float,
    predicted_right: float,
    left_peak: int,
    right_peak: int,
    expected_width_px: float,
) -> tuple[float, float, float]:
    """Return V1.7 dark/bright references for a 50% edge threshold.

    This function intentionally changes only the gray-level references used
    after V10 has selected its left/right gradient peaks. The V10 peak search,
    Peak Order=1, topology checks, tracked anchors and retry logic remain
    untouched.

    The definition is copied from ``sem_cd_measure_200k_batch_V1_7.py``:

    * dark reference: mean of the central 44% of the tracked dark band
      (28%--72% between the predicted edges);
    * left bright reference: mean of the outside band ending one pixel before
      the selected left peak;
    * right bright reference: mean of the outside band starting one pixel
      after the selected right peak;
    * outside-band width: max(3 px, round(0.18 * expected CD in pixels)).
    """
    width = float(predicted_right - predicted_left)
    if not np.isfinite(width) or width <= 0:
        return math.nan, math.nan, math.nan

    dark_start = predicted_left + 0.28 * width
    dark_end = predicted_right - 0.28 * width
    dark_mean = safe_mean_segment(profile, dark_start, dark_end)

    band = max(3, int(round(0.18 * max(float(expected_width_px), 4.0))))
    left_bright_mean = safe_mean_segment(
        profile,
        float(left_peak - band),
        float(left_peak - 1),
    )
    right_bright_mean = safe_mean_segment(
        profile,
        float(right_peak + 1),
        float(right_peak + band),
    )
    return dark_mean, left_bright_mean, right_bright_mean


def edge_topology_levels(
    profile: np.ndarray,
    edge_x: float,
    side: str,
    band_px: int,
    gap_px: int,
) -> tuple[float, float]:
    """Return (outside_line_level, inside_trench_level) around an edge.

    The bands intentionally avoid the exact gradient maximum.  For a valid SPACE
    boundary, the line side must be brighter on average than the trench side.
    This rejects small dark/bright marks inside a trench that happen to create a
    local gradient peak.
    """
    x = int(round(edge_x))
    band = max(2, int(band_px))
    gap = max(0, int(gap_px))
    if side == "left":
        outside = safe_mean_segment(profile, x - gap - band, x - gap)
        inside = safe_mean_segment(profile, x + gap, x + gap + band)
    elif side == "right":
        inside = safe_mean_segment(profile, x - gap - band, x - gap)
        outside = safe_mean_segment(profile, x + gap, x + gap + band)
    else:
        raise ValueError(f"unknown side: {side}")
    return outside, inside


def adaptive_peak_indices(
    profile: np.ndarray,
    gradient: np.ndarray,
    start: int,
    end: int,
    sign: str,
    peak_order: int,
    center_x: float,
    mad_multiplier: float,
    dark_level: float,
    min_contrast: float,
    anchor_x: float,
    max_anchor_distance: float,
    topology_band_px: int,
    topology_gap_px: int,
    topology_min_contrast: float,
) -> tuple[Optional[int], float, bool]:
    """Find a direction-correct peak near the expected edge anchor.

    V1 selected the first gradient peak encountered from the SPACE centre.  That
    can lock onto a speckle or a false dark patch inside the structure.  V2 keeps
    the gradient rule but additionally requires:
      1. proximity to the global/tracked edge anchor;
      2. a bright-line -> dark-trench average-intensity topology;
      3. a peak level clearly above the trench dark reference.
    """
    start = max(0, start)
    end = min(gradient.size, end)
    if end - start < 3:
        return None, math.nan, False

    local_g = gradient[start:end]
    transformed = -local_g if sign == "negative" else local_g
    abs_local = np.abs(local_g)
    baseline = float(np.median(abs_local))
    mad = robust_mad(abs_local)
    threshold = max(baseline + mad_multiplier * mad, 1e-6)

    peaks, _ = find_peaks(transformed, height=threshold)
    global_peaks = [int(start + p) for p in peaks if transformed[p] > 0]

    def peak_is_valid(p: int) -> bool:
        if sign == "negative":
            if p >= center_x:
                return False
            side = "left"
        else:
            if p <= center_x:
                return False
            side = "right"
        if abs(float(p) - anchor_x) > max_anchor_distance:
            return False
        if local_peak_level(profile, p) < dark_level + min_contrast:
            return False
        outside, inside = edge_topology_levels(
            profile, p, side, topology_band_px, topology_gap_px
        )
        if not np.isfinite(outside) or not np.isfinite(inside):
            return False
        return (outside - inside) >= topology_min_contrast

    valid_peaks = [p for p in global_peaks if peak_is_valid(p)]
    # The expected edge location is a stronger discriminator than "first from
    # centre" when a false internal dark spot exists.  Peak order is still
    # honoured after sorting by anchor distance.
    valid_peaks.sort(key=lambda p: (abs(float(p) - anchor_x), abs(float(p) - center_x)))
    if len(valid_peaks) >= peak_order:
        selected = valid_peaks[peak_order - 1]
        return selected, float(abs(gradient[selected])), False

    # An explicitly requested higher rank must not silently become rank 1.
    if peak_order > 1:
        return None, math.nan, False

    # Endpoint / plateau fallback, subject to the same topology checks.
    eligible: list[int] = []
    for idx in range(start, end):
        correct_sign = (
            sign == "negative" and idx < center_x and gradient[idx] < 0
        ) or (sign == "positive" and idx > center_x and gradient[idx] > 0)
        if correct_sign and peak_is_valid(idx):
            eligible.append(idx)
    if not eligible:
        return None, math.nan, False
    selected = max(
        eligible,
        key=lambda idx: (
            abs(float(gradient[idx])),
            -abs(float(idx) - anchor_x),
        ),
    )
    strength = float(abs(gradient[selected]))
    if strength >= threshold:
        return int(selected), strength, True
    return None, math.nan, False


def interpolate_threshold_crossing(
    profile: np.ndarray,
    threshold: float,
    start_peak_x: int,
    center_x: float,
    side: str,
) -> Optional[float]:
    """Search the first 50% crossing from the bright line toward the SPACE."""
    if not np.isfinite(threshold):
        return None

    if side == "left":
        stop = min(profile.size - 1, int(math.ceil(center_x)))
        for x in range(max(0, start_peak_x), stop):
            p0 = float(profile[x])
            p1 = float(profile[x + 1])
            if p0 >= threshold and p1 < threshold:
                denom = p1 - p0
                if abs(denom) < EPS:
                    return None
                return float(x + (threshold - p0) / denom)
    elif side == "right":
        stop = max(0, int(math.floor(center_x)))
        for x in range(min(profile.size - 1, start_peak_x - 1), stop - 1, -1):
            p0 = float(profile[x])
            p1 = float(profile[x + 1])
            if p0 < threshold and p1 >= threshold:
                denom = p1 - p0
                if abs(denom) < EPS:
                    return None
                return float(x + (threshold - p0) / denom)
    else:
        raise ValueError(f"unknown side: {side}")
    return None


def bounded_threshold_crossing(profile, threshold, peak, center, side, lo, hi, anchor):
    """Find a polarity-correct crossing on either side of the selected peak.

    Stay inside the existing anchor search window and on the proper side of
    the feature. Rank by selected-peak proximity; all pair topology/width/jump
    checks in detect_edge_pair still apply. No edge positions are fabricated.
    """
    if not np.isfinite(threshold):
        return None
    lo = max(0, int(math.floor(lo)))
    hi = min(len(profile) - 1, int(math.ceil(hi)))
    picks = []
    for x in range(lo, hi):
        p0, p1 = float(profile[x]), float(profile[x + 1])
        correct = (p0 >= threshold > p1) if side == "left" else (p0 < threshold <= p1)
        if not correct:
            continue
        value = x + (threshold - p0) / (p1 - p0)
        if (side == "left" and value < center) or (side == "right" and value > center):
            picks.append(value)
    return min(picks, key=lambda x: (abs(x - peak), abs(x - anchor))) if picks else None


def detect_edge_pair(
    smoothed_profile: np.ndarray,
    candidate_center_x: float,
    params: MeasurementParams,
    *,
    anchor_left_x: Optional[float] = None,
    anchor_right_x: Optional[float] = None,
    search_radius_px: Optional[float] = None,
    strict_tracking: bool = False,
) -> EdgeDetection:
    """Detect one SPACE edge pair, optionally anchored to tracked positions."""
    n = smoothed_profile.size
    target_cd_px = params.target_cd_nm / params.pixel_size_nm
    predicted_left = (
        float(anchor_left_x)
        if anchor_left_x is not None and np.isfinite(anchor_left_x)
        else candidate_center_x - target_cd_px / 2.0
    )
    predicted_right = (
        float(anchor_right_x)
        if anchor_right_x is not None and np.isfinite(anchor_right_x)
        else candidate_center_x + target_cd_px / 2.0
    )
    if predicted_right <= predicted_left:
        return EdgeDetection(valid=False, failure_reason="invalid_edge_anchors")

    if search_radius_px is None:
        left_start, left_end = clamp_int_window(
            predicted_left - params.search_range_out_px,
            predicted_left + params.search_range_in_px,
            n,
        )
        right_start, right_end = clamp_int_window(
            predicted_right - params.search_range_in_px,
            predicted_right + params.search_range_out_px,
            n,
        )
        left_anchor_limit = max(params.search_range_in_px, params.search_range_out_px)
        right_anchor_limit = left_anchor_limit
    else:
        radius = float(search_radius_px)
        left_start, left_end = clamp_int_window(
            predicted_left - radius, predicted_left + radius, n
        )
        right_start, right_end = clamp_int_window(
            predicted_right - radius, predicted_right + radius, n
        )
        left_anchor_limit = radius
        right_anchor_limit = radius

    width = predicted_right - predicted_left
    dark_start = predicted_left + 0.30 * width
    dark_end = predicted_right - 0.30 * width
    dark_level = safe_median_segment(smoothed_profile, dark_start, dark_end)
    if not np.isfinite(dark_level):
        return EdgeDetection(valid=False, failure_reason="invalid_dark_level")

    gradient = central_difference(smoothed_profile)
    left_peak, left_strength, left_fallback = adaptive_peak_indices(
        smoothed_profile,
        gradient,
        left_start,
        left_end,
        sign="negative",
        peak_order=params.peak_order,
        center_x=candidate_center_x,
        mad_multiplier=params.peak_mad_multiplier,
        dark_level=dark_level,
        min_contrast=params.peak_min_contrast_8bit,
        anchor_x=predicted_left,
        max_anchor_distance=left_anchor_limit,
        topology_band_px=params.topology_band_px,
        topology_gap_px=params.topology_gap_px,
        topology_min_contrast=params.topology_min_contrast_8bit,
    )
    if left_peak is None:
        return EdgeDetection(
            valid=False,
            failure_reason="no_valid_left_peak_near_anchor",
            anchor_left_x=predicted_left,
            anchor_right_x=predicted_right,
        )

    right_peak, right_strength, right_fallback = adaptive_peak_indices(
        smoothed_profile,
        gradient,
        right_start,
        right_end,
        sign="positive",
        peak_order=params.peak_order,
        center_x=candidate_center_x,
        mad_multiplier=params.peak_mad_multiplier,
        dark_level=dark_level,
        min_contrast=params.peak_min_contrast_8bit,
        anchor_x=predicted_right,
        max_anchor_distance=right_anchor_limit,
        topology_band_px=params.topology_band_px,
        topology_gap_px=params.topology_gap_px,
        topology_min_contrast=params.topology_min_contrast_8bit,
    )
    if right_peak is None:
        return EdgeDetection(
            valid=False,
            failure_reason="no_valid_right_peak_near_anchor",
            left_peak_x=float(left_peak),
            left_peak_strength=left_strength,
            left_peak_fallback=left_fallback,
            anchor_left_x=predicted_left,
            anchor_right_x=predicted_right,
        )

    # Preserve the V10 peak validation exactly. These local peak levels are
    # deliberately not used as V13's threshold bright references.
    v10_left_peak_level = local_peak_level(smoothed_profile, left_peak)
    v10_right_peak_level = local_peak_level(smoothed_profile, right_peak)
    if not all(
        np.isfinite(v) for v in (dark_level, v10_left_peak_level, v10_right_peak_level)
    ):
        return EdgeDetection(valid=False, failure_reason="invalid_dark_or_peak_level")
    if v10_left_peak_level <= dark_level or v10_right_peak_level <= dark_level:
        return EdgeDetection(
            valid=False, failure_reason="peak_not_brighter_than_trench"
        )

    # V13 is identical to V10 except for this threshold-reference definition.
    # The selected peaks are still the V10 peaks. Only the dark and bright
    # gray-level references used to place the 50% crossings come from V1.7.
    threshold_dark, left_bright, right_bright = v17_threshold_reference_levels(
        smoothed_profile,
        predicted_left,
        predicted_right,
        left_peak,
        right_peak,
        target_cd_px,
    )
    if params.engine_kind == "V10":
        threshold_dark, left_bright, right_bright = (
            dark_level,
            v10_left_peak_level,
            v10_right_peak_level,
        )
    threshold_mode = (
        "v17_side_band_mean50"
        if params.engine_kind == "V13"
        else "v10_center40_peak3_median"
    )
    if not all(np.isfinite(v) for v in (threshold_dark, left_bright, right_bright)):
        return EdgeDetection(
            valid=False,
            failure_reason="invalid_v17_threshold_reference_levels",
            left_peak_x=float(left_peak),
            right_peak_x=float(right_peak),
            left_peak_strength=left_strength,
            right_peak_strength=right_strength,
            dark_level=dark_level,
            left_peak_level=v10_left_peak_level,
            right_peak_level=v10_right_peak_level,
            threshold_reference_mode=threshold_mode,
            threshold_dark_level=threshold_dark,
            left_threshold_bright_level=left_bright,
            right_threshold_bright_level=right_bright,
            anchor_left_x=predicted_left,
            anchor_right_x=predicted_right,
        )
    if left_bright <= threshold_dark or right_bright <= threshold_dark:
        return EdgeDetection(
            valid=False,
            failure_reason="v17_bright_reference_not_above_dark_reference",
            left_peak_x=float(left_peak),
            right_peak_x=float(right_peak),
            left_peak_strength=left_strength,
            right_peak_strength=right_strength,
            dark_level=dark_level,
            left_peak_level=v10_left_peak_level,
            right_peak_level=v10_right_peak_level,
            threshold_reference_mode=threshold_mode,
            threshold_dark_level=threshold_dark,
            left_threshold_bright_level=left_bright,
            right_threshold_bright_level=right_bright,
            anchor_left_x=predicted_left,
            anchor_right_x=predicted_right,
        )

    left_threshold = threshold_dark + (params.edge_threshold_left_pct / 100.0) * (
        left_bright - threshold_dark
    )
    right_threshold = threshold_dark + (params.edge_threshold_right_pct / 100.0) * (
        right_bright - threshold_dark
    )

    left_edge = interpolate_threshold_crossing(
        smoothed_profile, left_threshold, left_peak, candidate_center_x, "left"
    )
    if params.threshold_search == "bounded":
        left_edge = bounded_threshold_crossing(
            smoothed_profile,
            left_threshold,
            left_peak,
            candidate_center_x,
            "left",
            left_start,
            left_end - 1,
            predicted_left,
        )
    if left_edge is None:
        return EdgeDetection(
            valid=False,
            failure_reason="left_threshold_crossing_not_found",
            left_peak_x=float(left_peak),
            right_peak_x=float(right_peak),
            left_peak_strength=left_strength,
            right_peak_strength=right_strength,
            left_threshold=left_threshold,
            right_threshold=right_threshold,
            dark_level=dark_level,
            left_peak_level=v10_left_peak_level,
            right_peak_level=v10_right_peak_level,
            threshold_reference_mode=threshold_mode,
            threshold_dark_level=threshold_dark,
            left_threshold_bright_level=left_bright,
            right_threshold_bright_level=right_bright,
            left_peak_fallback=left_fallback,
            right_peak_fallback=right_fallback,
            anchor_left_x=predicted_left,
            anchor_right_x=predicted_right,
        )

    right_edge = interpolate_threshold_crossing(
        smoothed_profile, right_threshold, right_peak, candidate_center_x, "right"
    )
    if params.threshold_search == "bounded":
        right_edge = bounded_threshold_crossing(
            smoothed_profile,
            right_threshold,
            right_peak,
            candidate_center_x,
            "right",
            right_start,
            right_end - 1,
            predicted_right,
        )
    if right_edge is None:
        return EdgeDetection(
            valid=False,
            failure_reason="right_threshold_crossing_not_found",
            left_peak_x=float(left_peak),
            right_peak_x=float(right_peak),
            left_peak_strength=left_strength,
            right_peak_strength=right_strength,
            left_threshold=left_threshold,
            right_threshold=right_threshold,
            dark_level=dark_level,
            left_peak_level=v10_left_peak_level,
            right_peak_level=v10_right_peak_level,
            threshold_reference_mode=threshold_mode,
            threshold_dark_level=threshold_dark,
            left_threshold_bright_level=left_bright,
            right_threshold_bright_level=right_bright,
            left_edge_x=float(left_edge),
            left_peak_fallback=left_fallback,
            right_peak_fallback=right_fallback,
            anchor_left_x=predicted_left,
            anchor_right_x=predicted_right,
        )

    local_cd_px = float(right_edge - left_edge)
    if not np.isfinite(local_cd_px) or local_cd_px <= 0:
        return EdgeDetection(valid=False, failure_reason="nonpositive_local_cd")

    left_outside, left_inside = edge_topology_levels(
        smoothed_profile,
        left_edge,
        "left",
        params.topology_band_px,
        params.topology_gap_px,
    )
    right_outside, right_inside = edge_topology_levels(
        smoothed_profile,
        right_edge,
        "right",
        params.topology_band_px,
        params.topology_gap_px,
    )
    topology_values = (left_outside, left_inside, right_inside, right_outside)
    if not all(np.isfinite(v) for v in topology_values):
        return EdgeDetection(valid=False, failure_reason="invalid_edge_topology_levels")
    if left_outside - left_inside < params.topology_min_contrast_8bit:
        return EdgeDetection(
            valid=False, failure_reason="left_edge_wrong_brightness_topology"
        )
    if right_outside - right_inside < params.topology_min_contrast_8bit:
        return EdgeDetection(
            valid=False, failure_reason="right_edge_wrong_brightness_topology"
        )

    if not (left_edge < candidate_center_x < right_edge):
        return EdgeDetection(
            valid=False, failure_reason="candidate_center_not_inside_trench"
        )

    min_cd = params.valid_cd_min_factor * target_cd_px
    max_cd = params.valid_cd_max_factor * target_cd_px
    if not (min_cd <= local_cd_px <= max_cd):
        return EdgeDetection(valid=False, failure_reason="local_cd_outside_hard_range")

    if strict_tracking:
        anchor_cd = predicted_right - predicted_left
        max_cd_delta = params.tracking_max_cd_change_fraction * max(
            anchor_cd, target_cd_px
        )
        if abs(local_cd_px - anchor_cd) > max_cd_delta:
            return EdgeDetection(valid=False, failure_reason="tracking_cd_jump")
        if abs(left_edge - predicted_left) > params.tracking_max_edge_jump_px:
            return EdgeDetection(valid=False, failure_reason="tracking_left_edge_jump")
        if abs(right_edge - predicted_right) > params.tracking_max_edge_jump_px:
            return EdgeDetection(valid=False, failure_reason="tracking_right_edge_jump")
        midpoint = 0.5 * (left_edge + right_edge)
        anchor_midpoint = 0.5 * (predicted_left + predicted_right)
        if abs(midpoint - anchor_midpoint) > params.tracking_max_center_shift_px:
            return EdgeDetection(valid=False, failure_reason="tracking_center_jump")

    return EdgeDetection(
        valid=True,
        left_peak_x=float(left_peak),
        right_peak_x=float(right_peak),
        left_peak_strength=left_strength,
        right_peak_strength=right_strength,
        left_threshold=left_threshold,
        right_threshold=right_threshold,
        dark_level=dark_level,
        left_peak_level=v10_left_peak_level,
        right_peak_level=v10_right_peak_level,
        threshold_reference_mode=threshold_mode,
        threshold_dark_level=threshold_dark,
        left_threshold_bright_level=left_bright,
        right_threshold_bright_level=right_bright,
        left_edge_x=float(left_edge),
        right_edge_x=float(right_edge),
        local_cd_px=local_cd_px,
        left_peak_fallback=left_fallback,
        right_peak_fallback=right_fallback,
        left_outside_level=left_outside,
        left_inside_level=left_inside,
        right_inside_level=right_inside,
        right_outside_level=right_outside,
        anchor_left_x=predicted_left,
        anchor_right_x=predicted_right,
    )


def measurement_area_bounds(
    image_shape: tuple[int, int], params: MeasurementParams
) -> tuple[int, int, int, int, bool]:
    height, width = image_shape
    raw_x0 = int(round(params.center_x_px - params.meas_area_width_px / 2.0))
    raw_x1 = raw_x0 + params.meas_area_width_px
    raw_y0 = int(round(params.center_y_px - params.meas_area_height_px / 2.0))
    raw_y1 = raw_y0 + params.meas_area_height_px

    x0 = max(0, raw_x0)
    x1 = min(width, raw_x1)
    y0 = max(0, raw_y0)
    y1 = min(height, raw_y1)
    truncated = (x0, x1, y0, y1) != (raw_x0, raw_x1, raw_y0, raw_y1)
    return x0, x1, y0, y1, truncated


def global_profile(
    image: np.ndarray, params: MeasurementParams
) -> tuple[np.ndarray, int, int]:
    y0 = int(round(params.center_y_px - params.extend_length_px / 2.0))
    y1 = y0 + params.extend_length_px
    y0 = max(0, y0)
    y1 = min(image.shape[0], y1)
    if y1 <= y0:
        raise ValueError("Extend Length does not overlap the image")
    profile = image[y0:y1, :].mean(axis=0)
    return profile, y0, y1


def candidate_basin_metrics(
    profile: np.ndarray,
    left_edge: float,
    right_edge: float,
) -> tuple[float, float, float, float]:
    """Average brightness metrics for discriminating true trenches from line stains."""
    width = max(float(right_edge - left_edge), 1.0)
    basin_mean = safe_mean_segment(
        profile, left_edge + 0.15 * width, right_edge - 0.15 * width
    )
    shoulder_width = max(5.0, 0.22 * width)
    left_shoulder = safe_mean_segment(
        profile, left_edge - shoulder_width, left_edge - 1.0
    )
    right_shoulder = safe_mean_segment(
        profile, right_edge + 1.0, right_edge + shoulder_width
    )
    contrast = min(left_shoulder, right_shoulder) - basin_mean
    return basin_mean, left_shoulder, right_shoulder, contrast


def candidate_score(
    center_x: float,
    edge: EdgeDetection,
    profile: np.ndarray,
    x0: int,
    x1: int,
    params: MeasurementParams,
) -> tuple[float, float, float, float, float]:
    target_px = params.target_cd_nm / params.pixel_size_nm
    meas_half = max(params.meas_area_width_px / 2.0, 1.0)
    center_score = math.exp(-abs(center_x - params.center_x_px) / meas_half)
    width_score = math.exp(-abs(edge.local_cd_px - target_px) / max(target_px, EPS))
    basin_mean, left_shoulder, right_shoulder, contrast = candidate_basin_metrics(
        profile, edge.left_edge_x, edge.right_edge_x
    )
    roi_dynamic = float(
        np.percentile(profile[x0:x1], 95) - np.percentile(profile[x0:x1], 5)
    )
    gray_score = float(np.clip(contrast / max(roi_dynamic, EPS), 0.0, 1.0))
    roi_grad = np.abs(central_difference(profile[x0:x1]))
    grad_reference = float(np.percentile(roi_grad, 95)) if roi_grad.size else 1.0
    mean_strength = (edge.left_peak_strength + edge.right_peak_strength) / 2.0
    grad_score = float(np.clip(mean_strength / max(grad_reference, EPS), 0.0, 1.0))
    # Darkness contributes explicitly in V2.  A fake dark patch in the middle of
    # a line can have a low single-pixel value, but its full-basin mean is higher.
    roi_low = float(np.percentile(profile[x0:x1], 5))
    roi_high = float(np.percentile(profile[x0:x1], 95))
    darkness_score = 1.0 - float(
        np.clip((basin_mean - roi_low) / max(roi_high - roi_low, EPS), 0.0, 1.0)
    )
    score = (
        0.20 * center_score
        + 0.22 * width_score
        + 0.23 * gray_score
        + 0.15 * grad_score
        + 0.20 * darkness_score
    )
    return (
        float(score),
        float(contrast),
        float(basin_mean),
        float(left_shoulder),
        float(right_shoulder),
    )


@dataclass
class V115FixedParams:
    """Minimal interface required by the validated V1.15 Basin-First module."""

    pixel_size_nm: float
    center_x_px: float
    center_y_px: float
    extend_length_nm: float
    preferred_trenches: int
    min_candidate_trenches: int
    adaptive_min_width_nm: float
    adaptive_max_width_nm: float
    random_seed: int
    brightness_samples: int
    brightness_block_height_px: int
    brightness_core_fraction: float
    brightness_trim_fraction: float
    brightness_candidate_pool: int
    brightness_kmeans_restarts: int
    width_profile_samples: int
    width_profile_block_height_px: int
    width_profile_trim_fraction: float
    width_profile_sigma_px: float
    separator_kmeans_restarts: int
    separator_min_run_fraction: float
    separator_close_gap_fraction: float
    separator_min_run_px: int
    separator_close_gap_px: int
    pattern_preferred_basin_gap: int
    basin_first_enabled: bool
    basin_width_min_factor: float
    basin_width_max_factor: float
    basin_reject_edge_basins: bool
    basin_validation_samples: int
    basin_validation_block_height_px: int
    basin_validation_sigma_px: float
    basin_validation_min_iou: float
    basin_min_vertical_presence: float
    basin_min_separator_presence: float
    basin_max_center_jitter_px: float
    basin_max_center_jitter_fraction: float
    basin_max_width_cv: float


def v115_block_profile(image: np.ndarray, y: float, block_height: int) -> np.ndarray:
    y0 = max(0, int(round(y - block_height / 2.0)))
    y1 = min(image.shape[0], y0 + max(1, int(block_height)))
    if y1 <= y0:
        return np.full(image.shape[1], np.nan, dtype=np.float64)
    return np.mean(image[y0:y1], axis=0)


def build_v115_fixed(
    params: MeasurementParams, local_center_x: float, local_center_y: float
) -> V115FixedParams:
    return V115FixedParams(
        pixel_size_nm=params.pixel_size_nm,
        center_x_px=local_center_x,
        center_y_px=local_center_y,
        extend_length_nm=params.extend_length_px * params.pixel_size_nm,
        preferred_trenches=min(3, params.max_number),
        min_candidate_trenches=min(3, params.min_number),
        adaptive_min_width_nm=params.locator_adaptive_min_width_nm,
        adaptive_max_width_nm=params.locator_adaptive_max_width_nm,
        random_seed=params.random_seed,
        brightness_samples=params.locator_brightness_samples,
        brightness_block_height_px=params.locator_brightness_block_height_px,
        brightness_core_fraction=params.locator_brightness_core_fraction,
        brightness_trim_fraction=params.locator_brightness_trim_fraction,
        brightness_candidate_pool=params.locator_brightness_candidate_pool,
        brightness_kmeans_restarts=params.locator_brightness_kmeans_restarts,
        width_profile_samples=params.locator_width_profile_samples,
        width_profile_block_height_px=params.locator_width_profile_block_height_px,
        width_profile_trim_fraction=params.locator_width_profile_trim_fraction,
        width_profile_sigma_px=params.locator_width_profile_sigma_px,
        separator_kmeans_restarts=params.locator_separator_kmeans_restarts,
        separator_min_run_fraction=params.locator_separator_min_run_fraction,
        separator_close_gap_fraction=params.locator_separator_close_gap_fraction,
        separator_min_run_px=params.locator_separator_min_run_px,
        separator_close_gap_px=params.locator_separator_close_gap_px,
        pattern_preferred_basin_gap=params.locator_pattern_preferred_basin_gap,
        basin_first_enabled=True,
        basin_width_min_factor=params.locator_basin_width_min_factor,
        basin_width_max_factor=params.locator_basin_width_max_factor,
        basin_reject_edge_basins=params.locator_basin_reject_edge_basins,
        basin_validation_samples=params.locator_basin_validation_samples,
        basin_validation_block_height_px=params.locator_basin_validation_block_height_px,
        basin_validation_sigma_px=params.locator_basin_validation_sigma_px,
        basin_validation_min_iou=params.locator_basin_validation_min_iou,
        basin_min_vertical_presence=params.locator_basin_min_vertical_presence,
        basin_min_separator_presence=params.locator_basin_min_separator_presence,
        basin_max_center_jitter_px=params.locator_basin_max_center_jitter_px,
        basin_max_center_jitter_fraction=params.locator_basin_max_center_jitter_fraction,
        basin_max_width_cv=params.locator_basin_max_width_cv,
    )


def _row_float(row: dict[str, Any], key: str, default: float = math.nan) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return float(default)


def _row_bool(row: dict[str, Any], key: str, default: bool = False) -> bool:
    value = row.get(key, default)
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return bool(value)


def _row_center(row: dict[str, Any]) -> float:
    return _row_float(
        row, "candidate_center_x_px", _row_float(row, "basin_center_x_px")
    )


def _basin_contains(row: dict[str, Any], x: float) -> bool:
    left = _row_float(row, "basin_left_x_px")
    right = _row_float(row, "basin_right_x_px")
    return np.isfinite(left) and np.isfinite(right) and left <= x <= right


def choose_center_left_right_v115(
    selected_rows: list[dict[str, Any]],
    all_rows: list[dict[str, Any]],
    center_x_local: float,
    params: MeasurementParams,
) -> list[tuple[str, dict[str, Any]]]:
    """Use V1.15 basin identity, but force the user's known center trench and nearest true left/right trench."""
    finite_rows = [r for r in all_rows if np.isfinite(_row_center(r))]
    if not finite_rows:
        return []
    valid_rows = [r for r in finite_rows if _row_bool(r, "basin_valid")]
    pool = valid_rows if valid_rows else finite_rows

    containing = [r for r in pool if _basin_contains(r, center_x_local)]
    if containing:
        center_row = min(containing, key=lambda r: abs(_row_center(r) - center_x_local))
    else:
        center_row = min(pool, key=lambda r: abs(_row_center(r) - center_x_local))

    center_pos = _row_center(center_row)
    center_index = int(
        center_row.get("basin_sequence_index", center_row.get("candidate_id", 0))
    )
    selected_keys = {
        int(r.get("basin_sequence_index", r.get("candidate_id", -99999)))
        for r in selected_rows
    }

    # Prefer rows V1.15 classified as true trench; if too few, fall back to valid basins.
    true_pool = [r for r in pool if _row_bool(r, "brightness_true_trench")]
    side_pool = true_pool if len(true_pool) >= 3 else pool
    preferred_gap = max(1, int(params.locator_pattern_preferred_basin_gap))

    def side_key(r: dict[str, Any]) -> tuple[float, ...]:
        idx = int(r.get("basin_sequence_index", r.get("candidate_id", 0)))
        gap_penalty = abs(abs(idx - center_index) - preferred_gap)
        selected_penalty = 0.0 if idx in selected_keys else 1.0
        true_penalty = 0.0 if _row_bool(r, "brightness_true_trench") else 1.0
        quality = _row_float(r, "basin_quality_score", 0.0)
        return (
            gap_penalty,
            selected_penalty,
            true_penalty,
            abs(_row_center(r) - center_pos),
            -quality,
        )

    left_pool = [r for r in side_pool if _row_center(r) < center_pos]
    right_pool = [r for r in side_pool if _row_center(r) > center_pos]
    chosen: list[tuple[str, dict[str, Any]]] = [("center", center_row)]
    if left_pool:
        chosen.append(("left", min(left_pool, key=side_key)))
    if right_pool:
        chosen.append(("right", min(right_pool, key=side_key)))
    order = {"center": 0, "left": 1, "right": 2}
    chosen.sort(key=lambda item: order[item[0]])
    return chosen[: params.max_number]


def detect_space_candidates(
    image: np.ndarray, params: MeasurementParams, region_mask=None
) -> tuple[list[Candidate], dict[str, Any]]:
    """V1.15 coarse basin locator -> V2 global edge refinement -> V2 128-point measurement downstream."""
    x0, x1, y0, y1, roi_truncated = measurement_area_bounds(image.shape, params)
    if (
        roi_truncated
        or (x1 - x0) != params.meas_area_width_px
        or (y1 - y0) != params.meas_area_height_px
    ):
        raise ValueError(
            f"定位ROI必须完整位于图像内；请求宽高={params.meas_area_width_px}x{params.meas_area_height_px}px，实际={x1 - x0}x{y1 - y0}px。调整center-x/y或meas-area-width/height"
        )

    roi = image[y0:y1, x0:x1]
    from cdsem_localization import (
        analyze_mode,
        masked_block_profile,
        select_bright_line_trenches,
        mode_params,
        supported_profile,
        observation_mask,
    )

    profile_mask = (
        observation_mask(region_mask, params.target_cd_nm / params.pixel_size_nm)
        if region_mask is not None
        else None
    )
    mask = profile_mask[y0:y1, x0:x1] if profile_mask is not None else None
    mode_diag = analyze_mode(roi, mask, params.locator_mode, params.locator_majority)
    if params.locator_mode == "auto":
        from dataclasses import replace

        params = replace(params, locator_mode=mode_diag["locator_mode_selected"])
    params = mode_params(params)
    block_profile_fn = masked_block_profile(roi, mask, v115_block_profile)
    local_center_x = params.center_x_px - x0
    local_center_y = params.center_y_px - y0
    fixed = build_v115_fixed(params, local_center_x, local_center_y)
    initial_summary = {
        "y_start": int(round(local_center_y - params.extend_length_px / 2.0)),
        "y_end": int(round(local_center_y + params.extend_length_px / 2.0)),
    }

    try:
        if params.locator_mode == "bright-line":
            selected_rows, all_rows, v115_summary = select_bright_line_trenches(
                roi,
                initial_summary,
                fixed,
                block_profile_fn,
                params.locator_majority,
                mask,
            )
        else:
            selected_rows, all_rows, v115_summary = (
                basin_v115.select_basin_first_trenches(
                    roi, initial_summary, fixed, block_profile_fn
                )
            )
        locator_exception = ""
    except Exception as exc:
        selected_rows, all_rows = [], []
        v115_summary = {
            "detection_success": False,
            "reason": f"V1.15 exception: {type(exc).__name__}: {exc}",
        }
        locator_exception = f"{type(exc).__name__}: {exc}"

    from cdsem_refinement import choose_objects

    chosen = choose_objects(
        selected_rows, all_rows, local_center_x, params, choose_center_left_right_v115
    )
    if params.locator_mode == "bright-line":
        chosen = [
            (role, row)
            for role, row in chosen
            if row.get("basin_valid") and row.get("brightness_true_trench")
        ]
    raw_global, gy0, gy1 = global_profile(image, params)
    if region_mask is not None:
        raw_global = supported_profile(image[gy0:gy1], profile_mask[gy0:gy1])
    profile = moving_average_reflect(raw_global, params.smoothing_pixel)
    target_px = params.target_cd_nm / params.pixel_size_nm

    selected_keys = {
        int(r.get("basin_sequence_index", r.get("candidate_id", -99999)))
        for r in selected_rows
    }
    candidates: list[Candidate] = []
    chosen_key_to_role: dict[int, str] = {}

    for role, row in chosen:
        idx = int(row.get("basin_sequence_index", row.get("candidate_id", -1)))
        chosen_key_to_role[idx] = role
        rough_center_local = _row_center(row)
        rough_width = _row_float(row, "estimated_width_px")
        if not np.isfinite(rough_width) or rough_width <= 0:
            rough_width = max(
                3.0,
                _row_float(row, "basin_right_x_px")
                - _row_float(row, "basin_left_x_px"),
            )
        rough_center = rough_center_local + x0
        rough_left = rough_center - rough_width / 2.0
        rough_right = rough_center + rough_width / 2.0

        # IMPORTANT: from here onward use the original V2 exact edge detector.
        edge = detect_edge_pair(
            profile,
            rough_center,
            params,
            anchor_left_x=rough_left,
            anchor_right_x=rough_right,
            search_radius_px=params.tracking_retry_radius_px,
            strict_tracking=False,
        )
        if not edge.valid:
            # Same V2 detector, just retry using Target-CD-based prediction around the V1.15 center.
            edge = detect_edge_pair(profile, rough_center, params)

        if edge.valid:
            center = 0.5 * (edge.left_edge_x + edge.right_edge_x)
            global_left = edge.left_edge_x
            global_right = edge.right_edge_x
            global_cd = edge.local_cd_px
            dark_level = edge.dark_level
            left_strength = edge.left_peak_strength
            right_strength = edge.right_peak_strength
            basin_mean, left_shoulder, right_shoulder, contrast = (
                candidate_basin_metrics(profile, global_left, global_right)
            )
            score, _, _, _, _ = candidate_score(center, edge, profile, x0, x1, params)
            source = "v115_basin+v2_global_refine"
        else:
            # Preserve V1.15 rough geometry as the V2 tracking anchor. The 32-row local profiles
            # in measure_single_space may still recover an edge that is weak in the 360-row profile.
            center = rough_center
            global_left = rough_left
            global_right = rough_right
            global_cd = rough_width
            dark_level = safe_median_segment(
                profile,
                rough_left + 0.30 * rough_width,
                rough_right - 0.30 * rough_width,
            )
            left_strength = math.nan
            right_strength = math.nan
            basin_mean, left_shoulder, right_shoulder, contrast = (
                candidate_basin_metrics(profile, rough_left, rough_right)
            )
            score = _row_float(row, "basin_quality_score", 0.0)
            source = "v115_basin_rough_anchor_v2_global_refine_failed"

        if params.candidate_min_basin_contrast_8bit is not None and (
            not np.isfinite(contrast)
            or contrast < params.candidate_min_basin_contrast_8bit
        ):
            row["v113_candidate_rejected"] = "candidate_min_contrast"
            continue
        candidates.append(
            Candidate(
                center_x_px=float(center),
                predicted_left_px=float(center - target_px / 2.0),
                predicted_right_px=float(center + target_px / 2.0),
                global_left_edge_px=float(global_left),
                global_right_edge_px=float(global_right),
                global_cd_px=float(global_cd),
                dark_level=float(dark_level),
                left_peak_strength=float(left_strength),
                right_peak_strength=float(right_strength),
                mean_grad_strength=float(np.nanmean([left_strength, right_strength]))
                if np.isfinite(left_strength) or np.isfinite(right_strength)
                else math.nan,
                contrast=float(contrast),
                score=float(score),
                is_center_reference=(role == "center"),
                source=source,
                role=role,
                basin_mean=float(basin_mean),
                left_shoulder_mean=float(left_shoulder),
                right_shoulder_mean=float(right_shoulder),
                v115_basin_sequence_index=idx,
                v115_basin_left_x_px=_row_float(row, "basin_left_x_px") + x0,
                v115_basin_right_x_px=_row_float(row, "basin_right_x_px") + x0,
                v115_estimated_width_px=float(rough_width),
                v115_quality_score=_row_float(row, "basin_quality_score"),
                v115_brightness_core_mean=_row_float(row, "brightness_core_mean_raw"),
                v115_original_selected=(idx in selected_keys),
            )
        )

    role_order = {"center": 0, "left": 1, "right": 2}
    candidates.sort(key=lambda c: role_order.get(c.role, 9))
    triplet_complete = {"center", "left", "right"}.issubset(
        {c.role for c in candidates}
    )

    debug_rows: list[dict[str, Any]] = []
    for row in all_rows:
        idx = int(row.get("basin_sequence_index", row.get("candidate_id", -1)))
        out = dict(row)
        # convert important V1.15 local-X fields to original-image X while keeping local copies
        for key in (
            "basin_left_x_px",
            "basin_right_x_px",
            "basin_center_x_px",
            "candidate_center_x_px",
            "adaptive_left_x_px",
            "adaptive_right_x_px",
            "adaptive_center_x_px",
            "basin_center_median_px",
        ):
            val = _row_float(out, key)
            if np.isfinite(val):
                out[key + "_local"] = val
                out[key] = val + x0
        out["v115_original_selected"] = idx in selected_keys
        out["selected"] = idx in chosen_key_to_role
        out["selected_role"] = chosen_key_to_role.get(idx, "")
        debug_rows.append(out)

    diagnostics = {
        **mode_diag,
        "locator_threshold_raw": v115_summary.get("separator_bright_threshold_raw"),
        "global_profile_y0": gy0,
        "global_profile_y1": gy1,
        "roi_truncated": roi_truncated,
        "candidate_seed_count": len(all_rows),
        "candidate_valid_count": sum(_row_bool(r, "basin_valid") for r in all_rows),
        "candidate_filtered_count": len(chosen),
        "candidate_selected_count": len(candidates),
        "candidate_detection_rejected_count": sum(
            "refine_failed" in c.source for c in candidates
        ),
        "candidate_brightness_rejected_count": sum(
            _row_bool(r, "brightness_valid")
            and not _row_bool(r, "brightness_true_trench")
            for r in all_rows
        ),
        "candidate_other_rejected_count": 0,
        "triplet_complete": triplet_complete,
        "selected_roles": ",".join(c.role for c in candidates),
        "candidate_debug_rows": debug_rows,
        "v115_detection_success": bool(v115_summary.get("detection_success", False)),
        "v115_reason": str(v115_summary.get("reason", "")),
        "v115_edge_mode": str(v115_summary.get("edge_mode", "")),
        "v115_locator_exception": locator_exception,
    }
    return candidates, diagnostics


def sample_y_positions(params: MeasurementParams) -> np.ndarray:
    y_start = params.center_y_px - params.extend_length_px / 2.0
    if params.integer_sampling:
        start = int(math.ceil(y_start))
        return np.arange(start, start + params.sample_number, dtype=float)
    y_end = y_start + params.extend_length_px - 1.0
    return np.linspace(y_start, y_end, params.sample_number)


def average_rows_profile(
    image: np.ndarray, y_center: float, average_range_px: int
) -> tuple[Optional[np.ndarray], int, int, int]:
    half = average_range_px // 2
    y0_raw = int(round(y_center)) - half
    y1_raw = y0_raw + average_range_px
    y0 = max(0, y0_raw)
    y1 = min(image.shape[0], y1_raw)
    actual = max(0, y1 - y0)
    if actual <= 0:
        return None, y0, y1, actual
    return image[y0:y1, :].mean(axis=0), y0, y1, actual


def rolling_median_reference(values: np.ndarray, window: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    n = values.size
    reference = np.full(n, np.nan, dtype=np.float64)
    finite_global = values[np.isfinite(values)]
    global_median = float(np.median(finite_global)) if finite_global.size else math.nan
    radius = max(0, window // 2)
    for i in range(n):
        start = max(0, i - radius)
        end = min(n, i + (window - radius))
        local = values[start:end]
        local = local[np.isfinite(local)]
        reference[i] = float(np.median(local)) if local.size else global_median
    return reference


def mark_flyers(
    values: np.ndarray, params: MeasurementParams
) -> tuple[np.ndarray, np.ndarray]:
    reference = rolling_median_reference(values, params.flyer_local_reference_window)
    delta = values - reference
    flyers = (
        np.isfinite(values)
        & np.isfinite(reference)
        & (
            (delta < -params.flyer_left_offset_px)
            | (delta > params.flyer_right_offset_px)
        )
    )
    return flyers, reference


def stats_with_flyer_mode(
    left_px: np.ndarray,
    right_px: np.ndarray,
    local_cd_nm: np.ndarray,
    left_flyer: np.ndarray,
    right_flyer: np.ndarray,
    left_reference: np.ndarray,
    right_reference: np.ndarray,
    params: MeasurementParams,
) -> tuple[float, float, float, float, int]:
    mode = params.flyer_mode.lower()
    left_work = left_px.copy()
    right_work = right_px.copy()
    cd_work = local_cd_nm.copy()

    if mode == "remove":
        keep = ~(left_flyer | right_flyer)
        left_work = left_work[keep]
        right_work = right_work[keep]
        cd_work = cd_work[keep]
    elif mode == "replace":
        left_work[left_flyer] = left_reference[left_flyer]
        right_work[right_flyer] = right_reference[right_flyer]
        cd_work = (right_work - left_work) * params.pixel_size_nm
    # reserve: unchanged

    finite = np.isfinite(left_work) & np.isfinite(right_work) & np.isfinite(cd_work)
    left_work = left_work[finite]
    right_work = right_work[finite]
    cd_work = cd_work[finite]
    n = cd_work.size
    if n == 0:
        return math.nan, math.nan, math.nan, math.nan, int(n)
    if n < max(2, params.standard_deviation_ddof + 1):
        return float(np.mean(cd_work)), math.nan, math.nan, math.nan, int(n)

    mean_cd = float(np.mean(cd_work))
    left_nm = left_work * params.pixel_size_nm
    right_nm = right_work * params.pixel_size_nm
    ler_left = params.roughness_multiplier * float(
        np.std(left_nm - np.mean(left_nm), ddof=params.standard_deviation_ddof)
    )
    ler_right = params.roughness_multiplier * float(
        np.std(right_nm - np.mean(right_nm), ddof=params.standard_deviation_ddof)
    )
    lwr = params.roughness_multiplier * float(
        np.std(cd_work - mean_cd, ddof=params.standard_deviation_ddof)
    )
    return mean_cd, ler_left, ler_right, lwr, int(n)


def robust_series_outliers(
    values: np.ndarray,
    window: int,
    mad_multiplier: float,
    minimum_tolerance: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return outlier mask, rolling reference and a robust global tolerance."""
    reference = rolling_median_reference(values, window)
    residual = values - reference
    finite = residual[np.isfinite(residual)]
    sigma = 1.4826 * robust_mad(finite) if finite.size else 0.0
    tolerance = max(float(minimum_tolerance), float(mad_multiplier) * sigma)
    outliers = (
        np.isfinite(values) & np.isfinite(reference) & (np.abs(residual) > tolerance)
    )
    return outliers, reference, tolerance


def nearest_measured_anchor(
    sample_index: int,
    y_values: np.ndarray,
    left_series: np.ndarray,
    right_series: np.ndarray,
    candidate: Candidate,
) -> tuple[float, float, str]:
    valid = np.flatnonzero(np.isfinite(left_series) & np.isfinite(right_series))
    if valid.size:
        nearest = int(
            valid[np.argmin(np.abs(y_values[valid] - y_values[sample_index]))]
        )
        return (
            float(left_series[nearest]),
            float(right_series[nearest]),
            f"sample_{nearest}",
        )
    return (
        float(candidate.global_left_edge_px),
        float(candidate.global_right_edge_px),
        "global_profile",
    )


def update_sample_row_from_edge(
    row: dict[str, Any], edge: EdgeDetection, local_cd_nm: float
) -> None:
    row.update(
        {
            "left_peak_x": edge.left_peak_x,
            "right_peak_x": edge.right_peak_x,
            "left_peak_strength": edge.left_peak_strength,
            "right_peak_strength": edge.right_peak_strength,
            "left_threshold": edge.left_threshold,
            "right_threshold": edge.right_threshold,
            "dark_level": edge.dark_level,
            "left_peak_level": edge.left_peak_level,
            "right_peak_level": edge.right_peak_level,
            "threshold_reference_mode": edge.threshold_reference_mode,
            "threshold_dark_level": edge.threshold_dark_level,
            "left_threshold_bright_level": edge.left_threshold_bright_level,
            "right_threshold_bright_level": edge.right_threshold_bright_level,
            "left_edge_x": edge.left_edge_x,
            "right_edge_x": edge.right_edge_x,
            "local_cd_nm": local_cd_nm,
            "left_peak_fallback": edge.left_peak_fallback,
            "right_peak_fallback": edge.right_peak_fallback,
            "left_outside_level": edge.left_outside_level,
            "left_inside_level": edge.left_inside_level,
            "right_inside_level": edge.right_inside_level,
            "right_outside_level": edge.right_outside_level,
            "anchor_left_x": edge.anchor_left_x,
            "anchor_right_x": edge.anchor_right_x,
        }
    )


def fit_centerline_and_width_factor(
    y_values: np.ndarray,
    left_px: np.ndarray,
    right_px: np.ndarray,
) -> dict[str, float]:
    """Orthogonal/PCA fit of centerline, reported as both x=k*y+c and y=a*x+b.

    For the near-vertical trench, x=k*y+c is numerically stable.  The equivalent
    y=a*x+b coefficients are also returned as requested.  Perpendicular width is
    horizontal Δx multiplied by |v_y| = 1/sqrt(1+k^2).
    """
    y = np.asarray(y_values, dtype=np.float64)
    left = np.asarray(left_px, dtype=np.float64)
    right = np.asarray(right_px, dtype=np.float64)
    mask = np.isfinite(y) & np.isfinite(left) & np.isfinite(right)
    if int(mask.sum()) < 2:
        return {
            "k": 0.0,
            "c": math.nan,
            "a": math.inf,
            "b": math.nan,
            "angle_from_vertical_deg": 0.0,
            "normal_width_factor": 1.0,
        }

    x = 0.5 * (left[mask] + right[mask])
    yy = y[mask]
    x_mean = float(np.mean(x))
    y_mean = float(np.mean(yy))
    points = np.column_stack([x - x_mean, yy - y_mean])
    cov = points.T @ points / max(points.shape[0], 1)
    values, vectors = np.linalg.eigh(cov)
    v = vectors[:, int(np.argmax(values))]
    vx, vy = float(v[0]), float(v[1])
    if vy < 0:
        vx, vy = -vx, -vy
    norm = math.hypot(vx, vy)
    if norm <= EPS:
        vx, vy = 0.0, 1.0
    else:
        vx, vy = vx / norm, vy / norm

    # near-vertical centerline: x = k*y + c
    if abs(vy) > EPS:
        k = vx / vy
        c = x_mean - k * y_mean
    else:
        k = math.copysign(math.inf, vx if abs(vx) > EPS else 1.0)
        c = math.nan

    # requested equivalent y = a*x + b
    if abs(vx) > EPS:
        a = vy / vx
        b = y_mean - a * x_mean
    else:
        a = math.copysign(math.inf, vy if abs(vy) > EPS else 1.0)
        b = math.nan

    angle = math.degrees(math.atan2(vx, vy))
    factor = abs(vy)  # projection of horizontal Δx onto the line normal
    return {
        "k": float(k),
        "c": float(c),
        "a": float(a),
        "b": float(b),
        "angle_from_vertical_deg": float(angle),
        "normal_width_factor": float(factor),
    }


def measure_single_space(
    image_name: str,
    image: np.ndarray,
    candidate: Candidate,
    space_index: int,
    params: MeasurementParams,
    region_mask=None,
) -> tuple[SpaceResult, list[dict[str, Any]], dict[str, Any]]:
    y_values = sample_y_positions(params)
    from cdsem_regions import sample_in_region

    eligible = np.array(
        [
            sample_in_region(
                region_mask, candidate.center_x_px, y, params.average_range_px
            )
            for y in y_values
        ],
        dtype=bool,
    )
    sample_rows: list[dict[str, Any]] = []

    left_series = np.full(params.sample_number, np.nan, dtype=np.float64)
    right_series = np.full(params.sample_number, np.nan, dtype=np.float64)
    cd_series_nm = np.full(params.sample_number, np.nan, dtype=np.float64)
    peak_fallback = np.zeros(params.sample_number, dtype=bool)
    recovered = np.zeros(params.sample_number, dtype=bool)
    postfilter_rejected = np.zeros(params.sample_number, dtype=bool)

    min_actual_rows = params.minimum_average_range_fraction * params.average_range_px
    min_cd_nm = params.valid_cd_min_factor * params.target_cd_nm
    max_cd_nm = params.valid_cd_max_factor * params.target_cd_nm

    for sample_index, y in enumerate(y_values):
        _, y0, y1, actual_rows = average_rows_profile(image, y, params.average_range_px)
        sample_rows.append(
            {
                "image_name": image_name,
                "space_index": space_index,
                "space_role": candidate.role,
                "space_center_x": candidate.center_x_px,
                "sample_index": sample_index,
                "sample_y": float(y),
                "average_y0": y0,
                "average_y1": y1,
                "actual_average_range": actual_rows,
                "tracking_anchor_source": "",
                "tracking_retry_used": False,
                "recovered": False,
                "postfilter_rejected": False,
                "left_peak_x": math.nan,
                "right_peak_x": math.nan,
                "left_peak_strength": math.nan,
                "right_peak_strength": math.nan,
                "left_threshold": math.nan,
                "right_threshold": math.nan,
                "dark_level": math.nan,
                "left_peak_level": math.nan,
                "right_peak_level": math.nan,
                "threshold_reference_mode": "",
                "threshold_dark_level": math.nan,
                "left_threshold_bright_level": math.nan,
                "right_threshold_bright_level": math.nan,
                "left_outside_level": math.nan,
                "left_inside_level": math.nan,
                "right_inside_level": math.nan,
                "right_outside_level": math.nan,
                "anchor_left_x": math.nan,
                "anchor_right_x": math.nan,
                "left_edge_x": math.nan,
                "right_edge_x": math.nan,
                "local_cd_nm": math.nan,
                "local_cd_x_nm": math.nan,
                "local_cd_slant_nm": math.nan,
                "left_is_flyer": False,
                "right_is_flyer": False,
                "cd_out_of_recommended_range": False,
                "left_peak_fallback": False,
                "right_peak_fallback": False,
                "valid": False,
                "failure_reason": "",
            }
        )

    # Start at Center Y and expand outwards.  Each new point is anchored to the
    # nearest already-valid point, which prevents a single row from jumping to
    # another contrast feature.
    processing_order = sorted(
        range(params.sample_number),
        key=lambda idx: (abs(float(y_values[idx]) - params.center_y_px), idx),
    )

    for sample_index in processing_order:
        y = float(y_values[sample_index])
        row = sample_rows[sample_index]
        row["background_excluded"] = not bool(eligible[sample_index])
        if not eligible[sample_index]:
            row["failure_reason"] = "background_excluded"
            continue
        raw_profile, _, _, actual_rows = average_rows_profile(
            image, y, params.average_range_px
        )
        if raw_profile is None or actual_rows < min_actual_rows:
            row["failure_reason"] = "insufficient_average_rows"
            continue

        anchor_left, anchor_right, anchor_source = nearest_measured_anchor(
            sample_index, y_values, left_series, right_series, candidate
        )
        row["tracking_anchor_source"] = anchor_source
        smoothed = moving_average_reflect(raw_profile, params.smoothing_pixel)
        edge = detect_edge_pair(
            smoothed,
            candidate.center_x_px,
            params,
            anchor_left_x=anchor_left,
            anchor_right_x=anchor_right,
            search_radius_px=params.tracking_search_radius_px,
            strict_tracking=True,
        )

        # If propagation from a neighbouring row failed, retry once from the
        # global edge pair.  This avoids a bad local point poisoning the track.
        if not edge.valid and anchor_source != "global_profile":
            row["tracking_retry_used"] = True
            retry = detect_edge_pair(
                smoothed,
                candidate.center_x_px,
                params,
                anchor_left_x=candidate.global_left_edge_px,
                anchor_right_x=candidate.global_right_edge_px,
                search_radius_px=params.tracking_retry_radius_px,
                strict_tracking=True,
            )
            if retry.valid:
                edge = retry
                row["tracking_anchor_source"] = "global_retry"

        if not edge.valid:
            row["failure_reason"] = edge.failure_reason
            # Preserve whatever partial diagnostic values exist.
            update_sample_row_from_edge(row, edge, math.nan)
            continue

        local_cd_nm = edge.local_cd_px * params.pixel_size_nm
        left_series[sample_index] = edge.left_edge_x
        right_series[sample_index] = edge.right_edge_x
        cd_series_nm[sample_index] = local_cd_nm
        peak_fallback[sample_index] = (
            edge.left_peak_fallback or edge.right_peak_fallback
        )
        update_sample_row_from_edge(row, edge, local_cd_nm)
        row["cd_out_of_recommended_range"] = not (min_cd_nm <= local_cd_nm <= max_cd_nm)
        row["valid"] = True
        row["failure_reason"] = ""

    initial_left = left_series.copy()
    initial_right = right_series.copy()

    # First robust pass identifies suspicious points and supplies rolling anchors
    # for a second, tighter measurement attempt.  Missing points are also retried
    # when neighbouring valid points provide a reliable expected location.
    left_bad, left_ref, _ = robust_series_outliers(
        left_series,
        params.postfilter_window,
        params.postfilter_mad_multiplier,
        params.postfilter_min_edge_tolerance_px,
    )
    right_bad, right_ref, _ = robust_series_outliers(
        right_series,
        params.postfilter_window,
        params.postfilter_mad_multiplier,
        params.postfilter_min_edge_tolerance_px,
    )
    cd_px = right_series - left_series
    cd_bad, cd_ref, _ = robust_series_outliers(
        cd_px,
        params.postfilter_window,
        params.postfilter_mad_multiplier,
        params.postfilter_min_cd_tolerance_px,
    )
    suspicious = left_bad | right_bad | cd_bad

    if params.recover_rejected_points:
        recovery_targets = np.flatnonzero(
            suspicious | ~(np.isfinite(left_series) & np.isfinite(right_series))
        )
        for sample_index in recovery_targets:
            row = sample_rows[int(sample_index)]
            if not eligible[sample_index]:
                continue
            if row["actual_average_range"] < min_actual_rows:
                continue
            anchor_left = (
                float(left_ref[sample_index])
                if np.isfinite(left_ref[sample_index])
                else float(candidate.global_left_edge_px)
            )
            anchor_right = (
                float(right_ref[sample_index])
                if np.isfinite(right_ref[sample_index])
                else float(candidate.global_right_edge_px)
            )
            if anchor_right <= anchor_left:
                continue
            raw_profile, _, _, _ = average_rows_profile(
                image, float(y_values[sample_index]), params.average_range_px
            )
            if raw_profile is None:
                continue
            smoothed = moving_average_reflect(raw_profile, params.smoothing_pixel)
            edge = detect_edge_pair(
                smoothed,
                candidate.center_x_px,
                params,
                anchor_left_x=anchor_left,
                anchor_right_x=anchor_right,
                search_radius_px=params.tracking_retry_radius_px,
                strict_tracking=True,
            )
            if not edge.valid:
                continue
            local_cd_nm = edge.local_cd_px * params.pixel_size_nm
            left_series[sample_index] = edge.left_edge_x
            right_series[sample_index] = edge.right_edge_x
            cd_series_nm[sample_index] = local_cd_nm
            peak_fallback[sample_index] = (
                edge.left_peak_fallback or edge.right_peak_fallback
            )
            recovered[sample_index] = True
            update_sample_row_from_edge(row, edge, local_cd_nm)
            row["tracking_anchor_source"] = "rolling_median_recovery"
            row["tracking_retry_used"] = True
            row["recovered"] = True
            row["valid"] = True
            row["failure_reason"] = ""

    # Final post-filter. Gross jumps are not treated as ordinary Flyer=Reserve
    # points because they are geometric misdetections rather than real roughness.
    left_bad, left_ref, left_tol = robust_series_outliers(
        left_series,
        params.postfilter_window,
        params.postfilter_mad_multiplier,
        params.postfilter_min_edge_tolerance_px,
    )
    right_bad, right_ref, right_tol = robust_series_outliers(
        right_series,
        params.postfilter_window,
        params.postfilter_mad_multiplier,
        params.postfilter_min_edge_tolerance_px,
    )
    cd_px = right_series - left_series
    cd_bad, cd_ref, cd_tol = robust_series_outliers(
        cd_px,
        params.postfilter_window,
        params.postfilter_mad_multiplier,
        params.postfilter_min_cd_tolerance_px,
    )
    final_bad = left_bad | right_bad | cd_bad
    postfilter_rejected[final_bad] = True
    for sample_index in np.flatnonzero(final_bad):
        row = sample_rows[int(sample_index)]
        row["postfilter_rejected"] = True
        row["valid"] = False
        row["failure_reason"] = "robust_postfilter_rejected"
    left_series[final_bad] = np.nan
    right_series[final_bad] = np.nan
    cd_series_nm[final_bad] = np.nan

    # V1.12: recover only final invalid slots, before all statistics.
    import sys as _sys
    from cdsem_refinement import restore_edges
    from cdsem_refinement import refine_path_and_erf, apply_final_flyers

    refine_path_and_erf(
        _sys.modules[__name__],
        image,
        candidate,
        params,
        y_values,
        left_series,
        right_series,
        cd_series_nm,
        sample_rows,
        postfilter_rejected,
    )
    restore_edges(
        _sys.modules[__name__],
        image,
        candidate,
        params,
        y_values,
        left_series,
        right_series,
        cd_series_nm,
        sample_rows,
        postfilter_rejected,
        left_tol,
        right_tol,
        cd_tol,
    )
    apply_final_flyers(
        left_series,
        right_series,
        cd_series_nm,
        sample_rows,
        params,
        _sys.modules[__name__],
    )
    # A forced continuity fill must NEVER turn excluded background into data.
    left_series[~eligible] = np.nan
    right_series[~eligible] = np.nan
    cd_series_nm[~eligible] = np.nan
    for index in np.flatnonzero(~eligible):
        sample_rows[index].update(
            valid=False,
            failure_reason="background_excluded",
            left_edge_x=math.nan,
            right_edge_x=math.nan,
            local_cd_nm=math.nan,
            continuity_synthetic=False,
            continuity_source="background_excluded",
            continuity_original_valid=False,
        )
    synthetic_mask = np.array(
        [r["continuity_synthetic"] for r in sample_rows], dtype=bool
    )

    valid_mask = (
        np.isfinite(left_series) & np.isfinite(right_series) & np.isfinite(cd_series_nm)
    )
    valid_indices = np.flatnonzero(valid_mask)
    left_valid = left_series[valid_mask]
    right_valid = right_series[valid_mask]
    cd_valid_nm = cd_series_nm[valid_mask]

    left_flyer_valid, left_flyer_ref = mark_flyers(left_valid, params)
    right_flyer_valid, right_flyer_ref = mark_flyers(right_valid, params)
    for local_idx, original_idx in enumerate(valid_indices):
        sample_rows[int(original_idx)]["left_is_flyer"] = bool(
            sample_rows[int(original_idx)].get("left_is_flyer", False)
            or left_flyer_valid[local_idx]
        )
        sample_rows[int(original_idx)]["right_is_flyer"] = bool(
            sample_rows[int(original_idx)].get("right_is_flyer", False)
            or right_flyer_valid[local_idx]
        )

    mean_cd, ler_left, ler_right, lwr, stats_count = stats_with_flyer_mode(
        left_valid,
        right_valid,
        cd_valid_nm,
        left_flyer_valid,
        right_flyer_valid,
        left_flyer_ref,
        right_flyer_ref,
        __import__("dataclasses").replace(params, flyer_mode="reserve"),
    )

    # After the original V2 statistics, fit trench direction and report both widths.
    fit = fit_centerline_and_width_factor(y_values, left_series, right_series)
    normal_factor = float(fit["normal_width_factor"])
    cd_slant_valid_nm = cd_valid_nm * normal_factor
    mean_cd_slant = (
        float(np.mean(cd_slant_valid_nm)) if cd_slant_valid_nm.size else math.nan
    )
    lwr_slant = (
        params.roughness_multiplier
        * float(
            np.std(
                cd_slant_valid_nm - np.mean(cd_slant_valid_nm),
                ddof=params.standard_deviation_ddof,
            )
        )
        if cd_slant_valid_nm.size >= max(2, params.standard_deviation_ddof + 1)
        else math.nan
    )
    for original_idx in valid_indices:
        sample_rows[int(original_idx)]["local_cd_x_nm"] = float(
            cd_series_nm[int(original_idx)]
        )
        sample_rows[int(original_idx)]["local_cd_slant_nm"] = float(
            cd_series_nm[int(original_idx)] * normal_factor
        )
        sample_rows[int(original_idx)]["centerline_angle_from_vertical_deg"] = float(
            fit["angle_from_vertical_deg"]
        )
        sample_rows[int(original_idx)]["normal_width_factor"] = normal_factor

    valid_count = int(valid_mask.sum())
    eligible_count = max(1, int(eligible.sum()))
    valid_fraction = valid_count / float(eligible_count)
    stable = float(
        np.sum(valid_mask & ~synthetic_mask)
    ) / eligible_count >= params.minimum_valid_fraction and stats_count >= max(
        2, params.standard_deviation_ddof + 1
    )
    warnings: list[str] = []
    if not stable:
        if synthetic_mask.any():
            detected_fraction = (
                float(np.sum(valid_mask & ~synthetic_mask)) / params.sample_number
            )
            warnings.append(
                f"detected_fraction={detected_fraction:.3f}; filled_fraction={valid_fraction:.3f}; required detected fraction={params.minimum_valid_fraction:.3f}"
            )
        else:
            warnings.append(
                f"valid_fraction={valid_fraction:.3f} below {params.minimum_valid_fraction:.3f}"
            )
    fallback_count = int(peak_fallback[valid_mask].sum())
    hard_rejected_count = int(postfilter_rejected.sum())
    recovered_count = int(recovered.sum())
    out_count = int(
        sum(
            bool(row["cd_out_of_recommended_range"])
            for row in sample_rows
            if row["valid"]
        )
    )
    if synthetic_mask.any():
        warnings.append(
            f"{int(synthetic_mask.sum())} synthetic fills included in metrics; not counted as detected for quality"
        )
    redetected_count = sum(
        r["continuity_source"] == "relaxed_redetect" for r in sample_rows
    )
    if redetected_count:
        warnings.append(f"{redetected_count} samples accepted by continuity relaxation")
    if hard_rejected_count:
        warnings.append(
            f"{hard_rejected_count} gross edge jumps removed by robust postfilter"
        )
    if recovered_count:
        warnings.append(f"{recovered_count} samples recovered by rolling-anchor retry")
    if fallback_count:
        warnings.append(f"{fallback_count} valid samples used gradient-peak fallback")
    if params.flyer_mode.lower() == "reserve" and (
        left_flyer_valid.any() or right_flyer_valid.any()
    ):
        warnings.append("remaining flyers marked and retained in statistics")

    result = SpaceResult(
        image_name=image_name,
        space_index=space_index,
        space_role=candidate.role,
        space_center_x=candidate.center_x_px,
        valid_sample_count=valid_count,
        valid_fraction=valid_fraction,
        stable=stable,
        flyer_left_count=int(left_flyer_valid.sum()),
        flyer_right_count=int(right_flyer_valid.sum()),
        cd_out_of_range_count=out_count,
        peak_fallback_count=fallback_count,
        hard_rejected_count=hard_rejected_count,
        recovered_count=recovered_count,
        mean_cd_nm=mean_cd,
        ler_left_nm=ler_left,
        ler_right_nm=ler_right,
        lwr_nm=lwr,
        mean_cd_x_nm=mean_cd,
        lwr_x_nm=lwr,
        mean_cd_slant_nm=mean_cd_slant,
        lwr_slant_nm=lwr_slant,
        centerline_angle_from_vertical_deg=float(fit["angle_from_vertical_deg"]),
        normal_width_factor=normal_factor,
        fit_x_eq_ky_plus_c_k=float(fit["k"]),
        fit_x_eq_ky_plus_c_c_px=float(fit["c"]),
        fit_y_eq_ax_plus_b_a=float(fit["a"]),
        fit_y_eq_ax_plus_b_b_px=float(fit["b"]),
        warning="; ".join(warnings),
    )

    annotation_payload = {
        "background_excluded": ~eligible,
        "space_index": space_index,
        "candidate": candidate,
        "y_values": y_values,
        "left_edges": left_series,
        "right_edges": right_series,
        "initial_left_edges": initial_left,
        "initial_right_edges": initial_right,
        "postfilter_rejected": postfilter_rejected,
        "recovered": recovered,
        "left_flyers": np.array(
            [bool(row["left_is_flyer"]) for row in sample_rows], dtype=bool
        ),
        "right_flyers": np.array(
            [bool(row["right_is_flyer"]) for row in sample_rows], dtype=bool
        ),
        "stable": stable,
        "centerline_fit": fit,
        "continuity_sources": [r["continuity_source"] for r in sample_rows],
        "viterbi_selected": np.array(
            [r.get("viterbi_selected", False) for r in sample_rows]
        ),
        "erf_fitted": np.array([r.get("erf_fitted", False) for r in sample_rows]),
        "continuity_synthetic": synthetic_mask,
        "continuity_original_valid": np.array(
            [r["continuity_original_valid"] for r in sample_rows]
        ),
        "postfilter_tolerances": {
            "left_px": left_tol,
            "right_px": right_tol,
            "cd_px": cd_tol,
        },
    }
    return result, sample_rows, annotation_payload


def finite_mean(values: Iterable[float]) -> float:
    arr = np.asarray(list(values), dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr)) if arr.size else math.nan


def measure_image(
    path: Path, params: MeasurementParams, region_mask=None
) -> ImageMeasurement:
    image, load_warnings = read_image(path)
    display_image = image.copy()
    if params.pattern_kind == "line":
        image = 255.0 - image
    height, width = image.shape
    warnings = list(load_warnings)

    if (width, height) != (params.expected_width_px, params.expected_height_px):
        raise ValueError(
            f"image size is {width}x{height}; expected "
            f"{params.expected_width_px}x{params.expected_height_px}"
        )

    x0, x1, y0, y1, roi_truncated = measurement_area_bounds(image.shape, params)
    if roi_truncated:
        warnings.append("measurement area was truncated at an image boundary")

    candidates, candidate_diag = detect_space_candidates(image, params, region_mask)
    for candidate_row in candidate_diag.get("candidate_debug_rows", []):
        candidate_row["image_name"] = path.name
    if not candidates:
        image_row = {
            "image_name": path.name,
            "image_width": width,
            "image_height": height,
            "pixel_size_nm": params.pixel_size_nm,
            "center_x": params.center_x_px,
            "center_y": params.center_y_px,
            "candidate_seed_count": candidate_diag.get("candidate_seed_count", 0),
            "triplet_complete": False,
            "selected_roles": "",
            "selected_space_count": 0,
            "stable_space_count": 0,
            "ACD_nm": math.nan,
            "LER_left_nm": math.nan,
            "LER_right_nm": math.nan,
            "LWR_nm": math.nan,
            "ACD_x_nm": math.nan,
            "LWR_x_nm": math.nan,
            "ACD_slant_nm": math.nan,
            "LWR_slant_nm": math.nan,
            "mean_centerline_angle_from_vertical_deg": math.nan,
            "valid": False,
            "warning": "; ".join(warnings + ["no valid SPACE candidate found"]),
        }
        return ImageMeasurement(
            image_row=image_row,
            space_rows=[],
            sample_rows=[],
            annotation_payload={
                "image": image,
                "display_image": display_image,
                "bounds": (x0, x1, y0, y1),
                "spaces": [],
                "candidate_diag": candidate_diag,
            },
        )

    space_results: list[SpaceResult] = []
    sample_rows: list[dict[str, Any]] = []
    space_annotations: list[dict[str, Any]] = []
    for space_index, candidate in enumerate(candidates, start=1):
        result, rows, annotation = measure_single_space(
            path.name, image, candidate, space_index, params, region_mask
        )
        space_results.append(result)
        sample_rows.extend(rows)
        space_annotations.append(annotation)

    stable_results = [result for result in space_results if result.stable]
    aggregate_results = (
        stable_results
        if len(stable_results) >= params.min_number
        else [
            result
            for result in space_results
            if np.isfinite(result.mean_cd_nm)
            and np.isfinite(result.ler_left_nm)
            and np.isfinite(result.ler_right_nm)
            and np.isfinite(result.lwr_nm)
        ]
    )

    triplet_complete = bool(candidate_diag.get("triplet_complete", False))
    image_valid = len(stable_results) >= params.min_number
    if params.require_center_left_right:
        image_valid = image_valid and triplet_complete
    if not triplet_complete and params.require_center_left_right:
        warnings.append(
            "centre-left-right triplet incomplete; at least one adjacent trench was not accepted"
        )
    if len(stable_results) < params.min_number:
        warnings.append(
            f"stable SPACE count {len(stable_results)} below Min Number={params.min_number}; "
            "aggregate uses available unstable SPACE results when possible"
        )
    for result in space_results:
        if result.warning:
            warnings.append(f"SPACE {result.space_index}: {result.warning}")

    image_row = {
        "image_name": path.name,
        "image_width": width,
        "image_height": height,
        "pixel_size_nm": params.pixel_size_nm,
        "center_x": params.center_x_px,
        "center_y": params.center_y_px,
        "candidate_seed_count": candidate_diag.get("candidate_seed_count", 0),
        "triplet_complete": triplet_complete,
        "selected_roles": candidate_diag.get("selected_roles", ""),
        "candidate_valid_count": candidate_diag.get("candidate_valid_count", 0),
        "candidate_brightness_rejected_count": candidate_diag.get(
            "candidate_brightness_rejected_count", 0
        ),
        "selected_space_count": len(space_results),
        "stable_space_count": len(stable_results),
        # Compatibility columns remain the original V2 x-axis results.
        "ACD_nm": finite_mean(result.mean_cd_nm for result in aggregate_results),
        "LER_left_nm": finite_mean(result.ler_left_nm for result in aggregate_results),
        "LER_right_nm": finite_mean(
            result.ler_right_nm for result in aggregate_results
        ),
        "LWR_nm": finite_mean(result.lwr_nm for result in aggregate_results),
        "ACD_x_nm": finite_mean(result.mean_cd_x_nm for result in aggregate_results),
        "LWR_x_nm": finite_mean(result.lwr_x_nm for result in aggregate_results),
        "ACD_slant_nm": finite_mean(
            result.mean_cd_slant_nm for result in aggregate_results
        ),
        "LWR_slant_nm": finite_mean(
            result.lwr_slant_nm for result in aggregate_results
        ),
        "mean_centerline_angle_from_vertical_deg": finite_mean(
            result.centerline_angle_from_vertical_deg for result in aggregate_results
        ),
        "valid": image_valid,
        "warning": "; ".join(warnings),
    }

    return ImageMeasurement(
        image_row=image_row,
        space_rows=[asdict(result) for result in space_results],
        sample_rows=sample_rows,
        annotation_payload={
            "image": image,
            "display_image": display_image,
            "bounds": (x0, x1, y0, y1),
            "spaces": space_annotations,
            "candidate_diag": candidate_diag,
        },
    )


def save_annotated_image(
    source_name: str,
    payload: dict[str, Any],
    output_path: Path,
    params: MeasurementParams,
) -> None:
    image = payload.get("display_image", payload["image"])
    x0, x1, y0, y1 = payload["bounds"]
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    ax.imshow(image, cmap="gray", vmin=0, vmax=255, origin="upper")
    ax.add_patch(
        Rectangle(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            fill=False,
            linewidth=1.2,
            linestyle="--",
            label="Meas Area",
        )
    )
    ax.axvline(params.center_x_px, linewidth=0.8, linestyle=":")
    ax.axhline(params.center_y_px, linewidth=0.8, linestyle=":")
    y_start = params.center_y_px - params.extend_length_px / 2.0
    y_end = y_start + params.extend_length_px - 1.0
    ax.axhline(y_start, linewidth=0.7, linestyle="--")
    ax.axhline(y_end, linewidth=0.7, linestyle="--")

    for annotation in payload.get("spaces", []):
        candidate: Candidate = annotation["candidate"]
        y_values = annotation["y_values"]
        left_edges = annotation["left_edges"]
        right_edges = annotation["right_edges"]
        left_flyers = annotation["left_flyers"]
        right_flyers = annotation["right_flyers"]
        initial_left = annotation.get("initial_left_edges", left_edges)
        initial_right = annotation.get("initial_right_edges", right_edges)
        rejected = annotation.get(
            "postfilter_rejected", np.zeros_like(y_values, dtype=bool)
        )
        recovered = annotation.get("recovered", np.zeros_like(y_values, dtype=bool))

        ax.axvline(
            candidate.center_x_px,
            linewidth=0.9,
            alpha=0.8,
            label=f"{candidate.role} {params.pattern_kind} center",
        )
        ax.axvline(
            candidate.global_left_edge_px, linewidth=0.55, linestyle="--", alpha=0.6
        )
        ax.axvline(
            candidate.global_right_edge_px, linewidth=0.55, linestyle="--", alpha=0.6
        )
        from cdsem_refinement import draw_overlay

        draw_overlay(ax, annotation, params.edge_continuity)
        if params.viterbi_enabled or params.erf_enabled:
            for values in (left_edges, right_edges):
                if params.viterbi_enabled:
                    ax.plot(
                        values,
                        y_values,
                        "-",
                        linewidth=0.65,
                        color="deepskyblue",
                        label="Viterbi path",
                    )
                fitted = annotation.get(
                    "erf_fitted", np.zeros_like(y_values, dtype=bool)
                ) & np.isfinite(values)
                if fitted.any():
                    ax.plot(
                        values[fitted],
                        y_values[fitted],
                        "s",
                        markersize=3,
                        fillstyle="none",
                        color="purple",
                        label="ERF fitted",
                    )
        fit = annotation.get("centerline_fit", {})
        k = float(fit.get("k", math.nan))
        c = float(fit.get("c", math.nan))
        if np.isfinite(k) and np.isfinite(c):
            fit_y = np.array([float(np.nanmin(y_values)), float(np.nanmax(y_values))])
            fit_x = k * fit_y + c
            ax.plot(
                fit_x,
                fit_y,
                "-",
                linewidth=0.7,
                alpha=0.7,
                label=f"{candidate.role} center fit",
            )

        valid_left = np.isfinite(left_edges)
        valid_right = np.isfinite(right_edges)
        ax.plot(
            left_edges[valid_left],
            y_values[valid_left],
            ".",
            markersize=2.5,
            label=f"{candidate.role} accepted left",
        )
        ax.plot(
            right_edges[valid_right],
            y_values[valid_right],
            ".",
            markersize=2.5,
            label=f"{candidate.role} accepted right",
        )
        rejected_left = rejected & np.isfinite(initial_left)
        rejected_right = rejected & np.isfinite(initial_right)
        if np.any(rejected_left):
            ax.plot(
                initial_left[rejected_left],
                y_values[rejected_left],
                "x",
                markersize=4.5,
                label=f"{candidate.role} rejected",
            )
        if np.any(rejected_right):
            ax.plot(
                initial_right[rejected_right],
                y_values[rejected_right],
                "x",
                markersize=4.5,
            )
        recovered_left = recovered & valid_left
        recovered_right = recovered & valid_right
        if np.any(recovered_left):
            ax.plot(
                left_edges[recovered_left],
                y_values[recovered_left],
                "o",
                fillstyle="none",
                markersize=4,
                label=f"{candidate.role} recovered",
            )
        if np.any(recovered_right):
            ax.plot(
                right_edges[recovered_right],
                y_values[recovered_right],
                "o",
                fillstyle="none",
                markersize=4,
            )
        if np.any(left_flyers & valid_left):
            mask = left_flyers & valid_left
            ax.plot(left_edges[mask], y_values[mask], "+", markersize=4)
        if np.any(right_flyers & valid_right):
            mask = right_flyers & valid_right
            ax.plot(right_edges[mask], y_values[mask], "+", markersize=4)

    ax.set_xlim(0, image.shape[1] - 1)
    ax.set_ylim(image.shape[0] - 1, 0)
    ax.set_xlabel("X (pixel)")
    ax.set_ylabel("Y (pixel)")
    ax.set_title(source_name)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        # Deduplicate labels to keep the legend readable.
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc="upper right", fontsize=7)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
