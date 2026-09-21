"""V1_16: consolidated active implementation, derived from V13_modified."""

from __future__ import annotations
import math
import numpy as np

import re
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional
import cv2
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths


DEFAULT_ROOT = Path(".")


NUMBER_PATTERN = r"\d+(?:\.\d+)?"


@dataclass
class GeneralConfig:
    root_dir: Path = DEFAULT_ROOT
    output_dir: Optional[Path] = None
    pixel_size_nm: float = 1.3181
    magnification_token: str = "200K"
    image_suffix: str = ".png"
    force_pattern: Optional[str] = None
    via_reference_nm: Optional[float] = None
    trench_reference_nm: Optional[float] = None
    space_reference_nm: Optional[float] = None
    save_debug_masks: bool = False
    overwrite_output: bool = True
    continue_on_error: bool = True
    # 排除非常靠近图像边界、容易被截断的结构。相对设计尺寸再取更大值。
    absolute_border_margin_px: int = 3

    def resolved_output_dir(self) -> Path:
        if self.output_dir is not None:
            return self.output_dir
        return self.root_dir / "CD_measure_output_200K_V1_16"

    def validate(self) -> None:
        if self.pixel_size_nm <= 0:
            raise ValueError("pixel_size_nm 必须大于 0")
        if not self.root_dir.exists():
            raise FileNotFoundError(f"输入根目录不存在：{self.root_dir}")
        if not self.root_dir.is_dir():
            raise NotADirectoryError(f"输入路径不是文件夹：{self.root_dir}")
        if (self.force_pattern or "").casefold() not in {"trench", "line"}:
            raise ValueError("必须通过 pattern 指定整个输入目录为 trench 或 line")
        for name, value in (
            ("via_reference_nm", self.via_reference_nm),
            ("trench_reference_nm", self.trench_reference_nm),
            ("space_reference_nm", self.space_reference_nm),
        ):
            if value is not None and value <= 0:
                raise ValueError(f"{name} 必须大于 0")


@dataclass
class ParsedImage:
    path: Path
    pattern: str
    folder_name: str
    folder_nominal_nm: float
    pretreatment_nm: float
    dose_mj: float
    sample_id: str
    design_via_nm: float = math.nan
    design_trench_nm: float = math.nan
    design_space_nm: float = math.nan
    metadata_source: str = "image_estimated"
    # 为兼容 V1.6 的输出表结构保留；V1.7 不做逐图自动分类，因此恒为 NaN。
    auto_pattern_score: float = math.nan
    auto_gradient_x_to_y_ratio: float = math.nan


def moving_average_reflect(profile: np.ndarray, window: int) -> np.ndarray:
    p = np.asarray(profile, dtype=np.float64)
    w = max(1, int(window))
    if w <= 1:
        return p.copy()
    if w % 2 == 0:
        w += 1
    pad = w // 2
    padded = np.pad(p, (pad, pad), mode="reflect")
    kernel = np.ones(w, dtype=np.float64) / float(w)
    return np.convolve(padded, kernel, mode="valid")


def unicode_imread(path: Path) -> np.ndarray:
    try:
        raw = np.fromfile(str(path), dtype=np.uint8)
        img = cv2.imdecode(raw, cv2.IMREAD_UNCHANGED)
    except Exception as exc:
        raise RuntimeError(
            f"读取图像失败：{path}；{type(exc).__name__}: {exc}"
        ) from exc
    if img is None:
        raise RuntimeError(f"OpenCV 无法解码图像：{path}")
    return img


def normalize_gray_uint8(image: np.ndarray) -> np.ndarray:
    img = np.asarray(image)
    if img.ndim == 3:
        if img.shape[2] == 4:
            img = cv2.cvtColor(img, cv2.COLOR_BGRA2GRAY)
        else:
            img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    if img.ndim != 2:
        raise ValueError(f"期望二维灰度图，实际 shape={img.shape}")
    if img.dtype == np.uint8:
        return img.copy()
    f = img.astype(np.float32)
    finite = np.isfinite(f)
    if not np.any(finite):
        raise ValueError("图像中没有有限灰度值")
    if not np.all(finite):
        f = np.where(finite, f, float(np.median(f[finite])))
    lo, hi = np.percentile(f, [1.0, 99.0])
    if hi <= lo + 1e-6:
        lo, hi = float(np.min(f)), float(np.max(f))
    if hi <= lo + 1e-6:
        return np.full(img.shape, 127, dtype=np.uint8)
    return np.clip((f - lo) * 255.0 / (hi - lo), 0, 255).astype(np.uint8)


def read_gray_png(path: Path) -> np.ndarray:
    return normalize_gray_uint8(unicode_imread(path))


def unicode_imwrite(path: Path, image: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower() or ".png"
    ok, buf = cv2.imencode(suffix, image)
    if not ok:
        raise RuntimeError(f"无法编码输出图像：{path}")
    buf.tofile(str(path))


def parse_folder_nominal_nm(name: str) -> float:
    m = re.search(rf"(?P<n>{NUMBER_PATTERN})\s*NM", name, flags=re.IGNORECASE)
    return float(m.group("n")) if m else math.nan


def is_png(path: Path) -> bool:
    """只检查文件类型；不检查文件名、倍率字段或所在文件夹名称。"""
    return path.is_file() and path.suffix.casefold() == ".png"


def _binary_run_lengths(mask: np.ndarray) -> list[tuple[bool, int, int, int]]:
    values = np.asarray(mask, dtype=bool).ravel()
    if values.size == 0:
        return []
    changes = np.flatnonzero(values[1:] != values[:-1]) + 1
    bounds = np.concatenate(([0], changes, [values.size]))
    return [
        (bool(values[start]), int(end - start), int(start), int(end))
        for start, end in zip(bounds[:-1], bounds[1:])
    ]


def estimate_trench_references_nm(
    gray: np.ndarray, pixel_size_nm: float
) -> tuple[float, float]:
    """从全高列均值的暗/亮游程估计 Trench 宽度和 Space 宽度。"""
    W = gray.shape[1]
    profile = gray.mean(axis=0).astype(np.float32)
    smooth_window = max(3, min(21, int(round(W / 120))))
    if smooth_window % 2 == 0:
        smooth_window += 1
    smooth = moving_average_reflect(profile, smooth_window).astype(np.float32)
    normalized = cv2.normalize(
        smooth.reshape(1, -1), None, 0, 255, cv2.NORM_MINMAX
    ).astype(np.uint8)
    threshold, _ = cv2.threshold(
        normalized, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
    )
    from cdsem_localization import interior_band_threshold

    interior_threshold = interior_band_threshold(smooth)
    dark_mask = (
        smooth < interior_threshold
        if interior_threshold is not None
        else normalized.ravel() <= threshold
    )
    runs = _binary_run_lengths(dark_mask)
    dark_widths = [
        length
        for is_dark, length, start, end in runs
        if is_dark and start > 0 and end < W and 3 <= length <= 0.30 * W
    ]
    bright_widths = [
        length
        for is_dark, length, start, end in runs
        if (not is_dark) and start > 0 and end < W and 3 <= length <= 0.45 * W
    ]
    trench_px = float(np.median(dark_widths)) if dark_widths else 30.0 / pixel_size_nm
    space_px = (
        float(np.median(bright_widths)) if bright_widths else 70.0 / pixel_size_nm
    )
    trench_px = float(np.clip(trench_px, 4.0, max(6.0, 0.25 * W)))
    space_px = float(np.clip(space_px, 3.0, max(5.0, 0.45 * W)))
    return trench_px * pixel_size_nm, space_px * pixel_size_nm


def parse_image_from_content(
    path: Path,
    folder_name: str,
    gray: np.ndarray,
    config: GeneralConfig,
) -> ParsedImage:
    # V1.7：图案类型只由调用参数决定，并统一应用于 root 及全部子文件夹。
    # 不再依据单张图内容、文件名或文件夹名改变 Via/Trench 类型。
    pattern = (config.force_pattern or "").casefold()
    folder_nm = parse_folder_nominal_nm(folder_name)
    common = dict(
        path=path,
        pattern=pattern,
        folder_name=folder_name,
        folder_nominal_nm=folder_nm,
        pretreatment_nm=math.nan,
        dose_mj=math.nan,
        sample_id=path.stem,
        auto_pattern_score=math.nan,
        auto_gradient_x_to_y_ratio=math.nan,
    )
    work_gray = 255.0 - gray if pattern == "line" else gray
    trench_nm, space_nm = estimate_trench_references_nm(work_gray, config.pixel_size_nm)
    if config.trench_reference_nm is not None:
        trench_nm = float(config.trench_reference_nm)
    if config.space_reference_nm is not None:
        space_nm = float(config.space_reference_nm)
    if config.trench_reference_nm is not None and config.space_reference_nm is not None:
        source = "cli_reference"
    elif (
        config.trench_reference_nm is not None or config.space_reference_nm is not None
    ):
        source = "image_estimated_or_partial_cli"
    else:
        source = "image_estimated"
    return ParsedImage(
        **common,
        design_trench_nm=trench_nm,
        design_space_nm=space_nm,
        metadata_source=source,
    )


def discover_images(
    config: GeneralConfig,
) -> tuple[list[ParsedImage], list[dict[str, Any]]]:
    records: list[ParsedImage] = []
    errors: list[dict[str, Any]] = []
    root = config.root_dir
    out_dir = config.resolved_output_dir().resolve()
    output_markers = {}

    def is_previous_output(directory):
        if directory not in output_markers:
            try:
                settings = json.loads(
                    (directory / "settings.json").read_text(encoding="utf-8")
                )
                output_markers[directory] = settings.get("script_version") in {
                    "V13_modified",
                    "V1_14",
                    "V1_15",
                    "V1_16",
                }
            except (OSError, ValueError, AttributeError):
                output_markers[directory] = False
        return output_markers[directory]

    paths = sorted(
        (p for p in root.rglob("*") if is_png(p)), key=lambda p: str(p).casefold()
    )
    for path in paths:
        try:
            resolved = path.resolve()
            if resolved == out_dir or out_dir in resolved.parents:
                continue
            parents = [
                parent
                for parent in path.parents
                if parent != root and root in parent.parents
            ]
            if any(is_previous_output(parent) for parent in parents):
                continue
        except Exception:
            pass
        try:
            relative_parent = path.parent.relative_to(root)
            # 保存原始相对目录；根目录中的图片使用 "."，便于输出镜像输入目录树。
            folder_name = relative_parent.as_posix() if relative_parent.parts else "."
            gray = read_gray_png(path)
            records.append(parse_image_from_content(path, folder_name, gray, config))
        except Exception as exc:
            errors.append(
                {
                    "image": path.name,
                    "path": str(path),
                    "pattern": config.force_pattern,
                    "stage": "content_scale_estimation",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                }
            )
    return records, errors


class NoPatternRegion(ValueError):
    """No supported line/trench region; this is a skip, not zero roughness."""


@dataclass
class PatternRegion:
    mask: np.ndarray
    bounds: tuple[int, int, int, int]
    diagnostics: dict


def trim_region_ends(region, fraction=0.05):
    """One approximate length estimator: mean foreground y-span per column.

    For near-vertical equal-width structures this approximates mean trench length.
    Trim the SAME ceil(fraction * mean_length) at both ends of each supported
    column, preserving all interior background gaps. Coordinates remain original.
    """
    if not np.isfinite(fraction) or not 0 <= fraction < 0.5:
        raise ValueError("end-trim-fraction必须在[0,0.5)之间")
    mask = region.mask
    supported = mask.any(axis=0)
    if not supported.any():
        raise NoPatternRegion("没有可裁剪的trench区域")
    starts = mask.argmax(axis=0)
    ends = mask.shape[0] - mask[::-1].argmax(axis=0)
    mean_length = float(np.mean((ends - starts)[supported]))
    trim = int(np.ceil(fraction * mean_length))
    yy = np.arange(mask.shape[0])[:, None]
    kept = mask & (yy >= starts[None, :] + trim) & (yy < ends[None, :] - trim)
    ys, xs = np.nonzero(kept)
    if not len(ys):
        raise NoPatternRegion("端部裁剪后无有效trench区域")
    bounds = (int(xs.min()), int(xs.max() + 1), int(ys.min()), int(ys.max() + 1))
    diag = dict(region.diagnostics)
    diag.update(
        roi_length_estimator="mean_foreground_column_y_span",
        roi_mean_trench_length_px=mean_length,
        roi_end_trim_fraction=float(fraction),
        roi_end_trim_px=trim,
        roi_end_excluded_pixels=int(mask.sum() - kept.sum()),
        roi_before_trim_x0=region.bounds[0],
        roi_before_trim_x1=region.bounds[1],
        roi_before_trim_y0=region.bounds[2],
        roi_before_trim_y1=region.bounds[3],
        roi_x0=bounds[0],
        roi_x1=bounds[1],
        roi_y0=bounds[2],
        roi_y1=bounds[3],
        roi_foreground_fraction=float(kept.mean()),
    )
    return PatternRegion(kept, bounds, diag)


def manual_pattern_region(gray, args):
    """Explicit manual observation rectangle, still subject to end exclusion."""
    height, width = gray.shape
    cx = args.center_x if args.center_x is not None else (width - 1) / 2
    cy = args.center_y if args.center_y is not None else (height - 1) / 2
    w = args.meas_area_width if args.meas_area_width is not None else width
    h = args.meas_area_height if args.meas_area_height is not None else height
    x0, y0 = max(0, int(round(cx - w / 2))), max(0, int(round(cy - h / 2)))
    x1, y1 = min(width, x0 + w), min(height, y0 + h)
    mask = np.zeros(gray.shape, bool)
    mask[y0:y1, x0:x1] = True
    return PatternRegion(mask, (x0, x1, y0, y1), {"roi_mode": "manual"})


def find_pattern_region(
    gray, target_px, pattern="trench", min_contrast=6.0, minimum_length=24
):
    """Select persistent, width-plausible paired boundaries, not plain texture.

    Block profiles are used ONLY for locating the region. Measurement always
    uses the unchanged input pixels. Per-row inside/outside contrast rejects
    solid backgrounds, horizontal bands, and nonpersistent texture. This is a
    geometric heuristic, not a semantic classifier of arbitrary SEM patterns.
    """
    from cdsem_localization import interior_band_threshold
    from cdsem_locator import _false_runs

    raw = np.asarray(gray, dtype=float)
    if raw.ndim != 2 or not np.isfinite(raw).all():
        raise ValueError("背景识别要求有限的二维灰度图")
    work = 255.0 - raw if pattern == "line" else raw
    height, width = work.shape
    target_px = max(4.0, float(target_px))
    mask = np.zeros((height, width), dtype=bool)
    smooth = gaussian_filter1d(work, sigma=1, axis=1, mode="reflect")
    block = max(8, min(24, int(target_px / 3)))
    step = max(4, block // 2)
    checked, accepted = 0, 0
    for y0 in range(0, height, step):
        y1 = min(height, y0 + block)
        if y1 - y0 < 4:
            continue
        profile = np.median(smooth[y0:y1], axis=0)
        # Horizontal high-frequency residual is an empirical noise floor.
        residual = profile - gaussian_filter1d(profile, 2.0, mode="reflect")
        noise = 1.4826 * np.median(np.abs(residual - np.median(residual)))
        threshold = max(float(min_contrast), 5.0 * float(noise))
        peaks, properties = find_peaks(
            -profile, prominence=threshold, distance=max(3, int(0.6 * target_px))
        )
        widths, _, lefts, rights = peak_widths(-profile, peaks, rel_height=0.5)
        envelope_threshold = interior_band_threshold(profile, target_px)
        if envelope_threshold is not None:
            intervals = [
                (a, b)
                for a, b in _false_runs(profile >= envelope_threshold)
                if a > 0 and b < width and 0.65 * target_px <= b - a <= 1.35 * target_px
            ]
            if intervals:
                lefts = np.array([a - 0.5 for a, b in intervals])
                rights = np.array([b - 0.5 for a, b in intervals])
                widths = rights - lefts
                peaks = (lefts + rights) / 2
        for peak, size, left, right in zip(peaks, widths, lefts, rights):
            checked += 1
            if not 0.45 * target_px <= size <= 1.8 * target_px:
                continue
            band = max(3, int(0.18 * size))
            l, r = int(np.floor(left)), int(np.ceil(right))
            if l - band < 0 or r + band >= width:
                continue
            a, b = int(left + 0.3 * size), int(right - 0.3 * size) + 1
            if b <= a:
                continue
            inside = np.median(smooth[y0:y1, a:b], axis=1)
            outside_l = np.median(smooth[y0:y1, l - band : l], axis=1)
            outside_r = np.median(smooth[y0:y1, r + 1 : r + band + 1], axis=1)
            present = (
                np.minimum(outside_l - inside, outside_r - inside) >= threshold * 0.6
            )
            if present.mean() < 0.6:
                continue
            accepted += 1
            margin = max(4, int(0.30 * target_px))
            xx0, xx1 = max(0, l - margin), min(width, r + margin + 1)
            mask[y0:y1, xx0:xx1] |= present[:, None]
    # Remove short disconnected features; retain separate valid regions/gaps.
    count, labels, stats, _ = cv2.connectedComponentsWithStats(mask.astype(np.uint8), 8)
    kept = np.zeros_like(mask)
    for label in range(1, count):
        if stats[label, cv2.CC_STAT_HEIGHT] >= minimum_length:
            kept[labels == label] = True
    ys, xs = np.nonzero(kept)
    if not len(ys):
        raise NoPatternRegion(
            "未找到持续的 trench/line 成对边界；背景已跳过。可检查 ROI 图或使用 --no-auto-roi"
        )
    bounds = (int(xs.min()), int(xs.max() + 1), int(ys.min()), int(ys.max() + 1))
    return PatternRegion(
        kept,
        bounds,
        dict(
            roi_mode="auto",
            roi_x0=bounds[0],
            roi_x1=bounds[1],
            roi_y0=bounds[2],
            roi_y1=bounds[3],
            roi_foreground_fraction=float(kept.mean()),
            roi_candidate_checks=checked,
            roi_accepted_blocks=accepted,
            roi_min_contrast_8bit=min_contrast,
            roi_min_length_px=minimum_length,
        ),
    )


def sample_in_region(mask, center_x, y, average_range):
    """The full along-line averaging window must stay in supported pixels."""
    if mask is None:
        return True
    x = int(round(center_x))
    start = int(round(y)) - average_range // 2
    end = start + average_range
    return (
        0 <= x < mask.shape[1]
        and 0 <= start < end <= mask.shape[0]
        and bool(mask[start:end, x].all())
    )


def save_region_overlay(gray, region, path):
    color = cv2.cvtColor(np.clip(gray, 0, 255).astype(np.uint8), cv2.COLOR_GRAY2BGR)
    outside = ~region.mask
    color[outside] = (color[outside].astype(float) * 0.35).astype(np.uint8)
    x0, x1, y0, y1 = region.bounds
    cv2.rectangle(color, (x0, y0), (x1 - 1, y1 - 1), (0, 220, 0), 1)
    unicode_imwrite(path, color)
