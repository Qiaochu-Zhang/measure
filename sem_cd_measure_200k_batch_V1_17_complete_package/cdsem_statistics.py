"""V1_17: consolidated active implementation, derived from V13_modified."""

from __future__ import annotations
import math
import numpy as np

from pathlib import Path
from typing import Any, Iterable
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cdsem_engine as edges


EPS = 1e-12


def finite_array(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=np.float64)
    return arr[np.isfinite(arr)]


def finite_mean(values: Iterable[float]) -> float:
    arr = finite_array(values)
    return float(np.mean(arr)) if arr.size else math.nan


def three_sigma(values: np.ndarray, ddof: int = 0) -> float:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size < max(2, ddof + 1):
        return math.nan
    return 3.0 * float(np.std(arr - np.mean(arr), ddof=ddof))


def group_mean(values: np.ndarray, group_size: int) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    output: list[float] = []
    for start in range(0, data.size, group_size):
        block = data[start : start + group_size]
        finite = block[np.isfinite(block)]
        output.append(float(np.mean(finite)) if finite.size else math.nan)
    return np.asarray(output, dtype=np.float64)


def detrend_preserve_mean(y: np.ndarray, values: np.ndarray, degree: int) -> np.ndarray:
    yy = np.asarray(y, dtype=np.float64)
    vv = np.asarray(values, dtype=np.float64)
    result = vv.copy()
    mask = np.isfinite(yy) & np.isfinite(vv)
    if int(mask.sum()) <= degree:
        return result
    coefficients = np.polyfit(yy[mask], vv[mask], degree)
    fit = np.polyval(coefficients, yy[mask])
    result[mask] = vv[mask] - fit + float(np.mean(vv[mask]))
    return result


def apply_stat_mode(
    y_values: np.ndarray,
    left_px: np.ndarray,
    right_px: np.ndarray,
    mode: str,
    v5_params: edges.MeasurementParams,
    group_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = np.asarray(y_values, dtype=np.float64)
    left = np.asarray(left_px, dtype=np.float64).copy()
    right = np.asarray(right_px, dtype=np.float64).copy()

    if mode == "all_submean":
        return y, left, right
    if mode == "group4_mean":
        return (
            group_mean(y, group_size),
            group_mean(left, group_size),
            group_mean(right, group_size),
        )
    if mode == "every4":
        indices = np.arange(0, y.size, group_size, dtype=int)
        return y[indices], left[indices], right[indices]
    if mode == "every2":
        indices = np.arange(0, y.size, 2, dtype=int)
        return y[indices], left[indices], right[indices]
    if mode in {"linear_detrend", "quadratic_detrend"}:
        degree = 1 if mode == "linear_detrend" else 2
        return (
            y,
            detrend_preserve_mean(y, left, degree),
            detrend_preserve_mean(y, right, degree),
        )

    valid = np.isfinite(left) & np.isfinite(right)
    valid_indices = np.flatnonzero(valid)
    if valid_indices.size:
        left_valid = left[valid]
        right_valid = right[valid]
        left_flyer, left_ref = edges.mark_flyers(left_valid, v5_params)
        right_flyer, right_ref = edges.mark_flyers(right_valid, v5_params)
        if mode == "flyer_remove":
            remove = left_flyer | right_flyer
            left[valid_indices[remove]] = np.nan
            right[valid_indices[remove]] = np.nan
            return y, left, right
        if mode == "flyer_replace":
            left_valid[left_flyer] = left_ref[left_flyer]
            right_valid[right_flyer] = right_ref[right_flyer]
            left[valid_indices] = left_valid
            right[valid_indices] = right_valid
            return y, left, right
    if mode in {"flyer_remove", "flyer_replace"}:
        return y, left, right
    raise ValueError(f"unknown stat mode: {mode}")


def sequence_metrics(
    y_values: np.ndarray,
    left_px: np.ndarray,
    right_px: np.ndarray,
    pixel_size_nm: float,
    normal_factor: float,
    ddof: int,
) -> dict[str, Any]:
    y = np.asarray(y_values, dtype=np.float64)
    left = np.asarray(left_px, dtype=np.float64)
    right = np.asarray(right_px, dtype=np.float64)
    valid = np.isfinite(y) & np.isfinite(left) & np.isfinite(right) & (right > left)
    cd_x = np.full(left.shape, np.nan, dtype=np.float64)
    cd_x[valid] = (right[valid] - left[valid]) * pixel_size_nm
    cd_slant = cd_x * normal_factor
    left_nm = left * pixel_size_nm
    right_nm = right * pixel_size_nm
    center_nm = 0.5 * (left_nm + right_nm)

    rho = math.nan
    sigma_left = math.nan
    sigma_right = math.nan
    covariance = math.nan
    lwr_from_components = math.nan
    if int(valid.sum()) >= 3:
        l_res = left_nm[valid] - np.mean(left_nm[valid])
        r_res = right_nm[valid] - np.mean(right_nm[valid])
        sigma_left = float(np.std(l_res, ddof=ddof))
        sigma_right = float(np.std(r_res, ddof=ddof))
        covariance = (
            float(np.mean(l_res * r_res))
            if ddof == 0
            else float(np.cov(l_res, r_res, ddof=ddof)[0, 1])
        )
        if sigma_left > EPS and sigma_right > EPS:
            rho = float(covariance / (sigma_left * sigma_right))
        width_variance = max(0.0, sigma_left**2 + sigma_right**2 - 2.0 * covariance)
        lwr_from_components = 3.0 * math.sqrt(width_variance)

    return {
        "y": y,
        "left": left,
        "right": right,
        "valid": valid,
        "cd_x": cd_x,
        "cd_slant": cd_slant,
        "mean_cd_px": finite_mean((right - left)[valid]),
        "mean_cd_x_nm": finite_mean(cd_x),
        "mean_cd_slant_nm": finite_mean(cd_slant),
        "ler_left_nm": three_sigma(left_nm, ddof=ddof),
        "ler_right_nm": three_sigma(right_nm, ddof=ddof),
        "lwr_x_nm": three_sigma(cd_x, ddof=ddof),
        "lwr_slant_nm": three_sigma(cd_slant, ddof=ddof),
        "centerline_3sigma_nm": three_sigma(center_nm, ddof=ddof),
        "rho_left_right": rho,
        "sigma_left_nm": sigma_left,
        "sigma_right_nm": sigma_right,
        "cov_left_right_nm2": covariance,
        "lwr_from_variance_components_nm": lwr_from_components,
        "valid_count": int(valid.sum()),
        "valid_fraction": float(valid.sum() / max(1, left.size)),
    }


def choose_aggregate_space_indices(
    space_rows: list[dict[str, Any]], min_number: int
) -> list[int]:
    stable = [row for row in space_rows if bool(row.get("stable", False))]
    selected = (
        stable
        if len(stable) >= min_number
        else [
            row
            for row in space_rows
            if np.isfinite(float(row.get("mean_cd_nm", math.nan)))
            and np.isfinite(float(row.get("lwr_nm", math.nan)))
        ]
    )
    return [int(row["space_index"]) for row in selected]


EPS = 1e-12


GEOMETRIES = {
    "x_axis": {
        "folder": "x_axis",
        "title": "X-axis distance",
        "metrics": {
            "CD": {
                "machine_col": "machine_CD_nm",
                "calc_col": "ACD_x_nm",
                "signed_pct_col": "CD_x_signed_diff_pct",
                "abs_pct_col": "CD_x_abs_diff_pct",
                "unit": "nm",
                "note": "Horizontal same-Y width: right_x - left_x.",
            },
            "LER_left": {
                "machine_col": "machine_LER_left_nm",
                "calc_col": "LER_left_nm",
                "signed_pct_col": "LER_left_signed_diff_pct",
                "abs_pct_col": "LER_left_abs_diff_pct",
                "unit": "nm",
                "note": "Left-edge x-position roughness.",
            },
            "LWR": {
                "machine_col": "machine_LWR_nm",
                "calc_col": "LWR_x_nm",
                "signed_pct_col": "LWR_x_signed_diff_pct",
                "abs_pct_col": "LWR_x_abs_diff_pct",
                "unit": "nm",
                "note": "3sigma of grouped local CD using x-axis width.",
            },
        },
    },
    "trench_normal": {
        "folder": "slant",
        "title": "Trench-normal / slant",
        "metrics": {
            "CD": {
                "machine_col": "machine_CD_nm",
                "calc_col": "ACD_slant_nm",
                "signed_pct_col": "CD_slant_signed_diff_pct",
                "abs_pct_col": "CD_slant_abs_diff_pct",
                "unit": "nm",
                "note": "Width projected onto the fitted trench normal.",
            },
            "LER_left": {
                "machine_col": "machine_LER_left_nm",
                "calc_col": "LER_left_nm",
                "signed_pct_col": "LER_left_signed_diff_pct",
                "abs_pct_col": "LER_left_abs_diff_pct",
                "unit": "nm",
                "note": (
                    "Same V2 left-edge x roughness as x-axis mode. "
                    "V9 does not invent a new projected LER definition."
                ),
            },
            "LWR": {
                "machine_col": "machine_LWR_nm",
                "calc_col": "LWR_slant_nm",
                "signed_pct_col": "LWR_slant_signed_diff_pct",
                "abs_pct_col": "LWR_slant_abs_diff_pct",
                "unit": "nm",
                "note": "3sigma of grouped local CD after trench-normal width projection.",
            },
        },
    },
}


def _to_numeric_array(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)


def _safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3 or np.std(x) <= EPS or np.std(y) <= EPS:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def metric_statistics(
    df: pd.DataFrame,
    geometry: str,
    metric: str,
) -> dict[str, Any]:
    spec = GEOMETRIES[geometry]["metrics"][metric]
    machine = _to_numeric_array(df[spec["machine_col"]])
    calc = _to_numeric_array(df[spec["calc_col"]])

    valid = np.isfinite(machine) & np.isfinite(calc) & (np.abs(machine) > EPS)
    m = machine[valid]
    c = calc[valid]
    if m.size == 0:
        raise RuntimeError(f"No valid values for geometry={geometry}, metric={metric}")

    diff = c - m
    signed_pct = diff / np.abs(m) * 100.0
    abs_pct = np.abs(signed_pct)

    if m.size > 1:
        diff_sd = float(np.std(diff, ddof=1))
        pct_sd = float(np.std(signed_pct, ddof=1))
    else:
        # One statistics pair cannot estimate spread or agreement limits.
        diff_sd = math.nan
        pct_sd = math.nan

    ba_bias = float(np.mean(diff))
    ba_low = ba_bias - 1.96 * diff_sd
    ba_high = ba_bias + 1.96 * diff_sd

    return {
        "geometry": geometry,
        "geometry_label": GEOMETRIES[geometry]["title"],
        "metric": metric,
        "n": int(m.size),
        "machine_mean_nm": float(np.mean(m)),
        "python_mean_nm": float(np.mean(c)),
        "mean_difference_nm": float(np.mean(diff)),
        "mae_nm": float(np.mean(np.abs(diff))),
        "rmse_nm": float(np.sqrt(np.mean(diff**2))),
        "mean_signed_diff_pct": float(np.mean(signed_pct)),
        "mape_pct": float(np.mean(abs_pct)),
        "median_abs_diff_pct": float(np.median(abs_pct)),
        "std_signed_diff_pct": pct_sd,
        "pearson_r": _safe_pearson(m, c),
        "bland_altman_bias_nm": ba_bias,
        "bland_altman_lower_95_nm": ba_low,
        "bland_altman_upper_95_nm": ba_high,
        "metric_note": spec["note"],
    }


def _plot_title(metric: str, geometry: str, suffix: str) -> str:
    return f"{metric} — {GEOMETRIES[geometry]['title']} — {suffix}"


def _save_machine_vs_python(
    df: pd.DataFrame,
    geometry: str,
    metric: str,
    out_path: Path,
    dpi: int,
) -> None:
    spec = GEOMETRIES[geometry]["metrics"][metric]
    machine = _to_numeric_array(df[spec["machine_col"]])
    calc = _to_numeric_array(df[spec["calc_col"]])
    valid = np.isfinite(machine) & np.isfinite(calc)
    m, c = machine[valid], calc[valid]

    fig = plt.figure(figsize=(7.2, 6.0))
    ax = fig.add_axes([0.12, 0.12, 0.82, 0.80])
    ax.scatter(m, c, alpha=0.75)

    if m.size:
        lo = float(min(np.min(m), np.min(c)))
        hi = float(max(np.max(m), np.max(c)))
        pad = max((hi - lo) * 0.06, 1e-6)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--")

    stats = metric_statistics(df, geometry, metric)
    ax.set_xlabel(f"Machine {metric} (nm)")
    ax.set_ylabel(f"Python {metric} (nm)")
    ax.set_title(_plot_title(metric, geometry, "Machine vs Python"))
    ax.grid(True, alpha=0.25)
    ax.text(
        0.02,
        0.98,
        (
            f"n={stats['n']}\n"
            f"MAPE={stats['mape_pct']:.2f}%\n"
            f"Bias={stats['mean_signed_diff_pct']:+.2f}%\n"
            f"r={stats['pearson_r']:.3f}"
        ),
        transform=ax.transAxes,
        va="top",
    )
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _save_error_by_image(
    df: pd.DataFrame,
    geometry: str,
    metric: str,
    out_path: Path,
    dpi: int,
) -> None:
    spec = GEOMETRIES[geometry]["metrics"][metric]
    machine = _to_numeric_array(df[spec["machine_col"]])
    calc = _to_numeric_array(df[spec["calc_col"]])
    valid = np.isfinite(machine) & np.isfinite(calc) & (np.abs(machine) > EPS)
    rows = df.loc[valid].copy()
    m = machine[valid]
    c = calc[valid]
    err = (c - m) / np.abs(m) * 100.0

    if "machine_excel_row" in rows.columns:
        x = pd.to_numeric(rows["machine_excel_row"], errors="coerce").to_numpy(
            dtype=float
        )
        if not np.all(np.isfinite(x)):
            x = np.arange(1, len(rows) + 1, dtype=float)
            xlabel = "Compared image index"
        else:
            xlabel = "Data.xlsx row"
    else:
        x = np.arange(1, len(rows) + 1, dtype=float)
        xlabel = "Compared image index"

    fig = plt.figure(figsize=(9.2, 5.2))
    ax = fig.add_axes([0.10, 0.16, 0.86, 0.76])
    ax.scatter(x, err, alpha=0.78)
    ax.axhline(0.0, linestyle="--")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Signed difference vs Machine (%)")
    ax.set_title(_plot_title(metric, geometry, "Per-image percent error"))
    ax.grid(True, alpha=0.25)

    # Label the five largest absolute percentage errors.
    if len(err):
        worst = np.argsort(np.abs(err))[-min(5, len(err)) :]
        for idx in worst:
            image_name = str(rows.iloc[idx].get("image_name", int(idx) + 1))
            label = Path(image_name).stem
            ax.annotate(
                label,
                (x[idx], err[idx]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
            )

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _save_error_histogram(
    df: pd.DataFrame,
    geometry: str,
    metric: str,
    out_path: Path,
    dpi: int,
) -> None:
    spec = GEOMETRIES[geometry]["metrics"][metric]
    machine = _to_numeric_array(df[spec["machine_col"]])
    calc = _to_numeric_array(df[spec["calc_col"]])
    valid = np.isfinite(machine) & np.isfinite(calc) & (np.abs(machine) > EPS)
    err = (calc[valid] - machine[valid]) / np.abs(machine[valid]) * 100.0

    fig = plt.figure(figsize=(7.6, 5.2))
    ax = fig.add_axes([0.12, 0.15, 0.83, 0.76])
    ax.hist(err, bins="auto", alpha=0.8)
    ax.axvline(0.0, linestyle="--")
    if err.size:
        ax.axvline(float(np.mean(err)), linestyle=":")
    ax.set_xlabel("Signed difference vs Machine (%)")
    ax.set_ylabel("Image count")
    ax.set_title(_plot_title(metric, geometry, "Percent-error distribution"))
    ax.grid(True, axis="y", alpha=0.25)

    stats = metric_statistics(df, geometry, metric)
    ax.text(
        0.02,
        0.97,
        (
            f"Bias={stats['mean_signed_diff_pct']:+.2f}%\n"
            f"MAPE={stats['mape_pct']:.2f}%\n"
            f"SD={stats['std_signed_diff_pct']:.2f}%"
        ),
        transform=ax.transAxes,
        va="top",
    )
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _save_bland_altman(
    df: pd.DataFrame,
    geometry: str,
    metric: str,
    out_path: Path,
    dpi: int,
) -> None:
    spec = GEOMETRIES[geometry]["metrics"][metric]
    machine = _to_numeric_array(df[spec["machine_col"]])
    calc = _to_numeric_array(df[spec["calc_col"]])
    valid = np.isfinite(machine) & np.isfinite(calc)
    m, c = machine[valid], calc[valid]

    means = 0.5 * (m + c)
    diff = c - m
    bias = float(np.mean(diff)) if diff.size else math.nan
    sd = float(np.std(diff, ddof=1)) if diff.size > 1 else 0.0
    lower = bias - 1.96 * sd
    upper = bias + 1.96 * sd

    fig = plt.figure(figsize=(7.8, 5.4))
    ax = fig.add_axes([0.12, 0.15, 0.83, 0.76])
    ax.scatter(means, diff, alpha=0.75)
    ax.axhline(bias, linestyle="-")
    ax.axhline(lower, linestyle="--")
    ax.axhline(upper, linestyle="--")
    ax.set_xlabel(f"Mean of Machine and Python {metric} (nm)")
    ax.set_ylabel("Python - Machine (nm)")
    ax.set_title(_plot_title(metric, geometry, "Bland–Altman"))
    ax.grid(True, alpha=0.25)
    ax.text(
        0.02,
        0.98,
        (f"Bias={bias:+.4f} nm\n95% limits=[{lower:+.4f}, {upper:+.4f}] nm"),
        transform=ax.transAxes,
        va="top",
    )
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _save_geometry_summary_bars(
    summary: pd.DataFrame,
    geometry: str,
    out_dir: Path,
    dpi: int,
) -> list[dict[str, Any]]:
    part = summary[summary["geometry"] == geometry].copy()
    order = ["CD", "LER_left", "LWR"]
    part["metric_order"] = part["metric"].map({name: i for i, name in enumerate(order)})
    part = part.sort_values("metric_order")

    records = []

    fig = plt.figure(figsize=(7.2, 5.0))
    ax = fig.add_axes([0.13, 0.16, 0.82, 0.75])
    ax.bar(part["metric"], part["mape_pct"])
    ax.set_ylabel("MAPE vs Machine (%)")
    ax.set_title(f"{GEOMETRIES[geometry]['title']} — MAPE summary")
    ax.grid(True, axis="y", alpha=0.25)
    path = out_dir / "summary_MAPE.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    records.append(
        {
            "geometry": geometry,
            "metric": "ALL",
            "plot_type": "summary_MAPE",
            "relative_path": str(path),
            "description": "MAPE statistics for CD, LER-left and LWR.",
        }
    )

    fig = plt.figure(figsize=(7.2, 5.0))
    ax = fig.add_axes([0.13, 0.16, 0.82, 0.75])
    ax.bar(part["metric"], part["mean_signed_diff_pct"])
    ax.axhline(0.0, linestyle="--")
    ax.set_ylabel("Mean signed difference vs Machine (%)")
    ax.set_title(f"{GEOMETRIES[geometry]['title']} — Bias summary")
    ax.grid(True, axis="y", alpha=0.25)
    path = out_dir / "summary_Bias.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    records.append(
        {
            "geometry": geometry,
            "metric": "ALL",
            "plot_type": "summary_Bias",
            "relative_path": str(path),
            "description": "Mean signed percent bias for CD, LER-left and LWR.",
        }
    )
    return records
