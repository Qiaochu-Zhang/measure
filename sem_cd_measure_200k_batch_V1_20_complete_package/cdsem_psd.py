"""V1_20: rotated, within-image ensemble PSD for LER and LWR.

Compute each contiguous edge segment before averaging *powers*, never edge
coordinates. V10 and V13, and different images, remain separate populations.
Native periodograms remain available for support and normalization audits.
"""

from pathlib import Path
from cdsem_regions import output_image_path
import json
import math
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from scipy import signal
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def runs(mask):
    edges = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def segment_spectra(values, spacing_nm, args, hard_gaps=None):
    """Return every usable segment separately, preserving physical spacing."""
    x = np.asarray(values, float).copy()
    filled = np.zeros(len(x), bool)
    hard_gaps = (
        np.zeros(len(x), bool) if hard_gaps is None else np.asarray(hard_gaps, bool)
    )
    if not np.isfinite(spacing_nm) or spacing_nm <= 0:
        return []
    if args.psd_gap_mode == "interpolate":
        for a, b in runs(~np.isfinite(x)):
            if (
                a > 0
                and b < len(x)
                and b - a <= args.psd_max_gap
                and not hard_gaps[a:b].any()
            ):
                x[a:b] = np.linspace(x[a - 1], x[b], b - a + 2)[1:-1]
                filled[a:b] = True
    spectra = []
    for a, b in runs(np.isfinite(x)):
        if b - a < args.psd_min_segment:
            continue
        detrend = False if args.psd_detrend == "none" else args.psd_detrend
        f, power = signal.periodogram(
            x[a:b],
            fs=1 / spacing_nm,
            window=args.psd_window,
            detrend=detrend,
            scaling="density",
            return_onesided=True,
        )
        spectra.append(
            (
                f,
                power,
                dict(
                    segment_index=len(spectra) + 1,
                    start_slot=int(a),
                    stop_slot_exclusive=int(b),
                    used_slots=int(b - a),
                    interpolated_slots=int(filled[a:b].sum()),
                    spacing_nm=float(spacing_nm),
                    span_nm=float((b - a - 1) * spacing_nm),
                    frequency_step_per_nm=float(f[1] - f[0]),
                ),
            )
        )
    return spectra


def scalar_metrics(f, power, split_wavelength):
    df = float(f[1] - f[0])
    variance = float(power.sum() * df)
    low = float(power[f <= 1 / split_wavelength].sum() * df)
    return dict(
        variance_psd_nm2=variance,
        roughness_psd_3sigma_nm=3 * math.sqrt(max(0, variance)),
        low_band_variance_nm2=low,
        high_band_variance_nm2=variance - low,
    )


def atomic_csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(
            mode="w",
            encoding="utf-8-sig",
            newline="",
            dir=path.parent,
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            frame.to_csv(stream, index=False)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


IMAGE_SUMMARY_COLUMNS = [
    "image", "status", "method", "image_key", "engine", "axis", "signal", "reason",
    "contributing_structures", "contributing_edges", "segment_count",
    "spectrum_averaging", "group_size", "frequency_bin_count", "positive_psd_bin_count",
    "native_min_frequency_step_per_nm", "native_max_frequency_step_per_nm",
    "min_frequency_per_nm", "max_frequency_per_nm", "weight",
]
IMAGE_CURVE_COLUMNS = [
    "image", "image_key", "engine", "axis", "signal", "frequency_per_nm", "wavelength_nm",
    "psd_nm3", "contributing_segments", "contributing_structures", "total_weight_nm",
    "native_resolution_best_per_nm", "native_resolution_worst_per_nm", "spectrum_averaging",
]
BAND_COLUMNS = [
    "image", "image_key", "engine", "axis", "band", "status",
    "lower_frequency_per_nm", "upper_frequency_per_nm", "supported_segments",
    "supported_structures", "native_frequency_bin_count", "min_observed_frequency_per_nm",
    "max_observed_frequency_per_nm", "LER_variance_nm2", "LWR_variance_nm2",
    "center_variance_nm2", "LER_3sigma_nm", "LWR_3sigma_nm", "center_3sigma_nm",
    "LER_to_LWR_power_ratio", "edge_correlation_from_psd", "integration",
]


def paired_segments(curves):
    """Pair left/right powers on the same physical structure and valid segment."""
    grouped = {}
    for meta, f, power in curves:
        key = (meta["space_index"], meta["segment_index"])
        entry = grouped.setdefault(key, dict(meta=meta, f=f, powers={}))
        if not np.array_equal(entry["f"], f):
            raise ValueError("Paired edge spectra must have identical native frequencies")
        entry["powers"][meta["signal"]] = power
    segments = []
    for entry in grouped.values():
        powers = entry["powers"]
        if not {"left", "right", "width", "center"}.issubset(powers):
            raise ValueError("Missing paired left/right/width/center spectrum")
        segments.append(dict(
            meta=entry["meta"], f=entry["f"],
            LER=(powers["left"] + powers["right"]) / 2,
            LWR=powers["width"], center=powers["center"],
            weight=entry["meta"]["used_slots"] * entry["meta"]["spacing_nm"],
        ))
    return segments


def ensemble_spectrum(segments):
    """Align densities without extrapolating short segments into low frequencies.

    Interpolation aligns different rotated spacings/lengths; it does not improve
    their independent resolution. Band integrals use the original FFT bins.
    """
    grid = min(segments, key=lambda s: s["f"][1])["f"][1:].copy() if segments else np.array([])
    n = len(grid)
    sums = {name: np.zeros(n) for name in ("LER", "LWR", "center")}
    weights, counts = np.zeros(n), np.zeros(n, dtype=int)
    best, worst = np.full(n, np.inf), np.zeros(n)
    structures = [set() for _ in grid]
    for segment in segments:
        f, weight = segment["f"], segment["weight"]
        supported = (grid >= f[1] * (1 - 1e-12)) & (grid <= f[-1] * (1 + 1e-12))
        for name in sums:
            sums[name][supported] += weight * np.interp(grid[supported], f, segment[name])
        weights[supported] += weight
        counts[supported] += 1
        best[supported] = np.minimum(best[supported], f[1])
        worst[supported] = np.maximum(worst[supported], f[1])
        for i in np.flatnonzero(supported):
            structures[i].add(segment["meta"]["space_index"])
    return dict(
        f=grid, weight=weights, segments=counts,
        structures=np.array([len(ids) for ids in structures]),
        best_df=np.where(weights > 0, best, np.nan),
        worst_df=np.where(weights > 0, worst, np.nan),
        **{name: np.divide(values, weights, out=np.full(n, np.nan), where=weights > 0)
           for name, values in sums.items()},
    )


def native_band_metrics(segments, args):
    """Integrate observed positive native FFT bins; missing support is NaN."""
    low = 1 / args.psd_split_wavelength
    high = 1 / args.psd_high_split_wavelength
    rows = []
    for band, lower, upper in (("LF", 0.0, low), ("MF", low, high), ("HF", high, np.inf)):
        supported, frequencies = [], []
        for segment in segments:
            f = segment["f"]
            mask = (f > lower) & (f <= upper)
            if mask.any():
                supported.append((segment, mask))
                frequencies.extend(f[mask])
        weight = sum(s["weight"] for s, _ in supported)
        values = {
            name: sum(s["weight"] * float(s[name][mask].sum()) * s["f"][1]
                      for s, mask in supported) / weight if weight else math.nan
            for name in ("LER", "LWR", "center")
        }
        ler, lwr = values["LER"], values["LWR"]
        rows.append(dict(
            band=band, status="OK" if supported else "NO_FREQUENCY_SUPPORT",
            lower_frequency_per_nm=lower,
            upper_frequency_per_nm=upper if np.isfinite(upper) else math.nan,
            supported_segments=len(supported),
            supported_structures=len({s["meta"]["space_index"] for s, _ in supported}),
            native_frequency_bin_count=len(frequencies),
            min_observed_frequency_per_nm=min(frequencies, default=math.nan),
            max_observed_frequency_per_nm=max(frequencies, default=math.nan),
            **{f"{name}_variance_nm2": value for name, value in values.items()},
            **{f"{name}_3sigma_nm": 3 * math.sqrt(max(0, value)) if np.isfinite(value) else math.nan
               for name, value in values.items()},
            LER_to_LWR_power_ratio=ler / lwr if lwr > 0 else math.nan,
            edge_correlation_from_psd=1 - lwr / (2 * ler) if ler > 0 else math.nan,
            integration="native positive FFT bins; length-weighted supported segments; no DC",
        ))
    return rows


def plot_image_psd(ensemble, key, engine, args, segment_count):
    """Exactly one two-curve panel per image/engine, including empty placeholders."""
    fig, ax = plt.subplots(figsize=(10, 7))
    frequency = ensemble["f"]
    for kind, color in (("LER", "#d62728"), ("LWR", "#1f77b4")):
        values = ensemble[kind]
        ax.loglog(frequency, np.where(values > 0, values, np.nan),
                  color=color, linewidth=1.2, marker=".", markersize=2.5, label=kind)
    low, high = 1 / args.psd_split_wavelength, 1 / args.psd_high_split_wavelength
    xmin = min(float(frequency[0]), low / 3) if len(frequency) else low / 10
    xmax = max(float(frequency[-1]), high * 3) if len(frequency) else high * 10
    ax.set_xlim(xmin, xmax)
    for boundary, wavelength, label in (
        (low, args.psd_split_wavelength, "LF / MF"),
        (high, args.psd_high_split_wavelength, "MF / HF"),
    ):
        ax.axvline(boundary, color="#666666", linestyle="--", linewidth=1)
        ax.text(boundary, 0.88, f" {label}\n {boundary:.4g} 1/nm\n ({wavelength:g} nm)",
                transform=ax.get_xaxis_transform(), va="top", ha="left", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75))
    for a, b, label in ((xmin, low, "LF"), (low, high, "MF"), (high, xmax, "HF")):
        ax.text(math.sqrt(a * b), 0.96, label, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=10, fontweight="bold")
    ax.set(xlabel="Spatial frequency along rotated centerline (1/nm)",
           ylabel="PSD (nm$^3$)")
    ax.set_title(f"{Path(key).name} / {engine} / rotated ensemble PSD", pad=12)
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(loc="lower left")
    if not segment_count:
        ax.text(0.5, 0.5, "No qualified contiguous segment supports PSD.\nNo zero spectrum is inferred.",
                transform=ax.transAxes, ha="center", va="center")
    elif not any(np.any(ensemble[name] > 0) for name in ("LER", "LWR")):
        ax.text(0.5, 0.5, "All measured PSD bins are zero; log axis cannot display zero.",
                transform=ax.transAxes, ha="center", va="center")
    count = int(max(ensemble["structures"], default=0))
    fig.text(0.1, 0.04,
             f"{count} structures; {segment_count} paired segments; {len(frequency)} frequency bins. "
             "No coordinate concatenation.\n"
             "LER = mean left/right edge power; LWR = width power. Dashed boundaries are configurable.",
             fontsize=8)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    return fig


class PSDBatch:
    def __init__(self, args, output_dir):
        if args.psd_axis != "normal":
            raise ValueError("V1_20 PSD supports rotated/normal coordinates only")
        self.args = args
        self.images = []
        self.output = Path(output_dir) / "PSD"
        self.objects, self.curves, self.audit, self.manifest = [], [], [], []

    def register_image(self, key):
        if key not in self.images:
            self.images.append(key)

    def add_measurement(self, measurement, key, engine, params):
        if engine not in ("V10", "V13"):
            raise ValueError("PSD engine must be V10 or V13")
        self.register_image(key)
        annotations = measurement.annotation_payload.get(
            "psd_all_spaces", measurement.annotation_payload.get("spaces", [])
        )
        annotations = [ann for ann in annotations if ann.get("stable", False)]
        self.manifest.append(
            dict(
                image_key=key,
                engine=engine,
                detected_structures=len(annotations),
                main_selected_structures=len(measurement.space_rows),
                population="all_qualified_structures_in_ROI",
                axis="normal",
                average_range_px=params.average_range_px,
                smoothing_pixel=params.smoothing_pixel,
                edge_source=self.args.psd_edge_source,
                group_size=1,
                spectrum_averaging=True,
                measurement_valid=measurement.image_row.get("valid", False),
                warning=measurement.image_row.get("warning", ""),
            )
        )
        for row in measurement.annotation_payload.get(
            "psd_all_sample_rows", measurement.sample_rows
        ):
            self.audit.append(dict(image_key=key, engine=engine, **row))
        for ann in annotations:
            y = np.asarray(ann["y_values"], float)
            left = np.asarray(ann["left_edges"], float).copy()
            right = np.asarray(ann["right_edges"], float).copy()
            if len(y) < 2 or not np.allclose(
                np.diff(y), y[1] - y[0], atol=1e-9, rtol=1e-7
            ):
                raise ValueError(
                    "PSD requires uniformly spaced original sampling slots"
                )
            synthetic = np.asarray(
                ann.get("continuity_synthetic", np.zeros(len(y))), bool
            )
            background = np.asarray(
                ann.get("background_excluded", np.zeros(len(y))), bool
            )
            valid = (
                np.isfinite(left) & np.isfinite(right) & (right > left) & ~background
            )
            if not self.args.psd_include_synthetic:
                valid &= ~synthetic
            left[~valid] = np.nan
            right[~valid] = np.nan
            theta = np.deg2rad(
                ann.get("centerline_fit", {}).get("angle_from_vertical_deg", 0)
            )
            axes = ("normal",)
            for axis in axes:
                c, s = (
                    (float(np.cos(theta)), float(np.sin(theta)))
                    if axis == "normal"
                    else (1.0, 0.0)
                )
                if abs(c) < 1e-6:
                    continue
                l = (left * c - y * s) * params.pixel_size_nm
                r = (right * c - y * s) * params.pixel_size_nm
                spacing = (y[1] - y[0]) * params.pixel_size_nm / abs(c)
                for kind, values in dict(
                    left=l, right=r, width=r - l, center=(l + r) / 2
                ).items():
                    meta = dict(
                        image_key=key,
                        engine=engine,
                        method=engine,
                        axis=axis,
                        signal=kind,
                        space_index=ann["space_index"],
                        space_role=ann["candidate"].role,
                        group_size=1,
                        spectrum_averaging=False,
                        average_range_px=params.average_range_px,
                        smoothing_pixel=params.smoothing_pixel,
                        input_slots=len(y),
                        input_valid_count=int(valid.sum()),
                        excluded_background_count=int(background.sum()),
                        synthetic_count=int(synthetic.sum()),
                        measurement_stable=bool(ann["stable"]),
                        angle_from_vertical_deg=float(np.rad2deg(theta)),
                    )
                    spectra = segment_spectra(
                        values, spacing, self.args, hard_gaps=background
                    )
                    if not spectra:
                        self.objects.append(
                            {
                                **meta,
                                "status": "SKIPPED",
                                "reason": "no_long_enough_contiguous_segment",
                                "used_slots": 0,
                            }
                        )
                    for f, power, info in spectra:
                        combined = {**meta, **info}
                        self.objects.append(
                            {
                                **combined,
                                "status": "OK",
                                "reason": "",
                                **scalar_metrics(
                                    f, power, self.args.psd_split_wavelength
                                ),
                            }
                        )
                        self.curves.append((combined, f, power))

    def finalize(self):
        self.output.mkdir(parents=True, exist_ok=True)
        summary = pd.DataFrame(self.objects)
        if summary.empty:
            summary = pd.DataFrame(
                columns=["image_key", "engine", "status", "signal", "axis"]
            )
        summary["image"] = summary["image_key"].map(lambda key: Path(str(key)).name)
        if "method" not in summary:
            summary["method"] = summary["engine"]
        leading = ["image", "status", "method", "image_key", "engine"]
        summary = summary[leading + [c for c in summary if c not in leading]]
        records = [
            dict(meta, frequency_per_nm=float(frequency),
                 wavelength_nm=1 / float(frequency) if frequency > 0 else math.nan,
                 psd_nm3=float(value))
            for meta, f, power in self.curves for frequency, value in zip(f, power)
        ]
        curves = pd.DataFrame(records) if records else pd.DataFrame(columns=[
            "image_key", "engine", "signal", "frequency_per_nm", "psd_nm3"
        ])
        image_rows, image_curves, band_rows = [], [], []
        for key in self.images:
            for engine in ("V10", "V13"):
                raw = [c for c in self.curves
                       if c[0]["image_key"] == key and c[0]["engine"] == engine]
                segments = paired_segments(raw)
                ensemble = ensemble_spectrum(segments)
                for kind in ("LER", "LWR"):
                    finite = np.isfinite(ensemble[kind])
                    image_rows.append(dict(
                        image=Path(key).name, image_key=key, engine=engine,
                        method=engine, axis="normal", signal=kind,
                        status="OK" if finite.any() else "SKIPPED",
                        reason="" if finite.any() else "no_long_enough_qualified_segment",
                        contributing_structures=len({s["meta"]["space_index"] for s in segments}),
                        contributing_edges=(2 if kind == "LER" else 1) * len(
                            {s["meta"]["space_index"] for s in segments}),
                        segment_count=len(segments), spectrum_averaging=True, group_size=1,
                        frequency_bin_count=int(finite.sum()),
                        positive_psd_bin_count=int((finite & (ensemble[kind] > 0)).sum()),
                        native_min_frequency_step_per_nm=min(
                            (s["f"][1] for s in segments), default=math.nan),
                        native_max_frequency_step_per_nm=max(
                            (s["f"][1] for s in segments), default=math.nan),
                        min_frequency_per_nm=float(ensemble["f"][finite][0]) if finite.any() else math.nan,
                        max_frequency_per_nm=float(ensemble["f"][finite][-1]) if finite.any() else math.nan,
                        weight="N_slots * rotated_spacing_nm; per-frequency support only",
                    ))
                    for i, frequency in enumerate(ensemble["f"]):
                        image_curves.append(dict(
                            image=Path(key).name, image_key=key, engine=engine,
                            axis="normal", signal=kind, frequency_per_nm=float(frequency),
                            wavelength_nm=float(1 / frequency), psd_nm3=float(ensemble[kind][i]),
                            contributing_segments=int(ensemble["segments"][i]),
                            contributing_structures=int(ensemble["structures"][i]),
                            total_weight_nm=float(ensemble["weight"][i]),
                            native_resolution_best_per_nm=float(ensemble["best_df"][i]),
                            native_resolution_worst_per_nm=float(ensemble["worst_df"][i]),
                            spectrum_averaging=True,
                        ))
                for row in native_band_metrics(segments, self.args):
                    band_rows.append(dict(image=Path(key).name, image_key=key,
                                          engine=engine, axis="normal", **row))
                if not self.args.skip_statistics_plots:
                    fig = plot_image_psd(ensemble, key, engine, self.args, len(segments))
                    path = self.output / engine / output_image_path(key, ".psd.png")
                    path.parent.mkdir(parents=True, exist_ok=True)
                    try:
                        fig.savefig(path, dpi=self.args.plot_dpi)
                    finally:
                        plt.close(fig)
        image_summary = pd.DataFrame(image_rows, columns=IMAGE_SUMMARY_COLUMNS)
        aggregate = pd.DataFrame(image_curves, columns=IMAGE_CURVE_COLUMNS)
        bands = pd.DataFrame(band_rows, columns=BAND_COLUMNS)
        audit, manifest = pd.DataFrame(self.audit), pd.DataFrame(self.manifest)
        for frame, name in (
            (summary, "per_structure_psd_summary.csv"),
            (curves, "per_structure_psd_curves.csv"),
            (image_summary, "image_psd_summary.csv"),
            (aggregate, "image_psd_curves.csv"),
            (bands, "image_psd_band_summary.csv"),
            (audit, "edge_coordinates.csv"), (manifest, "measurement_manifest.csv"),
        ):
            atomic_csv(frame, self.output / name)
        for engine in ("V10", "V13"):
            subset = image_summary[image_summary.engine == engine]
            atomic_csv(subset, self.output / f"{engine}_psd_summary.csv")
            atomic_csv(aggregate[aggregate.engine == engine], self.output / f"{engine}_psd_curves.csv")
            with pd.ExcelWriter(self.output / f"PSD_{engine}.xlsx", engine="openpyxl") as writer:
                for frame, sheet in (
                    (subset, "image_summary"),
                    (bands[bands.engine == engine], "frequency_bands"),
                    (summary[summary.engine == engine], "per_structure_segments"),
                    (manifest[manifest.engine == engine] if not manifest.empty else manifest, "manifest"),
                ):
                    frame.to_excel(writer, sheet_name=sheet, index=False)
                for ws in writer.book.worksheets:
                    ws.freeze_panes = "D2"
                    ws.auto_filter.ref = ws.dimensions
        settings = {k: v for k, v in vars(self.args).items() if k.startswith("psd_")}
        settings.update(
            version="V1_20", axis="normal", group_size=1, spectrum_averaging=True,
            population="all qualified structures within ROI; not main max-number subset",
            method="native segment periodograms then within-image power ensemble",
            LER="(PSD_left + PSD_right) / 2; never PSD of averaged coordinates",
            LWR="PSD(right - left), paired valid slots, no group4",
            weighting="N_slots * rotated_spacing_nm per segment; equal left/right weights",
            frequency_alignment="longest native physical-period grid; linear density interpolation only inside each segment's positive frequency support; no extrapolation; no padding",
            band_integration="native positive FFT bins, sum(PSD * native_df), then length-weighted mean over supported segments; no DC; not integration of interpolated plot",
            segmentation="each contiguous valid segment; never concatenate gaps or trenches",
            normal_frequency_coordinate="nominal distance along fitted mean centerline: y/cos(theta)",
            low_mid_boundary_per_nm=1 / self.args.psd_split_wavelength,
            mid_high_boundary_per_nm=1 / self.args.psd_high_split_wavelength,
            boundary_source="configurable analysis recipe, not universal paper thresholds",
            units={"frequency": "1/nm", "density": "nm^3"}, noise_debiased=False,
            palasantzas_fit=False,
            reference_doi="10.1109/CSTIC.2017.7919761",
        )
        (self.output / "PSD_settings.json").write_text(
            json.dumps(settings, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        return self.output
