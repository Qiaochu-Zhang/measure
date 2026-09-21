"""V1_17 raw-edge PSD: separate engines, structures and contiguous segments.

No group averaging, Welch averaging, cross-structure or cross-image averaging.
An FFT is never taken across a missing/background gap. Mean/linear detrending
is a declared trend subtraction, not spatial smoothing or block averaging.
"""

from pathlib import Path
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


class PSDBatch:
    def __init__(self, args, output_dir):
        self.args = args
        self.output = Path(output_dir) / "PSD"
        self.objects, self.curves, self.audit, self.manifest = [], [], [], []

    def add_measurement(self, measurement, key, engine, params):
        if engine not in ("V10", "V13"):
            raise ValueError("PSD engine must be V10 or V13")
        self.manifest.append(
            dict(
                image_key=key,
                engine=engine,
                detected_structures=len(measurement.space_rows),
                average_range_px=params.average_range_px,
                smoothing_pixel=params.smoothing_pixel,
                edge_source=self.args.psd_edge_source,
                group_size=1,
                spectrum_averaging=False,
                measurement_valid=measurement.image_row.get("valid", False),
                warning=measurement.image_row.get("warning", ""),
            )
        )
        for row in measurement.sample_rows:
            self.audit.append(dict(image_key=key, engine=engine, **row))
        for ann in measurement.annotation_payload.get("spaces", []):
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
            axes = (
                ("normal", "unrotated")
                if self.args.psd_axis == "both"
                else (self.args.psd_axis,)
            )
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
            {
                **meta,
                "frequency_per_nm": float(frequency),
                "wavelength_nm": 1 / float(frequency) if frequency > 0 else math.nan,
                "psd_nm3": float(value),
            }
            for meta, f, power in self.curves
            for frequency, value in zip(f, power)
        ]
        curves = pd.DataFrame(records)
        if curves.empty:
            curves = pd.DataFrame(
                columns=["image_key", "engine", "signal", "frequency_per_nm", "psd_nm3"]
            )
        audit = pd.DataFrame(self.audit)
        manifest = pd.DataFrame(self.manifest)
        atomic_csv(summary, self.output / "per_structure_psd_summary.csv")
        atomic_csv(curves, self.output / "per_structure_psd_curves.csv")
        atomic_csv(audit, self.output / "edge_coordinates.csv")
        atomic_csv(manifest, self.output / "measurement_manifest.csv")
        for engine in ("V10", "V13"):
            subset = summary[summary.engine == engine]
            atomic_csv(subset, self.output / f"{engine}_psd_summary.csv")
            atomic_csv(
                curves[curves.engine == engine],
                self.output / f"{engine}_psd_curves.csv",
            )
            with pd.ExcelWriter(
                self.output / f"PSD_{engine}.xlsx", engine="openpyxl"
            ) as writer:
                subset.to_excel(
                    writer, sheet_name="per_structure_segments", index=False
                )
                manifest[manifest.engine == engine].to_excel(
                    writer, sheet_name="manifest", index=False
                ) if not manifest.empty else pd.DataFrame().to_excel(
                    writer, sheet_name="manifest", index=False
                )
                for ws in writer.book.worksheets:
                    ws.freeze_panes = "D2"
                    ws.auto_filter.ref = ws.dimensions
            if not self.args.skip_statistics_plots:
                for key in dict.fromkeys(
                    m["image_key"] for m, _, _ in self.curves if m["engine"] == engine
                ):
                    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
                    for ax, kind in zip(
                        axes.flat, ("left", "right", "width", "center")
                    ):
                        for meta, f, power in self.curves:
                            if (
                                meta["engine"] != engine
                                or meta["image_key"] != key
                                or meta["signal"] != kind
                            ):
                                continue
                            positive = (f > 0) & (power > 0)
                            ax.loglog(
                                f[positive],
                                power[positive],
                                alpha=0.7,
                                label=f"{meta['axis']} S{meta['space_index']} seg{meta['segment_index']}",
                            )
                        ax.set(
                            title=kind,
                            xlabel="Spatial frequency (1/nm)",
                            ylabel="PSD (nm^3)",
                        )
                        if ax.lines:
                            ax.legend(fontsize=6)
                    fig.suptitle(
                        f"{key} / {engine} / separate raw segments, NO averaging"
                    )
                    fig.tight_layout()
                    path = self.output / engine / Path(key).with_suffix(".psd.png")
                    path.parent.mkdir(parents=True, exist_ok=True)
                    try:
                        fig.savefig(path, dpi=self.args.plot_dpi)
                    finally:
                        plt.close(fig)
        settings = {k: v for k, v in vars(self.args).items() if k.startswith("psd_")}
        settings.update(
            group_size=1,
            spectrum_averaging=False,
            method="periodogram",
            segmentation="each contiguous valid segment separately; no concatenation",
            normal_frequency_coordinate="nominal distance along fitted mean centerline: y/cos(theta)",
            units={"frequency": "1/nm", "density": "nm^3"},
            noise_debiased=False,
        )
        (self.output / "PSD_settings.json").write_text(
            json.dumps(settings, ensure_ascii=False, indent=2) + "\n"
        )
        return self.output
