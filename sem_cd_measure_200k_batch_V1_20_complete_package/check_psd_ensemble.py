"""V1_20 PSD physics, support, population and two-panel-export regressions."""

from argparse import Namespace
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from cdsem_psd import (
    PSDBatch, ensemble_spectrum, native_band_metrics, paired_segments, plot_image_psd,
    segment_spectra,
)
from check_line_pitch import stripe_image


def wiggling_image():
    """Deterministic synthetic SEM-like trenches; no real instrument data."""
    from scipy.ndimage import gaussian_filter

    y, x = np.indices((512, 1024), dtype=float)
    gray = np.full(y.shape, 210.0)
    center_shift = 0.025 * (y - 256) + 4 * np.sin(2 * np.pi * y / 160)
    half_width = 30 + 0.5 * np.sin(2 * np.pi * y / 32)
    for center in (100, 260, 420, 580, 740, 900):
        gray[np.abs(x - center - center_shift) < half_width] = 40
    gray = gaussian_filter(gray, sigma=0.8)
    gray += np.random.default_rng(120).normal(0, 0.8, gray.shape)
    return np.clip(np.rint(gray), 0, 255).astype(np.uint8)


def analytical_tests(require, app):
    args = app.build_parser().parse_args(["--pattern", "trench"])
    require(args.psd_axis == "normal" and args.psd_detrend == "constant", "rotated low-frequency defaults")
    args.psd_window, args.plot_dpi = "boxcar", 60
    params = Namespace(average_range_px=1, smoothing_pixel=1, pixel_size_nm=1)
    y = np.arange(512.0)
    wave = 2 * np.sin(2 * np.pi * y / 128)

    def measurement(mode="wiggle", tilt=0.0):
        annotations = []
        for index, sign in ((1, 1), (2, -1)):
            left = 100 * index + np.tan(tilt) * y + sign * wave
            right = 100 * index + 60 + np.tan(tilt) * y + sign * wave * (1 if mode == "wiggle" else -1)
            annotations.append(dict(
                y_values=y, left_edges=left, right_edges=right, stable=True,
                candidate=SimpleNamespace(role="candidate"), space_index=index,
                centerline_fit={"angle_from_vertical_deg": np.rad2deg(tilt)},
            ))
        return SimpleNamespace(space_rows=[{}, {}], sample_rows=[], image_row={"valid": True},
                               annotation_payload={"spaces": annotations})

    with TemporaryDirectory(prefix="v120-psd-physics-") as directory:
        batch = PSDBatch(args, directory)
        batch.add_measurement(measurement(), "中文/sample.tif", "V10", params)
        segments = paired_segments(batch.curves)
        ensemble = ensemble_spectrum(segments)
        require(len(ensemble["f"]) == 256, "512 spatial slots produce 256 positive FFT bins")
        np.testing.assert_allclose(ensemble["LER"], segments[0]["LER"][1:], rtol=1e-10, atol=1e-20)
        require(abs(ensemble["f"][np.argmax(ensemble["LER"])] - 1 / 128) < 1e-12, "known low-frequency peak")
        require(np.max(ensemble["LWR"]) < 1e-20, "same-direction wiggling leaves width unchanged")
        low = native_band_metrics(segments, args)[0]
        require(abs(low["LER_variance_nm2"] - 2) < 1e-12, "native LF energy normalization")
        require(abs(low["edge_correlation_from_psd"] - 1) < 1e-12, "correlated edge interpretation")
        require(abs(low["center_variance_nm2"] - 2) < 1e-12, "centerline wiggling retained")

        breathing = PSDBatch(args, directory)
        breathing.add_measurement(measurement("breathing"), "width.png", "V13", params)
        anti = native_band_metrics(paired_segments(breathing.curves), args)[0]
        require(abs(anti["LWR_variance_nm2"] - 8) < 1e-12, "opposite edges generate width power")
        require(abs(anti["edge_correlation_from_psd"] + 1) < 1e-12, "anticorrelation interpretation")
        require(anti["center_variance_nm2"] < 1e-20, "no false center wiggling")

        # A short segment contributes neither invented LF points nor zeros there.
        short = dict(segments[0], f=segments[0]["f"] * 16, weight=32.0,
                     meta=dict(segments[0]["meta"], space_index=99))
        pooled = ensemble_spectrum([segments[0], short])
        require(pooled["structures"][0] == 1 and pooled["structures"][-1] == 2,
                "frequency-dependent contributor support")
        np.testing.assert_allclose(pooled["LER"][:15], segments[0]["LER"][1:16], atol=1e-20)
        no_low = native_band_metrics([short], args)[0]
        require(no_low["status"] == "NO_FREQUENCY_SUPPORT" and np.isnan(no_low["LER_variance_nm2"]),
                "unsupported LF is missing rather than zero")

        tilted = PSDBatch(args, directory)
        tilted.add_measurement(measurement(tilt=np.deg2rad(10)), "tilt.png", "V10", params)
        require({m["axis"] for m, _, _ in tilted.curves} == {"normal"}, "never compute unrotated PSD")
        require(abs(tilted.curves[0][0]["spacing_nm"] - 1 / np.cos(np.deg2rad(10))) < 1e-12,
                "rotated physical sampling interval")
        # V13 must never borrow a usable V10 spectrum, even for empty/failing images.
        batch.register_image("empty.png")
        batch.finalize()
        folder = Path(directory) / "PSD"
        require(len(list(folder.rglob("*.psd.png"))) == 4, "two images produce exactly four PSD files")
        summary = pd.read_csv(folder / "image_psd_summary.csv")
        require(summary[summary.engine == "V13"].status.eq("SKIPPED").all(), "engines remain isolated")
        curves = pd.read_csv(folder / "image_psd_curves.csv")
        require(set(curves.signal) == {"LER", "LWR"} and set(curves.engine) == {"V10"}, "two aggregate signals only")
        for key in ("中文/sample.tif", "empty.png"):
            for engine in ("V10", "V13"):
                require((folder / engine / app.regions.output_image_path(key, ".psd.png")).is_file(), "placeholder path mirrors original extension")
        fig = plot_image_psd(ensemble, "sample.tif", "V10", args, 2)
        try:
            lines = fig.axes[0].lines
            require(len(fig.axes) == 1 and len(lines) == 4, "single panel, two PSD curves and two vertical boundaries")
            require([line.get_label() for line in lines[:2]] == ["LER", "LWR"], "legend and signal order")
            require(lines[0].get_color() != lines[1].get_color(), "distinct LER/LWR colors")
            require(all(line.get_linestyle() == "--" for line in lines[2:]), "dashed band boundaries")
            require(np.allclose([line.get_xdata()[0] for line in lines[2:]], [0.01, 0.05]), "physical boundary locations")
        finally:
            plt.close(fig)
        skipped = PSDBatch(replace_namespace(args, skip_statistics_plots=True), Path(directory) / "skip")
        skipped.register_image("empty.png")
        skipped.finalize()
        require(not list((Path(directory) / "skip").rglob("*.png")), "explicit plot suppression retained")

    gap = wave.copy()
    gap[100:102] = np.nan
    interpolate = replace_namespace(args, psd_gap_mode="interpolate")
    hard = np.zeros(len(gap), bool)
    hard[100:102] = True
    require(len(segment_spectra(gap, 1, interpolate, hard_gaps=hard)) == 2, "never interpolate background")
    for extra in (["--psd-high-split-wavelength", "100"], ["--psd-high-split-wavelength", "0"],
                  ["--psd-high-split-wavelength", "nan"]):
        parsed = app.build_parser().parse_args(["--pattern", "trench", *extra])
        try:
            app.validate_args(parsed, extra)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid band recipe accepted")
    return dict(status="PASS", known_LF_LER_variance_nm2=2.0, known_breathing_LWR_variance_nm2=8.0,
                analytic_positive_frequency_bins=256, engines_isolated=True, unsupported_band_is_null=True,
                exact_plots_per_image=2, no_unrotated_psd=True)


def replace_namespace(args, **kwargs):
    return Namespace(**dict(vars(args), **kwargs))


def integration_tests(temp, execute, require, app):
    root = temp / "v120_ensemble_input"
    app.regions.unicode_imwrite(root / "中文 子目录/sample.tif", stripe_image())
    app.regions.unicode_imwrite(root / "中文 子目录/sample.png", stripe_image())
    captured = []
    original = app.edges.measure_image

    def inspect(*positional, **kwargs):
        measured = original(*positional, **kwargs)
        params = positional[1]
        if params.average_range_px == 1:
            payload = measured.annotation_payload
            qualified = {r["source_trench_id"] for r in payload["trench_selection"] if r["qualified"]}
            require({a["space_index"] for a in payload["psd_all_spaces"]} == qualified,
                    "all qualified PSD source IDs, without quota truncation or duplicate selected IDs")
            require(len(qualified) > 1 and len(payload["spaces"]) == 1, "PSD population exceeds max-number=1")
            captured.append(len(qualified))
        return measured

    with patch.object(app.edges, "measure_image", side_effect=inspect):
        out, frame = execute("v120_ensemble_full_population", root=root, plots=True,
                             extra=["--max-number", "1", "--no-annotated-images", "--plot-dpi", "60"])
    require(len(captured) == 4, "two images times two independent single-row engines")
    folder = out / "PSD"
    require(len(list(folder.rglob("*.psd.png"))) == 4, "mixed extensions do not overwrite PSD figures")
    summary = pd.read_csv(folder / "image_psd_summary.csv")
    require(len(summary) == 8 and summary.status.eq("OK").all(), "two signals times two engines times two images")
    require(summary.contributing_structures.gt(1).all() and summary.axis.eq("normal").all(), "aggregate full population, rotated only")
    require(not pd.read_csv(out / "processing_errors.csv").size, "PSD exports introduce no processing error")
    settings = json.loads((folder / "PSD_settings.json").read_text())
    require(settings["spectrum_averaging"] and not settings["noise_debiased"], "record ensemble method and noise limits")
    bands = pd.read_csv(folder / "image_psd_band_summary.csv")
    require(len(bands) == 12 and set(bands.band) == {"LF", "MF", "HF"}, "all band records exported")
    # Reused main edges also include all sources but retain their actual smoothing.
    out_reuse, _ = execute("v120_ensemble_measurement_source", root=root,
                          extra=["--max-number", "1", "--no-annotated-images", "--psd-edge-source", "measurement"])
    raw = pd.read_csv(out_reuse / "PSD/per_structure_psd_summary.csv")
    require(raw.average_range_px.eq(32).all(), "measurement source preserves recipe")
    require(pd.read_csv(out_reuse / "PSD/image_psd_summary.csv").contributing_structures.gt(1).all(), "measurement source includes unselected qualified trenches")
    return dict(status="PASS", image_count=2, PSD_png_count=4, max_number=1,
                all_qualified_counts=captured, signals=["LER", "LWR"], axes=["normal"],
                actual_native_positive_bins=sorted(set(summary.frequency_bin_count.astype(int))))
