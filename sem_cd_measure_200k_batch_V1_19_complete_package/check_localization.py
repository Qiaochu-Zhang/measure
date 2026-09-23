"""Synthetic V1_19 regression fixtures and assertions (no network required)."""

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

import cdsem_localization as loc
import cdsem_regions as regions


def synthetic_pattern(dark_line=False, noise=2.0, seed=22, stain_width=40):
    image = np.full((440, 840), 210.0)
    for center in [100, 260, 420, 580, 740]:
        image[:, center - 30 : center + 30] = 40
    if dark_line:
        for center in [180, 340, 500, 660]:
            image[:, center - stain_width // 2 : center + stain_width // 2] = 105
    image = gaussian_filter(image, 1) + np.random.default_rng(seed).normal(
        0, noise, image.shape
    )
    padded = np.full((620, 1080), 155, np.uint8)
    padded[80:520, 120:960] = np.clip(image, 0, 255).astype(np.uint8)
    return padded


def analytical_tests(require, app):
    result = {}
    for dark_line in (False, True):
        for noise in (0, 2, 8):
            image = synthetic_pattern(dark_line, noise)
            region = regions.trim_region_ends(regions.find_pattern_region(image, 60))
            x0, x1, y0, y1 = region.bounds
            observed = loc.observation_mask(region.mask, 60)[y0:y1, x0:x1]
            diag = loc.analyze_mode(image[y0:y1, x0:x1], observed)
            expected = "dark-line" if dark_line else "bright-line"
            require(
                diag["locator_mode_selected"] == expected, "auto morphology " + expected
            )
            require(region.diagnostics["roi_end_trim_px"] == 22, "5% of 440px")
            require(
                not region.mask[:102].any() and not region.mask[498:].any(),
                "both 22px end regions excluded",
            )
            require(
                not region.mask[:, :150].any() and not region.mask[:, 930:].any(),
                "lateral background excluded",
            )
            for mode in ("dark-line", "bright-line"):
                explicit = loc.analyze_mode(image[y0:y1, x0:x1], observed, mode)
                require(
                    explicit["locator_mode_selected"] == mode,
                    "explicit mode precedence",
                )
            # Positive affine brightness change should preserve the auto decision.
            shifted = loc.analyze_mode(0.65 * image[y0:y1, x0:x1] + 20, observed)
            require(
                shifted["locator_mode_selected"] == expected,
                "adaptive intensity threshold",
            )
            result[f"{expected}_noise_{noise}"] = diag
    narrow = synthetic_pattern(True, stain_width=14)
    region = regions.trim_region_ends(regions.find_pattern_region(narrow, 60))
    observed = loc.observation_mask(region.mask, 60)
    diag = loc.analyze_mode(narrow, observed)
    require(
        diag["locator_mode_selected"] == "dark-line",
        "auto sees narrow fake stripe outside trench edge-context mask",
    )
    result["narrow_fake_line"] = diag
    mask = np.zeros((500, 100), bool)
    mask[100:300, 10:30] = True
    mask[50:450, 60:80] = True
    mask[240:260, 60:80] = False
    base = regions.PatternRegion(mask, (10, 80, 50, 450), {})
    trimmed = regions.trim_region_ends(base)
    require(
        trimmed.diagnostics["roi_mean_trench_length_px"] == 300,
        "one mean length estimate",
    )
    require(trimmed.diagnostics["roi_end_trim_px"] == 15, "common end exclusion in px")
    require(
        trimmed.mask[115:285, 20].all() and not trimmed.mask[:115, 20].any(),
        "short trench ends trimmed locally",
    )
    require(not trimmed.mask[240:260, 60:80].any(), "interior gap remains excluded")
    require(
        np.array_equal(regions.trim_region_ends(base, 0).mask, mask),
        "zero trim compatibility",
    )
    result["unequal_lengths_and_internal_gap"] = "PASS"

    # Background replacement cannot affect the coarse profile or mode choice.
    image = synthetic_pattern(True)
    region = regions.trim_region_ends(regions.find_pattern_region(image, 60))
    modified = image.copy()
    modified[~region.mask] = np.random.default_rng(15).integers(
        0, 256, (~region.mask).sum()
    )
    require(
        np.array_equal(
            loc.supported_profile(image, region.mask),
            loc.supported_profile(modified, region.mask),
        ),
        "excluded pixels must not enter locator profile",
    )
    result["background_profile_invariance"] = "PASS"

    for args in (
        ["--end-trim-fraction", "-0.1"],
        ["--end-trim-fraction", "0.5"],
        ["--locator-majority", "0.5"],
        ["--locator-majority", "1.1"],
        ["--end-trim-fraction", "nan"],
    ):
        parsed = app.build_parser().parse_args(["--pattern", "trench", *args])
        try:
            app.validate_args(parsed, args)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid V1_19 parameter accepted: " + str(args))
    result["parameter_validation"] = "PASS"
    return result


def integration_tests(temp, execute, require):
    result = {}
    for dark_line in (False, True):
        expected = "dark-line" if dark_line else "bright-line"
        root = temp / (expected + "_input")
        regions.unicode_imwrite(
            root / "中文 子目录/offset.png", synthetic_pattern(dark_line)
        )
        runs = {}
        for mode in ("auto", expected):
            output, frame = execute(
                expected + "_" + mode,
                root=root,
                extra=[
                    "--locator-mode",
                    mode,
                    "--max-number",
                    "3",
                    "--skip-psd",
                    "--no-annotated-images",
                ],
            )
            require(frame.status.eq("OK").all(), expected + " integration status")
            require(
                frame.locator_mode_selected.eq(expected).all(), "per-image chosen mode"
            )
            candidates = pd.read_csv(output / "locator_candidates.csv")
            selected = candidates[candidates.selected]
            for engine in ("V10", "V13"):
                centers = sorted(
                    selected.loc[selected.engine == engine, "candidate_center_x_px"]
                )
                require(
                    len(centers) == 3 and np.allclose(centers, [380, 540, 700], atol=2),
                    "correct adjacent true trenches, original-image coordinates",
                )
            runs[mode] = pd.read_csv(output / "per_sample_results.csv")
            result[expected + "_" + mode] = dict(status="PASS", centers_px=centers)
        for column in ("未旋转_left_edge_x_px", "未旋转_right_edge_x_px", "valid"):
            require(
                np.allclose(
                    runs["auto"][column], runs[expected][column], equal_nan=True
                ),
                "auto and explicit selected mode must use identical edge pipeline",
            )

    root = temp / "end_gap_input"
    image = synthetic_pattern(False)
    image[285:300] = 155
    regions.unicode_imwrite(root / "gap.png", image)
    output, frame = execute(
        "end_gap_psd_continuity",
        root=root,
        extra=[
            "--max-number",
            "3",
            "--extend-length",
            "440",
            "--average-range",
            "8",
            "--edge-continuity",
            "100",
            "--psd-gap-mode",
            "interpolate",
            "--psd-max-gap",
            "64",
            "--no-annotated-images",
        ],
    )
    require(
        frame.status.iloc[0] in {"OK", "REVIEW"}, "end/gap measurement retains support"
    )
    samples = pd.read_csv(output / "per_sample_results.csv")
    require(
        samples.background_excluded.any(), "explicit long span includes excluded ends"
    )
    require(
        not samples.loc[samples.background_excluded, "valid"].any(),
        "no end/gap synthetic filling",
    )
    valid = samples[samples.valid]
    require(
        (valid.average_y0 >= 102).all() and (valid.average_y1 <= 498).all(),
        "complete averaging windows lie inside trimmed ends",
    )
    require(
        ((valid.average_y1 <= 285) | (valid.average_y0 >= 300)).all(),
        "no averaging across internal gap",
    )
    psd = pd.read_csv(output / "PSD/edge_coordinates.csv")
    require(
        not psd.loc[psd.background_excluded, "valid"].any(),
        "PSD excludes ends even with interpolation",
    )
    require(
        psd.loc[
            (psd.sample_y < 102) | (psd.sample_y >= 498), "background_excluded"
        ].all(),
        "PSD raw rows honor end masks",
    )
    result["ends_gap_psd_continuity"] = dict(
        status="PASS", valid_main_samples=len(valid)
    )
    return result
