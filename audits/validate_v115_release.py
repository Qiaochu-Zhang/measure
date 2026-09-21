#!/usr/bin/env python3
"""Check V1_15 docs, unchanged numerical code, and the independently unpacked ZIP."""

import ast
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import zipfile

import pandas as pd
from openpyxl import Workbook

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "sem_cd_measure_200k_batch_V1_15_complete_package"
BASELINE = ROOT / "sem_cd_measure_200k_batch_V1_14_complete_package"


def run(args, cwd):
    result = subprocess.run(
        [sys.executable, *map(str, args)],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return result.stdout


def main():
    sys.path.insert(0, str(PACKAGE))
    import sem_cd_measure_200k_batch_V1_15 as app

    report = {}
    guide = (ROOT / "V1_15_RUN_GUIDE.md").read_text()
    actions = app.build_parser()._actions
    for action in actions:
        assert any(f"`{name}`" in guide for name in action.option_strings), (
            action.option_strings
        )
    examples = []
    for line in guide.splitlines():
        if (
            line.startswith("python sem_cd_measure_200k_batch_V1_15.py ")
            and "--help" not in line
        ):
            args = shlex.split(line)[2:]
            app.validate_args(app.build_parser().parse_args(args), args)
            examples.append(line)
    report["documentation"] = dict(
        cli_options=len(actions), validated_examples=len(examples)
    )
    unchanged = []
    for name in (
        "cdsem_locator.py",
        "cdsem_refinement.py",
        "cdsem_statistics.py",
        "cdsem_psd.py",
    ):
        old = (BASELINE / name).read_text().replace("V1_14", "V1_15")
        new = (PACKAGE / name).read_text()
        assert ast.dump(ast.parse(old)) == ast.dump(ast.parse(new)), name
        unchanged.append(name)
    report["unchanged_algorithm_modules_after_version_rename"] = unchanged
    functions = []
    for package in (BASELINE, PACKAGE):
        tree = ast.parse((package / "cdsem_engine.py").read_text())
        functions.append(
            {
                node.name: ast.dump(node)
                for node in tree.body
                if isinstance(node, ast.FunctionDef)
            }
        )
    changed = [
        name for name in functions[0] if functions[0][name] != functions[1].get(name)
    ]
    assert changed == ["detect_space_candidates", "measure_image"], changed
    report["edge_engine_changed_functions"] = changed
    archive = ROOT / (PACKAGE.name + ".zip")
    with tempfile.TemporaryDirectory(prefix="v115-release-") as temp:
        temp = Path(temp)
        with zipfile.ZipFile(archive) as zipped:
            assert zipped.testzip() is None
            zipped.extractall(temp)
        unpacked = temp / PACKAGE.name
        manifest = (unpacked / "PACKAGE_SHA256.txt").read_text().splitlines()
        for line in manifest:
            digest, path = line.split("  ", 1)
            assert (
                hashlib.sha256((unpacked / path).read_bytes()).hexdigest() == digest
            ), path
        report["archive"] = dict(
            files_verified=len(manifest),
            sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
        )
        report["unpacked_quick_self_check"] = run(
            ["self_check_V1_15.py", "--quick"], unpacked
        ).strip()
        book = Workbook()
        book.active.title = "v2"
        book.active.append(["A", "name", "C", "CD", "LER", "LWR"])
        book.active.append([None, "test.png", None, 60, 2, 1])
        book.save(unpacked / "reference.xlsx")
        output = unpacked / "full_export_smoke"
        run(
            [
                "sem_cd_measure_200k_batch_V1_15.py",
                "--root",
                "examples/input_dark_line",
                "--output",
                output,
                "--pattern",
                "trench",
                "--pixel-size",
                "1",
                "--trench-reference-nm",
                "60",
                "--machine-excel",
                "reference.xlsx",
                "--save-debug-masks",
            ],
            unpacked,
        )
        images = pd.read_csv(output / "image_summary.csv")
        errors = pd.read_csv(output / "processing_errors.csv")
        machine = pd.read_csv(output / "machine_comparison_detailed.csv")
        assert images.status.eq("OK").all() and errors.empty
        assert images.locator_mode_selected.eq("dark-line").all()
        assert len(machine) > 0
        psd = pd.read_csv(output / "PSD/per_structure_psd_summary.csv")
        assert set(psd.engine) == {"V10", "V13"}
        assert psd.status.eq("OK").all()
        png_count = len(list(output.rglob("*.png")))
        assert png_count > 10
        report["unpacked_full_export"] = dict(
            status="PASS",
            image_statuses=list(images.status),
            locator_modes=list(images.locator_mode_selected),
            processing_errors=len(errors),
            machine_rows=len(machine),
            psd_curves=len(psd),
            png_exports=png_count,
            mixed_cd_nm=float(images["旋转_mixed_CD_nm"].iloc[0]),
        )
    evidence = ROOT / "audits/v1_15_release_evidence.json"
    evidence.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
