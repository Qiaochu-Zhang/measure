#!/usr/bin/env python3
"""Build the self-contained V1_20 archive; exclude caches and run outputs."""

from argparse import ArgumentParser
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import zipfile


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "sem_cd_measure_200k_batch_V1_20_complete_package"


def package_files():
    root_suffixes = {".py", ".md", ".txt", ".json", ".ps1", ".sh"}
    for path in sorted(PACKAGE.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        relative = path.relative_to(PACKAGE)
        if len(relative.parts) == 1 and (
            path.suffix in root_suffixes or path.name == ".gitignore"
        ):
            yield path
        elif (
            relative.parts[0] == "examples"
            and path.suffix == ".png"
            and relative.parts[1]
            in {
                "input_trench",
                "input_line",
                "input_bright_line",
                "input_dark_line",
                "input_wiggling",
            }
        ):
            yield path


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--replace", action="store_true", help="replace an existing V1_20 archive"
    )
    args = parser.parse_args()
    archive = ROOT / (PACKAGE.name + ".zip")
    if archive.exists() and not args.replace:
        raise SystemExit(f"Archive exists; use --replace to rebuild: {archive}")
    report = PACKAGE / "validation_report.json"
    if not report.is_file():
        raise SystemExit("Run self_check_V1_20.py --output validation_report.json first")
    manifest = PACKAGE / "PACKAGE_SHA256.txt"
    lines = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  "
        f"{path.relative_to(PACKAGE).as_posix()}"
        for path in package_files()
        if path != manifest
    ]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with TemporaryDirectory(prefix=".v120-package-", dir=ROOT) as directory:
        temporary = Path(directory) / archive.name
        with zipfile.ZipFile(
            temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
        ) as output:
            for path in package_files():
                info = zipfile.ZipInfo(
                    (Path(PACKAGE.name) / path.relative_to(PACKAGE)).as_posix(),
                    date_time=(2026, 10, 28, 0, 0, 0),
                )
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = (
                    0o100755 if path.suffix == ".sh" else 0o100644
                ) << 16
                output.writestr(info, path.read_bytes())
        with zipfile.ZipFile(temporary) as result:
            assert result.testzip() is None
            print(
                f"Created {archive.name}: {len(result.infolist())} files, "
                f"{temporary.stat().st_size} bytes"
            )
        temporary.replace(archive)
    print("SHA256:", hashlib.sha256(archive.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
