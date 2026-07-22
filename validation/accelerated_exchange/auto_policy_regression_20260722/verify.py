#!/usr/bin/env python3
"""Verify the exact-build periodic g-xTB automatic-policy regression."""

from __future__ import annotations

import csv
import hashlib
import io
import re
import tarfile
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CP2K_REVISION = "f44008823d3319547f34ef335561256816a1a031"
SAVE_REVISION = "718629fbc86e0b362491cf70dd4198d0d82082b5"
CP2K_BINARY = "8850e8a39c14fbd172ab89a7992cee69e492b3d6ab039451985c312711e3e0aa"
CH4_ENERGY = -40.468866070692435
O2_ENERGY = -150.340821038026832
O2_CLI_DELTA = 5.911715561524e-12
O2_GRADIENT_DELTA = 3.300347844518e-06
FLOAT = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest() -> None:
    listed: set[Path] = set()
    for line in (ROOT / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        relative = Path(name.removeprefix("*"))
        path = ROOT / relative
        require(path.is_file(), f"missing manifest file: {relative}")
        require(sha256(path) == expected, f"checksum mismatch: {relative}")
        listed.add(relative)
    actual = {
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if path.is_file()
        and path != ROOT / "SHA256SUMS"
        and "__pycache__" not in path.parts
    }
    require(listed == actual, f"manifest coverage mismatch: {listed ^ actual}")


def verify_summary() -> None:
    rows = list(csv.DictReader((ROOT / "summary.tsv").open(), delimiter="\t"))
    require(len(rows) == 1, "summary must contain exactly one suite")
    row = rows[0]
    expected = {
        "suite": "xTB/regtest-tblite-gxtb",
        "calculations": "48",
        "matcher_results": "196",
        "correct": "196",
        "wrong": "0",
        "failed": "0",
        "normally_terminated": "48",
        "cp2k_revision": CP2K_REVISION,
        "save_tblite_revision": SAVE_REVISION,
        "status": "PASS",
    }
    require(row == expected, f"unexpected regression summary: {row}")
    identity = (ROOT / "build_identity.txt").read_text()
    require(f"cp2k_revision={CP2K_REVISION}" in identity, "wrong CP2K source")
    require(f"save_tblite_revision={SAVE_REVISION}" in identity, "wrong provider source")
    require(f"cp2k_executable_sha256={CP2K_BINARY}" in identity, "wrong binary hash")


def energy(text: str) -> float:
    matches = re.findall(
        rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})",
        text,
    )
    require(bool(matches), "missing CP2K total energy")
    return float(matches[-1])


def verify_selected_outputs() -> None:
    ch4 = (ROOT / "outputs/CH4_gxtb_kp_auto_default.inp.out").read_text()
    o2 = (ROOT / "outputs/O2_gxtb_uks_reference_cli.inp.out").read_text()
    for name, text in (("CH4", ch4), ("O2", o2)):
        require(text.count("PROGRAM ENDED AT") == 1, f"{name} did not terminate once")
        require(
            re.search(r"CP2K\| source code revision number:\s+f44008823", text),
            f"{name} carries the wrong CP2K revision",
        )
    require(energy(ch4) == CH4_ENERGY, "CH4 energy changed")
    require(
        "GXTB-ACCELERATION mode=AUTO exchange=SYMMETRY_FUSED gradient=STREAMED transform=MIXED_RADIX_FFT"
        in ch4,
        "missing effective AUTO selector",
    )
    require("GXTB-ACP-MESH STREAMED nFull=27 batch=8 fullStorage=0" in ch4, "missing ACP stream")
    require("GXTB-ACP-MESH SPARSE-REVERSE" in ch4, "missing sparse ACP reverse")
    require(energy(o2) == O2_ENERGY, "O2 energy changed")
    match = re.search(rf"Energy CP2K/CLI/absdiff:\s+{FLOAT}\s+{FLOAT}\s+({FLOAT})", o2)
    require(match and float(match.group(1)) == O2_CLI_DELTA, "O2 CLI delta changed")
    match = re.search(rf"Gradient diff sum/max:\s+{FLOAT}\s+({FLOAT})", o2)
    require(match and float(match.group(1)) == O2_GRADIENT_DELTA, "O2 gradient delta changed")


def verify_raw_suite() -> None:
    archive = ROOT / "raw/focused_gxtb_regression_raw.tar.gz"
    prefix = "xTB/regtest-tblite-gxtb/"
    with tarfile.open(archive, "r:gz") as handle:
        members = {member.name: member for member in handle.getmembers() if member.isfile()}
        inputs = sorted(
            name
            for name in members
            if name.startswith(prefix)
            and "/" not in name[len(prefix) :]
            and name.endswith(".inp")
        )
        outputs = sorted(
            name
            for name in members
            if name.startswith(prefix)
            and "/" not in name[len(prefix) :]
            and name.endswith(".inp.out")
        )
        require(len(inputs) == 48, f"expected 48 inputs, found {len(inputs)}")
        require(len(outputs) == 48, f"expected 48 outputs, found {len(outputs)}")
        for name in outputs:
            stream = handle.extractfile(members[name])
            require(stream is not None, f"cannot read {name}")
            text = stream.read().decode(errors="replace")
            require(text.count("PROGRAM ENDED AT") == 1, f"incomplete raw output: {name}")
            require(
                re.search(r"CP2K\| source code revision number:\s+f44008823", text),
                f"wrong raw source revision: {name}",
            )
        manifest_name = prefix + "TEST_FILES.toml"
        stream = handle.extractfile(members[manifest_name])
        require(stream is not None, "missing raw matcher manifest")
        tests = tomllib.load(io.BytesIO(stream.read()))
        require(len(tests) == 48, f"expected 48 matcher cases, found {len(tests)}")
        require(sum(len(matchers) for matchers in tests.values()) == 196, "matcher count changed")


if __name__ == "__main__":
    verify_manifest()
    verify_summary()
    verify_selected_outputs()
    verify_raw_suite()
    print("automatic-policy exact-build regression: PASS (196/196)")
