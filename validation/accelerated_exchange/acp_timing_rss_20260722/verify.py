#!/usr/bin/env python3
"""Verify the repeated ACP timing/RSS archive and its acceptance gates."""

from __future__ import annotations

import csv
import hashlib
import math
import re
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
ACCEPTED = RAW / "gxtb-acp-timing-20260722"
FLOAT = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"
EXPECTED_BINARY = "df6466552a495e94e710174dbf468ec1765f1a4690ef955f9129ef983f35790b"
EXPECTED_PROVIDER = "fe210c64a4c4fa6897668a8657dd234046143c55bda2f7c9279d24108f2f152a"
EXPECTED_INPUT = "7e1d23a4d9e0d66df81caa3744ca005342eb802599750f6888e8de38d2dd7f4f"
EXPECTED_ENERGY = -40.473748967057013
EXPECTED_STALE_CAMPAIGN_HASH = "5ebb02ec5199458145da0ed202ff6805251800b0bb448139cbca9c61d2f7ba45"
EXPECTED_FINAL_CAMPAIGN_HASH = "c8e30c828868c82561208512c9bbdc2c2f80dec78e6383e8fcd5f6c98a17fad2"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_portable_manifest() -> None:
    entries: set[Path] = set()
    for line in (ROOT / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        relative = Path(name.removeprefix("*"))
        path = ROOT / relative
        require(path.is_file(), f"missing manifest file: {relative}")
        require(sha256(path) == expected, f"portable checksum mismatch: {relative}")
        entries.add(relative)
    actual = {
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if path.is_file() and path != ROOT / "SHA256SUMS" and "__pycache__" not in path.parts
    }
    require(entries == actual, f"portable manifest coverage mismatch: {entries ^ actual}")


def verify_original_manifest() -> None:
    mismatches: list[tuple[str, str, str]] = []
    prefix = "/home/kuehne88/work/gxtb-acp-timing-20260722/"
    for line in (ACCEPTED / "SHA256SUMS").read_text().splitlines():
        expected, absolute = line.split(maxsplit=1)
        require(absolute.startswith(prefix), f"unexpected raw manifest path: {absolute}")
        relative = absolute[len(prefix) :]
        actual = sha256(ACCEPTED / relative)
        if actual != expected:
            mismatches.append((relative, expected, actual))
    require(
        mismatches
        == [
            (
                "provenance/campaign.txt",
                EXPECTED_STALE_CAMPAIGN_HASH,
                EXPECTED_FINAL_CAMPAIGN_HASH,
            )
        ],
        f"unexpected original-manifest mismatches: {mismatches}",
    )


def parse_observables(path: Path):
    text = path.read_text(errors="replace")
    require(text.count("PROGRAM ENDED AT") == 1, f"incomplete CP2K output: {path}")
    energies = re.findall(
        rf"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+({FLOAT})", text
    )
    forces = [
        tuple(map(float, row))
        for row in re.findall(
            rf"^ FORCES\|\s+\d+\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+{FLOAT}\s*$",
            text,
            re.MULTILINE,
        )
    ]
    stress_blocks = re.findall(
        r"STRESS\| Analytical stress tensor \[bar\](.*?)(?:STRESS\| 1/3 Trace)",
        text,
        re.DOTALL,
    )
    require(energies and forces and stress_blocks, f"missing observables: {path}")
    stress = [
        tuple(map(float, row))
        for row in re.findall(
            rf"^ STRESS\|\s+[xyz]\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
            stress_blocks[-1],
            re.MULTILINE,
        )
    ]
    require(len(stress) == 3, f"malformed stress tensor: {path}")
    return float(energies[-1]), forces, stress, text


def maximum_delta(left, right) -> float:
    require(len(left) == len(right), "paired observable lengths differ")
    return max(abs(a - b) for row_a, row_b in zip(left, right) for a, b in zip(row_a, row_b))


def median_mad(values):
    median = statistics.median(values)
    return median, statistics.median(abs(value - median) for value in values)


def verify_accepted_campaign() -> None:
    require((ACCEPTED / "campaign-exit-status.txt").read_text().strip() == "0", "campaign failed")
    hashes = (ACCEPTED / "provenance/binary-provider-input.sha256").read_text()
    for expected in (EXPECTED_BINARY, EXPECTED_PROVIDER, EXPECTED_INPUT):
        require(expected in hashes, f"missing provenance hash {expected}")
    campaign = (ACCEPTED / "provenance/campaign.txt").read_text()
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        require(f"{name}=1" in campaign, f"missing singleton thread proof for {name}")
    require("cpu=90" in campaign and "campaign_end=" in campaign, "incomplete campaign metadata")

    rows = list(csv.DictReader((ACCEPTED / "summary.tsv").open(), delimiter="\t"))
    require(len(rows) == 12, f"expected 12 runs, found {len(rows)}")
    require(sum(row["measured"] == "yes" for row in rows) == 10, "wrong measured-run count")
    require(sum(row["measured"] == "no" for row in rows) == 2, "wrong warm-up count")
    require(all(row["exit_code"] == "0" for row in rows), "nonzero run exit")
    require(all(float(row["energy_eh"]) == EXPECTED_ENERGY for row in rows), "energy mismatch")

    by_rep: dict[int, dict[str, tuple]] = {}
    for row in rows:
        sequence = int(row["sequence"])
        repetition = int(row["repetition"])
        mode = row["mode"]
        run_dir = ACCEPTED / "results" / f"{sequence:02d}-r{repetition}-{mode.lower()}"
        require((run_dir / "exit-status.txt").read_text().strip() == "0", f"bad exit: {run_dir}")
        affinity = (run_dir / "affinity-preexec.txt").read_text()
        require("expected_cpu=90" in affinity, f"wrong expected CPU: {run_dir}")
        require(re.search(r"^Cpus_allowed_list:\s*90$", affinity, re.MULTILINE), f"bad affinity: {run_dir}")
        parsed = parse_observables(run_dir / "cp2k.out")
        require(parsed[0] == EXPECTED_ENERGY, f"output energy mismatch: {run_dir}")
        if mode == "DENSE":
            require("GXTB-ACP-MESH STREAMED" not in parsed[3], f"dense selector leaked: {run_dir}")
        else:
            require("GXTB-ACP-MESH STREAMED nFull=" in parsed[3], f"missing streamed marker: {run_dir}")
            require("GXTB-ACP-MESH SPARSE-REVERSE projectorImages=" in parsed[3], f"missing sparse marker: {run_dir}")
        by_rep.setdefault(repetition, {})[mode] = parsed[:3]

    max_de = max_df = max_ds = 0.0
    for repetition, pair in by_rep.items():
        require(set(pair) == {"DENSE", "STREAMED"}, f"unpaired repetition {repetition}")
        dense, streamed = pair["DENSE"], pair["STREAMED"]
        max_de = max(max_de, abs(dense[0] - streamed[0]))
        max_df = max(max_df, maximum_delta(dense[1], streamed[1]))
        max_ds = max(max_ds, maximum_delta(dense[2], streamed[2]))
    require(max_de == 0.0, f"paired energy delta {max_de}")
    require(max_df <= 7.0e-17, f"paired force delta {max_df}")
    require(max_ds <= 4.0e-12, f"paired stress delta {max_ds}")

    minimum_margin = math.inf
    for path in sorted((ACCEPTED / "provenance").glob("*-budget.txt")):
        values = dict(line.split("=", 1) for line in path.read_text().splitlines())
        margin = int(values["computed_margin_kib"])
        floor = int(values["minimum_margin_kib"])
        require(margin >= floor, f"unsafe memory margin: {path}")
        minimum_margin = min(minimum_margin, margin)
    require(minimum_margin == 308890012, f"unexpected minimum margin {minimum_margin}")

    archived_stats = {
        row["mode"]: row
        for row in csv.DictReader((ACCEPTED / "statistics.tsv").open(), delimiter="\t")
    }
    comparison = {
        row["mode"]: row
        for row in csv.DictReader((ROOT / "comparison.tsv").open(), delimiter="\t")
    }
    measured = [row for row in rows if row["measured"] == "yes"]
    for mode in ("DENSE", "STREAMED"):
        selected = [row for row in measured if row["mode"] == mode]
        walls = [float(row["wall_s"]) for row in selected]
        rss = [int(row["peak_rss_kib"]) for row in selected]
        wall_median, wall_mad = median_mad(walls)
        rss_median, rss_mad = median_mad(rss)
        expected = archived_stats[mode]
        require(len(selected) == int(expected["n"]) == 5, f"wrong n for {mode}")
        require(abs(wall_median - float(expected["median_wall_s"])) < 5e-7, f"wall median {mode}")
        require(abs(wall_mad - float(expected["mad_wall_s"])) < 5e-7, f"wall MAD {mode}")
        require(rss_median == int(expected["median_peak_rss_kib"]), f"RSS median {mode}")
        require(rss_mad == int(expected["mad_peak_rss_kib"]), f"RSS MAD {mode}")
        table = comparison[mode]
        require(float(table["min_wall_s"]) == min(walls), f"wall minimum {mode}")
        require(float(table["max_wall_s"]) == max(walls), f"wall maximum {mode}")
        require(int(table["min_peak_rss_kib"]) == min(rss), f"RSS minimum {mode}")
        require(int(table["max_peak_rss_kib"]) == max(rss), f"RSS maximum {mode}")

    print(f"accepted runs: {len(rows)}; measured pairs: 5")
    print(f"max paired dE={max_de:.16e} Ha")
    print(f"max paired dF={max_df:.16e} Ha/bohr")
    print(f"max paired dStress={max_ds:.16e} bar")
    print(f"minimum launch margin={minimum_margin} KiB")


def verify_failed_attempts() -> None:
    attempt1 = (RAW / "gxtb-acp-timing-20260722-attempt1" / "campaign.log").read_text()
    attempt2 = (RAW / "gxtb-acp-timing-20260722-attempt2" / "campaign.log").read_text()
    attempt3 = RAW / "gxtb-acp-timing-20260722-attempt3"
    require("No such file or directory" in attempt1, "attempt1 reason changed")
    require("/usr/bin/time" in attempt2 and "No such file or directory" in attempt2, "attempt2 reason changed")
    require(not (attempt3 / "campaign-exit-status.txt").exists(), "attempt3 unexpectedly qualified")
    require(not (attempt3 / "statistics.tsv").exists(), "attempt3 unexpectedly has statistics")
    require(len(list((attempt3 / "results").glob("*"))) == 2, "attempt3 scope changed")


if __name__ == "__main__":
    verify_portable_manifest()
    verify_original_manifest()
    verify_accepted_campaign()
    verify_failed_attempts()
    print("ACP timing/RSS archive: PASS")
