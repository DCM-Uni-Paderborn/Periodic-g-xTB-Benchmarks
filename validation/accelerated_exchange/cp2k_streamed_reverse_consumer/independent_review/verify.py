#!/usr/bin/env python3
"""Independent verifier for the frozen CP2K streamed-reverse Linux evidence."""

from __future__ import annotations

import hashlib
import json
import math
import re
import tarfile
import tempfile
from pathlib import Path


ARCHIVE_ROOT = Path(__file__).resolve().parents[1]
RAW_ARCHIVE = (
    ARCHIVE_ROOT
    / "raw_archive"
    / "cp2k_gxtb_streamed_reverse_consumer_evidence_20260716.tar.gz"
)
_EXTRACTED = tempfile.TemporaryDirectory(prefix="gxtb-streamed-reverse-")
_EXTRACTED_ROOT = Path(_EXTRACTED.name).resolve()
with tarfile.open(RAW_ARCHIVE, "r:gz") as archive:
    members = archive.getmembers()
    for member in members:
        if member.issym() or member.islnk():
            raise AssertionError(f"archive link is not permitted: {member.name}")
        destination = (_EXTRACTED_ROOT / member.name).resolve()
        if destination != _EXTRACTED_ROOT and _EXTRACTED_ROOT not in destination.parents:
            raise AssertionError(f"archive path escapes extraction root: {member.name}")
    archive.extractall(_EXTRACTED_ROOT, members=members)

EVIDENCE_ROOT = _EXTRACTED_ROOT / "cp2k_gxtb_streamed_reverse_consumer_evidence"
ROOT = EVIDENCE_ROOT / "linux_matrix_terok"
MODE_ROOT = EVIDENCE_ROOT / "linux_mode_rss_terok"

CASES = {
    "k290_rks_3d_fd": (2, "709388624dda15c066e7b9b5d3e97fe57665a4d53e6f58851685dae298ca4e54"),
    "shifted_spglib_rks_3d": (1, "a93fbc314c6d9feeb71d24045d6a70e3c762de11f6eed9cd2a932156bea6db68"),
    "tr_rks_3d_fd": (2, "1b03853f044ccb1e58917dfd0443b846ffd967bfbad817380af553388f91677c"),
    "fullmesh_rks_3d_fd": (2, "184ac5aa81c88563177f4c67274dcfed7a6db12e52fcf3fa39b4ccb7af8951e2"),
    "fullmesh_uks_3d_fd": (2, "5875d3c725f5e939d738c9fcc7bf1b237915a0ca3d177e3cb85a838c19f17bc8"),
    "rks_1d_fd": (2, "539ac4fe8c8b22532bafc2551844834d487f88a89dcfa85619b85a3dece83963"),
    "spglib_rks_2d_fd": (2, "e2337936594342504351ebf472039461a49649bdbaa9a3b1f33fa9ae76b622a6"),
}
RANKS = (1, 2, 4)
QUAL_RE = re.compile(
    r"STREAMED-REVERSE dOverlap=\s*(\S+) dForce=\s*(\S+) "
    r"dStress=\s*(\S+) peak=(\d+) expected=(\d+)"
)
ENERGY_RE = re.compile(r"ENERGY\| Total FORCE_EVAL.*?([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)$")
FORCE_RE = re.compile(
    r"^\s*FORCES\|\s+\d+\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+[-+\d.Ee]+\s*$"
)
STRESS_RE = re.compile(
    r"^\s*STRESS\|\s+[xyz]\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*$"
)
FD_STRESS_RE = re.compile(r"DEBUG\| Sum of differences\s+([-+\d.Ee]+)\s*$")
FD_FORCE_RE = re.compile(r"DEBUG\| Sum of differences:\s+([-+\d.Ee]+)\s+[-+\d.Ee]+\s*$")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_manifest(run: Path) -> None:
    lines = (run / "SHA256SUMS").read_text().splitlines()
    assert lines, f"empty manifest: {run}"
    for line in lines:
        expected, name = line.split(None, 1)
        name = name.lstrip(" *")
        assert digest(run / name) == expected, f"hash mismatch: {run / name}"


def parse_last_observables(text: str) -> tuple[float, list[float], list[float]]:
    lines = text.splitlines()
    energies = [float(m.group(1)) for m in map(ENERGY_RE.search, lines) if m]
    forces = [tuple(float(x) for x in m.groups()) for m in map(FORCE_RE.match, lines) if m]
    stress_blocks: list[list[tuple[float, float, float]]] = []
    current_stress: list[tuple[float, float, float]] | None = None
    for line in lines:
        if "STRESS| Analytical stress tensor [bar]" in line:
            current_stress = []
            stress_blocks.append(current_stress)
            continue
        if current_stress is not None and len(current_stress) < 3:
            match = STRESS_RE.match(line)
            if match:
                current_stress.append(tuple(float(x) for x in match.groups()))
    assert energies
    # DEBUG finite-difference runs suppress the ordinary force/stress tables;
    # ENERGY_FORCE runs retain them.  Whenever present, every test system has
    # two atoms and the last two atomic rows plus last three analytical-stress
    # rows are the terminal evaluation used for rank checks.
    force_values = [x for row in forces[-2:] for x in row] if forces else []
    complete_stress = [block for block in stress_blocks if len(block) == 3]
    stress_values = [x for row in complete_stress[-1] for x in row] if complete_stress else []
    return energies[-1], force_values, stress_values


def max_delta(left: list[float], right: list[float]) -> float:
    assert len(left) == len(right)
    return max(abs(a - b) for a, b in zip(left, right))


max_overlap = 0.0
max_force = 0.0
max_stress = 0.0
max_rank_energy = 0.0
max_rank_force = 0.0
max_rank_stress = 0.0
max_fd_force_sum = 0.0
max_fd_stress_sum = 0.0
qualifications = 0
runs = 0
max_single_rss = 0
case_summary: dict[str, dict[str, float]] = {}

for case, (expected_qual, expected_input_hash) in CASES.items():
    reference = None
    case_wall: dict[str, float] = {}
    case_tree_rss: dict[str, int] = {}
    for ranks in RANKS:
        run = ROOT / f"{case}_p{ranks}"
        assert (run / "PASS").is_file(), f"missing PASS: {run}"
        assert (run / "returncode.txt").read_text().strip() == "0"
        assert (run / "input.sha256").read_text().split()[0] == expected_input_hash
        verify_manifest(run)
        text = (run / "run.out").read_text(errors="replace")
        assert text.count("PROGRAM ENDED AT") == 1
        matches = list(QUAL_RE.finditer(text))
        assert len(matches) == expected_qual, (run, len(matches), expected_qual)
        for match in matches:
            overlap, force, stress = (float(match.group(i)) for i in range(1, 4))
            peak, expected = int(match.group(4)), int(match.group(5))
            assert peak == expected
            assert max(overlap, force, stress) <= 1.0e-10
            assert all(math.isfinite(x) for x in (overlap, force, stress))
            max_overlap = max(max_overlap, overlap)
            max_force = max(max_force, force)
            max_stress = max(max_stress, stress)
        fd_force = [float(match.group(1)) for match in map(FD_FORCE_RE.search, text.splitlines()) if match]
        fd_stress = [float(match.group(1)) for match in map(FD_STRESS_RE.search, text.splitlines()) if match]
        if expected_qual == 2:
            assert len(fd_force) == 1 and len(fd_stress) == 1
            max_fd_force_sum = max(max_fd_force_sum, *map(abs, fd_force))
            max_fd_stress_sum = max(max_fd_stress_sum, *map(abs, fd_stress))
        qualifications += len(matches)
        runs += 1
        observables = parse_last_observables(text)
        if reference is None:
            reference = observables
        else:
            max_rank_energy = max(max_rank_energy, abs(observables[0] - reference[0]))
            if observables[1] and reference[1]:
                max_rank_force = max(max_rank_force, max_delta(observables[1], reference[1]))
            if observables[2] and reference[2]:
                max_rank_stress = max(max_rank_stress, max_delta(observables[2], reference[2]))
        usage = json.loads((run / "rusage.json").read_text())
        assert usage["returncode"] == 0
        case_wall[str(ranks)] = usage["elapsed_seconds"]
        case_tree_rss[str(ranks)] = usage["peak_tree_rss_kib"]
        max_single_rss = max(max_single_rss, usage["peak_single_process_rss_kib"])
    case_summary[case] = {
        **{f"wall_p{p}": case_wall[str(p)] for p in RANKS},
        **{f"tree_rss_kib_p{p}": case_tree_rss[str(p)] for p in RANKS},
    }

assert runs == 21
assert qualifications == 39

mode = {}
for name in ("dense", "streamed"):
    run = MODE_ROOT / f"{name}_p1"
    assert (run / "PASS").is_file()
    assert (run / "returncode.txt").read_text().strip() == "0"
    verify_manifest(run)
    text = (run / "run.out").read_text(errors="replace")
    assert text.count("PROGRAM ENDED AT") == 1
    obs = parse_last_observables(text)
    assert obs[1] and obs[2]
    usage = json.loads((run / "rusage.json").read_text())
    mode[name] = {"observables": obs, "usage": usage}

mode_energy_delta = abs(mode["dense"]["observables"][0] - mode["streamed"]["observables"][0])
mode_force_delta = max_delta(mode["dense"]["observables"][1], mode["streamed"]["observables"][1])
mode_stress_delta = max_delta(mode["dense"]["observables"][2], mode["streamed"]["observables"][2])

print(json.dumps({
    "status": "PASS",
    "runs": runs,
    "qualifications": qualifications,
    "max_dense_oracle_residual": {
        "overlap": max_overlap,
        "force": max_force,
        "stress": max_stress,
    },
    "max_p1_p2_p4_printed_delta": {
        "energy_hartree": max_rank_energy,
        "force_hartree_per_bohr": max_rank_force,
        "stress_bar": max_rank_stress,
    },
    "max_finite_difference_sum": {
        "force_hartree_per_bohr": max_fd_force_sum,
        "stress_atomic_units": max_fd_stress_sum,
    },
    "dense_streamed_p1_printed_delta": {
        "energy_hartree": mode_energy_delta,
        "force_hartree_per_bohr": mode_force_delta,
        "stress_bar": mode_stress_delta,
    },
    "dense_streamed_p1_resources": {
        key: {
            "wall_seconds": value["usage"]["elapsed_seconds"],
            "peak_tree_rss_kib": value["usage"]["peak_tree_rss_kib"],
            "peak_single_process_rss_kib": value["usage"]["peak_single_process_rss_kib"],
        }
        for key, value in mode.items()
    },
    "max_single_process_rss_kib": max_single_rss,
    "cases": case_summary,
}, indent=2, sort_keys=True))
