#!/usr/bin/env python3
"""Compare CP2K Gamma-point g-xTB ice energies with periodic save_tblite CLI."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
GXTB_PROTOCOL_ID = "dmc13-gxtb-spglib-reduced-v1"
CLI_PROTOCOL_ID = "dmc13-gxtb-spglib-gamma-cli-v1"
GXTB_INPUT_DIRECTORY = "gxtb_spglib_inputs"
GXTB_RUN_DIRECTORY = "runs_gxtb_spglib"
CLI_RUN_DIRECTORY = "runs_cli_gxtb_spglib"
ACCURACY = 0.1
DEFAULT_CAMPAIGN_MANIFEST = (
    ROOT.parent
    / "campaigns"
    / "gxtb-pbc-v1-20260714"
    / "build_manifest.json"
)


@dataclass(frozen=True)
class BuildIdentity:
    campaign_id: str
    cp2k: Path
    cp2k_sha256: str
    cp2k_library: Path
    cp2k_library_sha256: str
    tblite: Path
    tblite_sha256: str
    tblite_static_library: Path
    tblite_static_library_sha256: str
    cp2k_source_revision: str
    tblite_source_revision: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(command: list[str], cwd: Path | None = None) -> str:
    process = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    return (process.stdout + process.stderr).strip()


def cp2k_rpaths(cp2k: Path) -> list[Path]:
    process = subprocess.run(
        ["otool", "-l", str(cp2k)],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ValueError(f"cannot inspect CP2K RPATHs: {process.stderr.strip()}")
    rpaths: list[Path] = []
    waiting_for_path = False
    for raw_line in process.stdout.splitlines():
        line = raw_line.strip()
        if line == "cmd LC_RPATH":
            waiting_for_path = True
            continue
        if not waiting_for_path or not line.startswith("path "):
            continue
        value = line[5:].split(" (offset ", 1)[0]
        value = value.replace("@loader_path", str(cp2k.parent))
        value = value.replace("@executable_path", str(cp2k.parent))
        rpaths.append(Path(value).resolve())
        waiting_for_path = False
    return rpaths


def resolve_cp2k_library(cp2k: Path) -> Path:
    process = subprocess.run(
        ["otool", "-L", str(cp2k)],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ValueError(f"cannot inspect CP2K dependencies: {process.stderr.strip()}")
    dependencies = [
        line.strip().split(" (", 1)[0]
        for line in process.stdout.splitlines()[1:]
        if line.strip()
    ]
    matches = [
        dependency
        for dependency in dependencies
        if Path(dependency).name.startswith("libcp2k")
        and dependency.endswith(".dylib")
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one loaded libcp2k dylib, found {matches!r}"
        )
    dependency = matches[0]
    if dependency.startswith("@rpath/"):
        suffix = dependency.removeprefix("@rpath/")
        candidates = [rpath / suffix for rpath in cp2k_rpaths(cp2k)]
    elif dependency.startswith("@loader_path/"):
        candidates = [cp2k.parent / dependency.removeprefix("@loader_path/")]
    elif dependency.startswith("@executable_path/"):
        candidates = [cp2k.parent / dependency.removeprefix("@executable_path/")]
    else:
        candidates = [Path(dependency)]
    existing = [candidate.resolve() for candidate in candidates if candidate.is_file()]
    if len(existing) != 1:
        raise ValueError(
            "could not resolve the uniquely loaded libcp2k dylib; candidates: "
            + ", ".join(str(candidate) for candidate in candidates)
        )
    return existing[0]


def embedded_revision(text: str, label: str) -> str:
    match = re.search(rf"{re.escape(label)}\s*([0-9a-fA-F]{{7,40}})", text)
    if match is None:
        raise ValueError(f"cannot find {label.strip()} in build metadata")
    return match.group(1).lower()


def resolve_source_revision(source: Path, revision: str, label: str) -> str:
    source = source.resolve()
    head = command_output(["git", "rev-parse", "HEAD"], source)
    resolved = command_output(
        ["git", "rev-parse", f"{revision}^{{commit}}"],
        source,
    )
    if not re.fullmatch(r"[0-9a-f]{40}", head) or resolved != head:
        raise ValueError(
            f"{label} embedded revision {revision} does not equal source HEAD {head}"
        )
    return head


def build_identity(
    campaign_id: str,
    cp2k: Path,
    cp2k_library_expected: Path | None,
    tblite: Path,
    tblite_static_library: Path,
    cp2k_source: Path,
    tblite_source: Path,
    tblite_revision_expected: str,
) -> BuildIdentity:
    cp2k = cp2k.resolve()
    tblite = tblite.resolve()
    tblite_static_library = tblite_static_library.resolve()
    detected_library = resolve_cp2k_library(cp2k)
    if (
        cp2k_library_expected is not None
        and cp2k_library_expected.resolve() != detected_library
    ):
        raise ValueError(
            f"--cp2k-library resolves to {cp2k_library_expected.resolve()}, but "
            f"otool/RPATH selects {detected_library}"
        )
    if not tblite_static_library.is_file():
        raise ValueError(f"missing static save_tblite library: {tblite_static_library}")
    if tblite.parent.parent != tblite_static_library.parent.parent:
        raise ValueError(
            "tblite CLI and static library must come from the same install prefix"
        )
    cp2k_revision = embedded_revision(
        command_output([str(cp2k), "--version"]),
        "Source code revision",
    )
    cp2k_source_revision = resolve_source_revision(
        cp2k_source,
        cp2k_revision,
        "CP2K",
    )
    library_strings = command_output(["strings", str(detected_library)])
    revision_match = re.search(
        r"tblite source revision:\s*([0-9a-fA-F]{7,40}|unknown)",
        library_strings,
    )
    if revision_match is None:
        raise ValueError("cannot find tblite source revision in libcp2k metadata")
    tblite_embedded = revision_match.group(1).lower()
    if tblite_embedded != "unknown":
        embedded_full = command_output(
            ["git", "rev-parse", f"{tblite_embedded}^{{commit}}"],
            tblite_source,
        )
        if embedded_full != tblite_revision_expected:
            raise ValueError(
                "libcp2k embeds save_tblite revision "
                f"{tblite_embedded}, expected {tblite_revision_expected}"
            )
    tblite_source_revision = resolve_source_revision(
        tblite_source,
        tblite_revision_expected,
        "save_tblite",
    )
    return BuildIdentity(
        campaign_id=campaign_id,
        cp2k=cp2k,
        cp2k_sha256=sha256(cp2k),
        cp2k_library=detected_library,
        cp2k_library_sha256=sha256(detected_library),
        tblite=tblite,
        tblite_sha256=sha256(tblite),
        tblite_static_library=tblite_static_library,
        tblite_static_library_sha256=sha256(tblite_static_library),
        cp2k_source_revision=cp2k_source_revision,
        tblite_source_revision=tblite_source_revision,
    )


def read_campaign_manifest(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read campaign manifest {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"campaign manifest {path} is not a JSON object")
    return payload


def validate_campaign_identity(
    identity: BuildIdentity,
    manifest: dict[str, object],
) -> None:
    cp2k = manifest.get("cp2k")
    save_tblite = manifest.get("save_tblite")
    if not isinstance(cp2k, dict) or not isinstance(save_tblite, dict):
        raise ValueError("campaign manifest lacks cp2k/save_tblite identity blocks")
    checks = {
        "campaign_id": (identity.campaign_id, manifest.get("campaign_id")),
        "cp2k revision": (
            identity.cp2k_source_revision,
            cp2k.get("revision"),
        ),
        "CP2K launcher SHA256": (
            identity.cp2k_sha256,
            cp2k.get("binary_sha256"),
        ),
        "loaded libcp2k SHA256": (
            identity.cp2k_library_sha256,
            cp2k.get("loaded_library_sha256"),
        ),
        "save_tblite revision": (
            identity.tblite_source_revision,
            save_tblite.get("revision"),
        ),
        "save_tblite CLI SHA256": (
            identity.tblite_sha256,
            save_tblite.get("cli_sha256"),
        ),
        "static libtblite SHA256": (
            identity.tblite_static_library_sha256,
            save_tblite.get("static_library_sha256"),
        ),
    }
    mismatches = [
        f"{label}: actual {actual!r}, manifest {expected!r}"
        for label, (actual, expected) in checks.items()
        if actual != expected
    ]
    if mismatches:
        raise ValueError("campaign identity mismatch: " + "; ".join(mismatches))


def git_metadata(source: Path | None) -> dict[str, object] | None:
    if source is None:
        return None
    return {
        "path": str(source.resolve()),
        "revision": command_output(["git", "rev-parse", "HEAD"], source),
        "branch": command_output(["git", "branch", "--show-current"], source),
        "status": command_output(["git", "status", "--short"], source),
    }


def coordinate_labels(path: Path) -> list[str]:
    labels: list[str] = []
    in_coord = False
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        upper = line.upper()
        if upper == "&COORD":
            in_coord = True
            continue
        if in_coord and upper in {"&END COORD", "&END"}:
            break
        if not in_coord or not line or upper == "SCALED":
            continue
        fields = line.split()
        if len(fields) >= 4:
            labels.append(re.sub(r"[^A-Za-z].*$", "", fields[0]))
    return labels


def poscar_text(phase: str, geometry: dict[str, object], labels: list[str]) -> str:
    coords = geometry["coords"]
    if not isinstance(coords, list) or len(coords) != len(labels):
        raise ValueError(f"{phase}: coordinate/element count mismatch")
    grouped: OrderedDict[str, list[list[float]]] = OrderedDict()
    for element, coord in zip(labels, coords):
        grouped.setdefault(element, []).append(coord)
    cell = geometry["cell"]
    mode = str(geometry["mode"])
    if mode.lower() not in {"direct", "cartesian"}:
        raise ValueError(f"{phase}: unsupported coordinate mode {mode!r}")
    lines = [f"DMC-ICE13 ice {phase}", "1.0"]
    lines.extend("  " + " ".join(f"{float(value):.14f}" for value in vector) for vector in cell)
    lines.append("  " + " ".join(grouped))
    lines.append("  " + " ".join(str(len(values)) for values in grouped.values()))
    lines.append(mode)
    for values in grouped.values():
        lines.extend("  " + " ".join(f"{float(value):.14f}" for value in coord) for coord in values)
    return "\n".join(lines) + "\n"


def parse_cp2k_energy(path: Path) -> float | None:
    if not path.is_file():
        return None
    text = path.read_text(errors="ignore")
    if (
        "PROGRAM ENDED" not in text
        or "SCF run converged" not in text
        or "SCF run NOT converged" in text
        or "ABORT" in text
    ):
        return None
    values = re.findall(
        r"ENERGY\| Total FORCE_EVAL.*?([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*$",
        text,
        flags=re.MULTILINE,
    )
    return float(values[-1]) if values else None


def parse_cli_energy(path: Path) -> float | None:
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text()).get("energy")
        return float(value) if value is not None else None
    except (json.JSONDecodeError, TypeError, ValueError):
        return None


def parse_cli_atomic_energies(path: Path) -> list[float] | None:
    if not path.is_file():
        return None
    try:
        values = json.loads(path.read_text()).get("energies")
        if not isinstance(values, list):
            return None
        energies = [float(value) for value in values]
        return energies if all(math.isfinite(value) for value in energies) else None
    except (json.JSONDecodeError, TypeError, ValueError):
        return None


def archive_stale_file(path: Path) -> Path | None:
    """Move an incompatible artifact aside instead of silently adopting it."""
    if not path.is_file():
        return None
    digest = sha256(path)[:12]
    archived = path.with_name(f"{path.name}.stale-{digest}")
    serial = 1
    while archived.exists():
        archived = path.with_name(f"{path.name}.stale-{digest}.{serial}")
        serial += 1
    path.replace(archived)
    return archived


def cp2k_stamp_path(output: Path) -> Path:
    return output.with_suffix(".run.json")


def gamma_input_contract_errors(path: Path) -> list[str]:
    if not path.is_file():
        return ["production input is missing"]
    lines = {
        line.strip().upper()
        for line in path.read_text(errors="ignore").splitlines()
        if line.strip()
    }
    errors: list[str] = []
    for required in (
        f"# DMC13_GXTB_PROTOCOL {GXTB_PROTOCOL_ID}".upper(),
        "METHOD GXTB",
        "ACCURACY 0.1",
        "SCC_MIXER TBLITE",
        "ITERATIONS 300",
        "EPS_SCF 1.0E-9",
        "METHOD DIRECT_P_MIXING",
        "ALPHA 0.2",
        "CANONICALIZE TRUE",
    ):
        if required not in lines:
            errors.append(f"missing {required}")
    if "&KPOINTS" in lines:
        errors.append("Gamma production input must be implicit without &KPOINTS")
    return errors


def cp2k_validation_errors(
    phase: str,
    input_path: Path,
    output: Path,
    identity: BuildIdentity,
) -> list[str]:
    errors = gamma_input_contract_errors(input_path)
    stamp = cp2k_stamp_path(output)
    if parse_cp2k_energy(output) is None:
        errors.append("completed and converged CP2K energy is missing")
    if not stamp.is_file():
        errors.append("CP2K production stamp is missing")
        return errors
    try:
        payload = json.loads(stamp.read_text())
    except (json.JSONDecodeError, OSError):
        errors.append("CP2K production stamp is unreadable")
        return errors
    expected = {
        "campaign_id": identity.campaign_id,
        "method": "GXTB",
        "mesh": "gamma",
        "phase": phase,
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "cp2k_sha256": identity.cp2k_sha256,
        "cp2k_library_sha256": identity.cp2k_library_sha256,
        "tblite_static_library_sha256": identity.tblite_static_library_sha256,
        "cp2k_source_revision": identity.cp2k_source_revision,
        "tblite_source_revision": identity.tblite_source_revision,
    }
    for key, value in expected.items():
        if payload.get(key) != value:
            errors.append(f"CP2K stamp {key} mismatch")
    if input_path.is_file() and payload.get("input_sha256") != sha256(input_path):
        errors.append("CP2K stamp input_sha256 mismatch")
    if output.is_file() and payload.get("output_sha256") != sha256(output):
        errors.append("CP2K stamp output_sha256 mismatch")
    if payload.get("input_contract_valid") is not True:
        errors.append("CP2K stamp does not certify the input contract")
    if payload.get("adopted_existing_output") is not False:
        errors.append("CP2K output was adopted rather than produced by this executable")
    return errors


def cli_stamp_path(run_dir: Path) -> Path:
    return run_dir / "tblite.run.json"


def cli_stamp_valid(
    phase: str,
    run_dir: Path,
    expected_poscar_sha256: str,
    identity: BuildIdentity,
) -> bool:
    poscar = run_dir / "POSCAR"
    result_json = run_dir / "tblite.json"
    output = run_dir / "tblite.out"
    stamp = cli_stamp_path(run_dir)
    if (
        not poscar.is_file()
        or not result_json.is_file()
        or not output.is_file()
        or not stamp.is_file()
        or parse_cli_energy(result_json) is None
        or parse_cli_atomic_energies(result_json) is None
    ):
        return False
    try:
        payload = json.loads(stamp.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    return (
        payload.get("campaign_id") == identity.campaign_id
        and payload.get("cli_protocol_id") == CLI_PROTOCOL_ID
        and payload.get("gxtb_protocol_id") == GXTB_PROTOCOL_ID
        and payload.get("phase") == phase
        and payload.get("method") == "gxtb"
        and payload.get("accuracy") == ACCURACY
        and payload.get("tblite_sha256") == identity.tblite_sha256
        and payload.get("tblite_static_library_sha256")
        == identity.tblite_static_library_sha256
        and payload.get("tblite_source_revision") == identity.tblite_source_revision
        and payload.get("poscar_sha256") == expected_poscar_sha256
        and sha256(poscar) == expected_poscar_sha256
        and payload.get("result_sha256") == sha256(result_json)
        and payload.get("output_sha256") == sha256(output)
    )


def write_cli_stamp(
    phase: str,
    run_dir: Path,
    identity: BuildIdentity,
    command: list[str],
) -> None:
    poscar = run_dir / "POSCAR"
    result_json = run_dir / "tblite.json"
    output = run_dir / "tblite.out"
    payload = {
        "campaign_id": identity.campaign_id,
        "cli_protocol_id": CLI_PROTOCOL_ID,
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "phase": phase,
        "method": "gxtb",
        "accuracy": ACCURACY,
        "command": command,
        "tblite": str(identity.tblite),
        "tblite_sha256": identity.tblite_sha256,
        "tblite_static_library": str(identity.tblite_static_library),
        "tblite_static_library_sha256": identity.tblite_static_library_sha256,
        "tblite_source_revision": identity.tblite_source_revision,
        "poscar_sha256": sha256(poscar),
        "result_sha256": sha256(result_json),
        "output_sha256": sha256(output),
    }
    cli_stamp_path(run_dir).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--campaign-manifest",
        type=Path,
        default=DEFAULT_CAMPAIGN_MANIFEST,
    )
    parser.add_argument(
        "--cp2k",
        type=Path,
        required=True,
        help="current CP2K executable whose hash must match every production stamp",
    )
    parser.add_argument(
        "--cp2k-library",
        type=Path,
        help="optional expected libcp2k dylib; must equal the otool/RPATH result",
    )
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--tblite", type=Path, required=True)
    parser.add_argument("--tblite-static-library", type=Path, required=True)
    parser.add_argument("--tblite-source", type=Path, required=True)
    parser.add_argument(
        "--cp2k-input-root",
        type=Path,
        default=ROOT / GXTB_INPUT_DIRECTORY,
        help=f"validated input root (default: ROOT/{GXTB_INPUT_DIRECTORY})",
    )
    parser.add_argument(
        "--cp2k-run-root",
        type=Path,
        default=ROOT / GXTB_RUN_DIRECTORY,
        help=f"validated CP2K run root (default: ROOT/{GXTB_RUN_DIRECTORY})",
    )
    parser.add_argument(
        "--cli-run-root",
        type=Path,
        default=ROOT / CLI_RUN_DIRECTORY,
        help=f"isolated stamped CLI run root (default: ROOT/{CLI_RUN_DIRECTORY})",
    )
    parser.add_argument("--phase", action="append", choices=PHASES)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "data" / "dmc_ice13_gxtb_spglib_gamma_cli_check.csv",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=ROOT
        / "data"
        / "dmc_ice13_gxtb_spglib_gamma_cli_provenance.json",
    )
    parser.add_argument(
        "--atomic-csv",
        type=Path,
        default=ROOT / "data" / "dmc_ice13_gxtb_spglib_atomic_energies.csv",
    )
    args = parser.parse_args()
    phases = args.phase or list(PHASES)
    try:
        args.campaign_manifest = args.campaign_manifest.resolve()
        campaign_manifest = read_campaign_manifest(args.campaign_manifest)
        save_tblite_manifest = campaign_manifest.get("save_tblite")
        if not isinstance(save_tblite_manifest, dict):
            raise ValueError("campaign manifest lacks save_tblite identity block")
        identity = build_identity(
            str(campaign_manifest.get("campaign_id")),
            args.cp2k,
            args.cp2k_library,
            args.tblite,
            args.tblite_static_library,
            args.cp2k_source,
            args.tblite_source,
            str(save_tblite_manifest.get("revision")),
        )
        validate_campaign_identity(identity, campaign_manifest)
    except ValueError as error:
        parser.error(f"build identity gate failed: {error}")
    cp2k = identity.cp2k
    tblite = identity.tblite
    args.cp2k_input_root = args.cp2k_input_root.resolve()
    args.cp2k_run_root = args.cp2k_run_root.resolve()
    args.cli_run_root = args.cli_run_root.resolve()
    geometries = json.loads((ROOT / "data" / "geometries.json").read_text())
    rows: dict[str, dict[str, str]] = {}
    atomic_rows: dict[tuple[str, str], dict[str, str]] = {}
    failures: set[str] = set()

    # Refuse to start any save_tblite calculation unless every requested CP2K
    # reference is a current, hash-validated V1 Gamma production result.
    cp2k_records: dict[str, dict[str, Path]] = {}
    preflight_failures: list[str] = []
    for phase in phases:
        input_path = (
            args.cp2k_input_root / "gamma" / f"ice_{phase}_GXTB_gamma.inp"
        )
        output = (
            args.cp2k_run_root
            / "gamma"
            / phase
            / f"ice_{phase}_GXTB_gamma.out"
        )
        errors = cp2k_validation_errors(phase, input_path, output, identity)
        if errors:
            preflight_failures.append(f"{phase}: {'; '.join(errors)}")
        cp2k_records[phase] = {"input": input_path, "output": output}
    if preflight_failures:
        raise SystemExit(
            "Invalid CP2K Gamma production references; no CLI jobs were run:\n  "
            + "\n  ".join(preflight_failures)
        )

    cli_manifest: list[dict[str, object]] = []

    for phase in phases:
        run_dir = args.cli_run_root / phase
        run_dir.mkdir(parents=True, exist_ok=True)
        poscar = run_dir / "POSCAR"
        labels = coordinate_labels(ROOT / "inputs" / f"ice_{phase}_GFN2.inp")
        expected_poscar = poscar_text(phase, geometries[phase], labels)
        expected_poscar_sha256 = hashlib.sha256(expected_poscar.encode()).hexdigest()
        result_json = run_dir / "tblite.json"
        output = run_dir / "tblite.out"
        old_cli_is_valid = cli_stamp_valid(
            phase,
            run_dir,
            expected_poscar_sha256,
            identity,
        )
        archived: list[str] = []
        returncode = 0
        command = [
            str(tblite),
            "run",
            "--method",
            "gxtb",
            "--acc",
            str(ACCURACY),
            "--no-restart",
            "--json",
            result_json.name,
            poscar.name,
        ]
        if args.force or not old_cli_is_valid:
            for artifact in (poscar, result_json, output, cli_stamp_path(run_dir)):
                stale = archive_stale_file(artifact)
                if stale is not None:
                    archived.append(str(stale.resolve()))
            poscar.write_text(expected_poscar)
            env = os.environ.copy()
            env.update(
                {
                    "OMP_NUM_THREADS": "1",
                    "OPENBLAS_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1",
                    "VECLIB_MAXIMUM_THREADS": "1",
                    "OMP_WAIT_POLICY": "PASSIVE",
                }
            )
            with output.open("w") as handle:
                process = subprocess.run(
                    command,
                    cwd=run_dir,
                    env=env,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            returncode = process.returncode
            if (
                returncode == 0
                and parse_cli_energy(result_json) is not None
                and parse_cli_atomic_energies(result_json) is not None
            ):
                write_cli_stamp(
                    phase,
                    run_dir,
                    identity,
                    command,
                )
            elif cli_stamp_path(run_dir).is_file():
                cli_stamp_path(run_dir).unlink()

        cli_energy = parse_cli_energy(result_json)
        cli_is_valid = cli_stamp_valid(
            phase,
            run_dir,
            expected_poscar_sha256,
            identity,
        )

        cp2k_output = cp2k_records[phase]["output"]
        cp2k_energy = parse_cp2k_energy(cp2k_output)
        nwater = int(geometries[phase]["counts"]["O"])
        difference = (
            cp2k_energy - cli_energy
            if cp2k_energy is not None and cli_energy is not None
            else None
        )
        completed = (
            returncode == 0
            and cli_is_valid
            and cli_energy is not None
            and cp2k_energy is not None
            and difference is not None
            and math.isfinite(difference)
            and abs(difference) <= args.tolerance
        )
        rows[phase] = {
            "phase": phase,
            "gxtb_protocol_id": GXTB_PROTOCOL_ID,
            "cli_protocol_id": CLI_PROTOCOL_ID,
            "n_H2O": str(nwater),
            "cp2k_energy_hartree": (
                f"{cp2k_energy:.15f}" if cp2k_energy is not None else ""
            ),
            "save_tblite_cli_energy_hartree": (
                f"{cli_energy:.15f}" if cli_energy is not None else ""
            ),
            "cp2k_minus_cli_hartree": (
                f"{difference:.15e}" if difference is not None else ""
            ),
            "abs_difference_hartree": (
                f"{abs(difference):.15e}" if difference is not None else ""
            ),
            "abs_difference_per_H2O_hartree": (
                f"{abs(difference) / nwater:.15e}"
                if difference is not None
                else ""
            ),
            "completed": str(completed),
            "tolerance_hartree": f"{args.tolerance:.8e}",
            "tolerance_exceeded": str(
                difference is None
                or not math.isfinite(difference)
                or abs(difference) > args.tolerance
            ),
            "cli_returncode": str(returncode),
            "cp2k_hash_validated": "True",
            "cli_hash_validated": str(cli_is_valid),
        }
        print(
            f"{phase:4s} CP2K={cp2k_energy!s:>22s} CLI={cli_energy!s:>22s} "
            f"Delta={difference!s}",
            flush=True,
        )
        if not completed or difference is None or not math.isfinite(difference):
            failures.add(phase)

        per_atom = parse_cli_atomic_energies(result_json)
        if per_atom is None or len(per_atom) != len(labels):
            failures.add(f"{phase}/atomic_energies")
        else:
            for element in dict.fromkeys(labels):
                values = [energy for label, energy in zip(labels, per_atom) if label == element]
                atomic_rows[(phase, element)] = {
                    "phase": phase,
                    "element": element,
                    "N": str(len(values)),
                    "sum_hartree": f"{sum(values):.15f}",
                    "mean_hartree": f"{sum(values) / len(values):.15f}",
                    "min_hartree": f"{min(values):.15f}",
                    "max_hartree": f"{max(values):.15f}",
                }
        cli_manifest.append(
            {
                "phase": phase,
                "run_directory": str(run_dir.resolve()),
                "reused_hash_validated_result": old_cli_is_valid and not args.force,
                "archived_incompatible_artifacts": archived,
                "command": command,
                "returncode": returncode,
                "stamp": str(cli_stamp_path(run_dir).resolve()),
                "stamp_sha256": (
                    sha256(cli_stamp_path(run_dir))
                    if cli_stamp_path(run_dir).is_file()
                    else None
                ),
                "completed_and_hash_validated": cli_is_valid,
            }
        )

    fieldnames = [
        "phase",
        "gxtb_protocol_id",
        "cli_protocol_id",
        "n_H2O",
        "cp2k_energy_hartree",
        "save_tblite_cli_energy_hartree",
        "cp2k_minus_cli_hartree",
        "abs_difference_hartree",
        "abs_difference_per_H2O_hartree",
        "completed",
        "tolerance_hartree",
        "tolerance_exceeded",
        "cli_returncode",
        "cp2k_hash_validated",
        "cli_hash_validated",
    ]
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for phase in PHASES:
            if phase in rows:
                writer.writerow(rows[phase])

    atomic_fields = [
        "phase",
        "element",
        "N",
        "sum_hartree",
        "mean_hartree",
        "min_hartree",
        "max_hartree",
    ]
    args.atomic_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.atomic_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=atomic_fields, lineterminator="\n")
        writer.writeheader()
        for phase in PHASES:
            for element in ("H", "O"):
                if (phase, element) in atomic_rows:
                    writer.writerow(atomic_rows[(phase, element)])

    provenance = {
        "benchmark": "DMC-ICE13 Gamma CP2K/save_tblite validation",
        "campaign": {
            "id": identity.campaign_id,
            "manifest": str(args.campaign_manifest),
            "manifest_sha256": sha256(args.campaign_manifest),
        },
        "method": "GXTB",
        "gxtb_protocol_id": GXTB_PROTOCOL_ID,
        "cli_protocol_id": CLI_PROTOCOL_ID,
        "accuracy": ACCURACY,
        "tolerance_hartree": args.tolerance,
        "invocation": [sys.executable, *sys.argv],
        "coverage": {
            "requested_phases": phases,
            "expected_this_invocation": len(phases),
            "expected_full_benchmark": len(PHASES),
            "completed": sum(
                rows.get(phase, {}).get("completed") == "True" for phase in phases
            ),
        },
        "cp2k": {
            "executable": str(cp2k),
            "sha256": identity.cp2k_sha256,
            "loaded_library": str(identity.cp2k_library),
            "loaded_library_sha256": identity.cp2k_library_sha256,
            "source_revision_validated": identity.cp2k_source_revision,
            "version": command_output([str(cp2k), "--version"]),
            "source": git_metadata(args.cp2k_source),
            "input_root": str(args.cp2k_input_root),
            "run_root": str(args.cp2k_run_root),
            "references": [
                {
                    "phase": phase,
                    "input": str(cp2k_records[phase]["input"].resolve()),
                    "input_sha256": sha256(cp2k_records[phase]["input"]),
                    "output": str(cp2k_records[phase]["output"].resolve()),
                    "output_sha256": sha256(cp2k_records[phase]["output"]),
                    "stamp": str(
                        cp2k_stamp_path(cp2k_records[phase]["output"]).resolve()
                    ),
                    "stamp_sha256": sha256(
                        cp2k_stamp_path(cp2k_records[phase]["output"])
                    ),
                    "completed_and_hash_validated": True,
                }
                for phase in phases
            ],
        },
        "save_tblite": {
            "executable": str(tblite),
            "sha256": identity.tblite_sha256,
            "version": command_output([str(tblite), "--version"]),
            "source": git_metadata(args.tblite_source),
            "source_revision_validated": identity.tblite_source_revision,
            "static_library": str(identity.tblite_static_library),
            "static_library_sha256": identity.tblite_static_library_sha256,
            "run_root": str(args.cli_run_root),
            "run_manifest": cli_manifest,
        },
        "result_csv": {
            "path": str(args.csv.resolve()),
            "sha256": sha256(args.csv),
        },
        "atomic_energy_csv": {
            "path": str(args.atomic_csv.resolve()),
            "sha256": sha256(args.atomic_csv),
        },
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

    if failures:
        raise SystemExit(
            "Incomplete CLI comparison for: " + ", ".join(sorted(failures))
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
