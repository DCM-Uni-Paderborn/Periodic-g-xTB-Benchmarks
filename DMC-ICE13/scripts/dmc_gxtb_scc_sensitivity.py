#!/usr/bin/env python3
"""Check whether periodic g-xTB ice energies depend on guess or SCC mixer."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from pathlib import Path

import dmc_gxtb_gamma_cli_check as gamma_check


ROOT = gamma_check.ROOT
CONFIGURATIONS = (
    ("guess_sad", ("--guess", "sad")),
    ("guess_eeq", ("--guess", "eeq")),
    ("guess_eeqbc", ("--guess", "eeqbc")),
    ("guess_ceh", ("--guess", "ceh")),
    ("mixer_broyden", ("--guess", "eeqbc", "--broyden-start", "1", "--diis-start", "999")),
    ("mixer_diis4", ("--guess", "eeqbc", "--broyden-start", "999", "--diis-start", "4")),
    ("mixer_diis10", ("--guess", "eeqbc", "--broyden-start", "999", "--diis-start", "10")),
)


def environment() -> dict[str, str]:
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
    return env


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tblite", type=Path, required=True)
    parser.add_argument("--tblite-source", type=Path)
    parser.add_argument("--phase", action="append", choices=gamma_check.PHASES)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "data" / "dmc_ice13_gxtb_scc_sensitivity.csv",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=ROOT / "data" / "dmc_ice13_gxtb_scc_sensitivity_provenance.json",
    )
    args = parser.parse_args()
    phases = args.phase or ["VII"]
    tblite = args.tblite.resolve()
    geometries = json.loads((ROOT / "data" / "geometries.json").read_text())
    rows: list[dict[str, object]] = []
    failures: list[str] = []

    for phase in phases:
        run_dir = ROOT / "runs_cli" / "GXTB" / phase
        run_dir.mkdir(parents=True, exist_ok=True)
        poscar = run_dir / "POSCAR"
        labels = gamma_check.coordinate_labels(
            ROOT / "kpoint_inputs" / "gamma" / f"ice_{phase}_GXTB_gamma.inp"
        )
        poscar.write_text(
            gamma_check.poscar_text(phase, geometries[phase], labels),
            encoding="utf-8",
        )
        energies: dict[str, float] = {}
        returncodes: dict[str, int] = {}

        for name, options in CONFIGURATIONS:
            json_name = (
                f"root_{name.removeprefix('guess_')}.json"
                if name.startswith("guess_")
                else f"mixer_{name.removeprefix('mixer_')}.json"
            )
            result_json = run_dir / json_name
            energy = None if args.force else gamma_check.parse_cli_energy(result_json)
            returncode = 0
            if energy is None:
                command = [
                    str(tblite),
                    "run",
                    "--method",
                    "gxtb",
                    "--acc",
                    "0.1",
                    "--no-restart",
                    *options,
                    "--json",
                    json_name,
                    poscar.name,
                ]
                process = subprocess.run(
                    command,
                    cwd=run_dir,
                    env=environment(),
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
                returncode = process.returncode
                energy = gamma_check.parse_cli_energy(result_json)
            returncodes[name] = returncode
            if energy is not None:
                energies[name] = energy

        baseline = energies.get("guess_eeqbc")
        for name, options in CONFIGURATIONS:
            energy = energies.get(name)
            delta = energy - baseline if energy is not None and baseline is not None else None
            completed = returncodes[name] == 0 and energy is not None and delta is not None
            rows.append(
                {
                    "phase": phase,
                    "configuration": name,
                    "cli_options": " ".join(options),
                    "energy_hartree": f"{energy:.15f}" if energy is not None else "",
                    "delta_from_eeqbc_hartree": f"{delta:.15e}" if delta is not None else "",
                    "abs_delta_hartree": f"{abs(delta):.15e}" if delta is not None else "",
                    "completed": completed,
                    "returncode": returncodes[name],
                }
            )
            if not completed:
                failures.append(f"{phase}/{name}")
        spread = max(energies.values()) - min(energies.values()) if energies else float("inf")
        if spread > args.tolerance:
            failures.append(f"{phase}/energy_spread={spread:.6e}")
        print(f"{phase}: SCC guess/mixer energy spread = {spread:.12e} hartree", flush=True)

    columns = [
        "phase",
        "configuration",
        "cli_options",
        "energy_hartree",
        "delta_from_eeqbc_hartree",
        "abs_delta_hartree",
        "completed",
        "returncode",
    ]
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    try:
        result_path = args.csv.resolve().relative_to(ROOT.resolve())
    except ValueError:
        result_path = args.csv.resolve()
    provenance = {
        "benchmark": "DMC-ICE13 periodic g-xTB SCC-root sensitivity",
        "phases": phases,
        "accuracy": 0.1,
        "energy_spread_tolerance_hartree": args.tolerance,
        "configurations": {name: list(options) for name, options in CONFIGURATIONS},
        "save_tblite": {
            "executable": str(tblite),
            "sha256": gamma_check.sha256(tblite),
            "version": gamma_check.command_output([str(tblite), "--version"]),
            "source": gamma_check.git_metadata(args.tblite_source),
        },
        "result_csv": {
            "path": str(result_path),
            "sha256": gamma_check.sha256(args.csv),
        },
    }
    args.provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    if failures:
        raise SystemExit("SCC sensitivity validation failed: " + ", ".join(failures))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
