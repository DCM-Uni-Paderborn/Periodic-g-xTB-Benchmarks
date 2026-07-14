#!/usr/bin/env python3
"""Run method-owned CP2K density/Fock mixer diagnostics for periodic g-xTB."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import re
import subprocess
from pathlib import Path

import dmc_gxtb_gamma_cli_check as common


ROOT = common.ROOT
MESHES = ("gamma", "k111", "k222", "k333", "k444", "k555")
VARIANTS = ("cp2k_density", "cp2k_fock")


def diagnostic_input(source: str, phase: str, mesh: str, variant: str) -> str:
    if source.count("SCC_MIXER TBLITE") != 1:
        raise ValueError("Expected exactly one native g-xTB SCC_MIXER selection")
    text = source.replace("SCC_MIXER TBLITE", "SCC_MIXER CP2K", 1)
    project = f"ice_{phase}_GXTB_{mesh}"
    text = text.replace(f"PROJECT {project}", f"PROJECT {project}_{variant}", 1)
    scf_marker = "    &SCF\n"
    if text.count(scf_marker) != 1:
        raise ValueError("Cannot locate unique SCF section")
    if variant == "cp2k_density":
        text = text.replace(scf_marker, scf_marker + "      EPS_DIIS 0.0\n", 1)
    elif variant == "cp2k_fock":
        text = text.replace(
            scf_marker,
            scf_marker + "      EPS_DIIS 0.1\n      MAX_DIIS 7\n",
            1,
        )
        text = text.replace("        ALPHA 0.2\n", "        ALPHA 1.0\n        NMIXING 2\n", 1)
        if mesh != "gamma":
            kpoint_end = "    &END KPOINTS\n"
            if text.count(kpoint_end) != 1:
                raise ValueError("Cannot locate unique KPOINTS section")
            text = text.replace(kpoint_end, "      WAVEFUNCTIONS COMPLEX\n" + kpoint_end, 1)
    else:
        raise ValueError(f"Unknown mixer diagnostic {variant}")
    return text


def output_completed(path: Path) -> bool:
    if not path.is_file():
        return False
    text = path.read_text(errors="ignore")
    return (
        "PROGRAM ENDED" in text
        and "SCF run converged" in text
        and "SCF run NOT converged" not in text
        and "ABORT" not in text
        and "ENERGY| Total FORCE_EVAL" in text
    )


def parse_output(path: Path) -> tuple[float | None, int | None]:
    if not output_completed(path):
        return None, None
    text = path.read_text(errors="ignore")
    energies = re.findall(
        r"ENERGY\| Total FORCE_EVAL.*?([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*$",
        text,
        flags=re.MULTILINE,
    )
    steps = re.findall(r"SCF run converged in\s+(\d+) steps", text)
    return (float(energies[-1]) if energies else None, int(steps[-1]) if steps else None)


def stamp_valid(stamp: Path, inp: Path, output: Path, cp2k_sha: str) -> bool:
    if not output_completed(output) or not stamp.is_file():
        return False
    try:
        payload = json.loads(stamp.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    return (
        payload.get("input_sha256") == common.sha256(inp)
        and payload.get("cp2k_sha256") == cp2k_sha
        and payload.get("output_sha256") == common.sha256(output)
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


def run_case(
    cp2k: Path,
    cp2k_sha: str,
    phase: str,
    mesh: str,
    variant: str,
    force: bool,
) -> dict[str, object]:
    source = ROOT / "kpoint_inputs" / mesh / f"ice_{phase}_GXTB_{mesh}.inp"
    run_dir = ROOT / "runs_mixer_diagnostics" / mesh / variant / phase
    run_dir.mkdir(parents=True, exist_ok=True)
    inp = run_dir / f"ice_{phase}_GXTB_{mesh}_{variant}.inp"
    inp.write_text(diagnostic_input(source.read_text(), phase, mesh, variant))
    output = run_dir / f"ice_{phase}_GXTB_{mesh}_{variant}.out"
    stamp = run_dir / "run_stamp.json"
    if force:
        output.unlink(missing_ok=True)
        stamp.unlink(missing_ok=True)
    if stamp_valid(stamp, inp, output, cp2k_sha):
        returncode = 0
    else:
        process = subprocess.run(
            [str(cp2k), "-i", inp.name, "-o", output.name],
            cwd=run_dir,
            env=environment(),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            check=False,
        )
        returncode = process.returncode
        if returncode == 0 and output_completed(output):
            stamp.write_text(
                json.dumps(
                    {
                        "input_sha256": common.sha256(inp),
                        "cp2k_sha256": cp2k_sha,
                        "output_sha256": common.sha256(output),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
        else:
            stamp.unlink(missing_ok=True)
    energy, steps = parse_output(output)
    return {
        "phase": phase,
        "mesh": mesh,
        "variant": variant,
        "production_role": "diagnostic_only",
        "completed": output_completed(output),
        "returncode": returncode,
        "scf_steps": steps if steps is not None else "",
        "energy_hartree": f"{energy:.15f}" if energy is not None else "",
        "input_sha256": common.sha256(inp),
        "output_sha256": common.sha256(output) if output_completed(output) else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path)
    parser.add_argument("--phase", action="append", choices=common.PHASES)
    parser.add_argument("--mesh", choices=MESHES, default="k333")
    parser.add_argument("--variant", action="append", choices=VARIANTS)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "data" / "dmc_ice13_gxtb_mixer_diagnostics.csv",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=ROOT / "data" / "dmc_ice13_gxtb_mixer_diagnostics_provenance.json",
    )
    args = parser.parse_args()
    phases = args.phase or ["III"]
    variants = args.variant or list(VARIANTS)
    cp2k = args.cp2k.resolve()
    cp2k_sha = common.sha256(cp2k)
    cases = [(phase, args.mesh, variant) for phase in phases for variant in variants]
    rows: list[dict[str, object]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(run_case, cp2k, cp2k_sha, phase, mesh, variant, args.force): (phase, variant)
            for phase, mesh, variant in cases
        }
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            rows.append(row)
            print(
                f"{row['phase']} {row['mesh']} {row['variant']}: "
                f"completed={row['completed']} steps={row['scf_steps']} energy={row['energy_hartree']}",
                flush=True,
            )
    rows.sort(key=lambda row: (common.PHASES.index(str(row["phase"])), VARIANTS.index(str(row["variant"]))))
    columns = [
        "phase",
        "mesh",
        "variant",
        "production_role",
        "completed",
        "returncode",
        "scf_steps",
        "energy_hartree",
        "input_sha256",
        "output_sha256",
    ]
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance = {
        "benchmark": "DMC-ICE13 periodic g-xTB mixer diagnostics",
        "production_role": "diagnostic_only",
        "cp2k": {
            "executable": str(cp2k),
            "sha256": cp2k_sha,
            "version": common.command_output([str(cp2k), "--version"]),
            "source": common.git_metadata(args.cp2k_source),
        },
        "result_csv": {"path": str(args.csv.resolve()), "sha256": common.sha256(args.csv)},
    }
    args.provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    failed = [row for row in rows if not row["completed"]]
    if failed:
        raise SystemExit("Mixer diagnostics failed: " + ", ".join(f"{r['phase']}/{r['variant']}" for r in failed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
