#!/usr/bin/env python3
"""Run a diagnostic-only g-XTB EOS root-continuation pair.

The target calculation receives CP2K's converged Bloch-orbital restart from
the source cell.  This transfers only the orbital/density guess (and hence the
initial Mulliken populations and multipoles reconstructed from that density).
The internal save_tblite q/dipole/quadrupole state and FDIIS history are *not*
restartable in the frozen production build and are deliberately reinitialized.

Nothing produced here is eligible for promotion into the LC12 production EOS.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import run_goldzak12_benchmark as base


SCHEMA_VERSION = 1


def scale_tag(scale: float) -> str:
    return f"s{scale:.5f}".replace(".", "p")


def diagnostic_input(
    ref: base.Reference,
    mesh: str,
    scale: float,
    project: str,
    restart: Path | None,
) -> str:
    text = base.solid_input(ref, "GXTB", "ENERGY", mesh, ref.a_exp * scale, project)
    text = text.replace("        &RESTART OFF", "        &RESTART ON", 1)
    if restart is not None:
        text = text.replace(
            "  &DFT\n",
            f"  &DFT\n    WFN_RESTART_FILE_NAME {restart.resolve()}\n",
            1,
        )
        text = text.replace("      SCF_GUESS MOPAC", "      SCF_GUESS RESTART", 1)
    diagnostic_print = """    &PRINT
      &MO
        EIGENVALUES T
        EIGENVECTORS F
        OCCUPATION_NUMBERS T
        OCCUPATION_NUMBERS_STATS T 1.0E-8
        FILENAME __STD_OUT__
        &EACH
          QS_SCF 0
        &END EACH
      &END MO
      &MULLIKEN
        FILENAME __STD_OUT__
        &EACH
          QS_SCF 0
        &END EACH
      &END MULLIKEN
    &END PRINT
"""
    text = text.replace("    &QS\n", diagnostic_print + "    &QS\n", 1)
    base.validate_method_input(text, "GXTB")
    return text


def last_float(text: str, pattern: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.I | re.M)
    return float(matches[-1]) if matches else None


def last_int(text: str, pattern: str) -> int | None:
    matches = re.findall(pattern, text, flags=re.I | re.M)
    return int(matches[-1]) if matches else None


def observables(output: Path) -> dict[str, object]:
    text = output.read_text(errors="ignore")
    restart_read_explicit = any(
        marker in text
        for marker in (
            "WFN_RESTART| Reading restart file",
            "WFN_RESTART| Restart file",
        )
    )
    restart_requested = bool(
        re.search(r"SCF PARAMETERS\s+Density guess:\s+RESTART", text)
    )
    restart_missing_file_fallback = bool(
        re.search(
            r"(?is)User requested to restart the wavefunction.*?"
            r"file does not exist.*?Calculation continues using ATOMIC GUESS",
            text,
        )
    )
    restart_natom_mismatch_fallback = bool(
        re.search(r"(?i)READ RESTART\s*:\s*WARNING\s*:\s*DIFFERENT natom", text)
    )
    restart_read_error = bool(
        re.search(
            r"(?i)(restart file.{0,100}(?:not found|cannot|failed|error)|"
            r"(?:not found|cannot|failed|error).{0,100}restart file)",
            text,
        )
    )
    normal_completion = "PROGRAM ENDED" in text and "ABORT" not in text
    restart_read_implicit = (
        restart_requested
        and normal_completion
        and not restart_missing_file_fallback
        and not restart_natom_mismatch_fallback
        and not restart_read_error
    )
    total = last_float(text, r"^\s*Total energy:\s+([-+0-9.Ee]+)\s*$")
    entropic = last_float(text, r"^\s*Electronic entropic energy:\s+([-+0-9.Ee]+)\s*$")
    mo_records = [
        {
            "index": int(index),
            "eigenvalue_hartree": float(eigenvalue),
            "occupation": float(occupation),
        }
        for index, eigenvalue, occupation in re.findall(
            r"^\s*MO\|\s+(\d+)\s+([-+0-9.Ee]+)\s+[-+0-9.Ee]+\s+([-+0-9.Ee]+)\s*$",
            text,
            flags=re.M,
        )
    ]
    occupation_values = [float(record["occupation"]) for record in mo_records]
    maximum_occupation = max(occupation_values, default=0.0)
    fractional = [
        value
        for value in occupation_values
        if value > 1.0e-8 and value < maximum_occupation - 1.0e-8
    ]
    mulliken_matches = re.findall(
        r"^\s*(\d+)\s+([A-Za-z]{1,3})\s+(\d+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
        text,
        flags=re.M,
    )
    mulliken_by_atom: dict[int, dict[str, object]] = {}
    for atom, element, kind, population, charge in mulliken_matches:
        mulliken_by_atom[int(atom)] = {
            "atom": int(atom),
            "element": element,
            "kind": int(kind),
            "atomic_population": float(population),
            "net_charge": float(charge),
        }
    return {
        "total_energy_hartree_cp2k_label": total,
        "helmholtz_free_energy_hartree": total,
        "total_energy_extrapolated_t0_hartree": last_float(
            text,
            r"^\s*Total energy \(extrapolated to T->0\):\s+([-+0-9.Ee]+)\s*$",
        ),
        "force_eval_energy_hartree": last_float(
            text,
            r"^\s*ENERGY\| Total FORCE_EVAL .*?\[hartree\]\s+([-+0-9.Ee]+)\s*$",
        ),
        "electronic_entropic_energy_hartree": entropic,
        "total_minus_entropic_hartree": (
            total - entropic if total is not None and entropic is not None else None
        ),
        "internal_energy_from_free_minus_minus_ts_hartree": (
            total - entropic if total is not None and entropic is not None else None
        ),
        "fermi_energy_hartree": last_float(
            text, r"^\s*Fermi energy:\s+([-+0-9.Ee]+)\s*$"
        ),
        "electron_count": last_int(text, r"^\s*Number of electrons:\s+(\d+)\s*$"),
        "occupied_orbitals": last_int(
            text, r"^\s*Number of occupied orbitals:\s+(\d+)\s*$"
        ),
        "molecular_orbitals": last_int(
            text, r"^\s*Number of molecular orbitals:\s+(\d+)\s*$"
        ),
        "occupation_summary": {
            "records_across_irreducible_kpoints": len(mo_records),
            "minimum": min(occupation_values, default=None),
            "maximum": max(occupation_values, default=None),
            "fractional_count": len(fractional),
            "fractional_minimum": min(fractional, default=None),
            "fractional_maximum": max(fractional, default=None),
            "total_occupied_above_1e8_per_kpoint_max": last_int(
                text, r"^\s*MO\| Total occupied:\s+(\d+)\s*$"
            ),
        },
        "mulliken_atoms": [mulliken_by_atom[index] for index in sorted(mulliken_by_atom)],
        "scf_steps": last_int(text, r"SCF run converged in\s+(\d+)\s+steps"),
        "wfn_restart_read": restart_read_explicit or restart_read_implicit,
        "wfn_restart_read_explicit_log": restart_read_explicit,
        "wfn_restart_read_confirmation": (
            "explicit_log"
            if restart_read_explicit
            else (
                "normal_completion_with_restart_density_guess_and_no_restart_error"
                if restart_read_implicit
                else None
            )
        ),
        "wfn_restart_evidence": {
            "explicit_log_line": restart_read_explicit,
            "density_guess_restart_reported": restart_requested,
            "normal_program_completion": normal_completion,
            "missing_file_atomic_guess_fallback_detected": restart_missing_file_fallback,
            "natom_mismatch_fallback_detected": restart_natom_mismatch_fallback,
            "restart_read_error_detected": restart_read_error,
        },
        "mo_occupations_printed": bool(mo_records),
        "mulliken_printed": "MULLIKEN POPULATION ANALYSIS" in text.upper()
        and bool(mulliken_by_atom),
    }


def run_one(
    *,
    ref: base.Reference,
    mesh: str,
    scale: float,
    role: str,
    direction: str,
    restart: Path | None,
    cp2k: Path,
    campaign_identity: dict[str, object],
    threads: int,
    root: Path,
) -> tuple[Path, Path, Path, dict[str, object]]:
    tag = scale_tag(scale)
    project = f"{ref.solid}_GXTB_wfn_{direction}_{role}_{mesh}_{tag}"
    run_dir = root / role
    inp = run_dir / f"{project}.inp"
    out = run_dir / f"{project}.out"
    base.write_file(inp, diagnostic_input(ref, mesh, scale, project, restart))
    executable_identity = base.executable_fingerprint(cp2k)
    signature = base.job_signature(
        cp2k,
        inp,
        command_contract={
            "driver": "cp2k",
            "diagnostic": "gxtb_wfn_hysteresis",
            "role": role,
            "direction": direction,
            "omp_threads": threads,
            "production_eligible": False,
        },
        executable_identity=executable_identity,
        campaign_fingerprint=campaign_identity,
    )
    usable = base.output_ok(out)
    if not (usable and base.job_stamp_matches(out, signature)):
        code = base.run_cp2k(cp2k, inp, out, threads)
        usable = base.output_ok(out)
        base.write_job_stamp(out, signature, completed=usable, return_code=code)
        if code != 0 or not usable:
            raise RuntimeError(f"Diagnostic {role} failed; inspect {out}")
    restart_out = run_dir / f"{project}-RESTART.kp"
    if not restart_out.exists():
        raise RuntimeError(f"Diagnostic {role} did not write {restart_out}")
    obs = observables(out)
    if role == "target" and not obs["wfn_restart_read"]:
        raise RuntimeError(f"Target did not confirm reading the WFN restart: {out}")
    return inp, out, restart_out, obs


def artifact(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": base.sha256(path)}


def archive_failure(
    root: Path,
    *,
    stage: str,
    error: Exception,
    solid: str,
    mesh: str,
    source_scale: float,
    target_scale: float,
    campaign_identity: dict[str, object],
) -> Path:
    files = [artifact(path) for path in sorted(root.rglob("*")) if path.is_file()]
    manifest = root / "diagnostic_manifest.json"
    payload = {
        "schema_version": SCHEMA_VERSION,
        "diagnostic": "gxtb_wfn_hysteresis",
        "production_eligible": False,
        "completed": False,
        "failed_stage": stage,
        "error": str(error),
        "solid": solid,
        "mesh": mesh,
        "source_scale": source_scale,
        "target_scale": target_scale,
        "campaign_identity": campaign_identity,
        "transferred_state": "CP2K Bloch-orbital/density WFN guess; initial q/multipoles are reconstructed from that density",
        "not_transferred": "internal save_tblite q/dipole/quadrupole arrays and FDIIS history; these are reinitialized",
        "restart_confirmation_semantics": (
            "CP2K prints WFN_RESTART lines only when SCF/PRINT/RESTART LOG_PRINT_KEY is enabled. "
            "The unchanged diagnostic input therefore records an explicit line when available, otherwise "
            "normal PROGRAM ENDED completion with Density guess: RESTART and no restart-file error."
        ),
        "archived_artifacts": files,
    }
    base.write_file(manifest, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return manifest


def production_energy(solid: str, mesh: str, scale: float) -> tuple[Path, float | None]:
    project = f"{solid}_GXTB_eos_{mesh}_{scale_tag(scale)}"
    path = (
        base.ROOT
        / "runs"
        / "eos"
        / "GXTB"
        / solid
        / mesh
        / scale_tag(scale)
        / f"{project}.out"
    )
    return path, base.parse_energy(path) if path.exists() and base.output_ok(path) else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solid", required=True)
    parser.add_argument("--mesh", default="k444")
    parser.add_argument("--source-scale", type=float, required=True)
    parser.add_argument("--target-scale", type=float, required=True)
    parser.add_argument("--campaign-manifest", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    refs = {ref.solid: ref for ref in base.REFERENCES}
    if args.solid not in refs:
        parser.error(f"Unknown solid {args.solid!r}")
    if args.mesh != "k444":
        parser.error("This diagnostic is intentionally restricted to the LC12 k444 EOS mesh")
    try:
        campaign_identity, paths = base.validated_gxtb_campaign_from_manifest(
            args.campaign_manifest,
            args.cp2k_source,
            args.save_tblite_source,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    cp2k = paths["cp2k"]
    direction = f"{scale_tag(args.source_scale)}_to_{scale_tag(args.target_scale)}"
    root = (
        base.ROOT
        / "runs"
        / "gxtb_wfn_hysteresis"
        / args.solid
        / args.mesh
        / direction
    )
    try:
        seed_inp, seed_out, seed_restart, seed_obs = run_one(
            ref=refs[args.solid],
            mesh=args.mesh,
            scale=args.source_scale,
            role="seed",
            direction=direction,
            restart=None,
            cp2k=cp2k,
            campaign_identity=campaign_identity,
            threads=args.threads,
            root=root,
        )
    except RuntimeError as exc:
        manifest = archive_failure(
            root,
            stage="seed",
            error=exc,
            solid=args.solid,
            mesh=args.mesh,
            source_scale=args.source_scale,
            target_scale=args.target_scale,
            campaign_identity=campaign_identity,
        )
        raise RuntimeError(f"{exc}; archived in {manifest}") from exc
    try:
        target_inp, target_out, target_restart, target_obs = run_one(
            ref=refs[args.solid],
            mesh=args.mesh,
            scale=args.target_scale,
            role="target",
            direction=direction,
            restart=seed_restart,
            cp2k=cp2k,
            campaign_identity=campaign_identity,
            threads=args.threads,
            root=root,
        )
    except RuntimeError as exc:
        manifest = archive_failure(
            root,
            stage="target",
            error=exc,
            solid=args.solid,
            mesh=args.mesh,
            source_scale=args.source_scale,
            target_scale=args.target_scale,
            campaign_identity=campaign_identity,
        )
        raise RuntimeError(f"{exc}; archived in {manifest}") from exc
    source_prod, source_prod_energy = production_energy(args.solid, args.mesh, args.source_scale)
    target_prod, target_prod_energy = production_energy(args.solid, args.mesh, args.target_scale)
    target_energy = target_obs["total_energy_extrapolated_t0_hartree"]
    payload = {
        "schema_version": SCHEMA_VERSION,
        "diagnostic": "gxtb_wfn_hysteresis",
        "production_eligible": False,
        "completed": True,
        "solid": args.solid,
        "mesh": args.mesh,
        "source_scale": args.source_scale,
        "target_scale": args.target_scale,
        "lattice_source_A": refs[args.solid].a_exp * args.source_scale,
        "lattice_target_A": refs[args.solid].a_exp * args.target_scale,
        "temperature_K": 300.0,
        "eps_scf": 1.0e-9,
        "hamiltonian_contract": "METHOD GXTB; ACCURACY 0.05; native save_tblite FDIIS; MACDONALD k444; SPGLIB reduction",
        "transferred_state": "CP2K Bloch-orbital/density WFN guess; initial q/multipoles are reconstructed from that density",
        "not_transferred": "internal save_tblite q/dipole/quadrupole arrays and FDIIS history; these are reinitialized",
        "restart_confirmation_semantics": (
            "CP2K prints WFN_RESTART lines only when SCF/PRINT/RESTART LOG_PRINT_KEY is enabled. "
            "The unchanged diagnostic input therefore records an explicit line when available, otherwise "
            "normal PROGRAM ENDED completion with Density guess: RESTART and no restart-file error."
        ),
        "campaign_identity": campaign_identity,
        "seed_lineage": {
            "input": artifact(seed_inp),
            "output": artifact(seed_out),
            "wfn_restart": artifact(seed_restart),
            "observables": seed_obs,
            "independent_production_output": (
                artifact(source_prod) if source_prod.exists() else None
            ),
            "independent_production_energy_hartree": source_prod_energy,
        },
        "target": {
            "input": artifact(target_inp),
            "output": artifact(target_out),
            "wfn_restart": artifact(target_restart),
            "observables": target_obs,
            "independent_production_output": (
                artifact(target_prod) if target_prod.exists() else None
            ),
            "independent_production_energy_hartree": target_prod_energy,
            "wfn_minus_independent_energy_hartree": (
                target_energy - target_prod_energy
                if isinstance(target_energy, float) and target_prod_energy is not None
                else None
            ),
        },
    }
    manifest = root / "diagnostic_manifest.json"
    base.write_file(manifest, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"ok {args.solid} {args.source_scale:.5f}->{args.target_scale:.5f} "
        f"seed={seed_obs['total_energy_extrapolated_t0_hartree']:.12f} "
        f"target={target_obs['total_energy_extrapolated_t0_hartree']:.12f} "
        f"manifest={manifest}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
