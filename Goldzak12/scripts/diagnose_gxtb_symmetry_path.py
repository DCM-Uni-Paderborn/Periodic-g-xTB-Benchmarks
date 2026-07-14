#!/usr/bin/env python3
"""Archive an isolated g-XTB SPGLIB/full-grid symmetry diagnostic.

This is intentionally outside the LC12 production tree and cannot promote an
output.  ``spglib_reproduction`` repeats the production Hamiltonian and mesh;
``full_grid`` disables reduction while retaining the same shifted k444 mesh.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import diagnose_gxtb_wfn_hysteresis as wfn
import run_goldzak12_benchmark as base


def add_analysis_prints(text: str) -> str:
    block = """    &PRINT
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
    return text.replace("    &QS\n", block + "    &QS\n", 1)


def diagnostic_input(
    ref: base.Reference, mesh: str, scale: float, project: str, variant: str
) -> str:
    text = base.solid_input(ref, "GXTB", "ENERGY", mesh, ref.a_exp * scale, project)
    if variant == "spglib_reproduction":
        base.validate_method_input(text, "GXTB")
        return text
    if variant != "full_grid":
        raise ValueError(f"Unknown variant {variant}")
    text = text.replace("      SYMMETRY T", "      SYMMETRY F", 1)
    text = text.replace("      FULL_GRID F", "      FULL_GRID T", 1)
    text = add_analysis_prints(text)
    required = (
        "METHOD GXTB",
        "SCC_MIXER TBLITE",
        "SCHEME MACDONALD 4 4 4 0.125 0.125 0.125",
        "SYMMETRY F",
        "FULL_GRID T",
        "EPS_SCF 1.0E-9",
        "ACCURACY 0.05",
    )
    missing = [entry for entry in required if entry not in text]
    if missing:
        raise ValueError("Malformed full-grid diagnostic input: " + ", ".join(missing))
    return text


def abort_record(text: str) -> dict[str, object] | None:
    match = re.search(
        r"\[ABORT\]\s+(.*?)\n(?:.*\n){0,5}?.*?residual=\s*([-+0-9.Ee]+)",
        text,
        flags=re.S,
    )
    if match:
        return {"message": " ".join(match.group(1).split()), "residual": float(match.group(2))}
    match = re.search(r"\[ABORT\]\s+(.*)", text)
    return {"message": match.group(1).strip(), "residual": None} if match else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solid", default="MgO")
    parser.add_argument("--mesh", default="k444")
    parser.add_argument("--scale", type=float, default=0.85)
    parser.add_argument(
        "--variant", choices=("spglib_reproduction", "full_grid"), required=True
    )
    parser.add_argument("--campaign-manifest", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.mesh != "k444":
        parser.error("This diagnostic is restricted to the LC12 k444 EOS mesh")
    refs = {ref.solid: ref for ref in base.REFERENCES}
    if args.solid not in refs:
        parser.error(f"Unknown solid {args.solid!r}")
    try:
        campaign_identity, paths = base.validated_gxtb_campaign_from_manifest(
            args.campaign_manifest, args.cp2k_source, args.save_tblite_source
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    cp2k = paths["cp2k"]
    tag = wfn.scale_tag(args.scale)
    root = (
        base.ROOT
        / "runs"
        / "gxtb_symmetry_diagnostic"
        / args.solid
        / args.mesh
        / tag
        / args.variant
    )
    project = f"{args.solid}_GXTB_symdiag_{args.variant}_{args.mesh}_{tag}"
    inp = root / f"{project}.inp"
    out = root / f"{project}.out"
    base.write_file(
        inp, diagnostic_input(refs[args.solid], args.mesh, args.scale, project, args.variant)
    )
    signature = base.job_signature(
        cp2k,
        inp,
        command_contract={
            "driver": "cp2k",
            "diagnostic": "gxtb_symmetry_path",
            "variant": args.variant,
            "omp_threads": args.threads,
            "production_eligible": False,
        },
        campaign_fingerprint=campaign_identity,
    )
    if not base.job_stamp_matches(out, signature):
        code = base.run_cp2k(cp2k, inp, out, args.threads)
        completed = base.output_ok(out)
        base.write_job_stamp(out, signature, completed=completed, return_code=code)
    else:
        code = 0
        completed = True
    text = out.read_text(errors="ignore") if out.exists() else ""
    original = (
        base.ROOT
        / "runs"
        / "eos"
        / "GXTB"
        / args.solid
        / args.mesh
        / tag
        / f"{args.solid}_GXTB_eos_{args.mesh}_{tag}.out"
    )
    payload = {
        "schema_version": 1,
        "diagnostic": "gxtb_symmetry_path",
        "production_eligible": False,
        "solid": args.solid,
        "mesh": args.mesh,
        "scale": args.scale,
        "variant": args.variant,
        "campaign_identity": campaign_identity,
        "same_hamiltonian": True,
        "mesh_contract": (
            "MACDONALD shifted k444, SPGLIB reduction"
            if args.variant == "spglib_reproduction"
            else "same MACDONALD shifted k444 complete mesh, SYMMETRY F/FULL_GRID T"
        ),
        "completed": completed,
        "return_code": code,
        "abort": abort_record(text),
        "observables": wfn.observables(out) if out.exists() else None,
        "input": wfn.artifact(inp),
        "output": wfn.artifact(out) if out.exists() else None,
        "original_production_failure": wfn.artifact(original) if original.exists() else None,
        "fock_residual_stage": (
            "provider whole-mesh exchange Fock after cp2k_exchange_kmesh and Hermiticity, before foldback; "
            "the current SCC state already reflects native FDIIS, but FDIIS history is not separately observable"
        ),
    }
    manifest = root / "diagnostic_manifest.json"
    base.write_file(manifest, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"archived {args.variant}: completed={completed} rc={code} "
        f"abort={payload['abort']} manifest={manifest}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
