#!/usr/bin/env python3
"""Run every archived Part-I implementation gate as one reproducible audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parents[2]
GATES = (
    (
        "absolute_energy_parity",
        ROOT,
        ROOT / "tools/verify_absolute_energy_parity.py",
        (),
    ),
    (
        "accuracy_sensitivity",
        ROOT / "validation/accuracy_sensitivity_20260718",
        ROOT / "validation/accuracy_sensitivity_20260718/verify_accuracy.py",
        (),
    ),
    (
        "response_fix",
        ROOT / "validation/cp2k_response_fix_ab_20260719",
        ROOT / "validation/cp2k_response_fix_ab_20260719/verify_response_fix.py",
        (),
    ),
    (
        "energy_force_stress",
        ROOT / "validation/k222_XVII_derivatives",
        ROOT / "validation/k222_XVII_derivatives/verify_derivatives.py",
        (),
    ),
    (
        "kpoint_grid_bvk",
        ROOT / "validation/kpoint_grid_bvk_oracle_20260718",
        ROOT / "validation/kpoint_grid_bvk_oracle_20260718/verify_grid.py",
        (),
    ),
    (
        "model_revision",
        ROOT / "validation/model_revision_coarse_grid_ab_20260718",
        ROOT / "validation/model_revision_coarse_grid_ab_20260718/verify_comparison.py",
        (),
    ),
    (
        "native_derivative_hardening",
        ROOT / "validation/native_derivative_hardening_20260718",
        ROOT / "validation/native_derivative_hardening_20260718/verify_current_build.py",
        (),
    ),
    (
        "provider_revision",
        ROOT / "validation/provider_revision_bvk_ab_20260718",
        ROOT / "validation/provider_revision_bvk_ab_20260718/verify_comparison.py",
        (),
    ),
    (
        "qualified_energy_sentinels",
        ROOT / "validation/qualified_energy_sentinels_20260719",
        ROOT
        / "validation/qualified_energy_sentinels_20260719/verify_qualified_energy_sentinels.py",
        (),
    ),
    (
        "wigner_seitz_branch_diagnosis",
        ROOT / "validation/wigner_seitz_branch_diagnosis_20260718",
        ROOT / "validation/wigner_seitz_branch_diagnosis_20260718/verify_diagnosis.py",
        (),
    ),
    (
        "final_lowk_derivatives",
        REPOSITORY_ROOT / "validation/gxtb_final_lowk_derivatives_20260719",
        REPOSITORY_ROOT
        / "validation/gxtb_final_lowk_derivatives_20260719/verify_final_lowk_derivatives.py",
        (
            ".",
            "--legacy-partial-manifest",
            "legacy_partial_pbc_manifest.json",
            "--output",
            "verification.json",
        ),
    ),
    (
        "phase_viii_component_ablation",
        REPOSITORY_ROOT / "validation/dmc13_k222_viii_component_ablation_20260719",
        REPOSITORY_ROOT
        / "validation/dmc13_k222_viii_component_ablation_20260719/verify_component_ablation.py",
        (),
    ),
    (
        "phase_xvii_derivative_component_ablation",
        REPOSITORY_ROOT
        / "validation/dmc13_k222_xvii_derivative_component_ablation_20260719",
        REPOSITORY_ROOT
        / "validation/dmc13_k222_xvii_derivative_component_ablation_20260719/verify_derivative_component_ablation.py",
        (),
    ),
    (
        "provider_component_attribution",
        REPOSITORY_ROOT / "validation/provider_component_attribution_20260719",
        REPOSITORY_ROOT
        / "validation/provider_component_attribution_20260719/verify_attribution.py",
        (),
    ),
    (
        "pbc_h0_anisotropy_attribution",
        REPOSITORY_ROOT / "validation/pbc_h0_anisotropy_attribution_20260719",
        REPOSITORY_ROOT
        / "validation/pbc_h0_anisotropy_attribution_20260719/verify_h0_attribution.py",
        (),
    ),
    (
        "save_tblite_periodic_source_tests",
        REPOSITORY_ROOT / "validation/save_tblite_periodic_source_tests_20260719",
        REPOSITORY_ROOT
        / "validation/save_tblite_periodic_source_tests_20260719/verify_source_tests.py",
        (),
    ),
    (
        "three_route_k333_closure",
        REPOSITORY_ROOT / "validation/three_route_k333_closure_20260719",
        REPOSITORY_ROOT
        / "validation/three_route_k333_closure_20260719/generate_and_verify.py",
        (),
    ),
    (
        "seidler_recalculation_package",
        ROOT.parent / "seidler_dmc13_recalculation",
        ROOT.parent / "seidler_dmc13_recalculation/prepare_package.py",
        (),
    ),
    (
        "archive_sha256",
        ROOT,
        ROOT / "tools/verify_sha256_manifests.py",
        (str(ROOT),),
    ),
)


def digest_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def concise_tail(value: str) -> str:
    lines = [line.strip() for line in value.splitlines() if line.strip()]
    return lines[-1] if lines else ""


def run_gate(
    name: str, workdir: Path, script: Path, arguments: tuple[str, ...]
) -> dict[str, object]:
    if not script.is_file():
        raise AssertionError(f"missing gate script: {script}")
    completed = subprocess.run(
        [sys.executable, str(script), *arguments],
        cwd=workdir,
        text=True,
        capture_output=True,
        check=False,
    )
    return {
        "name": name,
        "status": "PASS" if completed.returncode == 0 else "FAIL",
        "returncode": completed.returncode,
        "script": str(script.relative_to(REPOSITORY_ROOT)),
        "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "stdout_sha256": digest_text(completed.stdout),
        "stderr_sha256": digest_text(completed.stderr),
        "stdout_tail": concise_tail(completed.stdout),
        "stderr_tail": concise_tail(completed.stderr),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    rows = [run_gate(*gate) for gate in GATES]
    failed = [row["name"] for row in rows if row["status"] != "PASS"]
    payload = {
        "status": "PASS" if not failed else "FAIL",
        "gate_count": len(rows),
        "passed_gate_count": len(rows) - len(failed),
        "failed_gates": failed,
        "gates": rows,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
