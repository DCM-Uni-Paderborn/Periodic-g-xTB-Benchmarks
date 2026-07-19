#!/usr/bin/env python3
"""Recompute the final CLI/native parity and CP2K response-fix A/B gates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
PHASES = (
    "Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI",
    "XIII", "XIV", "XV", "XVII",
)
N_WATER = {
    "Ih": 12, "II": 12, "III": 12, "IV": 16, "VI": 10, "VII": 12,
    "VIII": 8, "IX": 12, "XI": 8, "XIII": 28, "XIV": 12, "XV": 10,
    "XVII": 6,
}
REFERENCE = {
    "II": 0.31, "III": 1.25, "IV": 3.83, "VI": 1.78, "VII": 4.99,
    "VIII": 4.23, "IX": 0.60, "XI": 0.16, "XIII": 2.12, "XIV": 1.70,
    "XV": 1.74, "XVII": 1.75,
}
HARTREE_TO_KJMOL = 2625.4996394799
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(
            f"{label}: actual={actual:.15g}, expected={expected:.15g}, "
            f"tolerance={tolerance:.3g}"
        )


def cp2k_energy(path: Path) -> float:
    if not path.is_file():
        raise AssertionError(f"missing CP2K output: {path}")
    values: list[float] = []
    ended = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ENERGY_RE.match(line)
        if match:
            values.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if not ended or not values:
        raise AssertionError(f"incomplete CP2K output: {path}")
    status = path.with_name("exit_status")
    if status.is_file() and status.read_text().strip() != "0":
        raise AssertionError(f"nonzero exit status: {status}")
    return values[-1]


def cli_energy(path: Path) -> float:
    if not path.is_file():
        raise AssertionError(f"missing CLI result: {path}")
    value = float(json.loads(path.read_text())["energy"])
    if not math.isfinite(value):
        raise AssertionError(f"non-finite CLI energy: {path}")
    status = path.with_name("exit_status")
    if status.is_file() and status.read_text().strip() != "0":
        raise AssertionError(f"nonzero CLI exit status: {status}")
    return value / 8.0


def cp2k_set(directory: str) -> dict[str, float]:
    return {
        phase: cp2k_energy(RAW / directory / phase / "cp2k.out")
        for phase in PHASES
    }


def cli_set(directory: str) -> dict[str, float]:
    return {
        phase: cli_energy(RAW / directory / phase / "result.json")
        for phase in PHASES
    }


def relative(energies: dict[str, float], phase: str) -> float:
    return (
        energies[phase] / N_WATER[phase]
        - energies["Ih"] / N_WATER["Ih"]
    ) * HARTREE_TO_KJMOL


def first_manifest_hash(path: Path) -> str:
    value = path.read_text().splitlines()[0].split()[0]
    if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
        raise AssertionError(f"invalid SHA-256 manifest entry: {path}")
    return value


def verify_no_acp_parameter() -> dict[str, str]:
    parameter_root = RAW / "parameters"
    full_path = parameter_root / "gxtb_full.toml"
    no_acp_path = parameter_root / "gxtb_no_acp.toml"
    full = tomllib.loads(full_path.read_text(encoding="utf-8"))
    no_acp = tomllib.loads(no_acp_path.read_text(encoding="utf-8"))
    if "acp" not in full:
        raise AssertionError("full g-xTB parameter lacks the global ACP table")
    if "acp" in no_acp:
        raise AssertionError("No-ACP parameter still activates the global ACP table")
    full_without_acp = dict(full)
    del full_without_acp["acp"]
    if no_acp != full_without_acp:
        raise AssertionError(
            "No-ACP parameter changes content beyond the global ACP table"
        )
    ih_input = (RAW / "inputs/k222-no-acp/Ih/input.inp").read_text(
        encoding="utf-8", errors="replace"
    )
    if "PARAM " not in ih_input or "gxtb_no_acp.toml" not in ih_input:
        raise AssertionError("archived No-ACP input does not select the parameter file")
    return {
        "semantic_delta": "global ACP table removed only",
        "full_parameter_sha256": hashlib.sha256(full_path.read_bytes()).hexdigest(),
        "no_acp_parameter_sha256": hashlib.sha256(
            no_acp_path.read_bytes()
        ).hexdigest(),
    }


def compute() -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    no_acp_parameter = verify_no_acp_parameter()
    final = cp2k_set("cp2k_final_k222")
    old = cp2k_set("cp2k_old_d933_k222")
    cli = cli_set("cli_final_k222")
    tight_cli = {
        phase: cli_energy(RAW / "cli_tight_scc" / phase / "result.json")
        for phase in ("Ih", "VII")
    }
    tight_native_vii = cp2k_energy(
        RAW / "cp2k_native_tight_scc" / "VII" / "cp2k.out"
    )

    parity_rows: list[dict[str, object]] = []
    absolute_deltas = []
    relative_deltas = []
    for phase in PHASES:
        delta = final[phase] - cli[phase]
        absolute_deltas.append(delta)
        final_rel = 0.0 if phase == "Ih" else relative(final, phase)
        cli_rel = 0.0 if phase == "Ih" else relative(cli, phase)
        relative_delta = final_rel - cli_rel
        relative_deltas.append(relative_delta)
        parity_rows.append({
            "phase": phase,
            "cp2k_native_ha": final[phase],
            "cli_bvk_primitive_ha": cli[phase],
            "cp2k_minus_cli_ha": delta,
            "cp2k_relative_kj_mol": final_rel,
            "cli_relative_kj_mol": cli_rel,
            "relative_delta_kj_mol": relative_delta,
        })

    response_rows: list[dict[str, object]] = []
    old_errors = []
    final_errors = []
    for phase in PHASES[1:]:
        old_rel = relative(old, phase)
        final_rel = relative(final, phase)
        old_error = old_rel - REFERENCE[phase]
        final_error = final_rel - REFERENCE[phase]
        old_errors.append(old_error)
        final_errors.append(final_error)
        response_rows.append({
            "phase": phase,
            "reference_kj_mol": REFERENCE[phase],
            "pre_response_relative_kj_mol": old_rel,
            "final_relative_kj_mol": final_rel,
            "pre_response_error_kj_mol": old_error,
            "final_error_kj_mol": final_error,
            "absolute_error_improvement_kj_mol": abs(old_error) - abs(final_error),
        })

    old_sentinel = {
        phase: cp2k_energy(RAW / "cp2k_pre_response_sentinels" / phase / "cp2k.out")
        for phase in ("Ih", "XVII")
    }
    smear = {
        phase: cp2k_energy(RAW / "cp2k_final_smear300_sentinels" / phase / "cp2k.out")
        for phase in ("Ih", "XVII")
    }
    no_acp_before = {
        "Ih": cp2k_energy(RAW / "no_acp_ih/before/cp2k.out"),
        "XVII": cp2k_energy(RAW / "no_acp_xvii/before/cp2k.out"),
    }
    no_acp_final = {
        "Ih": cp2k_energy(RAW / "no_acp_ih/final/cp2k.out"),
        "XVII": cp2k_energy(RAW / "no_acp_xvii/final/cp2k.out"),
    }
    no_acp_before_rel = relative(no_acp_before, "XVII")
    no_acp_final_rel = relative(no_acp_final, "XVII")
    full_before_rel = relative(old, "XVII")
    full_final_rel = relative(final, "XVII")
    full_shift = full_final_rel - full_before_rel
    no_acp_shift = no_acp_final_rel - no_acp_before_rel

    ablation_rows = [
        {
            "gate": "exact_pre_response_vs_old_d933",
            "metric": "max_abs_energy_delta_ha_Ih_XVII",
            "value": max(abs(old_sentinel[p] - old[p]) for p in old_sentinel),
        },
        {
            "gate": "final_no_smear_vs_300K",
            "metric": "max_abs_energy_delta_ha_Ih_XVII",
            "value": max(abs(smear[p] - final[p]) for p in smear),
        },
        {
            "gate": "full_model_response_fix",
            "metric": "XVII_minus_Ih_relative_shift_kj_mol",
            "value": full_shift,
        },
        {
            "gate": "no_ACP_response_fix",
            "metric": "XVII_minus_Ih_relative_shift_kj_mol",
            "value": no_acp_shift,
        },
        {
            "gate": "no_ACP_fraction",
            "metric": "percent_of_full_relative_shift",
            "value": 100.0 * no_acp_shift / full_shift,
        },
    ]

    old_mae = sum(map(abs, old_errors)) / len(old_errors)
    final_mae = sum(map(abs, final_errors)) / len(final_errors)
    summary: dict[str, object] = {
        "mesh": "2x2x2",
        "phase_count_including_Ih": len(PHASES),
        "relative_phase_count": len(PHASES) - 1,
        "cli_native_max_abs_energy_delta_ha": max(map(abs, absolute_deltas)),
        "cli_native_rms_energy_delta_ha": math.sqrt(
            sum(value * value for value in absolute_deltas) / len(absolute_deltas)
        ),
        "cli_native_max_abs_relative_delta_kj_mol": max(map(abs, relative_deltas)),
        "tight_cli_max_abs_energy_change_ha": max(
            abs(tight_cli[phase] - cli[phase]) for phase in tight_cli
        ),
        "tight_native_VII_energy_change_ha": tight_native_vii - final["VII"],
        "tight_native_cli_VII_energy_delta_ha": (
            tight_native_vii - tight_cli["VII"]
        ),
        "tight_cli_native_VII_energy_delta_ha": final["VII"] - tight_cli["VII"],
        "tight_cli_native_VII_relative_delta_kj_mol": (
            relative(final, "VII") - relative(tight_cli, "VII")
        ),
        "pre_response_mae_kj_mol": old_mae,
        "final_mae_kj_mol": final_mae,
        "mae_improvement_kj_mol": old_mae - final_mae,
        "mae_improvement_percent": 100.0 * (old_mae - final_mae) / old_mae,
        "exact_pre_response_vs_old_d933_max_abs_energy_delta_ha": ablation_rows[0]["value"],
        "final_no_smear_vs_300K_max_abs_energy_delta_ha": ablation_rows[1]["value"],
        "full_model_XVII_Ih_response_shift_kj_mol": full_shift,
        "no_ACP_XVII_Ih_response_shift_kj_mol": no_acp_shift,
        "no_ACP_fraction_percent": 100.0 * no_acp_shift / full_shift,
        "no_ACP_parameter_semantic_delta": no_acp_parameter["semantic_delta"],
        "no_ACP_full_parameter_sha256": no_acp_parameter[
            "full_parameter_sha256"
        ],
        "no_ACP_parameter_sha256": no_acp_parameter["no_acp_parameter_sha256"],
        "final_cp2k_binary_sha256": first_manifest_hash(
            ROOT / "provenance/final/formal-binaries.sha256"
        ),
    }
    return parity_rows, response_rows, ablation_rows, summary


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_summary(summary: dict[str, object]) -> None:
    close(float(summary["cli_native_max_abs_energy_delta_ha"]),
          1.052083007380133e-7, 2.0e-13, "CLI/native maximum energy delta")
    if float(summary["cli_native_rms_energy_delta_ha"]) > 3.2e-8:
        raise AssertionError("CLI/native RMS energy delta exceeds 3.2e-8 Ha")
    if float(summary["cli_native_max_abs_relative_delta_kj_mol"]) > 3.0e-5:
        raise AssertionError("CLI/native relative-energy delta exceeds 3.0e-5 kJ/mol")
    if float(summary["tight_cli_max_abs_energy_change_ha"]) > 1.0e-10:
        raise AssertionError("tight CLI SCC changes a primitive-cell energy by more than 1e-10 Ha")
    if abs(float(summary["tight_native_VII_energy_change_ha"])) > 2.0e-13:
        raise AssertionError("tight native SCC changes the ice-VII energy by more than 2e-13 Ha")
    close(float(summary["tight_native_cli_VII_energy_delta_ha"]),
          -1.051994331646711e-7, 2.0e-13,
          "tight native/tight CLI ice-VII energy delta")
    close(float(summary["tight_cli_native_VII_energy_delta_ha"]),
          -1.051995468515088e-7, 2.0e-13, "tight CLI/native ice-VII energy delta")
    close(float(summary["tight_cli_native_VII_relative_delta_kj_mol"]),
          -2.143273491129e-5, 2.0e-10, "tight CLI/native ice-VII relative delta")
    close(float(summary["pre_response_mae_kj_mol"]),
          90.892218655178, 2.0e-9, "pre-response MAE")
    close(float(summary["final_mae_kj_mol"]),
          88.681375103723, 2.0e-9, "final MAE")
    close(float(summary["mae_improvement_kj_mol"]),
          2.210843551455, 3.0e-9, "MAE improvement")
    if float(summary["exact_pre_response_vs_old_d933_max_abs_energy_delta_ha"]) > 1.0e-12:
        raise AssertionError("exact pre-response build does not reproduce old d933 sentinels")
    if float(summary["final_no_smear_vs_300K_max_abs_energy_delta_ha"]) > 1.0e-12:
        raise AssertionError("explicit 300 K changes the final Ih/XVII energies")
    close(float(summary["full_model_XVII_Ih_response_shift_kj_mol"]),
          5.819640799, 3.0e-9, "full-model response shift")
    close(float(summary["no_ACP_XVII_Ih_response_shift_kj_mol"]),
          5.232844993338, 3.0e-9, "no-ACP response shift")
    if summary["final_cp2k_binary_sha256"] != (
        "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
    ):
        raise AssertionError("unexpected final CP2K binary hash")


def compare_written_summary(computed: dict[str, object]) -> None:
    recorded = json.loads((ROOT / "summary.json").read_text())
    if set(recorded) != set(computed):
        raise AssertionError("summary.json key set differs from recomputed summary")
    for key, value in computed.items():
        if isinstance(value, float):
            close(float(recorded[key]), value, 5.0e-13, f"summary.json {key}")
        elif recorded[key] != value:
            raise AssertionError(f"summary.json {key}: {recorded[key]!r} != {value!r}")


def compare_written_csv(path: Path, computed: list[dict[str, object]]) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        recorded = list(csv.DictReader(handle))
    if len(recorded) != len(computed):
        raise AssertionError(f"{path.name}: row count differs from recomputed data")
    for index, (actual_row, expected_row) in enumerate(zip(recorded, computed, strict=True)):
        if list(actual_row) != list(expected_row):
            raise AssertionError(f"{path.name}: columns differ")
        for key, expected in expected_row.items():
            actual = actual_row[key]
            if isinstance(expected, float):
                close(float(actual), expected, 5.0e-13,
                      f"{path.name} row {index + 1} column {key}")
            elif actual != str(expected):
                raise AssertionError(
                    f"{path.name} row {index + 1} column {key}: "
                    f"{actual!r} != {expected!r}"
                )


def verify_hardening() -> None:
    root = ROOT / "hardening_validation"

    def require_zero(path: Path) -> None:
        if path.read_text(encoding="utf-8").strip() != "0":
            raise AssertionError(f"nonzero validation status: {path}")

    require_zero(root / "unit_tests/kpoint_restart_transfer.rc")
    require_zero(root / "regtest_final/exit_status")
    for path in sorted((root / "focused_tests").glob("*/exit_status")):
        require_zero(path)

    focused = dict(
        line.split("=", 1)
        for line in (root / "focused_tests/summary.txt").read_text(
            encoding="utf-8"
        ).splitlines()
        if "=" in line
    )
    if focused.get("status") != "PASS":
        raise AssertionError("focused hardening gates did not pass")
    close(float(focused["si_general_shifted_delta_hartree"]), 0.0, 1.0e-14,
          "GENERAL/shifted-grid energy delta")
    close(float(focused["h2o_mixer_delta_hartree"]), 0.0, 1.0e-14,
          "density/Fock mixer energy delta")

    final_log = (root / "regtest_final/regtest.log").read_text(
        encoding="utf-8", errors="replace"
    )
    for marker in (
        "Number of FAILED  tests 0",
        "Number of WRONG   tests 0",
        "Number of CORRECT tests 78",
        "Total number of   tests 78",
        "Status: OK",
    ):
        if marker not in final_log:
            raise AssertionError(f"missing final regression marker: {marker}")

    pre_log = (root / "regtest_pre_refresh/regtest.log").read_text(
        encoding="utf-8", errors="replace"
    )
    for marker in (
        "Number of FAILED  tests 0",
        "Number of WRONG   tests 17",
        "Number of CORRECT tests 61",
        "Total number of   tests 78",
    ):
        if marker not in pre_log:
            raise AssertionError(f"missing pre-refresh regression marker: {marker}")

    for phase in ("Ih", "VII"):
        tight_output = (RAW / "cli_tight_scc" / phase / "process.out").read_text(
            encoding="utf-8", errors="replace"
        )
        for marker in (
            "energy convergence             1.0000000000000E-10 Eh",
            "density convergence            2.0000000000000E-09 e",
            "[Info] JSON dump of results written to 'result.json'",
        ):
            if marker not in tight_output:
                raise AssertionError(f"missing tight-CLI marker for {phase}: {marker}")

    tight_native_root = RAW / "cp2k_native_tight_scc" / "VII"
    tight_native_input = (tight_native_root / "input.inp").read_text(
        encoding="utf-8", errors="replace"
    )
    tight_native_output = (tight_native_root / "cp2k.out").read_text(
        encoding="utf-8", errors="replace"
    )
    for marker in ("ACCURACY 0.0001", "EPS_SCF 1.0E-12"):
        if marker not in tight_native_input:
            raise AssertionError(f"missing tight-native input marker: {marker}")
    for marker in ("PROGRAM ENDED AT", "eps_scf:                                        1.00E-12"):
        if marker not in tight_native_output:
            raise AssertionError(f"missing tight-native output marker: {marker}")
    if (tight_native_root / "exit_status").read_text().strip() != "0":
        raise AssertionError("tight-native ice-VII run has nonzero exit status")
    if "expected_cpu=93 allowed=93" not in (
        tight_native_root / "affinity_preexec.txt"
    ).read_text():
        raise AssertionError("tight-native ice-VII affinity proof is invalid")


def write_manifest() -> None:
    manifest = ROOT / "SHA256SUMS"
    paths = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file() and path != manifest
    )
    lines = []
    for path in paths:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(ROOT)}")
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def verify_manifest() -> None:
    manifest = ROOT / "SHA256SUMS"
    recorded: dict[Path, str] = {}
    for raw_line in manifest.read_text(encoding="utf-8").splitlines():
        digest, relative = raw_line.split(None, 1)
        path = ROOT / relative.strip()
        recorded[path] = digest
    actual = {
        path for path in ROOT.rglob("*")
        if path.is_file() and path != manifest
    }
    if set(recorded) != actual:
        missing = sorted(str(path.relative_to(ROOT)) for path in actual - set(recorded))
        stale = sorted(str(path.relative_to(ROOT)) for path in set(recorded) - actual)
        raise AssertionError(f"manifest coverage differs: missing={missing}, stale={stale}")
    for path, expected in recorded.items():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected:
            raise AssertionError(f"SHA-256 mismatch: {path.relative_to(ROOT)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="refresh derived tables")
    args = parser.parse_args()
    parity_rows, response_rows, ablation_rows, summary = compute()
    verify_summary(summary)
    if args.write:
        write_csv(ROOT / "cli_native_k222.csv", parity_rows)
        write_csv(ROOT / "response_fix_k222.csv", response_rows)
        write_csv(ROOT / "ablation_summary.csv", ablation_rows)
        (ROOT / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    compare_written_csv(ROOT / "cli_native_k222.csv", parity_rows)
    compare_written_csv(ROOT / "response_fix_k222.csv", response_rows)
    compare_written_csv(ROOT / "ablation_summary.csv", ablation_rows)
    compare_written_summary(summary)
    verify_hardening()
    if args.write:
        write_manifest()
    verify_manifest()
    print(json.dumps(summary, indent=2, sort_keys=True))
    print("CP2K response-fix, CLI/native, and hardening validation: pass")


if __name__ == "__main__":
    main()
