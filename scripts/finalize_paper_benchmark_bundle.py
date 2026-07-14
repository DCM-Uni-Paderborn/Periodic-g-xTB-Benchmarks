#!/usr/bin/env python3
"""Create one fail-closed, paper-facing DMC13/X23b/LC12 result bundle.

The three benchmark-specific finalizers remain the scientific publication
boundaries.  This script only accepts their completed, internally consistent
JSON/CSV pairs and translates them into a common table plus TeX number macros.
It deliberately emits nothing when any benchmark is incomplete.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Iterable, Mapping


METHODS = ("GFN1", "GFN2", "GXTB")
COMPARISON_METRICS = ("MAE", "RMSE", "MaxAE")
OUTPUT_STEM = "gxtb_periodic_benchmark_summary"
CSV_FIELDS = (
    "benchmark",
    "quantity",
    "scope",
    "method_id",
    "method_label",
    "N",
    "ME",
    "MAE",
    "RMSE",
    "MaxAE",
    "unit",
    "calculation",
    "mesh",
    "status",
    "source_json_sha256",
    "source_csv_sha256",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"required publication artifact is missing: {path}")
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read publication artifact {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"publication artifact is not a JSON object: {path}")
    return payload


def finite(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{label} is not numeric")
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{label} is not numeric") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite")
    return result


def integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} is not an integer")
    try:
        result = int(str(value))
    except ValueError as error:
        raise ValueError(f"{label} is not an integer") from error
    if str(result) != str(value) and not (
        isinstance(value, float) and value.is_integer()
    ):
        raise ValueError(f"{label} is not an exact integer")
    if result <= 0:
        raise ValueError(f"{label} must be positive")
    return result


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"required publication artifact is missing: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"publication CSV is empty: {path}")
    return rows


def close(left: float, right: float, label: str, tolerance: float = 5.0e-10) -> None:
    if not math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance):
        raise ValueError(f"{label} differs between JSON and CSV: {left} != {right}")


def source_record(path: Path, root: Path) -> dict[str, object]:
    try:
        relative = path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        relative = str(path.resolve())
    return {"path": relative, "sha256": sha256(path), "size_bytes": path.stat().st_size}


def metric_row(
    *,
    benchmark: str,
    quantity: str,
    scope: str,
    method_id: str,
    method_label: str,
    n: int,
    metrics: Mapping[str, object],
    unit: str,
    calculation: str,
    mesh: str,
    status: str,
    json_sha: str,
    csv_sha: str,
) -> dict[str, object]:
    return {
        "benchmark": benchmark,
        "quantity": quantity,
        "scope": scope,
        "method_id": method_id,
        "method_label": method_label,
        "N": n,
        **{name: finite(metrics.get(name), f"{benchmark}/{method_id}/{quantity}/{name}") for name in ("ME", "MAE", "RMSE", "MaxAE")},
        "unit": unit,
        "calculation": calculation,
        "mesh": mesh,
        "status": status,
        "source_json_sha256": json_sha,
        "source_csv_sha256": csv_sha,
    }


def validate_dmc(root: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    data = root / "DMC-ICE13" / "data"
    json_path = data / "dmc_ice13_gfn_gxtb_phasewise_summary.json"
    csv_path = data / "dmc_ice13_gfn_gxtb_phasewise_summary.csv"
    payload = read_json(json_path)
    rows = read_csv(csv_path)
    if payload.get("benchmark") != "DMC-ICE13" or payload.get("status") != "phasewise_kpoint_converged":
        raise ValueError("DMC-ICE13 publication artifact is not phase-wise converged")
    if payload.get("n_nonreference_phases") != 12:
        raise ValueError("DMC-ICE13 publication artifact does not contain 12 comparisons")
    fixed_contract = payload.get("fixed_k333_same_mesh_comparison")
    if (
        not isinstance(fixed_contract, dict)
        or fixed_contract.get("mesh") != "k333"
        or fixed_contract.get("not_a_phasewise_converged_result") is not True
    ):
        raise ValueError("DMC-ICE13 publication artifact lacks the fixed-k333 comparator")
    methods = payload.get("methods")
    if not isinstance(methods, dict) or set(methods) != set(METHODS):
        raise ValueError("DMC-ICE13 publication artifact lacks the three methods")
    by_method = {row.get("method_id", ""): row for row in rows}
    if set(by_method) != set(METHODS):
        raise ValueError("DMC-ICE13 publication CSV lacks the three methods")
    json_sha, csv_sha = sha256(json_path), sha256(csv_path)
    output: list[dict[str, object]] = []
    for method in METHODS:
        record = methods[method]
        if not isinstance(record, dict) or record.get("status") != "phasewise_kpoint_converged":
            raise ValueError(f"DMC-ICE13/{method} is not phase-wise converged")
        metrics = record.get("metrics_kjmol_per_h2o")
        if not isinstance(metrics, dict):
            raise ValueError(f"DMC-ICE13/{method} metrics are missing")
        csv_row = by_method[method]
        for metric in ("ME", "MAE", "RMSE", "MaxAE"):
            close(
                finite(metrics.get(metric), f"DMC-ICE13/{method}/{metric}"),
                finite(csv_row.get(f"{metric}_kJmol_per_H2O"), f"DMC CSV/{method}/{metric}"),
                f"DMC-ICE13/{method}/{metric}",
            )
        n = integer(record.get("n_nonreference_phases"), f"DMC-ICE13/{method}/N")
        if n != 12 or integer(csv_row.get("N_nonreference_phases"), f"DMC CSV/{method}/N") != n:
            raise ValueError(f"DMC-ICE13/{method} coverage mismatch")
        output.append(
            metric_row(
                benchmark="DMC-ICE13",
                quantity="relative_energy",
                scope="phasewise_kpoint_converged",
                method_id=method,
                method_label=str(record.get("method_label", method)),
                n=n,
                metrics=metrics,
                unit="kJ mol^-1 per H2O",
                calculation="same-mesh-Ih relative single points",
                mesh="phase-wise adaptive; |delta(N,N-1)| <= 0.05",
                status="publication_ready",
                json_sha=json_sha,
                csv_sha=csv_sha,
            )
        )
        fixed = record.get("fixed_k333_same_mesh_comparison")
        if not isinstance(fixed, dict) or fixed.get("mesh") != "k333":
            raise ValueError(f"DMC-ICE13/{method} fixed-k333 comparator is missing")
        fixed_status = str(fixed.get("status", ""))
        expected_status = (
            "numerically_unconverged_same_mesh_comparator"
            if method == "GXTB"
            else "same_mesh_comparator"
        )
        if fixed_status != expected_status or fixed.get("phasewise_kpoint_converged_value") is not False:
            raise ValueError(f"DMC-ICE13/{method} fixed-k333 status is invalid")
        fixed_metrics = fixed.get("metrics_kjmol_per_h2o")
        if not isinstance(fixed_metrics, dict):
            raise ValueError(f"DMC-ICE13/{method} fixed-k333 metrics are missing")
        if csv_row.get("fixed_k333_status") != fixed_status:
            raise ValueError(f"DMC-ICE13/{method} fixed-k333 status differs between JSON and CSV")
        for metric in ("ME", "MAE", "RMSE", "MaxAE"):
            close(
                finite(fixed_metrics.get(metric), f"DMC-ICE13/{method}/fixed-k333/{metric}"),
                finite(csv_row.get(f"fixed_k333_{metric}_kJmol_per_H2O"), f"DMC CSV/{method}/fixed-k333/{metric}"),
                f"DMC-ICE13/{method}/fixed-k333/{metric}",
            )
        output.append(
            metric_row(
                benchmark="DMC-ICE13",
                quantity="relative_energy",
                scope="fixed_k333_same_mesh_comparator",
                method_id=method,
                method_label=str(record.get("method_label", method)),
                n=n,
                metrics=fixed_metrics,
                unit="kJ mol^-1 per H2O",
                calculation="same-mesh-Ih relative single points",
                mesh="k333",
                status=fixed_status,
                json_sha=json_sha,
                csv_sha=csv_sha,
            )
        )
    return output, {"json": source_record(json_path, root), "csv": source_record(csv_path, root)}


def validate_x23b(root: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    data = root / "X23b" / "data"
    json_path = data / "x23b_gfn_gxtb_paper_summary.json"
    csv_path = data / "x23b_gfn_gxtb_paper_summary.csv"
    payload = read_json(json_path)
    csv_rows = read_csv(csv_path)
    if payload.get("benchmark") != "X23b" or payload.get("publication_status") != "publication_ready":
        raise ValueError("X23b publication artifact is not publication-ready")
    if payload.get("methods") != list(METHODS):
        raise ValueError("X23b publication artifact lacks the three methods")
    coverage = payload.get("coverage")
    if not isinstance(coverage, dict) or coverage.get("exact_common_coverage") is not True or coverage.get("common") != 23:
        raise ValueError("X23b publication artifact does not have exact 23/23 coverage")
    csv_sha = sha256(csv_path)
    if payload.get("publication_csv_sha256") != csv_sha:
        raise ValueError("X23b publication CSV hash differs from its JSON lineage")
    summary = payload.get("summary")
    if not isinstance(summary, list) or len(summary) != 6:
        raise ValueError("X23b publication summary must contain six rows")
    csv_by_key = {(row.get("method", ""), row.get("quantity", "")): row for row in csv_rows}
    expected_keys = {(method, quantity) for method in METHODS for quantity in ("lattice_energy_kJmol", "volume_error_percent")}
    if set(csv_by_key) != expected_keys:
        raise ValueError("X23b publication CSV has incomplete or duplicate method/quantity rows")
    json_sha = sha256(json_path)
    output: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for record in summary:
        if not isinstance(record, dict):
            raise ValueError("X23b summary contains a non-object row")
        method, quantity = str(record.get("method", "")), str(record.get("quantity", ""))
        key = (method, quantity)
        if key not in expected_keys or key in seen:
            raise ValueError("X23b summary has an unexpected or duplicate row")
        seen.add(key)
        csv_row = csv_by_key[key]
        for field in ("method_label", "calculation", "mesh"):
            if str(record.get(field, "")) != str(csv_row.get(field, "")):
                raise ValueError(
                    f"X23b/{method}/{quantity}/{field} differs between JSON and CSV"
                )
        for metric in ("ME", "MAE", "RMSE", "MaxAE"):
            close(finite(record.get(metric), f"X23b/{method}/{quantity}/{metric}"), finite(csv_row.get(metric), f"X23b CSV/{method}/{quantity}/{metric}"), f"X23b/{method}/{quantity}/{metric}")
        n = integer(record.get("N"), f"X23b/{method}/{quantity}/N")
        if n != 23 or integer(csv_row.get("N"), f"X23b CSV/{method}/{quantity}/N") != n:
            raise ValueError(f"X23b/{method}/{quantity} coverage mismatch")
        unit = "kJ mol^-1" if quantity == "lattice_energy_kJmol" else "percent"
        output.append(
            metric_row(
                benchmark="X23b",
                quantity=quantity,
                scope="full_common_set",
                method_id=method,
                method_label=str(record.get("method_label", method)),
                n=n,
                metrics=record,
                unit=unit,
                calculation=str(record.get("calculation", "")),
                mesh=str(record.get("mesh", "")),
                status="publication_ready",
                json_sha=json_sha,
                csv_sha=csv_sha,
            )
        )
    if seen != expected_keys:
        raise ValueError("X23b summary lacks required rows")
    return output, {"json": source_record(json_path, root), "csv": source_record(csv_path, root)}


def validate_lc12(root: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    data = root / "Goldzak12" / "data"
    json_path = data / "lc12_gfn_gxtb_paper_summary.json"
    csv_path = data / "lc12_gfn_gxtb_paper_summary.csv"
    payload = read_json(json_path)
    csv_rows = read_csv(csv_path)
    status = str(payload.get("status", ""))
    if payload.get("benchmark") != "LC12 (Goldzak12)" or status not in {"publication_ready", "publication_ready_reduced_coverage"}:
        raise ValueError("LC12 publication artifact is not publication-ready")
    methods = payload.get("methods")
    if not isinstance(methods, dict) or set(methods) != set(METHODS):
        raise ValueError("LC12 publication artifact lacks the three methods")
    protocol = payload.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("LC12 protocol is missing")
    common_n = integer(protocol.get("common_subset_count"), "LC12 common-subset N")
    summary = payload.get("summary_rows")
    if not isinstance(summary, list) or len(summary) != 6:
        raise ValueError("LC12 publication summary must contain six scope rows")
    csv_record = payload.get("paper_summary_csv")
    csv_sha = sha256(csv_path)
    if not isinstance(csv_record, dict) or csv_record.get("sha256") != csv_sha:
        raise ValueError("LC12 publication CSV hash differs from its JSON lineage")
    expected_keys = {(method, scope) for method in METHODS for scope in ("method_available_coverage", "three_method_common_subset")}
    csv_by_key = {(row.get("method_id", ""), row.get("scope", "")): row for row in csv_rows}
    if set(csv_by_key) != expected_keys:
        raise ValueError("LC12 publication CSV has incomplete or duplicate method/scope rows")
    json_sha = sha256(json_path)
    output: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    common_systems: str | None = None
    for record in summary:
        if not isinstance(record, dict):
            raise ValueError("LC12 summary contains a non-object row")
        method, scope = str(record.get("method_id", "")), str(record.get("scope", ""))
        key = (method, scope)
        if key not in expected_keys or key in seen:
            raise ValueError("LC12 summary has an unexpected or duplicate row")
        seen.add(key)
        csv_row = csv_by_key[key]
        for field in ("method_label", "systems", "eos_mesh", "result_mesh"):
            if str(record.get(field, "")) != str(csv_row.get(field, "")):
                raise ValueError(
                    f"LC12/{method}/{scope}/{field} differs between JSON and CSV"
                )
        n = integer(record.get("n_systems"), f"LC12/{method}/{scope}/N")
        if integer(csv_row.get("n_systems"), f"LC12 CSV/{method}/{scope}/N") != n:
            raise ValueError(f"LC12/{method}/{scope} coverage mismatch")
        if scope == "three_method_common_subset":
            if n != common_n:
                raise ValueError("LC12 common-subset coverage differs among methods")
            systems = str(record.get("systems", ""))
            if common_systems is None:
                common_systems = systems
            elif systems != common_systems:
                raise ValueError("LC12 common-subset systems differ among methods")
        for quantity, prefix, unit in (
            ("lattice_constant", "lattice", "angstrom"),
            ("cohesive_energy", "cohesive", "eV atom^-1"),
        ):
            metrics = {metric: record.get(f"{prefix}_{metric}_{'A' if prefix == 'lattice' else 'eV_per_atom'}") for metric in ("ME", "MAE", "RMSE", "MaxAE")}
            for metric in metrics:
                close(
                    finite(metrics[metric], f"LC12/{method}/{scope}/{quantity}/{metric}"),
                    finite(csv_row.get(f"{prefix}_{metric}_{'A' if prefix == 'lattice' else 'eV_per_atom'}"), f"LC12 CSV/{method}/{scope}/{quantity}/{metric}"),
                    f"LC12/{method}/{scope}/{quantity}/{metric}",
                )
            output.append(
                metric_row(
                    benchmark="LC12",
                    quantity=quantity,
                    scope=scope,
                    method_id=method,
                    method_label=str(record.get("method_label", method)),
                    n=n,
                    metrics=metrics,
                    unit=unit,
                    calculation="equation of state and cohesive-energy single point",
                    mesh=f"EOS {record.get('eos_mesh', '')}; result {record.get('result_mesh', '')}",
                    status=status,
                    json_sha=json_sha,
                    csv_sha=csv_sha,
                )
            )
    if seen != expected_keys:
        raise ValueError("LC12 summary lacks required rows")
    return output, {"json": source_record(json_path, root), "csv": source_record(csv_path, root)}


def csv_text(rows: Iterable[Mapping[str, object]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row[field] for field in CSV_FIELDS})
    return stream.getvalue()


def tex_token(value: str) -> str:
    pieces = re.findall(r"[A-Za-z0-9]+", value)
    return "".join(piece[0].upper() + piece[1:] for piece in pieces)


def build_gxtb_gfn2_comparisons(
    rows: Iterable[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Build only like-for-like g-XTB/GFN2 paper comparisons.

    LC12 method-available coverage can contain different solids and is
    intentionally excluded. Its three-method common subset is the only LC12
    scope for which a ratio or percentage change has an unambiguous meaning.
    """

    grouped: dict[tuple[str, str, str], dict[str, Mapping[str, object]]] = {}
    for row in rows:
        key = (str(row["benchmark"]), str(row["quantity"]), str(row["scope"]))
        method = str(row["method_id"])
        if method in grouped.setdefault(key, {}):
            raise ValueError(f"duplicate method row for comparison: {key}/{method}")
        grouped[key][method] = row

    comparisons: list[dict[str, object]] = []
    for (benchmark, quantity, scope), methods in sorted(grouped.items()):
        if benchmark == "LC12" and scope != "three_method_common_subset":
            continue
        if set(methods) != set(METHODS):
            raise ValueError(
                f"cannot build g-XTB/GFN2 comparison for {benchmark}/{quantity}/{scope}"
            )
        gxtb, gfn2 = methods["GXTB"], methods["GFN2"]
        n_gxtb = integer(gxtb["N"], "g-XTB comparison N")
        n_gfn2 = integer(gfn2["N"], "GFN2 comparison N")
        if n_gxtb != n_gfn2:
            raise ValueError(
                f"g-XTB/GFN2 coverage differs for {benchmark}/{quantity}/{scope}"
            )
        if str(gxtb["unit"]) != str(gfn2["unit"]):
            raise ValueError(
                f"g-XTB/GFN2 units differ for {benchmark}/{quantity}/{scope}"
            )
        record: dict[str, object] = {
            "benchmark": benchmark,
            "quantity": quantity,
            "scope": scope,
            "methods": ["GXTB", "GFN2"],
            "N": n_gxtb,
            "unit": str(gxtb["unit"]),
            "ME_delta_GXTB_minus_GFN2": finite(gxtb["ME"], "g-XTB ME")
            - finite(gfn2["ME"], "GFN2 ME"),
        }
        for metric in COMPARISON_METRICS:
            gxtb_value = finite(gxtb[metric], f"g-XTB {metric}")
            gfn2_value = finite(gfn2[metric], f"GFN2 {metric}")
            if gfn2_value <= 0.0:
                raise ValueError(
                    f"GFN2 {metric} must be positive for {benchmark}/{quantity}/{scope}"
                )
            record[f"{metric}_delta_GXTB_minus_GFN2"] = gxtb_value - gfn2_value
            record[f"{metric}_ratio_GXTB_over_GFN2"] = gxtb_value / gfn2_value
            record[f"{metric}_percent_change_GXTB_vs_GFN2"] = (
                100.0 * (gxtb_value - gfn2_value) / gfn2_value
            )
        comparisons.append(record)
    return comparisons


def tex_macros(
    rows: Iterable[Mapping[str, object]],
    comparisons: Iterable[Mapping[str, object]],
) -> str:
    lines = [
        "% AUTO-GENERATED by scripts/finalize_paper_benchmark_bundle.py",
        "% Do not edit numerical values by hand.",
    ]
    names: set[str] = set()
    for row in rows:
        prefix = "GXTB" + tex_token(str(row["benchmark"])) + tex_token(str(row["quantity"])) + tex_token(str(row["scope"])) + tex_token(str(row["method_id"]))
        for metric in ("N", "ME", "MAE", "RMSE", "MaxAE"):
            name = prefix + metric
            if name in names:
                raise ValueError(f"duplicate generated TeX macro: {name}")
            names.add(name)
            value = str(row[metric]) if metric == "N" else f"{float(row[metric]):.9f}"
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
    for comparison in comparisons:
        prefix = (
            "GXTBComparison"
            + tex_token(str(comparison["benchmark"]))
            + tex_token(str(comparison["quantity"]))
            + tex_token(str(comparison["scope"]))
            + "GXTBvsGFNtwo"
        )
        values = {
            "N": comparison["N"],
            "MEDelta": comparison["ME_delta_GXTB_minus_GFN2"],
        }
        for metric in COMPARISON_METRICS:
            values[f"{metric}Delta"] = comparison[
                f"{metric}_delta_GXTB_minus_GFN2"
            ]
            values[f"{metric}Ratio"] = comparison[
                f"{metric}_ratio_GXTB_over_GFN2"
            ]
            values[f"{metric}PercentChange"] = comparison[
                f"{metric}_percent_change_GXTB_vs_GFN2"
            ]
        for suffix, raw_value in values.items():
            name = prefix + suffix
            if name in names:
                raise ValueError(f"duplicate generated TeX macro: {name}")
            names.add(name)
            value = str(raw_value) if suffix == "N" else f"{float(raw_value):.9f}"
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
    return "\n".join(lines) + "\n"


def build_bundle(root: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    dmc_rows, dmc_sources = validate_dmc(root)
    x23b_rows, x23b_sources = validate_x23b(root)
    lc12_rows, lc12_sources = validate_lc12(root)
    rows = dmc_rows + x23b_rows + lc12_rows
    comparisons = build_gxtb_gfn2_comparisons(rows)
    return (
        {
            "schema_version": 2,
            "status": "publication_ready",
            "title": "Periodic g-xTB benchmark comparison with frozen GFN1/GFN2 baselines",
            "benchmarks": ["DMC-ICE13", "X23b", "LC12"],
            "methods": list(METHODS),
            "sources": {"DMC-ICE13": dmc_sources, "X23b": x23b_sources, "LC12": lc12_sources},
            "rows": rows,
            "gxtb_vs_gfn2_comparisons": comparisons,
        },
        rows,
    )


def finalize(root: Path, output_dir: Path | None = None) -> tuple[Path, Path, Path]:
    root = root.resolve()
    target = (output_dir or root / "paper").resolve()
    csv_path = target / f"{OUTPUT_STEM}.csv"
    json_path = target / f"{OUTPUT_STEM}.json"
    tex_path = target / f"{OUTPUT_STEM}.tex"
    outputs = (csv_path, json_path, tex_path)
    for path in outputs:
        path.unlink(missing_ok=True)
    target.mkdir(parents=True, exist_ok=True)
    temporaries = tuple(path.with_name(f".{path.name}.{os.getpid()}.tmp") for path in outputs)
    for path in temporaries:
        path.unlink(missing_ok=True)
    try:
        bundle, rows = build_bundle(root)
        csv_body = csv_text(rows)
        tex_body = tex_macros(rows, bundle["gxtb_vs_gfn2_comparisons"])
        bundle["generated_outputs"] = {
            "csv_sha256": hashlib.sha256(csv_body.encode()).hexdigest(),
            "tex_sha256": hashlib.sha256(tex_body.encode()).hexdigest(),
        }
        temporaries[0].write_text(csv_body)
        temporaries[1].write_text(json.dumps(bundle, indent=2, sort_keys=True) + "\n")
        temporaries[2].write_text(tex_body)
        for temporary, output in zip(temporaries, outputs, strict=True):
            os.replace(temporary, output)
    except BaseException:
        for path in (*temporaries, *outputs):
            path.unlink(missing_ok=True)
        raise
    return csv_path, json_path, tex_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        outputs = finalize(args.repository_root, args.output_dir)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print("\n".join(str(path) for path in outputs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
