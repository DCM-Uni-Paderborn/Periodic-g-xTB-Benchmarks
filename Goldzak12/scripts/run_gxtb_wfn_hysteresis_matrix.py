#!/usr/bin/env python3
"""Run the reviewed LC12 g-XTB WFN hysteresis matrix in parallel."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import subprocess
import sys
from pathlib import Path

import diagnose_gxtb_wfn_hysteresis as pair_driver
import run_goldzak12_benchmark as base


PAIRS = (
    ("LiH", 1.375, 1.30),
    ("LiH", 1.30, 1.375),
    ("LiH", 1.375, 1.45),
    ("LiH", 1.45, 1.375),
    ("LiH", 1.00, 1.375),
    ("MgO", 0.94, 0.96),
    ("MgO", 0.96, 0.94),
    ("MgO", 0.94, 0.98),
    ("MgO", 0.98, 0.94),
    ("MgO", 0.80, 0.82),
    ("MgO", 0.82, 0.80),
    ("MgO", 0.82, 0.85),
    ("MgO", 0.85, 0.82),
)


def label(solid: str, source: float, target: float) -> str:
    return f"{solid}_{pair_driver.scale_tag(source)}_to_{pair_driver.scale_tag(target)}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-manifest", type=Path, required=True)
    parser.add_argument("--cp2k-source", type=Path, required=True)
    parser.add_argument("--save-tblite-source", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=len(PAIRS))
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--run-label", required=True)
    args = parser.parse_args()
    log_root = base.ROOT / "logs" / "gxtb_wfn_hysteresis" / args.run_label
    log_root.mkdir(parents=True, exist_ok=True)
    pair_script = Path(__file__).with_name("diagnose_gxtb_wfn_hysteresis.py")

    def run(spec: tuple[str, float, float]) -> dict[str, object]:
        solid, source, target = spec
        name = label(solid, source, target)
        log = log_root / f"{name}.log"
        command = [
            sys.executable,
            str(pair_script),
            "--solid",
            solid,
            "--mesh",
            "k444",
            "--source-scale",
            str(source),
            "--target-scale",
            str(target),
            "--campaign-manifest",
            str(args.campaign_manifest),
            "--cp2k-source",
            str(args.cp2k_source),
            "--save-tblite-source",
            str(args.save_tblite_source),
            "--threads",
            str(args.threads),
        ]
        with log.open("w") as handle:
            proc = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
        manifest = (
            base.ROOT
            / "runs"
            / "gxtb_wfn_hysteresis"
            / solid
            / "k444"
            / f"{pair_driver.scale_tag(source)}_to_{pair_driver.scale_tag(target)}"
            / "diagnostic_manifest.json"
        )
        record = {
            "solid": solid,
            "source_scale": source,
            "target_scale": target,
            "return_code": proc.returncode,
            "log": pair_driver.artifact(log),
            "pair_manifest": pair_driver.artifact(manifest) if manifest.exists() else None,
        }
        print(f"done rc={proc.returncode} {name}", flush=True)
        return record

    records: list[dict[str, object]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run, spec) for spec in PAIRS]
        for future in concurrent.futures.as_completed(futures):
            records.append(future.result())
    matrix = {
        "schema_version": 1,
        "diagnostic": "gxtb_wfn_hysteresis_matrix",
        "production_eligible": False,
        "run_label": args.run_label,
        "threads_per_pair": args.threads,
        "parallel_pairs": args.jobs,
        "state_contract": {
            "transferred": "CP2K Bloch-orbital/density WFN guess; q/multipoles reconstructed from density",
            "not_transferred": "internal save_tblite q/dipole/quadrupole state and FDIIS history",
        },
        "pairs": sorted(
            records, key=lambda item: (str(item["solid"]), float(item["source_scale"]), float(item["target_scale"]))
        ),
    }
    manifest = log_root / "matrix_manifest.json"
    base.write_file(manifest, json.dumps(matrix, indent=2, sort_keys=True) + "\n")
    failed = [record for record in records if record["return_code"] != 0]
    print(f"matrix complete: {len(records) - len(failed)}/{len(records)} successful; {manifest}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
